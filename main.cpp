#include <iostream>
#include <boost/numeric/interval.hpp>
#include <string>
#include <functional>
#include <initializer_list>
#include <utility>
#include <vector>
#include <fstream>

// typedef boost::numeric::interval (with custom policy-ing) to Interval to work with regular functions (exp(...), ...)
typedef boost::numeric::interval<double, boost::numeric::interval_lib::policies<boost::numeric::interval_lib::save_state<
        boost::numeric::interval_lib::rounded_transc_std<double>>, boost::numeric::interval_lib::checking_base<double>>>
        Interval;

/**
 * converts a given interval to a string
 * @tparam T underlying numeric type of the argument
 * @tparam Policy policy of the argument
 * @param interval the interval to be stringified
 * @return a string object of the form "[(interval's lower bound), (interval's upper bound)]"
 */
template <class T, class Policy>
std::string to_string(const boost::numeric::interval<T, Policy>& interval) {
    return "[" + std::to_string(interval.lower()) + ", " + std::to_string(interval.upper()) + "]";
}

class IntervalTuple {
private:
    std::vector<Interval> tuple_;

public:
    // anti-pattern with <double>... not generic
    IntervalTuple(std::initializer_list<Interval> l) : tuple_(l) {}
    Interval operator[](std::size_t index) const { // doesn't check bounds...
        assert(index < tuple_.size() && index >= 0);
        return tuple_[index];
    }

    Interval operator[](std::size_t index) { // doesn't check bounds...
        assert(index < tuple_.size() && index >= 0);
        return tuple_[index];
    }
    [[nodiscard]] std::size_t size() const {return tuple_.size();}
    std::size_t size() {return tuple_.size();}

    /**
     * determines whether "this" interval tuple is a subset of a given interval tuple... this function assumes that the
     * interval tuple is a subset iff all the individual intervals in "this" are subsets of the corresponding intervals
     * in "set".
     * @param set the hypothesized interval tuple that is the super set of "this"
     * @return true if "this" is a subset of "set", false otherwise
     */
    bool is_subset(const IntervalTuple& set) {
        assert(set.tuple_.size() == this->tuple_.size());
        for (std::size_t i = 0; i < this->tuple_.size(); ++i) {
            if (!boost::numeric::subset(this->tuple_[i], set[i])) {
                // if any of the individual intervals in "this" are not subsets of the corresponding interval in
                // "super_set", then just return false (return early)
                return false;
            }
        }
        return true;
    }

    /**
     * determines whether "this" interval tuple intersected with another interval tuple is empty or non-empty... this
     * function assumes that two interval tuples have empty intersection iff the intersections of all the corresponding
     * interval elements are empty.
     * @param set the hypothesized interval tuple that is not intersected with "this"
     * @return true if "this" is a does not intersect with "set" at all, false otherwise
     */
    bool is_intersection_empty(const IntervalTuple& set) {
        assert(set.size() == tuple_.size());
        for (std::size_t i = 0; i < set.tuple_.size(); ++i) {
            if (!boost::numeric::empty(boost::numeric::intersect(tuple_[i], set[i]))) {
                return false;
            }
        }
        return true;
    }

    double largest_width() {
        double ret{0.0}; // smallest possible width anyhow, will always be caused to widen on the first interval
        for (auto & i : tuple_) {
            double w_i = boost::numeric::width(i);
            if (w_i > ret) {
                ret = w_i;
            }
        }
        return ret; // note: if the tuple is empty, this will return '0' not NAN or 'nullopt' or whatever...
    }

    /**
     * always bisect along the longest width
     * @return
     */
    std::pair<IntervalTuple, IntervalTuple> bisect() {
        // step 1. figure out which index has the largest width
        double ret{0.0}; // smallest possible width anyhow, will always be caused to widen on the first interval
        std::size_t index{0};
        for (std::size_t i = 0; i < tuple_.size(); ++i) {
            double w_i = boost::numeric::width(tuple_[i]);
            if (w_i > ret) {
                ret = w_i;
                index = i;
            }
        }

        // step 2. bisect along 'longest width direction'
        Interval lower_ith_interval = Interval(tuple_[index].lower(), boost::numeric::median(tuple_[index]));
        Interval upper_ith_interval = Interval(boost::numeric::median(tuple_[index]), tuple_[index].upper());

        // not sure why this is legal c++ code... isn't "tuple_" member var private... if so, why can we access in local var "lower_interval" and "upper_interval"?
        IntervalTuple lower_interval = *this; lower_interval.tuple_[index].set(tuple_[index].lower(), boost::numeric::median(tuple_[index]));
        IntervalTuple upper_interval = *this; upper_interval.tuple_[index].set(boost::numeric::median(tuple_[index]), tuple_[index].upper());

        return std::make_pair(lower_interval, upper_interval);
    }

    friend std::ostream& operator<<(std::ostream& os, const IntervalTuple& tuple) {
        std::string tuple_string;
        for (std::size_t i = 0; i < tuple.size(); ++i) {
            tuple_string += to_string(tuple[i]) + "\n";
        }
        return os << tuple_string;
    }
};


typedef std::function<IntervalTuple(IntervalTuple&)> InclusionFunction;

struct sivia_options {
    double epsilon;
};

// we need an inclusion function, an initial guess x0 (represented by an IntervalTuple in IR^n), and we need a set Y
// (represented by an IntervalTuple in IR^p)... we want to return X- and X+ s.t. X- < X < X+ (where '<' are is subset)
// and X- and X+ are represented by IntervalTuple in IR^n)
// what is an inclusion function ( notationally, [f]([x]) )? it maps [x] to [y] i.e. f: IR^n -> IR^p
std::pair<std::vector<IntervalTuple>, std::vector<IntervalTuple>> sivia(InclusionFunction& f_box, IntervalTuple& X_0, IntervalTuple& Y, sivia_options& options) {

    std::vector<IntervalTuple> X_minus{}; // initialize as empty-set
    std::vector<IntervalTuple> X_plus{}; // initialize as empty-set
    std::vector<IntervalTuple> L = {X_0}; // initialize with initial guess "X_0"

    while (!L.empty()) {
        // while our list is non-empty...
        IntervalTuple Xi = L.back(); L.pop_back();
        IntervalTuple Yi = f_box(Xi);

        if (Yi.is_subset(Y)) {
            X_minus.push_back(Xi);
            X_plus.push_back(Xi);
        } else if (Yi.is_intersection_empty(Y)) {
            continue; // discard Xi
        } else if (Xi.largest_width() < options.epsilon) {
            X_plus.push_back(Xi);
        } else {
            auto [X_lower, X_upper] = Xi.bisect();
            L.push_back(X_lower);
            L.push_back(X_upper);
        }
    }

    return std::make_pair(X_minus, X_plus);
}

void save_soln(std::vector<IntervalTuple> X_minus, std::vector<IntervalTuple> X_plus) {
    std::string X_minus_filename = "xminus.txt";
    std::string X_plus_filename = "xplus.txt";

    // create file obj and open file
    std::ofstream outfile;
    outfile.open(X_minus_filename, std::ios_base::out); // write over previous soln
    for (auto & X_minu : X_minus) {
        double x = X_minu[0].lower();
        double y = X_minu[1].lower();
        double w = boost::numeric::width(X_minu[0]);
        double h = boost::numeric::width(X_minu[1]);
        outfile << x << " " << y << " " << w << " " << h << "\n";
    }

    outfile.close();

    outfile.open(X_plus_filename, std::ios_base::out); // write over previous soln
    for (auto & X_plu : X_plus) {
        double x = X_plu[0].lower();
        double y = X_plu[1].lower();
        double w = boost::numeric::width(X_plu[0]);
        double h = boost::numeric::width(X_plu[1]);
        outfile << x << " " << y << " " << w << " " << h << "\n";
    }
    outfile.close();

    // save every element of each vector, where the Interval Tuple (members of IR^2) is [X_minus[i][0].lower(), X_minus[i][1].lower(), X_minus[i][0].width(), X_minus[i][1].width()]

}

int main() {
    double p1{5}; // unknown parameter 1
    double p2{1}; // unknown parameter 2
    double t1{0}; // time 1
    double t2{5}; // time 2
    double t3{10}; // time 3
    std::function<double(double)> my_system = [&](double t){ return p1 * std::exp(-p2 * t); }; // system dynamics

    // three measurements with fixed uncertainties, which are represented as tolerance bounds (i.e. intervals)
    Interval y1(my_system(0) - 1,my_system(0) + 1);
    Interval y2(my_system(5) - 1,my_system(5) + 1);
    Interval y3(my_system(10) - 1,my_system(10) + 1);

    // print the measurements (y1, y2, y3)
    std::cout << to_string(y1) << std::endl;
    std::cout << to_string(y2) << std::endl;
    std::cout << to_string(y3) << std::endl;
    IntervalTuple y_box = {y1, y2, y2};

    // determine the inclusion function of the system function
    std::function<Interval(double, Interval&, Interval&)> my_system2 = [](double t, Interval& p1, Interval& p2){
        return p1 * exp(p2 * t); // same thing as system function but p1 and p2 are "interval" types now
    };

    // initialize parameter guess
    Interval p1_initial_guess(-1, 7);
    Interval p2_initial_guess(-1, 2);
    auto resp1 = my_system2(0, p1_initial_guess, p2_initial_guess);
    auto resp2 = my_system2(5, p1_initial_guess, p2_initial_guess);
    auto resp3 = my_system2(10, p1_initial_guess, p2_initial_guess);
    std::cout << to_string(resp1) << std::endl;
    std::cout << to_string(resp2) << std::endl;
    std::cout << to_string(resp3) << std::endl;

    IntervalTuple resp_box = {resp1, resp2, resp3};

    std::cout << std::boolalpha << "resp_box is subset of y_box: " << resp_box.is_subset(y_box) << std::endl;
    std::cout << std::boolalpha << "resp_box intersect y_box is empty set: " << resp_box.is_intersection_empty(y_box) << std::endl;
    std::cout << std::boolalpha << "y_box is subset of resp_box: " << y_box.is_subset(resp_box) << std::endl;

    auto [lower, upper] = resp_box.bisect();
    std::cout << "lower: \n" << lower << std::endl;
    std::cout << "upper: \n" << upper << std::endl;


    // this inclusion function maps IR^2 (two parameters in parameter-space) to IR^3 (3 measurements in
    // measurement-space...
    InclusionFunction inclusion_function = [t1, t2, t3](IntervalTuple& x_box) {
        auto y_interval1 = x_box[0] * boost::numeric::exp(-x_box[1] * t1);
        auto y_interval2 = x_box[0] * boost::numeric::exp(-x_box[1] * t2);
        auto y_interval3 = x_box[0] * boost::numeric::exp(-x_box[1] * t3);
        IntervalTuple y_box = {y_interval1, y_interval2, y_interval3};
        return y_box;
    };

    IntervalTuple X_0 = {p1_initial_guess, p2_initial_guess};
    IntervalTuple Y = {y1, y2, y3};
    sivia_options options; options.epsilon = 0.1;
    auto [X_minus, X_plus] = sivia(inclusion_function, X_0, Y, options);

    std::cout << X_minus.size() << std::endl;
    std::cout << X_plus.size() << std::endl;

    save_soln(X_minus, X_plus);

    return 0;
}
