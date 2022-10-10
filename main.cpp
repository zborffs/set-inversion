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
    explicit IntervalTuple(std::vector<Interval>& v) : tuple_(v) {}
    IntervalTuple(std::initializer_list<Interval> l) : tuple_(l) {}

    Interval operator[](std::size_t index) const { // doesn't check bounds...
        assert(index < tuple_.size() && index >= 0);
        return tuple_[index];
    }

    Interval operator[](std::size_t index) { // doesn't check bounds...
        assert(index < tuple_.size() && index >= 0);
        return tuple_[index];
    }

    Interval& at(std::size_t index) {
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

/**
 * contractor for z = x + y
 */
void contractor_plus(Interval& z, Interval& x, Interval& y) {
    z = boost::numeric::intersect(z, x + y);
    x = boost::numeric::intersect(x, z - y);
    y = boost::numeric::intersect(y, z - x);
}

/**
 * contractor for z = exp(x)
 */
void contractor_exp(Interval& z, Interval& x) {
    z = boost::numeric::intersect(z, boost::numeric::exp(x));
    x = boost::numeric::intersect(x, boost::numeric::log(x));
}

/**
 * contractor for z = x * y
 */
void contractor_mult(Interval& z, Interval& x, Interval& y) {
    z = boost::numeric::intersect(z, x * y);
    x = boost::numeric::intersect(x, z / y);
    y = boost::numeric::intersect(y, z / x);
}

/**
 * contractor for z = x * y
 */
void contractor_negate(Interval& z, Interval& x) {
    z = -x;
    x = -z;
}

/**
 * contractor for z = x * y
 */
void contractor_sin(Interval& z, Interval& x) {
    z = boost::numeric::intersect(z, boost::numeric::sin(x));
    x = boost::numeric::intersect(x, boost::numeric::asin(z));
}

void contractor_custom(Interval& x, Interval& y) {
    Interval a(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Interval b(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Interval c(0, std::numeric_limits<double>::max());

    Interval x_old;
    Interval y_old;
    double my_epsilon{0.1};
    do {
        x_old = x;
        y_old = y;
        contractor_mult(a, x, y); // a = x * y
        contractor_sin(b, a); // b = sin(a) = sin(x * y)
        contractor_plus(c, x, b); // c = x + b = x + sin(x * y)
    } while (!(width(x_old) - width(x) < my_epsilon && width(y_old) - width(y) < my_epsilon));
}

// y is constant, because it is not contracting, it's just there to provide bounds on initial
void contractor_system_dynamics(IntervalTuple& x, const IntervalTuple& y, const std::vector<double>& time) {
    std::vector<Interval> a_vec(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        a_vec[i] = Interval(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    }
    IntervalTuple a(a_vec);

    std::vector<Interval> b_vec(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        b_vec[i] = Interval(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    }
    IntervalTuple b(b_vec);

    std::vector<Interval> c_vec(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        c_vec[i] = Interval(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    }
    IntervalTuple c(c_vec);

    std::vector<Interval> d_vec(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        d_vec[i] = Interval(y[i].lower(), y[i].upper());
    }
    IntervalTuple d(d_vec);

    std::vector<Interval> time_vec(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        time_vec[i] = Interval(time[i]);
    }
    IntervalTuple t(time_vec);

    Interval x0_old;
    Interval x1_old;
    double my_epsilon{0.1};
    do {
        x0_old = x[0];
        x1_old = x[1];

        for (std::size_t i = 0; i < a_vec.size(); ++i) {
            contractor_mult(a.at(i), x.at(1), t.at(i));
            contractor_negate(b.at(i), a.at(i));
            contractor_exp(c.at(i),  b.at(i));
            contractor_mult(d.at(i), x.at(0), c.at(i));
        }
    } while (!(width(x0_old) - width(x[0]) < my_epsilon && width(x1_old) - width(x[1]) < my_epsilon));
}


typedef std::function<IntervalTuple(IntervalTuple&)> InclusionFunction;

struct sivia_options {
    double epsilon;
    std::function<void(IntervalTuple&, const IntervalTuple&)> contractor;
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
            continue; // discard Xi... and all sub trees?
        } else if (Xi.largest_width() < options.epsilon) {
            X_plus.push_back(Xi);
        } else {
            options.contractor(Xi, Y);
            auto [X_lower, X_upper] = Xi.bisect();

            // don't add the X_lower or X_upper bounds IntervalTuples if any underlying Interval contains nan as bound
            bool add_lower{true};
            for (std::size_t i = 0; i < X_lower.size(); ++i) {
                if (isnan(X_lower[i].lower()) || isnan(X_lower[i].upper())) {
                    add_lower = false;
                }
            }

            bool add_upper{true};
            for (std::size_t i = 0; i < X_upper.size(); ++i) {
                if (isnan(X_upper[i].lower()) || isnan(X_upper[i].upper())) {
                    add_upper = false;
                }
            }

            if (add_lower) {
                L.push_back(X_lower);
            }

            if (add_upper) {
                L.push_back(X_upper);
            }
        }
    }

    return std::make_pair(X_minus, X_plus);
}

void save_soln(std::vector<IntervalTuple> X_minus, std::vector<IntervalTuple> X_plus) {
    std::string X_minus_filename = "../xminus.txt";
    std::string X_plus_filename = "../xplus.txt";

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
}

int main() {
    double p1{5}; // unknown parameter 1
    double p2{1}; // unknown parameter 2

    // create time variables
    double t_initial{0.1}; // first time
    double t_final{5.1}; // last time
    double sample_time{0.5}; // sample time
    std::size_t num_samples = (std::size_t)(t_final - t_initial) / sample_time;
    std::vector<double> time(num_samples);
    for (std::size_t i = 0; i < time.size(); ++i){
        time[i] = t_initial + (double)i * sample_time;
    }

    std::function<double(double)> f = [&](double t){ return p1 * std::exp(-p2 * t); }; // system dynamics

    // this inclusion function maps IR^2 (two parameters in parameter-space) to IR^p (p measurements in
    // measurement-space); it's the inclusion function analog to the system dynamics.
    InclusionFunction inclusion_function = [=](IntervalTuple& x_box) {
        std::vector<Interval> Y_vector(num_samples);
        for (std::size_t i = 0; i < num_samples; ++i) {
            Y_vector[i] = x_box[0] * boost::numeric::exp(-x_box[1] * time[i]);
        }
        IntervalTuple Y_box(Y_vector);
        return Y_box;
    };

    // compute all the measurements Y
    std::vector<Interval> Y_vector(num_samples);
    Interval tolerance_bounds(0.6, 1.4);
    for (std::size_t i = 0; i < num_samples; ++i) {
        Y_vector[i] = f(t_initial + (double)i * sample_time) * (tolerance_bounds); // additive noise like this ain't accurate
    }
    IntervalTuple Y(Y_vector); // copy...
    std::cout << Y << std::endl;

    // estimate unknown parameters
    Interval p1_0(0, 10.0); // cannot be < 0.0 otherwise unstable
    Interval p2_0(0, 10.0); // cannot be < 0.0 otherwise unstable
    IntervalTuple X_0 = {p1_0, p2_0};

    IntervalTuple Y_box = inclusion_function(X_0);
//    std::cout << Y_box << std::endl;

    // Is Y subset of [f]([x0])? If not, then the algorithm is not sound... so widen intervals of [x0].
    std::cout << std::boolalpha << "Y is subset of Y_box: " << Y.is_subset(Y_box) << std::endl;

    // Try contracting initial guess
    std::function<void(IntervalTuple&, const IntervalTuple&)> contractor = std::bind(&contractor_system_dynamics, std::placeholders::_1, std::placeholders::_2, time);
    std::cout << "Before Contraction:" << std::endl;
    std::cout << "X_0 = \n" << X_0;
    std::cout << "Y = \n" << Y << std::endl;
    contractor(X_0, Y);
    std::cout << "After Contraction:" << std::endl;
    std::cout << "X_0 = \n" << X_0;
    std::cout << "Y = \n" << Y << std::endl;

    // instantiate the SIVIA algorithm options
    sivia_options options{0.005, contractor};

    // run SIVIA algorithm
    auto [X_minus, X_plus] = sivia(inclusion_function, X_0, Y, options);

    // Count the number of elements in Xminus and Xplus (to check that soln btw runs are equal)...
    std::cout << X_minus.size() << std::endl;
    std::cout << X_plus.size() << std::endl;

    save_soln(X_minus, X_plus);

    return 0;
}
