cmake_minimum_required(VERSION 3.23)
project(SetInversion)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost 1.80.0 COMPONENTS program_options REQUIRED)

add_executable(SetInversion main.cpp)
target_link_libraries(
        SetInversion
    PRIVATE
        Boost::program_options)
