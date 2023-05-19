################### PREAMBLE ###############################

cmake_minimum_required(VERSION 3.15)

project(psr VERSION 0.1.0 LANGUAGES CXX)

################### PROJECT SETUP ##########################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # Release or Debug
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# ---- set C++ language standard ----

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

################### MAIN PROJECT TARGETS #####################

add_subdirectory(src)

