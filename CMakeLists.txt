cmake_minimum_required(VERSION 3.1)

project(psr LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(psr.x driver.cc psr.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()
find_package(Threads REQUIRED)
find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)

target_include_directories(psr.x PRIVATE "/usr/local/include")

target_link_libraries(psr.x ${CANTERA} ${CVODES} ${KINSOL} ${ACCELERATE} Threads::Threads)

install(TARGETS psr.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

