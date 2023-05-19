
add_executable(psr.x driver_s_curve.cc psr.cc psr.h kinsol_wrapper.h)

if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()
find_package(Threads REQUIRED)

target_include_directories(psr.x PRIVATE "/usr/local/include")
target_link_directories(   psr.x PRIVATE "/usr/local/lib")
target_link_libraries(     psr.x cantera yaml-cpp sundials_kinsol sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsoldense sundials_sunlinsolband ${ACCELERATE} Threads::Threads)

install(TARGETS psr.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

