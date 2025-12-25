cmake_minimum_required(VERSION 3.25)
project(sun-calc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "main"
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_library(sun-calc STATIC
        sun-calc/noaa_calc.h
        sun-calc/noaa_calc.cpp
        sun-calc/sun_calc.cpp
        sun-calc/sun_calc.h
)

# google tests
set(TEST_SOURCES
        test/main.cpp
)

enable_testing()
add_executable(sun-calc-test ${TEST_SOURCES} )
target_link_libraries(sun-calc-test GTest::gtest_main GTest::gmock_main sun-calc)
target_include_directories(sun-calc-test PRIVATE ${CMAKE_SOURCE_DIR}/sun-calc)
gtest_discover_tests(sun-calc-test)
