cmake_minimum_required(VERSION 3.25)

project(intersection LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    googleTest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        703bd9caab50b139428cea1aaff9974ebee5742e 
)

FetchContent_MakeAvailable(googleTest)


add_executable(intersection main.cpp)

message(STATUS ${gtest_SOURCE_DIR})
target_include_directories(intersection PRIVATE ${gtest_SOURCE_DIR}/include)
target_link_libraries(intersection PRIVATE gtest)

install(TARGETS intersection
        RUNTIME DESTINATION bin)