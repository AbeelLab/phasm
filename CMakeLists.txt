cmake_minimum_required(VERSION 3.7)
project(../phasm)

set(CMAKE_CXX_STANDARD 14)

# Load dependencies
find_package(ZLIB)
find_package(SeqAn REQUIRED)

set(SOURCE_FILES
        ../phasm/phasm.cpp
        ../phasm/phasm.h
)
add_library(../phasm SHARED ${SOURCE_FILES})