cmake_minimum_required(VERSION 3.7)
project(libphasm VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package(SeqAn REQUIRED)


message(STATUS "These dependencies where found:")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     ZLIB       ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     BZIP2      ${BZIP2_FOUND}      ${BZIP2_VERSION_STRING}")
message(   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")
message(STATUS "The requirements were met.")

# Add include directories.
include_directories(${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Main PHASM library configuration
add_library(phasm INTERFACE)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(phasm INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    )


# Depend on a library that we defined in the top-level file
target_link_libraries(phasm INTERFACE
    ${SEQAN_LIBRARIES}
) 

# 'make install' to the correct location
install(TARGETS phasm EXPORT PhasmConfig
    ARCHIVE  DESTINATION lib
    LIBRARY  DESTINATION lib
    RUNTIME  DESTINATION bin)  # This is for Windows
install(DIRECTORY include/ DESTINATION include)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT PhasmConfig DESTINATION share/phasm/cmake)

# This makes the project importable from the build directory
export(TARGETS phasm FILE PhasmConfig.cmake)

# Every library has unit tests, of course
add_executable(libphasm-test
        tests/libphasm-test.cpp
        tests/reads.cpp
        )

target_link_libraries(libphasm-test
    phasm catch ${SEQAN_LIBRARIES})

add_test(libphasm-test libphasm-test)
