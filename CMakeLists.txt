cmake_minimum_required(VERSION 3.7)
project(phasm VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Export compilation flags for vim autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Load dependencies
find_package(ZLIB)
find_package(OpenMP QUIET)
find_package(BZip2 QUIET)
find_package(SeqAn QUIET)
find_package(LEMON QUIET)

# SeqAn
if(NOT ${SeqAn_FOUND})
    message(STATUS "No system installation of SeqAn found.")

    find_path(SEQAN_SOURCE_ROOT_DIR CMakeLists.txt
        PATHS ${CMAKE_SOURCE_DIR}/src/external/seqan
        NO_DEFAULT_PATH
        DOC "Location of the SeqAn source directory")

    if(${SEQAN_SOURCE_ROOT_DIR} STREQUAL ${CMAKE_SOURCE_DIR}/src/external/seqan)
        message("Updating git submodules...")
        execute_process(COMMAND git submodule update --init -- src/external/seqan
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    if(EXISTS "${SEQAN_SOURCE_ROOT_DIR}/include/seqan/version.h")
        set(CMAKE_INCLUDE_PATH
            ${SEQAN_SOURCE_ROOT_DIR}/include
            ${CMAKE_INCLUDE_PATH})

        set(CMAKE_MODULE_PATH
            ${SEQAN_SOURCE_ROOT_DIR}/util/cmake
            ${CMAKE_MODULE_PATH})

        set(CMAKE_PREFIX_PATH
            ${SEQAN_SOURCE_ROOT_DIR}/util/cmake
            ${CMAKE_PREFIX_PATH})

        set(SEQAN_INCLUDE_PATH ${SEQAN_SOURCE_ROOT_DIR}/include)

        message(STATUS 
            "Found SeqAn source code in ${SEQAN_SOURCE_ROOT_DIR}")
        message(STATUS "This location will be preferred over a system "
                       "installation.")
    endif()
endif()


# Lemon Graph Library
if(NOT ${LEMON_FOUND})
    find_path(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
        PATHS ${CMAKE_SOURCE_DIR}/src/external/lemon
        NO_DEFAULT_PATH
        DOC "Location of LEMON source as a CMAKE subproject")

    if(EXISTS ${LEMON_SOURCE_ROOT_DIR})
        # Prepare LEMON library.
        # - make FindLEMON work
        set (LEMON_INCLUDE_PATH ${LEMON_SOURCE_ROOT_DIR})

        # - create config.h file
        configure_file(${LEMON_SOURCE_ROOT_DIR}/lemon/config.h.in
            ${LEMON_SOURCE_ROOT_DIR}/lemon/config.h)

        set(LEMON_INCLUDE_DIRS ${LEMON_SOURCE_ROOT_DIR})
        set(LEMON_LIBRARIES lemon)

        unset(LEMON_ROOT_DIR CACHE)
        unset(LEMON_DIR CACHE)
        unset(LEMON_INCLUDE_DIR CACHE)
        unset(LEMON_LIBRARY CACHE)
    endif()
endif()

# Catch unit testing library
find_path(CATCH_SOURCE_ROOT_DIR single_include/catch.hpp
        PATHS ${CMAKE_SOURCE_DIR}/src/external/catch
        NO_DEFAULT_PATH
        DOC "Location of the Catch complete source directory")
if(${CATCH_SOURCE_ROOT_DIR} STREQUAL ${CMAKE_SOURCE_DIR}/src/external/seqan)
    message("Updating git submodules...")
    execute_process(COMMAND git submodule update --init -- src/external/catch
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(EXISTS ${CATCH_SOURCE_ROOT_DIR})
    message(STATUS "Found Catch unit testing library at ${CATCH_SOURCE_ROOT_DIR}")
    set(CATCH_INCLUDE_DIRS ${CATCH_SOURCE_ROOT_DIR}/single_include)

    add_library(catch INTERFACE)
    target_include_directories(catch INTERFACE ${CATCH_INCLUDE_DIRS})

    enable_testing()
endif()


add_subdirectory(src/libphasm)

if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
  message (STATUS "CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.")
endif ()


