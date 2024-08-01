cmake_minimum_required(VERSION 3.22)
project(mirath_tcith C ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

# needed for nvim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
            ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

# define default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

set(ALLOWED_WARNINGS " -Wno-unused-function")
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG -O0 -Wall -fsanitize=address ${ALLOWED_WARNINGS}")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -march=native -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector ${ALLOWED_WARNINGS}")

set(COMMON_SOURCES
        ${PROJECT_SOURCE_DIR}/src/common/prng.c
        ${PROJECT_SOURCE_DIR}/src/common/fips202.c
)

add_executable(tests_arith src/tests/tests_arith.c ${COMMON_SOURCES})


