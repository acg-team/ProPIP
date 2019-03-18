# Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
# Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
# -------------------------------------------------------------------------
# This file is part of Castor
#
# Castor is a computer program whose purpose is to infer phylogentic trees
# under indel-aware and indel-non-aware substitution models for nucleotide,
# protein, and codon datasets
#
# This software is based and extends the following libraries:
# - the Bio++ libraries
#   developed by the Bio++ Development Team <http://biopp.univ-montp2.fr>
#
# - The Tree Search Heuristic Library (TSH-LIB)
#   developed by L. Gatti & M. Maiolo <http://bit.ly/tsh-lib>
#
# Castor is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.
#
# Castor is a free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Castor. If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.5)

set(CASTOR_SOFTWARENAME "Castor")
set(CASTOR_DESCRIPTION "Indel-Aware Phylogenetic Tree Inference")

project(${CASTOR_SOFTWARENAME}
        DESCRIPTION ${CASTOR_DESCRIPTION}
        LANGUAGES "CXX")

## Store the git hash of the current head
if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    file(READ "${PROJECT_SOURCE_DIR}/.git/HEAD"
            PROJECT_SOURCE_VERSION)
    if("${PROJECT_SOURCE_VERSION}" MATCHES "^ref:")
        string(REGEX REPLACE "^ref: *([^ \n\r]*).*" "\\1"
                PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
        file(READ "${PROJECT_SOURCE_DIR}/.git/${PROJECT_GIT_REF}"
                PROJECT_SOURCE_VERSION)
    endif()
    string(STRIP "${PROJECT_SOURCE_VERSION}"
            PROJECT_SOURCE_VERSION)

endif()

# Store the build date
if(WIN32)
    execute_process(COMMAND "cmd" " /c date /t"
            OUTPUT_VARIABLE DATE)
    string(REGEX REPLACE "[^0-9]*(..).*" "\\1" MONTH "${DATE}")
    set(MONTHS ""
            "Jan" "Feb" "Mar" "Apr" "May" "Jun"
            "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")
    list(GET MONTHS "${MONTH}" MONTH)
    string(REGEX REPLACE "[^/]*/(..)/(....).*" "\\1 ${MONTH} \\2"
            PROJECT_BUILD_DATE "${DATE}")
    execute_process(COMMAND "cmd" " /c echo %TIME%"
            OUTPUT_VARIABLE TIME)
    string(REGEX REPLACE "[^0-9]*(..:..:..).*" "\\1"
            PROJECT_BUILD_TIME "${TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
else()
    execute_process(COMMAND "date" "+%d %b %Y/%H:%M:%S"
            OUTPUT_VARIABLE DATE_TIME)
    string(REGEX REPLACE "([^/]*)/.*" "\\1"
            PROJECT_BUILD_DATE "${DATE_TIME}")
    string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
            PROJECT_BUILD_TIME "${DATE_TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
endif()

set(CASTOR_RELTIME ${PROJECT_BUILD_TIME})
set(CASTOR_RELDATE ${PROJECT_BUILD_DATE})
set(CASTOR_VERSION ${GIT_TAG_VERSION_STRIPPED})

project(${CASTOR_SOFTWARENAME} VERSION ${GIT_TAG_VERSION_STRIPPED})

add_definitions(-DPRJ_GITBRANCH=\"${PROJECT_GIT_REF}\"
        -DPRJ_VERSION=\"${CASTOR_VERSION}\"
        -DPRJ_GITREF=\"${PROJECT_SOURCE_VERSION}\"
        -DPRJ_DESC=\"${CASTOR_DESCRIPTION}\"
        -DPRJ_NAME=\"${CASTOR_SOFTWARENAME}\"
        -DPRJ_DATE=\"${CASTOR_RELDATE}\"
        -DPRJ_TIME=\"${CASTOR_RELTIME}\" )

set(CMAKE_CXX_STANDARD 14)

string(FIND "${CMAKE_BUILD_TYPE}" "static" STATIC)

# Using the icc/icpc/icl compilers
if (${USE_INTEL})
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Ofast -g -xSSE4.2 -axAVX,CORE-AVX2")
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-Ofast -g -xSSE4.2 -axAVX,CORE-AVX2")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -xSSE4.2 -axAVX,CORE-AVX2")
    SET(CMAKE_C_FLAGS_RELEASE "-O3 -g0 -xSSE4.2 -axAVX,CORE-AVX2")
    set(INTEL_LIB_ROOT /opt/intel/lib/intel64/)
else()
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 ")
    SET(CMAKE_C_FLAGS_RELEASE "-O3 -g0")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
endif()
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wshadow -Weffc++ -O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g")

message("-- Compilation will be performed with the following release of the software:
\tbranch   ${PROJECT_GIT_REF}
\tref      ${PROJECT_SOURCE_VERSION}
\ttime     ${CASTOR_RELDATE} ${CASTOR_RELTIME}
\tcurrent  ${CASTOR_VERSION} (latest version)")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
    message(STATUS "Compilation will be performed under Windows")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    message(STATUS "Compilation will be performed under Linux")
    find_package(glog 0.3.5 REQUIRED)
    if (${USE_INTEL})
        message(STATUS "Compilation will be performed using the Intel Compilers...")
        message(STATUS "Libraries will be searched in /opt/intel/lib/intel64")
        set(INTEL_LIBS_SHARED ${INTEL_LIB_ROOT}libimf.so ${INTEL_LIB_ROOT}libsvml.so ${INTEL_LIB_ROOT}libirng.so ${INTEL_LIB_ROOT}libintlc.so.5)
        set(INTEL_LIBS_STATIC ${INTEL_LIB_ROOT}libimf.a ${INTEL_LIB_ROOT}libsvml.a ${INTEL_LIB_ROOT}libirng.a)
    endif()
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    if (STATIC GREATER 0)
        message(FATAL_ERROR "You can not build this project statically on Mac OS. Ask Apple why! CMake will exit.")
    endif ()
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(glog 0.3.5 REQUIRED)
    if (${USE_INTEL})
        message(STATUS "Compilation will be performed using the Intel Compilers...")
        message(STATUS "Libraries will be searched in /opt/intel/lib/intel64")
        set(INTEL_LIBS_SHARED ${INTEL_LIB_ROOT}libimf.dylib ${INTEL_LIB_ROOT}libsvml.dylib ${INTEL_LIB_ROOT}libirng.dylib ${INTEL_LIB_ROOT}libintlc.dylib)
        set(INTEL_LIBS_STATIC ${INTEL_LIB_ROOT}libimf.a ${INTEL_LIB_ROOT}libsvml.a ${INTEL_LIB_ROOT}libirng.a)
    endif()
endif ()

if (${USE_INTEL})
    if(STATIC GREATER 0)
        foreach (INTEL_LIB IN LISTS INTEL_LIBS_STATIC)
            message(STATUS "   required: ${INTEL_LIB}")
        endforeach ()
    else()
        foreach (INTEL_LIB IN LISTS INTEL_LIBS_SHARED)
            message(STATUS "   required: ${INTEL_LIB}")
        endforeach ()
    endif()
endif()

# Dependencies not covered by find package should be found in the following directories
if (${CMAKE_PREFIX_PATH})
    include_directories("${CMAKE_PREFIX_PATH}/include")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}}/lib")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for headers in the following directory: ${CMAKE_PREFIX_PATH}/include")
endif ()

# Find Boost Libraries
if (STATIC GREATER 0)
set(Boost_USE_STATIC_LIBS        ON) # only find static libs
else ()
set(Boost_USE_STATIC_LIBS        OFF) # only find dynamic libs
endif()
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "  static lib:  ${Boost_LIBRARY_DIRS}" )
    message(STATUS "  include dir: ${Boost_INCLUDE_DIRS}" )
else ()

if (STATIC GREATER 0)
   message(STATUS "You have requested to link the static version of the Boost libraries.")
   message(STATUS "However, no static version of the Boost libraries has been found. " )
   message(STATUS "If you downloaded Boost source files, run the following command: " )
   message(STATUS "./b2 release --link=static install " )
endif()

endif()

# Find Tree-Search Heuristic Library
find_package(TshLib 3.0.0 REQUIRED)
if (TSHLIB_FOUND)
    include_directories(${TSHLIB_INCLUDE_DIRS})
endif ()

# Find BPP dependencies
find_package(bpp-core 4.0.0 REQUIRED)
find_package(bpp-seq  12.0.0 REQUIRED)
find_package(bpp-phyl 12.0.0 REQUIRED)

# Check for a parallel framework (openMP)
if(${USE_PARALLEL})
	message(STATUS "Compilation will check for a parallel framework" )

    find_package(OpenMP)

    if (OPENMP_FOUND)
        message(STATUS "${OpenMP_VERSION}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

        # Prepare static openMP libraries list
        foreach(loop_var ${OpenMP_CXX_LIBRARIES})
            string(REGEX MATCH "^(.*)\\.[^.]*$" dummy ${loop_var})
            set(MYFILE_WITHOUT_EXT ${CMAKE_MATCH_1}.a)
            list(APPEND OpenMP_CXX_LIBRARIES_STATIC ${MYFILE_WITHOUT_EXT})
        endforeach()

        message(STATUS "OpenMP Libraries Shared: ${OpenMP_CXX_LIBRARIES}")
        message(STATUS "OpenMP Libraries Static: ${OpenMP_CXX_LIBRARIES_STATIC}")

        add_definitions("-DUSE_OPENMP")
    endif()

	#endif()
endif()

# Add local include folder
include_directories("${PROJECT_SOURCE_DIR}/include")
# Files to compile
file(GLOB SOURCES src/*.cpp src/*.cc)
file(GLOB HEADERS include/*.hpp include/*.h)
add_executable(Castor ${SOURCES} ${HEADERS})
set_target_properties(Castor PROPERTIES PUBLIC_HEADER "${HEADERS}")


if (STATIC GREATER 0)
    message(STATUS "Compilation will produce a static executable")

    if(INTELCOMPILERVARS GREATER 0)
        message(STATUS "Adding static compilation flag to ICPC and ICC")
        SET(CMAKE_CXX_FLAGS_INTEL "${CMAKE_CXX_FLAGS_INTEL} -static-intel -static")
        SET(CMAKE_C_FLAGS_INTEL "${CMAKE_C_FLAGS_INTEL} -static-intel -static")
        SET(CMAKE_CXX_FLAGS_RELWITHINFOINTEL "${CMAKE_CXX_FLAGS_RELWITHINFOINTEL} -static-intel -static")
        SET(CMAKE_C_FLAGS_RELWITHINFOINTEL "${CMAKE_C_FLAGS_RELWITHINFOINTEL} -static-intel -static")
        SET(CMAKE_CXX_FLAGS_DEBUGINTEL "${CMAKE_CXX_FLAGS_DEBUGINTEL} -static-intel -static")
        SET(CMAKE_C_FLAGS_DEBUGINTEL "${CMAKE_C_FLAGS_DEBUGINTEL} -static-intel -static")
    endif()

    set_target_properties(Castor PROPERTIES LINK_FLAGS "-static")
    target_link_libraries(Castor tshlib-static)
    target_link_libraries(Castor glog::glog)
    target_link_libraries(Castor ${INTEL_LIBS_STATIC})
    target_link_libraries(Castor ${BPP_LIBS_STATIC})
    target_link_libraries(Castor ${Boost_LIBRARIES})
    target_link_libraries(Castor ${OpenMP_CXX_LIBRARIES_STATIC})
else ()

    target_link_libraries(Castor tshlib)
    target_link_libraries(Castor glog::glog)
    target_link_libraries(Castor ${INTEL_LIBS_SHARED})
    target_link_libraries(Castor ${BPP_LIBS_SHARED})
    target_link_libraries(Castor ${Boost_LIBRARIES})
    target_link_libraries(Castor ${OpenMP_CXX_LIBRARIES})
endif()

include (CTest)
add_test(example_test Castor)

# Test target
#enable_testing()
#add_subdirectory (tests)
