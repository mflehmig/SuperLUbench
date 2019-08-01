cmake_minimum_required (VERSION 2.6)
project (SuperLUbench)

# VERSION NUMBER
set (SuperLUbench_VERSION_MAJOR 1)
set (SuperLUbench_VERSION_MINOR 0)

# SET POLICIES
if(POLICY CMP0015)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0015 NEW)
endif()

# HARD-CODED PATHS
set(LIB_SEARCH_PATH1 "/usr/lib")
set(LIB_SEARCH_PATH2 "/usr/local")
set(LIB_SEARCH_PATH3 "/usr/share")
set(LIB_SEARCH_PATH4 "/usr/lib/x86_64-linux-gnu")
set(SUPERLU_LIB "")
set(SUPERLU_SRC "")
set(SUPERLU_MT_LIB "")
set(SUPERLU_MT_LIB "")
set(SUPERLU_DIST_LIB "/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib")
set(SUPERLU_DIST_SRC "/home/thomas/Code/libraries/SuperLU_DIST_6.1.0/SRC")
set(INCLUDE_DIR "/usr/include \
    /usr/src \
    /usr/src/linux-headers-4.15.0-29/include/linux/")

# FLAGS
set(CMAKE_BUILD_TYPE DEBUG) 
set(CMAKE_C_FLAGS "-std=c11")
set(CMAKE_C_FLAGS_DEBUG "-Wall -std=c11 -O0 -g")
set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")

# SEARCH FOR DEPENDENCIES
find_library(INTEL intel)
if(NOT INTEL)
    message("Intel library not found")
    
    find_library(LAPACK NAMES atlas lapack PATHS ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
    if(LAPACK)
        message(STATUS "LAPACK found")
        link_libraries(${LAPACK})
    else()
        message(SEND_ERROR "No LAPACK library found")
    endif()

    find_library(BLAS NAMES mkl openblas blas PATHS /home/thomas/Code/libraries/builds/CombBLAS_beta_16_2/lib ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
    link_libraries(${BLAS})
    find_library(COMBBLAS NAMES CombBLAS PATHS /home/thomas/Code/libraries/builds/CombBLAS_beta_16_2/lib)
    if(COMBBLAS)
        message(STATUS "COMBBLAS found")
        link_libraries(${COMBBLAS})
        message(${COMBBLAS})
    else()
        message(SEND_ERROR "No COMBBLAS library found")
    endif()
    
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        message(STATUS "MPI found")
    else()
        message(SEND_ERROR "No MPI library found")
    endif()
else()
    message(STATUS "Intel library found")
    link_libraries(${INTEL})
endif()

find_library(SUPERLU NAMES superlu PATHS ${SUPERLU_LIB} ${LIB_SEARCH_PATCH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
if(HAVE_SUPERLU)
    message(STATUS "Found SuperLU")
else()
    message("No SuperLU found. Target SuperLU_bench will fail.")
endif()

find_library(SUPERLU_DIST NAMES superlu_dist PATHS ${SUPERLU_DIST_LIB} ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
if(SUPERLU_DIST)
    message(STATUS "Found SuperLU_dist")
    link_libraries(${SUPERLU_DIST})
    #add_definitions(-D_DIST)
else()
    message("No SuperLU_dist found. Target SuperLU_Dist_bench will fail.")
endif()

find_library(SUPERLU_MT NAMES superlu_mt PATHS ${SUPERLU_MT_LIB} ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
if(HAVE_SUPERLU_MT)
    message(STATUS "Found SuperLU_MT")
else()
    message("No SuperLU_MT found. Target SuperLU_MT_bench will fail.")
endif()

find_library(LIBM NAMES m PATHS ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
if(LIBM)
    message(STATUS "Found libmath")
    link_libraries(${LIBM})
else()
    message("No libmath found.")
endif()

find_library(PARMETIS NAMES parmetis ${LIB_SEARCH_PATH1} ${LIB_SEARCH_PATH2} ${LIB_SEARCH_PATH3} ${LIB_SEARCH_PATH4})
if(PARMETIS)
    message(STATUS "Found parMETIS")
    link_libraries(${PARMETIS})
else()
    message("No parMETIS found")
endif()

#i think we don't need metis if we have parmetis i think we do
find_library(METIS NAMES metis ${LIB_SEARCH_PATH})
if(METIS)
    message(STATUS "Found metis")
    link_libraries(${METIS})
else()
    message("No metis found")
endif()

# FOR ALL TARGETS
#include_directories(${INCLUDE_DIR})
#include_directories(${MPI_C_INCLUDE_PATH})
#link_directories(${MPI_C_LINK_FLAGS})
set(CMAKE_COMPILE_FLAGS ${CMAKE_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
set(CMAKE_C_COMPILER ${MPI_CXX_COMPILER})

# TARGET: SuperLU_bench
add_executable (SuperLU_bench src/SuperLU/dlinsolx.c src/dreadMM src/Util)
target_include_directories(SuperLU_bench PRIVATE ${SUPERLU_INCLUDE} ${SUPERLU_SRC})
target_link_libraries(SuperLU_bench ${MPI_C_LIBRARIES})

# TARGET: SuperLU_MT_bench
add_executable (SuperLU_MT_bench src/SuperLU_MT/pdlinsolx.c src/dreadMM src/Util)
target_include_directories(SuperLU_MT_bench PRIVATE ${SUPERLU_MT_INCLUDE} ${SUPERLU_MT_SRC})
target_link_libraries(SuperLU_MT_bench ${MPI_C_LIBRARIES})

#TARGET: SuperLU_dist_bench
add_executable(SuperLU_Dist_bench src/SuperLU_Dist/pddrive_ABglobal.c src/dreadMM src/Util)
target_compile_definitions(SuperLU_Dist_bench PRIVATE _DIST=1)
target_include_directories(SuperLU_Dist_bench PRIVATE ${SUPERLU_DIST_LIB} ${SUPERLU_DIST_SRC} /home/thomas/Code/libraries/SuperLU_DIST_6.1.0/CBLAS /home/thomas/Code/libraries/SuperLU_DIST_6.1.0/SRC)
target_link_libraries(SuperLU_Dist_bench ${MPI_C_LIBRARIES})

# SET POLICIES
if(POLICY CMP0015)
    cmake_policy(POP)
endif()
