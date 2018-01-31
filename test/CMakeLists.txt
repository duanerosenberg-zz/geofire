cmake_minimum_required (VERSION 2.6)
project (GTest)


# The version number.
set(GTest_VERSION_MAJOR 1)
set(GTest_VERSION_MINOR 0)
 
# Following should be set in call to cmake as
# cmake -DCMAKE_BUILD_TYPE=DEBUG [or RELEASE]:
#set(CMAKE_BUILD_TYPE DEBUG) 

enable_language (Fortran)

# Following should _not_ be set in CMakeLists.txt:
#SET(CMAKE_C_COMPILER       gcc)
#SET(CMAKE_CXX_COMPILER     g++)
#SET(CMAKE_Fortran_COMPILER gfortran)

# Set additional compiler options, etc:
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} -Wall")

#add_subdirectory(../src/cdg)

# Set target soures:
set(CDG_EXE_SRC  gtest_class_acc.cpp)
set(CDG_BLAS_SRC ../src/cdg/blas/gtvector.cpp
                 ../src/cdg/blas/gtmatrix.cpp
   )

# Define libs/paths to link to (used based on target):
set(MPI_LIB_PATH 
    /scratch/duane.rosenberg/lib/mpich-3.2_gcc
   )
set(GPTL_LIB_PATH 
    /scratch/duane.rosenberg/lib/GPTL_src
   )
set(CDG_COMMON_LIBS  
    $(GPTL_LIB_PATH)/libgptl.a
#   $(MPI_LIB_PATH)/libmpich.a
   )

# Aggregate sources:
set(CDG_SRC  $(CDG_BLAS_SRC) 
             $(CDG_EXE_SRC)  
   )

add_definitions(
# Add the to the search path for include files --by target:
target_include_directories(gtest PUBLIC
                    ../src/cdg/include 
                    ../src/cdg/blas    
                    ../src/cdg/comm    
                    ../src/cdg/dd      
                    ../src/cdg/io      
                    ../src/cdg/sem     
                    ../src/cdg/utils   
                   )
 
# Fix compiler options here (don't let cmake do this for us): 
if (CMAKE_C_COMPILER MATCHES gcc)
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG   "-g")
endif (CMAKE_C_COMPILER MATCHES icc)
if (CMAKE_C_COMPILER MATCHES icc)
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG   "-g")
endif (CMAKE_C_COMPILER MATCHES gcc)
#if (CMAKE_Fortran_COMPILER MATCHES g77)
#  set(CMAKE_Fortran_FLAGS "-fno-backslash ${CMAKE_Fortran_FLAGS}")
#endif (CMAKE_Fortran_COMPILER MATCHES g77)

# Add executables here. If we want to build only one, XXX, 
# then invoke 'cmake --target XXX':
add_executable(gtest $(CDG_SRC))

target_link_libraries (gtest $(CDG_COMMON_LIBS))

# Install executables and scripts
install (TARGETS ${EXECUTABLES} RUNTIME DESTINATION "../bin")

