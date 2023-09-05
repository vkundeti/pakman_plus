project(pakman)
cmake_minimum_required(VERSION 3.13)

set(TARGET_NAME "pakman")
set(CMAKE_CXX_FLAGS "-std=c++14")

set(pakman_SRCS "")

file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)



foreach(cppfile ${cppfiles})
  # Exclude unit test executables and main #
  if (
       (NOT (cppfile MATCHES "_kernel.cpp" )) AND
       (NOT (cppfile MATCHES "_kernels.cpp")) AND 
       (NOT (cppfile MATCHES "_test.cpp")) AND 
       (NOT (cppfile MATCHES "_shad.cpp")) AND 
       (NOT (cppfile MATCHES "_tests.cpp" )) AND
       (NOT (cppfile MATCHES "incremental_pakgraph_builder.cpp")) AND
       (NOT (cppfile MATCHES "default_pakgraph_builder.cpp")) AND
       (NOT (cppfile MATCHES "concurrent_kmer_table_builder.cpp")) AND
       (NOT (cppfile MATCHES "default_kmer_table_builder.cpp")) 
     ) 
     list(APPEND pakman_SRCS ${cppfile})
  endif()
endforeach()

#message(STATUS "list of sources: ${pakman_SRCS}")


find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
add_executable(pakman-mpi ${pakman_SRCS})
#include
target_include_directories(pakman-mpi PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/../includes)
target_include_directories(pakman-mpi PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/../pakman_omp)


# Compile Definitions
target_compile_definitions(pakman-mpi PRIVATE LMER_LENGTH=8)
target_compile_definitions(pakman-mpi PRIVATE PERFORM_CONTIG_GENERATION=1)
target_compile_definitions(pakman-mpi PRIVATE DEBUG_P2TIME=1)
target_compile_definitions(pakman-mpi PRIVATE COMPACT_PGRAPH=1)
if (PAKMAN_WINDOW_SIZE)
  target_compile_definitions(pakman-mpi PRIVATE WINDW_SIZE=${PAKMAN_WINDOW_SIZE})
else()
  message(STATUS "USING DEFAULT WINDOW SIZE OF 31")
  target_compile_definitions(pakman-mpi PRIVATE WINDW_SIZE=31)
endif()

# Link #
target_link_libraries(pakman-mpi PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX 
  TBB::tbb TBB::tbbmalloc Threads::Threads)
target_include_directories(pakman-mpi PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/../pakman_omp/)

find_package(TBB REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


#set(BUILD_UNIT_TESTS "False")
if (BUILD_UNIT_TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  file(GLOB UNIT_TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
       ${CMAKE_CURRENT_SOURCE_DIR}/*_test.cpp )
  add_executable(pakman-unit-tests-exe ${UNIT_TEST_SRCS})
  target_include_directories(pakman-unit-tests-exe PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../includes)

  target_include_directories(pakman-unit-tests-exe PRIVATE 
    ${GTEST_INCLUDE_DIRS})
  target_include_directories(pakman-unit-tests-exe PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR})
  #includes from pakman_omp #
  target_include_directories(pakman-unit-tests-exe PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../pakman_omp/)

  target_link_libraries(pakman-unit-tests-exe PRIVATE ${GTEST_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} MPI::MPI_CXX OpenMP::OpenMP_CXX TBB::tbb 
      TBB::tbbmalloc Threads::Threads)
endif()
