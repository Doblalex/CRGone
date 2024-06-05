cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER g++)

project(Pace2024)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CPP_SOLVE_ILP src/solvescip.cpp src/myscipconshdlr.cpp src/myscipheuristic.cpp)

set(SOURCES1
  src/olcmexact.cpp
  src/fileops.cpp
  src/preprocess.cpp
  src/readdata.cpp
  src/separate.cpp
  src/scc.cpp
  ${CPP_SOLVE_ILP}

  # src/MySolverInterface.cpp
  src/decomposition.cpp
  src/partition.cpp

  # src/ModDecomp.cpp
  # src/Token.cpp
  # src/DyckWord.cpp
)

add_executable(olcmexact ${SOURCES1})

target_include_directories(olcmexact
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

set(IPOPT off)
set(WORHP off)
set(ZIMPL off)
set(READLINE off)
set(SOPLEX_DIR ${CMAKE_BINARY_DIR})
set(SCIP_DIR ${CMAKE_BINARY_DIR})
set(PAPILO off)
set(SHARED off)
add_subdirectory(${PROJECT_SOURCE_DIR}/soplex)
add_subdirectory(${PROJECT_SOURCE_DIR}/scip)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scip/cmake/Modules)

find_package(SCIP REQUIRED HINTS ${SCIP_DIR})

include_directories(${SCIP_INCLUDE_DIRS})

find_library(LIBM m)

if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(olcmexact ${SCIP_LIBRARIES} ${LIBM})