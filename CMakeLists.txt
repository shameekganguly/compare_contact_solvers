
cmake_minimum_required (VERSION 2.8.12)
project (COMPARE_COORDS)
SET(CMAKE_CXX_FLAGS "-std=c++11 -I/usr/include -I/usr/local/include")

# set source
include_directories(${PROJECT_SOURCE_DIR}/)

# set common dependencies
# - eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# - SAI2-URDF. TODO: this should not be needed. something fishy.
find_package(SAI2-URDF REQUIRED)
include_directories(${SAI2-URDF_INCLUDE_DIRS})
add_definitions(${SAI2-URDF_DEFINITIONS})

# - SAI2-MODEL. TODO: this should not be needed. something fishy.
find_package(SAI2-MODEL REQUIRED)
include_directories(${SAI2-MODEL_INCLUDE_DIRS})
add_definitions(${SAI2-MODEL_DEFINITIONS})

# - rbdl
find_package(RBDL QUIET)
include_directories(${RBDL_INCLUDE_DIR})

set(COMMON_LIBRARIES
	${SAI2-URDF_LIBRARIES}
  ${SAI2-MODEL_LIBRARIES}
)

# add apps
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY                ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG          ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL     ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE        ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PROJECT_SOURCE_DIR}/bin)
set (BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)

# create an executable
# ADD_EXECUTABLE (01-compare-kuka ../main.cpp)
ADD_EXECUTABLE (02-compare-toro ../main2.cpp)

# and link the library against the executable
# TARGET_LINK_LIBRARIES (01-compare-kuka ${COMMON_LIBRARIES})
TARGET_LINK_LIBRARIES (02-compare-toro ${COMMON_LIBRARIES})
