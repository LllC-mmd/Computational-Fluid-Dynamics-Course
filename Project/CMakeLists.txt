cmake_minimum_required(VERSION 3.16)
project(IUHM_vL)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")

file(GLOB IUHM_HEADERs "*.hpp")
file(GLOB IUHM_SRCs "*.cpp")
add_executable(IUHM_vL ${IUHM_SRC} ${IUHM_SRCs})

# set armadillo library location manually
set(armadillo_INCLUDE_DIRS /usr/local/opt/armadillo/include)
file(GLOB armadillo_LIBRARIES "/usr/local/opt/armadillo/lib/*.dylib")
message("Armadillo was found: ${armadillo_INCLUDE_DIRS}")
include_directories(${armadillo_INCLUDE_DIRS})

# set HDF5 library location manually
set(hdf5_INCLUDE_DIRS /usr/local/opt/hdf5/include)
file(GLOB hdf5_LIBRARIES "/usr/local/opt/hdf5/lib/*.dylib")

message("HDF5 was found: ${hdf5_INCLUDE_DIRS}")
include_directories(${hdf5_INCLUDE_DIRS})

target_link_libraries(IUHM_vL ${hdf5_LIBRARIES} ${armadillo_LIBRARIES})
