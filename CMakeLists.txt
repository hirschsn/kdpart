cmake_minimum_required(VERSION 3.0)
project(KDPart)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Include MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Include C++ file(s)
file(GLOB SOURCES
  util/*.cc
)

# Build shared library
add_library(kdpart SHARED kdpart.cc)
target_link_libraries(kdpart ${MPI_LIBRARIES} ${SOURCES})

# Build static library
add_library(kdpart_static STATIC kdpart.cc)
target_link_libraries(kdpart_static ${MPI_LIBRARIES} ${SOURCES})
set_target_properties(kdpart_static PROPERTIES OUTPUT_NAME kdpart)

# Build test
add_executable(kdpart_test_par kdpart_test_par.cc)
target_link_libraries(kdpart_test_par kdpart_static)

# Tests
enable_testing()
add_test(TestKDPart kdpart_test_par)

# Install
install(TARGETS kdpart DESTINATION lib)
install(TARGETS kdpart_static DESTINATION lib)
install(FILES kdpart.h DESTINATION include)
