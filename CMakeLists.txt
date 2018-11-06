cmake_minimum_required(VERSION 3.0)
project(KDPart)

# Include MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Build shared and static library
add_library(kdpart SHARED kdpart.cc)
target_link_libraries(kdpart ${MPI_LIBRARIES})
add_library(kdpart_static STATIC kdpart.cc)
target_link_libraries(kdpart_static ${MPI_LIBRARIES})
set_target_properties(kdpart_static PROPERTIES OUTPUT_NAME kdpart)

# Build test executable
add_executable(kdpart_test_par kdpart_test_par.cc)
target_link_libraries(kdpart_test_par kdpart)

# Tests
enable_testing()
add_test(TestKDPart kdpart_test_par)

# Install
install(TARGETS kdpart DESTINATION lib)
install(TARGETS kdpart_static DESTINATION lib)
install(FILES kdpart.h DESTINATION include)
