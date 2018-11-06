// See LICENSE for license details.

#include "mpi_global_vector.h"

namespace kdpart {
namespace util {

MPI_Datatype mpi_type<double>::value = MPI_DOUBLE;

} // namespace util
} // namespace kdpart
