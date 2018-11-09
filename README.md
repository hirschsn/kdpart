
# KDPart

This is a simple structure-of-arrays implementation of a k-d tree over a
discrete 3d domain for partitioning purposes (of a regular grid).


## Usage

You can use it in a sequential setting, see `make_parttree`, however
the parallel part depends on MPI, i.e. the code won't compile without.

For a parallel setting, generally see the testing code `kdpart_test_par.cc`.
Setup and repart example:

```c++
#include <array>
#include <vector>
#include <mpi.h>
#include "kdpart.h"

int size;
MPI_Comm_size(MPI_COMM_WORLD, &size);

std::array<int, 3> box = {{100, 100, 100}};
// Distributes the discrete domain [0, 100) x [0, 100) x [0, 100)
// evenly to all "size" processes.
// ("Auto" is deduced to "PartTreeStorage".)
auto t = kdpart::initial_part_par(size, box);

// ...

// Get the number of cells that a process owns under the partitioning
// given by "t"
int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
std::array<int, 3> lu, ro;
std::tie(lu, ro) = t.subdomain_bounds(rank);
int ncells = area(lu, ro); // (ro[0]-lu[0])*(ro[1]-lu[1])*(ro[2]-lu[2])

// Choose cell weights to be used for repartitioning
// (Must be ordered accorgin to "linearize" function if you run
//  over the local 3d box.)
std::vector<double> weights(ncells, 1.0);
// Repartition (collective operation).
// Process 0 collects "weights" from all processes, repartitions and
// redistributes the partitioned tree.
auto s = kdpart::repart_parttree_par(t, MPI_COMM_WORLD, weights)
```

Access to data:

```c++
PartTreeStorage t = /* ...(see above)... */;

// Get responsible process of cell (5, 11, 2):
std::array<int, 3> c = {{5, 11, 2}};
int rank = t.responsible_process(c);
// Is the same as:
int rank = t.node_of_cell(c).rank();

// Subdomain corners of process 12:
int rank  = 12;
std::array<int, 3> lu, ro;
std::tie(lu, ro) = t.subdomain_bounds(rank);
// Is the same as ("n" is an object of the template class NodeAcces,
// so better use "auto"):
auto n = t.node_of_rank(rank);
lu = n.lu();
ro = n.ro();
```

## License

See `LICENSE` file (MIT/Expat license).
