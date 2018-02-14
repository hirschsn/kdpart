// See LICENSE for license details.

// TODO: Allgatherv -> Gatherv

#ifndef KDPART_MPI_GLOBAL_VECTOR_H
#define KDPART_MPI_GLOBAL_VECTOR_H

#include <stdexcept>
#include <mpi.h>

namespace kdpart {
namespace util {

template <typename T>
struct mpi_type;

template <>
struct mpi_type<double> {
    static MPI_Datatype value;
};

MPI_Datatype mpi_type<double>::value = MPI_DOUBLE;

/** Allgather wrapper.
 */
template <typename T>
struct GlobalVector {
    GlobalVector(MPI_Comm comm, const std::vector<T>& local_data): local_data(local_data), comm(comm) {
        MPI_Comm_rank(comm, &my_rank);
        MPI_Comm_size(comm, &comm_size);
        sizes.resize(comm_size, 0);
        prefixes.resize(comm_size, 0);
        allgatherv();
    }

    T operator()(int rank, size_t index) {
        return global_data[prefixes[rank] + index];
    }

    int size(int rank) {
        return sizes[rank];
    }

private:
    const std::vector<T>& local_data;
    std::vector<T> global_data;
    // Store sizes as ints. Everything else is unreasonable when working with MPI.
    std::vector<int> sizes;
    std::vector<int> prefixes;
    int comm_size;
    int my_rank;
    MPI_Comm comm;

    void allgatherv() {
        if (local_data.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("Size of data too big for an int.");
        }
        int local_size = static_cast<int>(local_data.size());
        MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
        // C++17 or later: std::exclusive_scan
        prefixes[0] = 0;
        for (size_t i = 1; i < sizes.size(); ++i)
            prefixes[i] = prefixes[i - 1] + sizes[i - 1];
        int global_size = prefixes[comm_size - 1] + sizes[comm_size - 1];
        global_data.resize(global_size);
        MPI_Allgatherv(local_data.data(), local_size, mpi_type<T>::value, global_data.data(), sizes.data(), prefixes.data(), mpi_type<T>::value, comm);
    }
};

}
}

#endif