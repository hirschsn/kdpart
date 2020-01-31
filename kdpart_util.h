// See LICENSE for license details.

#pragma once
#ifndef KDPART_UTIL_H
#define KDPART_UTIL_H

#include <stdexcept>
#include <mpi.h>
#include <vector>
#include <limits>

namespace kdpart {
namespace util {

template <typename T>
inline MPI_Datatype _mpi_datatype();

template <>
inline MPI_Datatype _mpi_datatype<double>() {
    return MPI_DOUBLE;
};

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
        MPI_Allgatherv(local_data.data(), local_size, _mpi_datatype<T>(), global_data.data(), sizes.data(), prefixes.data(), _mpi_datatype<T>(), comm);
    }
};


/** Find the element numerically nearest to a given one in a range.
 * 
 * @param first Beginning of the range
 * @param last End of the range
 * @val Value to find.
 * @returns Iterator to nearest element.
 */
template <typename It, typename T>
It find_nearest(It first, It last, const T& val)
{
    It b = std::lower_bound(first, last, val);
    if (b == first || b == last)
        return b;
    // *(b-1) might be nearer than *b because lower_bound finds the first
    // element which is not less than val.
    else if (val - *(b - 1) <= *b - val)
        return b - 1;
    else
        return b;
}

/** Returns an iterator to the middle most minium element in a given range.
 * 
 * This function returns an iterator to the minimum element in the range [first, last).
 * If there are multiple mimimum elements, it returns an iterator to the one
 * closest to the "middle": (first + last) / 2.
 * 
 * E.g. invoking it on {1, 2, 3, 1, 1, 1} would return an iterator to the second "1" in the list.
 * 
 * @param first Beginning of the range
 * @param last End of the range
 * @param comp Compare function to use (signature: bool comp(const T& a, const T& b), where T = It::value_type)
 * @returns Iterator to the middle most minimum element.
 */
template <typename It, typename Comp>
It middle_most_min_element(It first, It last, Comp comp)
{
    typename It::difference_type radius = std::distance(first, last) / 2;
    It half = first + radius;

    // Start in the middle and work outwards
    It minel = half;
    for (int i = 1; i <= radius; ++i) {
        It e1 = half - i;
        It e2 = half + i;

        if (e1 >= first && comp(*e1, *minel))
            minel = e1;
        if (e2 < last && comp(*e2, *minel))
            minel = e2;
    }

    return minel;
}


/** Functor which calculates a 2d sums from a 3d functor.
 * Given a functor with an operator()({{i, j, k}}), call it v_{i, j, k}.
 * This class allows to calculate 2d sums of v_{i, j, k}.
 * Say, direction = 0, it can calculate sum_{j, k} of v_{i, j, k} for
 * a given i.
 * Also, it offers the possibility to return all such sums in a
 * std::vector.
 */
template <typename Func>
struct CodimSum {
    /** Constructor
     * 
     * @param load Functor (std::array<int, 3> -> double)
     */
    CodimSum(Func load): load(load) {}

    /** Calculates a 2d sub-sum of "load".
     * 
     * The summed up 2d plane is perpendicular to "dir", i.e. its codimension.
     * 
     * @param dir Codimension of the 2d area to sum over
     * @param coord Coordinate in said direction.
     * @param lu, ro Lower and upper summation bounds
     * @returns sum
     */ 
    double operator ()(size_t dir, int coord, const std::array<int, 3>& lu, const std::array<int, 3>& ro)
    {
        std::array<int, 3> i;
        double res = 0.0;
        size_t ndir = (dir + 1) % 3;
        size_t nndir = (dir + 2) % 3;
        i[dir] = coord;
        for (i[ndir] = lu[ndir]; i[ndir] < ro[ndir]; ++i[ndir]) {
            for (i[nndir] = lu[nndir]; i[nndir] < ro[nndir]; ++i[nndir]) {
                res += load(i);
            }
        }
        return res;
    }

    /** Returns all "ro[dir] - lu[dir]" elements calculated via
     * operator()(size_t, coord, const std::array<int, 3>&, const std::array<int, 3>&)
     * assembled in a std::vector of exactly this length.
     * @see operator()(size_t, coord, const std::array<int, 3>&, const std::array<int, 3>&)
     * @returns Vector in direction "dir" of all codimension sums
     */
    std::vector<double> operator()(size_t dir, const std::array<int, 3>& lu, const std::array<int, 3>& ro)
    {
        const size_t len = ro[dir] - lu[dir];
        std::vector<double> loads(len);
        for (size_t i = 0; i < len; ++i) {
            loads[i] = this->operator()(dir, lu[dir] + i, lu, ro);
        }
        return loads;
    }
private:
    Func load;
};

}
}

#endif