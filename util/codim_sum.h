// See LICENSE for license details.

#ifndef KDPART_CODIM_SUM_H
#define KDPART_CODIM_SUM_H

#include <array>
#include <vector>

namespace kdpart {
namespace util {

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