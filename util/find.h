// See LICENSE for license details.

#ifndef KDPART_UTIL_H
#define KDPART_UTIL_H

#include <iterator> // std::distance

namespace util {

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

}

#endif