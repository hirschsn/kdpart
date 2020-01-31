// See LICENSE for license details.

/** All functions are implemented in headers only.
 * This is intentional for NodeAccess and NodeTreeStorage, so a compiler
 * can optimize the accesses via the proxy object and the tree traversals
 * done in the storage class into direct iterations over arrays.
 * 
 * However, the rest of the functionality (tree creation, node splitting,
 * utils, etc.) is also implemented in headers since I was too lazy to
 * structure it into modules. Only include this header ("kdpart.h") in
 * a one compiled file only. Otherwise you will get linker errors.
 */
#ifndef KDPART_H_INCLUDED
#define KDPART_H_INCLUDED

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cstring> // memcpy for marshalling (type punning)
#include <cassert>
#include <mpi.h>

#include "util/mpi_global_vector.h"
#include "util/codim_sum.h"
#include "util/find.h"

namespace kdpart {

// Fwd Decl
struct PartTreeStorage;

// Fwd Decl for friend declarations
namespace marshall {
size_t marshall_size(const PartTreeStorage&);
size_t marshall_size_per_node(const PartTreeStorage&);
std::vector<char> marshall_parttree(const PartTreeStorage&);
PartTreeStorage unmarshall_parttree(std::vector<char>);
}

/** Represents a logic error which is thrown if data that is only available at
 * leaf nodes is tried to be accessed at inner nodes.
 */
struct NotALeafNodeError: public std::logic_error {
    NotALeafNodeError(const char *what): std::logic_error(what) {}
};

/** Proxy class for accessing the data of a single node inside a PartTreeStorage.
 */
template <typename TStorage, typename IntRet, typename PointRet>
struct NodeAccess {
    NodeAccess(TStorage& s, size_t i): s(s), i(i) {}

    inline IntRet inner() const { return s.inner[i]; }
    inline IntRet split_direction() const { return s.split_direction[i]; }
    inline IntRet split_coord() const { return s.split_coord[i]; }
    inline IntRet pstart() const { return s.pstart[i]; }
    inline IntRet pend() const { return s.pend[i]; }
    inline IntRet psplit() const { return s.psplit[i]; }

    inline IntRet rank() const {
        if (is_inner())
            throw NotALeafNodeError("Called rank() on inner node.");
        return pstart();
    }

    inline PointRet lu() const { return s.lu[i]; }
    inline PointRet ro() const { return s.ro[i]; }

    // Currently unused
    //inline NodeAccess parent() { return NodeAccess(s, (i - 1) / 2); }
    inline NodeAccess child1() { return NodeAccess(s, 2 * i + 1); }
    inline NodeAccess child2() { return NodeAccess(s, 2 * i + 2); }

    inline bool is_root() const { return i == 0; }

    inline size_t depth() const {
        size_t d = 0;
        for (size_t ii = i + 1; ii > 1; ii /= 2)
            d++;
        return d;
    }

    inline size_t nproc() const { return s.pend[i] - s.pstart[i] + 1; }

    inline bool is_inner() const { return inner(); }
    // Currently unused
    // !is_inner does not state anything about the existence of the node itself, therefore the check for parent().is_inner() is needed.
    //inline bool is_leaf() const { return !is_inner() && parent().is_inner(); }
private:
    TStorage& s;
    size_t i;
};


/** Struct of array storage for k-d trees holding a partitioning of processes
 * to a discrete 3d domain.
 */
struct PartTreeStorage {
    using point_type = std::array<int, 3>;
    using node_access_type = NodeAccess<PartTreeStorage&, int&, point_type&>;
    using const_node_access_type = NodeAccess<const PartTreeStorage&, int, const point_type&>;

    /** Constructs a new tree with a single uninitialized root node.
     */
    PartTreeStorage() { ensure_depth(0); }

    /** Returns the root node.
     */
    inline const_node_access_type root() const {
        return const_node_access_type(*this, 0);
    }

    /** Applies a function to each node in the tree with *no* guarantee of
     * traversal order.
     * 
     * @param f Void function with a single paramter of type const_node_access_type.
     */
    template <typename Func>
    void for_each(Func f) const {
        for (size_t i = 0; i < storage_size(); ++i) {
            if (is_existing_node(i))
                f(node(i));
        }
    }

    /** Returns the node corresponding to the subdomain of of rank "rank".
     * Precondition: Tree stores "rank", i.e. was created with "size" > "rank".
     * 
     * @param rank Subdomain number to find.
     */
    inline const_node_access_type node_of_rank(int rank) const {
        return find([rank](auto node) {
            return rank <= node.psplit() ? Descend::Left : Descend::Right;
        });
    }

    /** Returns the node responsible for point/cell "p".
     * Precondition: Tree was created with lu and ro s.t.:
     *               lu[i] <= p[i] < ro[i] for i = 0, 1, 2
     *               I.e. "p" is a coordinate inside its domain.
     * 
     * @param p Point/Cell to find.
     */
    inline const_node_access_type node_of_cell(const point_type& p) const {
        return find([&p](auto node){
            return p[node.split_direction()] < node.split_coord()? Descend::Left: Descend::Right;
        });
    }

    /** Retuirns the subdomain bounds of process "rank".
     * @see node_of_rank(int)
     */
    inline std::pair<point_type, point_type> subdomain_bounds(int rank) const {
        auto n = node_of_rank(rank);
        return std::make_pair(n.lu(), n.ro());
    }

    /** Finds the subdomain p lies in and returns its rank.
     * @see node_of_cell(const point_type&)
     */
    inline int responsible_process(const point_type& p) const {
        return node_of_cell(p).rank();
    }

    /** Applies a function to each node in the tree in a depth-first traversal.
     * Function f is being called on the parent node before any of its childs.
     * 
     * Changing nodes through "f" is explicitly supported. Especially: Children
     * which are added to a node through "f" are visited by the very same call
     * to "walk".
     * Used to initialize the tree (therefore, non-const).
     * 
     * @param f Void function with a single parameter of type const_node_access_type.
     */
    template <typename Func>
    void walk(Func f) {
        walk([](auto){ return true; }, f, 0);
    }

    /** Applies a function to each node in the tree in a depth-first traversal.
     * Descends into subtrees only if p(node) is true.
     * 
     * @see walk(Func f)
     * 
     * @param p Predicate that determines the descend into/pruning of a subtree
     * @param f Void function with a single parameter of type const_node_access_type.
     */
    template <typename Pred, typename Func>
    void walkp(Pred p, Func f) {
        walk(p, f, 0);
    }

private:
    // NodeAccess
    friend struct NodeAccess<PartTreeStorage&, int&, point_type&>;
    friend struct NodeAccess<const PartTreeStorage&, int, const point_type&>;
    // Marshalling
    friend size_t marshall::marshall_size(const PartTreeStorage&);
    friend size_t marshall::marshall_size_per_node(const PartTreeStorage&);
    friend std::vector<char> marshall::marshall_parttree(const PartTreeStorage&);
    friend PartTreeStorage marshall::unmarshall_parttree(std::vector<char>);
    // Tree creation
    template <typename LFunc, typename SFunc>
    friend PartTreeStorage make_parttree(int, std::array<int, 3>, std::array<int, 3>, LFunc, SFunc);
    // Comparison
    friend bool operator==(const PartTreeStorage& a, const PartTreeStorage& b);


    /** Returns a modifyable reference to the root node.
     */
    inline node_access_type root() {
        return node_access_type(*this, 0);
    }


    /** Ensures the tree can hold a fully refined binary tree of depth "new_depth".
     * Possibly reallocates data. Therefore, this method can throw std::bad_alloc.
     * 
     * @param new_depth Depth to support after this call.
     */
    inline void ensure_depth(size_t new_depth) {
        size_t new_size = (static_cast<size_t>(1) << (new_depth + 1)) - 1;
        if (storage_size() < new_size)
           apply_to_data_vectors([new_size](auto& v){ v.resize(new_size); });
    }

    /** Returns true if this tree holds a node with index i.
     * Precondition: "i" < storage_size(), otherwise UB.
     */
    inline bool is_existing_node(size_t i) const {
        // Node exists if parent is an inner node
        return inner[(i - 1) / 2];
    }

    /* Depth-first traversal the subtree of node i (including i itself)
     * calling function f on each node.
     * @param f Void function with a single parameter of type
     *          const_node_access_type
     */
    template <typename Pred, typename Func>
    void walk(Pred p, Func f, int i) {
        if (p(node(i))) {
            f(node(i));
            if (inner[i]) {
                walk(p, f, 2 * i + 1);
                walk(p, f, 2 * i + 2);
            }
        }
    }

    /// Descend direction for find(Chooser)
    enum class Descend {
        Left, Right
    };

    /** Finds a leaf node.
     * @param c Function. Applied to a value of type const_node_access_type it
     *          must return either Descend::Left or Descend::Right depending
     *          on in which subtree the desired leaf node is supposed to be.
     */
    template <typename Chooser>
    inline const_node_access_type find(Chooser c) const {
        return find(c, 0);
    }

    /** Finds a leaf node in the subtree of "i".
     * @see find(Chooser)
     */
    template <typename Chooser>
    const_node_access_type find(Chooser c, int i) const {
        if (!inner[i]) {
            return node(i);
        } else {
            if (c(node(i)) == Descend::Left)
                return find(c, 2 * i + 1);
            else
                return find(c, 2 * i + 2);
        }
    }

    /** Returns the number of nodes that can be hold in the current state.
     */
    inline size_t storage_size() const {
        return inner.size();
    }

    ///** Returns the number of nodes this tree holds.
    // */
    //size_t nnodes() const {
    //    return 2 * std::count(std::begin(inner), std::end(inner), 1) + inner[0];
    //}

    /** Returns a proxy object for easy access of all data associated with a single node.
     */
    inline const_node_access_type node(size_t i) const {
        return const_node_access_type(*this, i);
    }

    /** Returns a modifyable proxy object for easy read and write access of all data associated with a single node.
     * 
     * Settings of node properties is explicitly supported.
     */
    inline node_access_type node(size_t i) {
        return node_access_type(*this, i);
    }

    /** Apply a function to all raw std::vectors comprising this class.
     * The argument f will get called with
     * std::vectors of different types. Use auto lambas.
     */
    template <typename Func>
    void apply_to_data_vectors(Func f) {
        apply_to_data_vectors_impl<PartTreeStorage&>(*this, f);
    }

    /** Apply a function to all raw std::vectors comprising this class.
     * The argument f will get called with
     * std::vectors of different types. Use auto lambas.
     */
    template <typename Func>
    void apply_to_data_vectors(Func f) const {
        apply_to_data_vectors_impl<const PartTreeStorage&>(*this, f);
    }

    /** Templated implementation version of apply_to_data_vectors to
     * avoid having to implement it twice.
     */
    template <typename PTS, typename Func>
    static void apply_to_data_vectors_impl(PTS t, Func f) {
        f(t.inner);
        f(t.pstart);
        f(t.pend);
        f(t.lu);
        f(t.ro);
        f(t.split_direction);
        f(t.split_coord);
        f(t.psplit);
    }

    // The following data [i] are set when the parent, which is (i - 1) / 2,
    // is split. All values are explicitly 0 for non-existent nodes
    // (i.e. where inner[(i - 1) / 2] == 0).
    std::vector<int> inner; //< booleans: inner[i] is true if i is an inner node.
    std::vector<int> pstart, pend; //< Start and end of range of processes assigned to the subtree of node i
    std::vector<std::array<int, 3>> lu, ro; //< Lower left and upper right hand corner of subdomain i

    // The following data [i] are set by a split of node "i" itself.
    // All values are explicitly 0 for non-existent *and leaf nodes*.
    std::vector<int> split_direction; //< direction in which the node i is split if i is not a leaf node
    std::vector<int> split_coord; //< coordinate at which node i is split if i is not a leaf node
    std::vector<int> psplit; //< Splitting boundary of the process range. Processes <= last_left_process[i] are assigned to the left subtree of node i.
};

/** Equality comparison operator for the PartTreeStorage class.
 * Compares if all data is exactly the same. This is a sufficient condition
 * since it includes the size of the tree.
 */
bool operator==(const PartTreeStorage& a, const PartTreeStorage& b);

/** Inequality comparison operator for the PartTreeStorage class.
 * @see bool operator==(const PartTreeStorage& a, const PartTreeStorage& b)
 */
bool operator!=(const PartTreeStorage& a, const PartTreeStorage& b);

/** Utility functions
 */
namespace impl {

// int longest_edge(const std::array<int, 3>& lu, const std::array<int, 3>& ro)
// {
//     std::array<int, 3> local_box = {{ro[0] - lu[0], ro[1] - lu[1], ro[2] - lu[2] }};
//     return std::distance(std::begin(local_box), std::min_element(std::begin(local_box), std::end(local_box)));
// }

template <typename NodeAccess, typename SplitFunc> // too lazy to type the template params of NodeAccess here...
inline void split_node(NodeAccess& node, SplitFunc split)
{
    int splitdir = node.depth() % 3;
    // Split at the longest direction
    //int splitdir = longest_edge(node.lu(), node.ro());

    int split_at, nproc_left; // Split_at is a local coordinate
    std::tie(split_at, nproc_left) = split(splitdir, node.lu(), node.ro(), node.nproc());
    // Split in global coordinates, +1: split after "split_at".
    auto coord = node.lu()[splitdir] + split_at + 1;

    std::array<int, 3> middle_ro =  {{node.ro()[0], node.ro()[1], node.ro()[2]}};
    std::array<int, 3> middle_lu =  {{node.lu()[0], node.lu()[1], node.lu()[2]}};
    middle_lu[splitdir] = middle_ro[splitdir] = coord;

    // Set split state of node
    node.split_direction() = splitdir;
    node.split_coord() = coord;
    node.psplit() = node.pstart() + nproc_left - 1;
    node.inner() = 1;

    // Set state of child nodes
    // node.childX().inner() defaults to 0, so doesn't need to be set explicitly.
    node.child1().pstart() = node.pstart();
    node.child1().pend() = node.psplit();
    node.child1().lu() = node.lu();
    node.child1().ro() = middle_ro;

    node.child2().pstart() = node.psplit() + 1;
    node.child2().pend() = node.pend();
    node.child2().lu() = middle_lu;
    node.child2().ro() = node.ro();
}

} // kdpart::impl


/** Functions for representing an object of class PartTreeStorage
 * as a mere vector of chars to send via MPI.
 */
namespace marshall {

/** Returns the size that a number of chars a buffer has to have
 * in order to be able to hold a specific PartTreeStorage object.
 * 
 * @param t Object which size is to determine
 */
size_t marshall_size(const PartTreeStorage& t);

/** Returns the size that a single node inside a PartTreeStorage has.
 */
size_t marshall_size_per_node(const PartTreeStorage& t);

/** Returns a char vector representing a PartTreeStorage.
 * 
 * @param t PartTreeStorage to marshall
 */
std::vector<char> marshall_parttree(const PartTreeStorage& t);

/** Returns the PartTreeStorage encoded in a char buffer.
 * 
 * Warning: Using this method on a char buffer which has not been created
 *          with "marshall_parttree" will lead to silent faults and
 *          undefined behavior!
 * 
 * @param mdata Marshalled PartTreeStorage
 */
PartTreeStorage unmarshall_parttree(std::vector<char> mdata);

} // namespace "marshall"


/** Creates a PartTreeStorage.
 * 
 * Creates a tree with exactly "size" subdomains. The splitting of parent nodes
 * into child nodes (and hence, the quality of load distribution amongst subdomains)
 * is given by the function "split". It gets evaluated whenever a node is to be
 * split and determines the actual splitting.
 * 
 * Use this on a single node only. For parallel settings, use "make_parttee_par".
 * 
 * @param size Number of subdomains
 * @param lu Lower left hand corner of the domain
 * @param ro Upper right hand corner of the domain
 * @param load Function (std::array<int, 3> -> double) globally mapping cells to load values.
 * @param split Split function ("fast_splitting" or "quick_splitting")
 */
template <typename LFunc, typename SFunc>
PartTreeStorage make_parttree(int size, std::array<int, 3> lu, std::array<int, 3> ro, LFunc load, SFunc split)
{
    auto codimload = util::CodimSum<decltype(load)>(load);

    auto splitfunc = [&codimload, &split](int split_dir, std::array<int, 3> lu, std::array<int, 3> ro, int nproc) {
        return split(codimload(split_dir, lu, ro), nproc);
    };

    PartTreeStorage t;
    // Initialize the root node
    t.root().pstart() = 0;
    t.root().pend() = size - 1; // Upper bound is included
    t.root().lu() = lu;
    t.root().ro() = ro;

    t.walk([&t, &splitfunc](auto node) {
        // Need to split the node further?
        if (node.nproc() > 1) {
            t.ensure_depth(node.depth() + 1);
            impl::split_node(node, splitfunc);
        }
    });

    return t;
}

/** Make a PartTreeStorage in a parallel setting.
 * 
 * Collective function. Must be called by all ranks in the communicator "comm".
 * Builds the tree at the root node an broadcasts it to all nodes.
 * Don't call this function directly unless "load" can be evaluated globally on the root node!
 * 
 * @param comm Communicator to use
 * @param ro Box size
 * @param load Function (std::array<int, 3> -> double) evaluating the load of a cell. Must be evaluatable no only locally, but globally! Significant only at rank 0.
 * @param split Split function ("fast_splitting" or "quality_splitting")
 */
template <typename LFunc, typename SFunc>
PartTreeStorage make_parttree_par(MPI_Comm comm, std::array<int, 3> ro, LFunc load, SFunc split)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == 0) {
        PartTreeStorage t = make_parttree(size, {{0, 0, 0}}, ro, load, split);

        std::vector<char> data = marshall::marshall_parttree(t);
        int sz = static_cast<int>(data.size());
        MPI_Bcast(&sz, 1, MPI_INT, 0, comm);
        MPI_Bcast(data.data(), sz, MPI_BYTE, 0, comm);
        return t;
    } else {
        int sz;
        MPI_Bcast(&sz, 1, MPI_INT, 0, comm);
        std::vector<char> data(sz, static_cast<char>(0));
        MPI_Bcast(data.data(), sz, MPI_BYTE, 0, comm);

        return marshall::unmarshall_parttree(std::move(data));
    }
    // Not reached.
}

/********************/
/* PUBLIC INTERFACE */
/********************/

/** Determines a fast splitting of "loads".
 * 
 * Determines a fast splitting of "loads" into two parts based on a splitting
 * of "procs" into two parts of equal size. Or with a difference of 1 if
 * "nproc" is odd.
 * 
 * @param loads Vector of cost values
 * @param nproc Number of processes to split
 * 
 * @returns pair: index where to split "loads" at and number of processes assigned to the first subset
 */
std::pair<int, int> fast_splitting(std::vector<double> loads, int nproc);

/** Finds a more elaborate splitting of "loads" cost values and "procs" ranks.
 * 
 * Determines the best load splitting of loads into two parts.
 * The proc splitting is such that both subsets are guaranteed to have at least 1 element.
 * 
 * The function searches all possible "procs" splittings and their implied "loads" splittings
 * for the one that minimites:
 * max(prefix sum of load / nproc1, rest of load / rest of nproc).
 * 
 * Where "prefix sum of loads" is determined such that the fraction
 * "prefix sum of loads" / "rest of load" is closest to "nproc1" / "rest of nproc".
 * 
 * Also, the function uses a heuristic to avoid the splitting into "procs" subset of
 * very unequal cardinalities (see comment in function).
 *
 * @see get_quick_splitting
 */
std::pair<int, int> quality_splitting(std::vector<double> loads, int nproc);

/** Linearizes a 3d coordinate into a 1d index
 * 
 * Precondition: 0 <= c[i] < box[i] for i = 0, 1, 2
 * 
 * @param c 3d coordinate
 * @param box Box size.
 */
int linearize(const std::array<int, 3> c, const std::array<int, 3> box);

/** Repartitions a tree given weights for its cells.
 * 
 * Collective function. Must be called by all ranks in the communicator "comm".
 * 
 * @param s PartTreeStorage to be repartitioned
 * @param comm Communicator
 * @param cellweights Weights to use for repartitioning. Must be ordered according to "linearize".
 */
PartTreeStorage repart_parttree_par(const PartTreeStorage& s, MPI_Comm comm, const std::vector<double>& cellweights);

/** Returns a tree with "size" subdomain.
 * 
 * Suitable for an initial partitioning with more or less equally sized subdomains.
 * 
 * Non-collective. However, different processes calling with the same "size" and "ro"
 * will get the same tree.
 * 
 * @param size Number of subdomains (communicator size in MPI setting)
 * @param ro Box size
 */
PartTreeStorage initial_part_par(int size, std::array<int, 3> ro);

} // namespace "kdpart"

#endif
