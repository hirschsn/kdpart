// See LICENSE for license details.

#include "kdpart.h"

#include <random>
#include <chrono>

// Poor man's test framework
#undef NDEBUG
#include <cassert>
#define CHECK assert

#define UNUSED(x) ((void)(x))

void all_cells_assigned_check(const kdpart::PartTreeStorage& t, const std::array<int, 3>& box, int size)
{
    auto is_valid_rank = [&size](int r){ return r >= 0 && r <= size; };

    std::array<int, 3> i;
    for (i[0] = 0; i[0] < box[0]; ++i[0]) {
        for (i[1] = 0; i[1] < box[1]; ++i[1]) {
            for (i[2] = 0; i[2] < box[2]; ++i[2]) {
                int r = t.responsible_process(i);

                CHECK(is_valid_rank(r));
            }
        }
    }
}

void all_valid_subdomains_check(const kdpart::PartTreeStorage& t, const std::array<int, 3>& box)
{
    auto subdomain_has_valid_coords = [&box](const std::array<int, 3> lu, const std::array<int, 3>& ro){
        return lu[0] >= 0 && lu[0] < box[0]
                && lu[1] >= 0 && lu[1] < box[1]
                && lu[2] >= 0 && lu[2] < box[2]
                && ro[0] > 0 && ro[0] <= box[0]
                && ro[1] > 0 && ro[1] <= box[1]
                && ro[2] > 0 && ro[2] <= box[2];
    };

    auto subdomain_is_nonempty = [&box](const std::array<int, 3> lu, const std::array<int, 3>& ro){
        return lu[0] < ro[0] && lu[1] < ro[1] && lu[2] < ro[2];
    };

    t.for_each([&subdomain_has_valid_coords, &subdomain_is_nonempty](auto node){
        const auto& lu = node.lu();
        const auto& ro = node.ro();
        CHECK(subdomain_has_valid_coords(lu, ro));
        CHECK(subdomain_is_nonempty(lu, ro));
    });
}

int area(const std::array<int, 3>& lu, const std::array<int, 3>& ro)
{
    return (ro[0] - lu[0]) * (ro[1] - lu[1]) * (ro[2] - lu[2]);
}

void all_procs_assigned_check(const kdpart::PartTreeStorage& t, int size, const std::array<int, 3>& box)
{
    UNUSED(box);

    for (int i = 0; i < size; ++i) {
        auto n = t.node_of_rank(i);
        CHECK(n.pstart() == i);

        // The following is already checked in all_valid_subdomains_check
        /*
        std::array<int, 3> lu, ro;
        std::tie(lu, ro) = t.subdomain_bounds(i);
        CHECK(lu[0] >= 0 && lu[0] < box[0]
              && lu[1] >= 0 && lu[1] < box[1]
              && lu[2] >= 0 && lu[2] < box[2]
              && ro[0] > 0 && ro[0] <= box[0]
              && ro[1] > 0 && ro[1] <= box[1]
              && ro[2] > 0 && ro[2] <= box[2]);
        CHECK(lu[0] < ro[0] && lu[1] < ro[1] && lu[2] < ro[2]);
        CHECK(area(lu, ro) > 0); // essentially the same as above
        */
    }
}

void all_procs_have_cells(const kdpart::PartTreeStorage& t, int size, const std::array<int, 3>& box, bool check_equal_cell_dist=false)
{
    // Search for cells and check for procs instead of searching
    // for procs and checking for procs as in "all_procs_assigned_check".
    std::vector<size_t> ncells(size, 0);

    std::array<int, 3> i;
    for (i[0] = 0; i[0] < box[0]; ++i[0]) {
        for (i[1] = 0; i[1] < box[1]; ++i[1]) {
            for (i[2] = 0; i[2] < box[2]; ++i[2]) {
                int r = t.responsible_process(i);
                // "r" is guaranteed to be a number in [0, size). See "all_cells_assigned_check"
                ncells[r]++;
            }
        }
    }

    for (auto nc: ncells)
        CHECK(nc > 0);
    
    if (check_equal_cell_dist) {
        auto min = *std::min_element(std::begin(ncells), std::end(ncells));
        auto max = *std::max_element(std::begin(ncells), std::end(ncells));
        double quot = static_cast<double>(max) / min;
        CHECK(quot < 2.0);
    }
}

void all_procs_disjoint(const kdpart::PartTreeStorage& t, int size, const std::array<int, 3>& box)
{
    UNUSED(box);
    for (int r = 0; r < size; ++r) {
        auto n = t.node_of_rank(r);
        const auto& lu = n.lu();
        const auto& ro = n.ro();
        std::array<int, 3> i;
        for (i[0] = lu[0]; i[0] < ro[0]; ++i[0]) {
            for (i[1] = lu[1]; i[1] < ro[1]; ++i[1]) {
                for (i[2] = lu[2]; i[2] < ro[2]; ++i[2]) {
                    CHECK(t.responsible_process(i) == r);
                }
            }
        }

    }
}

void mashall_test(const kdpart::PartTreeStorage& t)
{
    std::vector<char> m = kdpart::marshall::marshall_parttree(t);
    kdpart::PartTreeStorage t2 = kdpart::marshall::unmarshall_parttree(m);

    CHECK(t == t2);
}

void test_for_equal_partitioning(const kdpart::PartTreeStorage& told, const kdpart::PartTreeStorage& tnew, int size, const std::array<int, 3>& box, const std::vector<double>& cellweights)
{
    kdpart::util::GlobalVector<double> global_load(MPI_COMM_WORLD, cellweights);

    // "Cellweight" vector is w.r.t. old tree (before repartitioning)!
    auto global_load_func = [&told, &global_load](const std::array<int, 3>& c){
        auto n = told.node_of_cell(c);

        // Transform c to process ("rank") local coordinates
        std::array<int, 3> loc_c, loc_box;
        for (auto i = 0; i < 3; ++i) {
            loc_c[i] = c[i] - n.lu()[i];
            loc_box[i] = n.ro()[i] - n.lu()[i];
        }

        auto i = kdpart::linearize(loc_c, loc_box);
        assert(global_load.size(n.rank()) > i);
        
        return global_load(n.rank(), i);
    };

    std::vector<double> load(size, 0.0);
    std::array<int, 3> i;
    for (i[0] = 0; i[0] < box[0]; ++i[0]) {
        for (i[1] = 0; i[1] < box[1]; ++i[1]) {
            for (i[2] = 0; i[2] < box[2]; ++i[2]) {
                int r = tnew.responsible_process(i);
                load[r] += global_load_func(i);
            }
        }
    }
    double max = *std::max_element(std::begin(load), std::end(load));
    double min = *std::min_element(std::begin(load), std::end(load));

    double quot = max / min;
    // Be conservative, here
    CHECK(quot < 2.0);
}

void check_it()
{
    std::random_device r;
    std::default_random_engine rng(r());
    //std::uniform_int_distribution<int> uniform_dist(50, 500);
    std::uniform_int_distribution<int> uniform_dist(50, 100);

    const int N = 100;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    if (rank == 0)
        std::cout << "Doing " << N << " checks with " << size << " procs." << std::endl;

    for (int i = 0; i < N; ++i) {
        std::array<int, 3> box = {{uniform_dist(rng), uniform_dist(rng), uniform_dist(rng)}};
        // Make sure every node has the same box
        MPI_Bcast(box.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);

        auto t = kdpart::initial_part_par(size, box);

        if (rank == 0)
            std::cout << "TEST " << i;

        // Check initial tree
        all_valid_subdomains_check(t, box);
        all_cells_assigned_check(t, box, size);
        all_procs_assigned_check(t, size, box);
        all_procs_have_cells(t, size, box, true);
        all_procs_disjoint(t, size, box);
        mashall_test(t);
        if (rank == 0)
            std::cout << " passed." << std::endl;

        if (rank == 0)
            std::cout << "REPART TEST " << i;

        // Get own subdomain size
        std::array<int, 3> lu, ro;
        std::tie(lu, ro) = t.subdomain_bounds(rank);
        int ncells = area(lu, ro);

        // Provide weights for own cells and repartition (reuse int rng).
        std::vector<double> weights(ncells);
        std::generate(std::begin(weights), std::end(weights), [&uniform_dist, &rng](){
            return static_cast<double>(uniform_dist(rng));
        });
        auto t2 = kdpart::repart_parttree_par(t, MPI_COMM_WORLD, weights);

        // Check repartitioned tree
        all_valid_subdomains_check(t2, box);
        all_cells_assigned_check(t2, box, size);
        all_procs_assigned_check(t2, size, box);
        all_procs_have_cells(t2, size, box, false); // not partitioned equally w.r.t. no. of cells, so no eq. dist check.
        all_procs_disjoint(t2, size, box);
        mashall_test(t2);
        test_for_equal_partitioning(t, t2, size, box, weights);

        if (rank == 0)
            std::cout << " passed." << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    check_it();
    MPI_Finalize();
}
