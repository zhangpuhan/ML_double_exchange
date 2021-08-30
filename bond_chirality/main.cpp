#include "lattice_base.h"
#include <iostream>
#include <vector>
#include <chrono>
#include "analysis.h"
#include "real_dynamics.h"


int main() {
    // time clock start:
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Square_Site> test_lattice;

    test_lattice.reserve(lattice_length * lattice_length);

    for (int i = 0; i < lattice_length * lattice_length; ++i) {
        test_lattice.emplace_back(i, lattice_length);
    }

    const double arr[] = {0.0,1.0,sqrt(2.0),2.0,sqrt(5.0),sqrt(8.0),3.0,sqrt(10.0),
        sqrt(13.0),4.0,sqrt(17.0),sqrt(18.0),sqrt(20.0),5.0,sqrt(26.0)};
//             sqrt(29.0),sqrt(32.0),sqrt(34.0),6.0,sqrt(37.0),sqrt(40.0),
//             sqrt(41.0),sqrt(45.0),7.0};
    const std::vector<double> distance_range(arr, arr + sizeof(arr)/ sizeof(arr[0]));
    initialize(test_lattice, distance_range, 3, 3);
    for (auto i = 0; i < 8; ++i) {
        symmetry_generator(test_lattice, i + 1);
    }
    print_bond(test_lattice);
    show_neighbors(test_lattice);
//    print_bond(test_lattice);
//    generate_random_spins(test_lattice);
////    generate_force(test_lattice);
//    generate_snapshots(300, test_lattice);
//    generate_force(test_lattice);
    // initialize all done.


//    for (auto &site : test_lattice) {
//        std::cout << site.get_index() << ", " << site.get_spin_x() << ", " << site.get_spin_y() << ", " << site.get_spin_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_force_x() << ", " << site.get_force_y() << ", " << site.get_force_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_t_x() << ", " << site.get_t_y() << ", " << site.get_t_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_c_x() << ", " << site.get_c_y() << ", " << site.get_c_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_i_force_x() << ", " << site.get_i_force_y() << ", " << site.get_i_force_z() << std::endl;
//
//    }

//    one_round(test_lattice);
//    for (auto &site : test_lattice) {
//        std::cout << site.get_index() << ", " << site.get_spin_x() << ", " << site.get_spin_y() << ", " << site.get_spin_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_force_x() << ", " << site.get_force_y() << ", " << site.get_force_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_t_x() << ", " << site.get_t_y() << ", " << site.get_t_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_c_x() << ", " << site.get_c_y() << ", " << site.get_c_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_i_force_x() << ", " << site.get_i_force_y() << ", " << site.get_i_force_z() << std::endl;
//
//    }
//    simulation(test_lattice, 3000, 0.01, 1.0);

//        for (auto &site : test_lattice) {
//        std::cout << site.get_index() << ", " << site.get_spin_x() << ", " << site.get_spin_y() << ", " << site.get_spin_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_force_x() << ", " << site.get_force_y() << ", " << site.get_force_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_t_x() << ", " << site.get_t_y() << ", " << site.get_t_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_c_x() << ", " << site.get_c_y() << ", " << site.get_c_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_i_force_x() << ", " << site.get_i_force_y() << ", " << site.get_i_force_z() << std::endl;
//        std::cout << site.get_index() << ", " << site.get_i_force_x() << ", " << site.get_i_force_y() << ", " << site.get_i_force_z() << std::endl;
//
//    }

//    generate_snapshots(50, test_lattice);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop - start);


    int count_bond = 0;
    for (const auto &bond_set : test_lattice[0].bond) {
        count_bond += int(bond_set.size());
    }

    std::cout << "Total number of bonds: " << count_bond << std::endl;
    std::cout << "Number of neighbors: " << test_lattice[0].neighbors.size() << std::endl;
    std::cout << "Number of bond types: " << test_lattice[0].bond.size() << std::endl;
    std::cout << "Time taken by process: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
