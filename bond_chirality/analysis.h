//
// Created by Puhan Zhang on 12/28/19.
//

#ifndef DE_C_TORCH_ANALYSIS_H
#define DE_C_TORCH_ANALYSIS_H

#include "lattice_base.h"
#include "help_function.h"
#include <random>
#include <ctime>
#include <chrono>
//#include <torch/torch.h>
#include <iostream>
#include <vector>


const double PI = 3.141592653589793238463;

//    double arr[] = {0.0,1.0,sqrt(2.0),2.0,sqrt(5.0),sqrt(8.0),3.0,sqrt(10.0),
//                    sqrt(13.0),4.0,sqrt(17.0),sqrt(18.0),sqrt(20.0),5.0};

//    double arr[] = {0.0,1.0,sqrt(2.0),2.0,sqrt(5.0),sqrt(8.0),3.0};


//arr[] = {0.0,1.0,sqrt(2.0),2.0,sqrt(5.0),sqrt(8.0),3.0,sqrt(10.0),
//sqrt(13.0),4.0,sqrt(17.0),sqrt(18.0),sqrt(20.0),5.0,sqrt(26.0),
//sqrt(29.0),sqrt(32.0),sqrt(34.0),6.0,sqrt(37.0),sqrt(40.0),
//sqrt(41.0),sqrt(45.0),7.0};


void initialize(std::vector<Square_Site>& test_lattice, const std::vector<double> &distance_range,
        double cut_off_distance, double cut_off_distance_addition=0.0) {
    true_distance(test_lattice);
    arrange_for_neighbors(test_lattice, distance_range);

    for (auto i = 0; i < 8; ++i) {
        symmetry_generator(test_lattice, i + 1);
    }

    flatten_neighbors(test_lattice);
    generate_bonds(test_lattice);
    cut_off_bonds(test_lattice, cut_off_distance);
    //cut_off_bonds_addition(test_lattice, cut_off_distance_addition);
}

void print_bond(std::vector<Square_Site>& test_lattice) {
    std::ofstream new_file_1("bond.csv");
    for (Square_Site site : test_lattice) {
        new_file_1 << site.get_index() << std::endl;
        for (const auto &bond_set : site.bond) {
            for (unsigned long i = 0; i < bond_set.size(); ++i) {
                new_file_1 << bond_set[i][0] << "," << bond_set[i][1];
                if (i != bond_set.size() - 1) {
                    new_file_1 << ",";
                }
            }
            new_file_1 << std::endl;
        }
        std::cout << "Site " << site.get_index() << " done." << std::endl;
    }
    new_file_1.close();
}

void generate_random_spins(std::vector<Square_Site>& test_lattice) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for (auto &site : test_lattice) {
        site.set_theta(2.0 * PI * distribution(generator));
        site.set_phi(std::acos(2.0 * distribution(generator) - 1));
        site.set_spins(std::sin(site.get_phi()) * std::cos(site.get_theta()),
                std::sin(site.get_phi()) * std::sin(site.get_theta()), std::cos(site.get_phi()));
    }
}

//void generate_force(std::vector<Square_Site>& test_lattice) {
////    auto options = torch::TensorOptions()
////            .dtype(torch::kFloat64);
////
////    double j1 = -3.0, j2 = -2.0, j3 = -1.0, j4 = -0.5;
////
////    for (auto &site : test_lattice) {
////        torch::Tensor local_tensor = torch::rand({1, 3}, options);
////        torch::Tensor neighbor_1 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_2 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_3 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_1r = torch::rand({4, 3}, options);
////
////        local_tensor[0][0] = site.get_spin_x();
////        local_tensor[0][1] = site.get_spin_y();
////        local_tensor[0][2] = site.get_spin_z();
////
////        neighbor_1r[0][0] = test_lattice[site.neighbors[1][1]].get_spin_x();
////        neighbor_1r[0][1] = test_lattice[site.neighbors[1][1]].get_spin_y();
////        neighbor_1r[0][2] = test_lattice[site.neighbors[1][1]].get_spin_z();
////
////        neighbor_1r[1][0] = test_lattice[site.neighbors[1][2]].get_spin_x();
////        neighbor_1r[1][1] = test_lattice[site.neighbors[1][2]].get_spin_y();
////        neighbor_1r[1][2] = test_lattice[site.neighbors[1][2]].get_spin_z();
////
////        neighbor_1r[2][0] = test_lattice[site.neighbors[1][3]].get_spin_x();
////        neighbor_1r[2][1] = test_lattice[site.neighbors[1][3]].get_spin_y();
////        neighbor_1r[2][2] = test_lattice[site.neighbors[1][3]].get_spin_z();
////
////        neighbor_1r[3][0] = test_lattice[site.neighbors[1][0]].get_spin_x();
////        neighbor_1r[3][1] = test_lattice[site.neighbors[1][0]].get_spin_y();
////        neighbor_1r[3][2] = test_lattice[site.neighbors[1][0]].get_spin_z();
////
////        for (int i = 0; i < 4; ++i) {
////            neighbor_1[i][0] = test_lattice[site.neighbors[1][i]].get_spin_x();
////            neighbor_1[i][1] = test_lattice[site.neighbors[1][i]].get_spin_y();
////            neighbor_1[i][2] = test_lattice[site.neighbors[1][i]].get_spin_z();
////
////            neighbor_2[i][0] = test_lattice[site.neighbors[2][i]].get_spin_x();
////            neighbor_2[i][1] = test_lattice[site.neighbors[2][i]].get_spin_y();
////            neighbor_2[i][2] = test_lattice[site.neighbors[2][i]].get_spin_z();
////
////            neighbor_3[i][0] = test_lattice[site.neighbors[3][i]].get_spin_x();
////            neighbor_3[i][1] = test_lattice[site.neighbors[3][i]].get_spin_y();
////            neighbor_3[i][2] = test_lattice[site.neighbors[3][i]].get_spin_z();
////        }
////
////        local_tensor.set_requires_grad(true);
////        auto energy = j1 * torch::sum(local_tensor * neighbor_1) +
////                j2 * torch::sum(local_tensor * neighbor_2) +
////                j3 * torch::sum(local_tensor * neighbor_3) +
////                j4 * torch::sum(torch::sum(local_tensor * neighbor_2, 1) * torch::sum(neighbor_1 * neighbor_1r, 1));
////
////        energy.backward();
////        auto force = -local_tensor.grad();
////
////        site.set_force(force[0][0].item().to<double>(),
////                force[0][1].item().to<double>(),
////                force[0][2].item().to<double>());
////    }
//    auto options = torch::TensorOptions()
//            .dtype(torch::kFloat64);
//
//    double j1 = -1.0, j2 = -2.0, j3 = -1.0, j4 = -0.5;
//    torch::Tensor magnetic_field = torch::rand({1, 3}, options);
//
//    magnetic_field[0][0] = 0.0;
//    magnetic_field[0][1] = 0.0;
//    magnetic_field[0][2] = 1.0;
//
//    for (auto &site : test_lattice) {
//        torch::Tensor local_tensor = torch::rand({1, 3}, options);
//        torch::Tensor neighbor_1 = torch::rand({4, 3}, options);
//        torch::Tensor neighbor_2 = torch::rand({4, 3}, options);
//        torch::Tensor neighbor_3 = torch::rand({4, 3}, options);
//        torch::Tensor neighbor_1r = torch::rand({4, 3}, options);
//
//        local_tensor[0][0] = site.get_spin_x();
//        local_tensor[0][1] = site.get_spin_y();
//        local_tensor[0][2] = site.get_spin_z();
//
//        neighbor_1r[0][0] = test_lattice[site.neighbors[1][1]].get_spin_x();
//        neighbor_1r[0][1] = test_lattice[site.neighbors[1][1]].get_spin_y();
//        neighbor_1r[0][2] = test_lattice[site.neighbors[1][1]].get_spin_z();
//
//        neighbor_1r[1][0] = test_lattice[site.neighbors[1][2]].get_spin_x();
//        neighbor_1r[1][1] = test_lattice[site.neighbors[1][2]].get_spin_y();
//        neighbor_1r[1][2] = test_lattice[site.neighbors[1][2]].get_spin_z();
//
//        neighbor_1r[2][0] = test_lattice[site.neighbors[1][3]].get_spin_x();
//        neighbor_1r[2][1] = test_lattice[site.neighbors[1][3]].get_spin_y();
//        neighbor_1r[2][2] = test_lattice[site.neighbors[1][3]].get_spin_z();
//
//        neighbor_1r[3][0] = test_lattice[site.neighbors[1][0]].get_spin_x();
//        neighbor_1r[3][1] = test_lattice[site.neighbors[1][0]].get_spin_y();
//        neighbor_1r[3][2] = test_lattice[site.neighbors[1][0]].get_spin_z();
//
//        for (int i = 0; i < 4; ++i) {
//            neighbor_1[i][0] = test_lattice[site.neighbors[1][i]].get_spin_x();
//            neighbor_1[i][1] = test_lattice[site.neighbors[1][i]].get_spin_y();
//            neighbor_1[i][2] = test_lattice[site.neighbors[1][i]].get_spin_z();
//
//            neighbor_2[i][0] = test_lattice[site.neighbors[2][i]].get_spin_x();
//            neighbor_2[i][1] = test_lattice[site.neighbors[2][i]].get_spin_y();
//            neighbor_2[i][2] = test_lattice[site.neighbors[2][i]].get_spin_z();
//
//            neighbor_3[i][0] = test_lattice[site.neighbors[3][i]].get_spin_x();
//            neighbor_3[i][1] = test_lattice[site.neighbors[3][i]].get_spin_y();
//            neighbor_3[i][2] = test_lattice[site.neighbors[3][i]].get_spin_z();
//        }
//
//        local_tensor.set_requires_grad(true);
//        auto energy = j1 * torch::sum(local_tensor * neighbor_1);
////                      j2 * torch::sum(local_tensor * neighbor_2) +
////                      j3 * torch::sum(local_tensor * neighbor_3);
//        // j4 * torch::sum(torch::sum(local_tensor * neighbor_2, 1) * torch::sum(neighbor_1 * neighbor_1r, 1));
//
////    for (auto &site : test_lattice) {
////        torch::Tensor local_tensor = torch::rand({1, 3}, options);
////
////        local_tensor[0][0] = site.get_i_spin_x();
////        local_tensor[0][1] = site.get_i_spin_y();
////        local_tensor[0][2] = site.get_i_spin_z();
////
////
////        local_tensor.set_requires_grad(true);
////        auto energy = j1 * torch::sum(local_tensor * magnetic_field);
//////                      j2 * torch::sum(local_tensor * neighbor_2) +
//////                      j3 * torch::sum(local_tensor * neighbor_3);
//        energy.backward();
//        auto force = -local_tensor.grad();
//        site.set_force(force[0][0].item().to<double>(),
//                       force[0][1].item().to<double>(),
//                       force[0][2].item().to<double>());
//    }
//}
//
//void generate_snapshots(const int number_of_snapshots, std::vector<Square_Site>& test_lattice) {
//    std::string file_dir = "snapshots/snap_";
//
//    for (int i = 0; i < number_of_snapshots; ++i) {
//        std::ofstream new_file_1(file_dir + std::to_string(i) + ".csv");
//        generate_random_spins(test_lattice);
//        generate_force(test_lattice);
//        for (auto &site : test_lattice) {
//            new_file_1 << site.get_x() << "," << site.get_y() << ",";
//            new_file_1 << 0.5 << ",";
//            new_file_1 << site.get_spin_x() << "," << site.get_spin_y() << "," << site.get_spin_z() << ",";
//            new_file_1 << site.get_force_x() << "," << site.get_force_y() << "," << site.get_force_z() << std::endl;
//        }
//        new_file_1.close();
//        std::cout << "file " + std::to_string(i) << " done." << std::endl;
//    }
//}

#endif //DE_C_TORCH_ANALYSIS_H
