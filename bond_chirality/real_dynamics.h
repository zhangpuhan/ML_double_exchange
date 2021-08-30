//
// Created by Puhan Zhang on 2/11/20.
//

#ifndef DE_C_TORCH_REAL_DYNAMICS_H
#define DE_C_TORCH_REAL_DYNAMICS_H
#include <armadillo>
#include "lattice_base.h"
#include "analysis.h"
#include <vector>

void generate_t_vector(Square_Site &site, double alpha = 0.01) {
    arma::vec B_vector(3);
    arma::vec X_vector(3);
    B_vector[0] = site.get_force_x(); B_vector[1] = site.get_force_y(); B_vector[2] = site.get_force_z();
    X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();
    arma::vec T_vector = B_vector + alpha * arma::cross(X_vector, B_vector);
    site.set_t_vector(T_vector[0], T_vector[1], T_vector[2]);
}

void generate_i_t_vector(Square_Site &site, double alpha = 0.01) {
    arma::vec B_vector(3);
    arma::vec X_vector(3);
    B_vector[0] = site.get_i_force_x(); B_vector[1] = site.get_i_force_y(); B_vector[2] = site.get_i_force_z();
    X_vector[0] = site.get_i_spin_x(); X_vector[1] = site.get_i_spin_y(); X_vector[2] = site.get_i_spin_z();
    arma::vec T_vector = B_vector + alpha * arma::cross(X_vector, B_vector);
    site.set_i_t_vector(T_vector[0], T_vector[1], T_vector[2]);
}

void generate_c_vector(Square_Site &site, double h = 0.01) {
    arma::mat k_matrix = arma::ones(3, 3);
    arma::vec X_vector(3);
    X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();
//    std::cout << arma::norm(X_vector, 2) << std::endl;
    k_matrix(0, 1) = 0.5 * h * site.get_t_z();
    k_matrix(0, 2) = -0.5 * h * site.get_t_y();
    k_matrix(1, 2) = 0.5 * h * site.get_t_x();
    k_matrix(1, 0) = -k_matrix(0, 1);
    k_matrix(2, 0) = -k_matrix(0, 2);
    k_matrix(2, 1) = -k_matrix(1, 2);
//    std::cout << k_matrix << ", " << k_matrix.i() << ", " <<  k_matrix.t() << std::endl;
    arma::vec c_vector = k_matrix.i() * k_matrix.t() * X_vector;
//    std::cout << arma::norm(c_vector, 2) << std::endl;
    site.set_c_vector(c_vector[0], c_vector[1], c_vector[2]);
}

void generate_new_spin(Square_Site &site, double h = 0.01) {
    arma::mat k_matrix = arma::ones(3, 3);
    arma::vec X_vector(3);
    X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();

    k_matrix(0, 1) = 0.5 * h * site.get_i_t_z();
    k_matrix(0, 2) = -0.5 * h * site.get_i_t_y();
    k_matrix(1, 2) = 0.5 * h * site.get_i_t_x();
    k_matrix(1, 0) = -k_matrix(0, 1);
    k_matrix(2, 0) = -k_matrix(0, 2);
    k_matrix(2, 1) = -k_matrix(1, 2);
    arma::vec c_vector = k_matrix.i() * k_matrix.t() * X_vector;
    site.set_spins(c_vector[0], c_vector[1], c_vector[2]);
}

void generate_half_spin_vector(Square_Site &site) {
    arma::vec X_vector(3);
    X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();
    arma::vec half_X_vector(3);
    half_X_vector[0] = site.get_c_x(); half_X_vector[1] = site.get_c_y(); half_X_vector[2] = site.get_c_z();
    arma::vec half_vector(3);
    half_vector = (X_vector + half_X_vector) * 0.5;
    site.set_i_spin_vector(half_vector[0], half_vector[1], half_vector[2]);
}

//void generate_intermediate_force(std::vector<Square_Site>& test_lattice) {
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
////        torch::Tensor local_tensor = torch::rand({1, 3}, options);
////        torch::Tensor neighbor_1 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_2 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_3 = torch::rand({4, 3}, options);
////        torch::Tensor neighbor_1r = torch::rand({4, 3}, options);
//
//        local_tensor[0][0] = site.get_i_spin_x();
//        local_tensor[0][1] = site.get_i_spin_y();
//        local_tensor[0][2] = site.get_i_spin_z();
//
//        neighbor_1r[0][0] = test_lattice[site.neighbors[1][1]].get_i_spin_x();
//        neighbor_1r[0][1] = test_lattice[site.neighbors[1][1]].get_i_spin_y();
//        neighbor_1r[0][2] = test_lattice[site.neighbors[1][1]].get_i_spin_z();
//
//        neighbor_1r[1][0] = test_lattice[site.neighbors[1][2]].get_i_spin_x();
//        neighbor_1r[1][1] = test_lattice[site.neighbors[1][2]].get_i_spin_y();
//        neighbor_1r[1][2] = test_lattice[site.neighbors[1][2]].get_i_spin_z();
//
//        neighbor_1r[2][0] = test_lattice[site.neighbors[1][3]].get_i_spin_x();
//        neighbor_1r[2][1] = test_lattice[site.neighbors[1][3]].get_i_spin_y();
//        neighbor_1r[2][2] = test_lattice[site.neighbors[1][3]].get_i_spin_z();
//
//        neighbor_1r[3][0] = test_lattice[site.neighbors[1][0]].get_i_spin_x();
//        neighbor_1r[3][1] = test_lattice[site.neighbors[1][0]].get_i_spin_y();
//        neighbor_1r[3][2] = test_lattice[site.neighbors[1][0]].get_i_spin_z();
//
//        for (int i = 0; i < 4; ++i) {
//            neighbor_1[i][0] = test_lattice[site.neighbors[1][i]].get_i_spin_x();
//            neighbor_1[i][1] = test_lattice[site.neighbors[1][i]].get_i_spin_y();
//            neighbor_1[i][2] = test_lattice[site.neighbors[1][i]].get_i_spin_z();
//
//            neighbor_2[i][0] = test_lattice[site.neighbors[2][i]].get_i_spin_x();
//            neighbor_2[i][1] = test_lattice[site.neighbors[2][i]].get_i_spin_y();
//            neighbor_2[i][2] = test_lattice[site.neighbors[2][i]].get_i_spin_z();
//
//            neighbor_3[i][0] = test_lattice[site.neighbors[3][i]].get_i_spin_x();
//            neighbor_3[i][1] = test_lattice[site.neighbors[3][i]].get_i_spin_y();
//            neighbor_3[i][2] = test_lattice[site.neighbors[3][i]].get_i_spin_z();
//        }
//
//        local_tensor.set_requires_grad(true);
//        auto energy = j1 * torch::sum(local_tensor * neighbor_1);
////                      j2 * torch::sum(local_tensor * neighbor_2) +
////                      j3 * torch::sum(local_tensor * neighbor_3);
//                      // j4 * torch::sum(torch::sum(local_tensor * neighbor_2, 1) * torch::sum(neighbor_1 * neighbor_1r, 1));
//
////        for (auto &site : test_lattice) {
////            torch::Tensor local_tensor = torch::rand({1, 3}, options);
////
////            local_tensor[0][0] = site.get_i_spin_x();
////            local_tensor[0][1] = site.get_i_spin_y();
////            local_tensor[0][2] = site.get_i_spin_z();
//
//
////            local_tensor.set_requires_grad(true);
////            auto energy = j1 * torch::sum(local_tensor * magnetic_field);
////                      j2 * torch::sum(local_tensor * neighbor_2) +
////                      j3 * torch::sum(local_tensor * neighbor_3);
//        energy.backward();
//        auto force = -local_tensor.grad();
//
//        site.set_i_force(force[0][0].item().to<double>(),
//                force[0][1].item().to<double>(),
//                        force[0][2].item().to<double>());
//    }
//}
//
//double calculate_energy(std::vector<Square_Site>& test_lattice) {
//    auto options = torch::TensorOptions()
//            .dtype(torch::kFloat64);
//
//    double j1 = -1.0, j2 = -2.0, j3 = -1.0, j4 = -0.5;
//    torch::Tensor magnetic_field = torch::rand({1, 3}, options);
//
//    magnetic_field[0][0] = 0.0;
//    magnetic_field[0][1] = 0.0;
//    magnetic_field[0][2] = 1.0;
//    double energy_total = 0.0;
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
//                      // j4 * torch::sum(torch::sum(local_tensor * neighbor_2, 1) * torch::sum(neighbor_1 * neighbor_1r, 1));
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
//        energy_total += energy.item().to<double>();
//        site.set_force(force[0][0].item().to<double>(),
//                       force[0][1].item().to<double>(),
//                       force[0][2].item().to<double>());
//    }
//    return energy_total / (lattice_length * lattice_length);
//}
//
//double calculate_mag(std::vector<Square_Site> &test_lattice) {
//    std::vector<double> mag(3, 0.0);
//    for (auto &site : test_lattice) {
//        mag[0] += site.get_spin_x();
//        mag[1] += site.get_spin_y();
//        mag[2] += site.get_spin_z();
//    }
//    mag[0] /= (lattice_length * lattice_length);
//    mag[1] /= (lattice_length * lattice_length);
//    mag[2] /= (lattice_length * lattice_length);
//    return std::sqrt(std::pow(mag[0], 2.0) + std::pow(mag[1], 2.0) + std::pow(mag[2], 2.0));
//}
//
//void one_round(std::vector<Square_Site> &test_lattice, double h = 0.01, double alpha = 0.01) {
//    for (auto &site : test_lattice) {
//        generate_t_vector(site, alpha);
//        generate_c_vector(site, h);
//        generate_half_spin_vector(site);
//    }
//    generate_intermediate_force(test_lattice);
////
//    for (auto &site : test_lattice) {
//        generate_i_t_vector(site, alpha);
//        generate_new_spin(site, h);
//    }
//}
//
//void simulation(std::vector<Square_Site> &test_lattice, int round = 100, double h = 0.01, double alpha = 0.01) {
//
//    std::ofstream new_file_1("energy.csv");
////    std::ofstream new_file_2("spin_sample_damp.csv");
//    std::ofstream new_file_3("mag.csv");
//    for (int i = 0; i < round; ++i) {
//        if (i == 0) {
//            generate_random_spins(test_lattice);
//        }
//        double energy = 0.0;
//        double magnetization = 0.0;
//        energy = calculate_energy(test_lattice);
//        magnetization = calculate_mag(test_lattice);
//
//        std::cout << "round: " << i << ", " << "energy = " << energy << ", magnetization = " << magnetization << std::endl;
//
//        new_file_1 << i << "," << energy << std::endl;
//        new_file_3 << i << "," << magnetization << std::endl;
//        one_round(test_lattice, h, alpha);
//    }
//    new_file_1.close();
//    new_file_3.close();
//}

#endif //DE_C_TORCH_REAL_DYNAMICS_H
