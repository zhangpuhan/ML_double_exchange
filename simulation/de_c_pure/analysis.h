//
// Created by Puhan Zhang on 2/18/20.
//

#ifndef DE_C_PURE_ANALYSIS_H
#define DE_C_PURE_ANALYSIS_H

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <random>
#include <ctime>
#include <chrono>
#include "lattice_base.h"
#include "armadillo"

void initialize(std::vector<Square_Site> &test_lattice) {
    std::ofstream share_file("share_file.csv", std::ios::out);
    generate_random_spins(test_lattice);
    arma::mat initial_set_up(test_lattice.size(), 12);
    initial_set_up.zeros();
    for (u_long i = 0; i < test_lattice.size(); ++i) {
        initial_set_up.at(i, 0) = test_lattice[i].get_spin_x();
        initial_set_up.at(i, 1) = test_lattice[i].get_spin_y();
        initial_set_up.at(i, 2) = test_lattice[i].get_spin_z();
    }
    initial_set_up.save(share_file, arma::csv_ascii);
    share_file.close();
}

void calculate_mag(std::vector<Square_Site> &test_lattice, int round = 0) {
    std::ofstream mag_file;
    if (std::filesystem::exists("mag.csv")) {
        mag_file.open("mag.csv", std::ios_base::app);
    }
    else {
        mag_file.open("mag.csv", std::ios_base::out);
    }
    std::vector<double> mag(3, 0.0);
    for (auto &site : test_lattice) {
        mag[0] += site.get_spin_x();
        mag[1] += site.get_spin_y();
        mag[2] += site.get_spin_z();
    }

    mag[0] /= (test_lattice.size());
    mag[1] /= (test_lattice.size());
    mag[2] /= (test_lattice.size());
    mag_file << round << "," << std::sqrt(std::pow(mag[0], 2.0) + std::pow(mag[1], 2.0) + std::pow(mag[2], 2.0)) << std::endl;
    mag_file.close();
}

arma::mat copy_from_ml(std::vector<Square_Site> &test_lattice) {
    arma::mat first_return_set_up(test_lattice.size(), 12);
    first_return_set_up.load("share_file.csv", arma::csv_ascii);
    for (u_long i = 0; i < test_lattice.size(); ++i) {
        test_lattice[i].set_spins(first_return_set_up.at(i, 0), first_return_set_up.at(i, 1), first_return_set_up.at(i, 2));
        test_lattice[i].set_force(first_return_set_up.at(i, 3), first_return_set_up.at(i, 4), first_return_set_up.at(i, 5));
    }
    return first_return_set_up;
}

void generate_t_vector(std::vector<Square_Site> &test_lattice, double alpha = 0.1) {
    for (auto &site : test_lattice) {
        arma::vec B_vector(3);
        arma::vec X_vector(3);
        B_vector[0] = site.get_force_x();
        B_vector[1] = site.get_force_y();
        B_vector[2] = site.get_force_z();
        X_vector[0] = site.get_spin_x();
        X_vector[1] = site.get_spin_y();
        X_vector[2] = site.get_spin_z();
        arma::vec T_vector = B_vector + alpha * arma::cross(X_vector, B_vector);
        site.set_t_vector(T_vector[0], T_vector[1], T_vector[2]);
    }
}

void generate_c_vector(std::vector<Square_Site> &test_lattice, double h = 0.01) {
    for (auto &site : test_lattice) {
        arma::mat k_matrix = arma::ones(3, 3);
        arma::vec X_vector(3);
        X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();
        k_matrix(0, 1) = 0.5 * h * site.get_t_z();
        k_matrix(0, 2) = -0.5 * h * site.get_t_y();
        k_matrix(1, 2) = 0.5 * h * site.get_t_x();
        k_matrix(1, 0) = -k_matrix(0, 1);
        k_matrix(2, 0) = -k_matrix(0, 2);
        k_matrix(2, 1) = -k_matrix(1, 2);
        arma::vec c_vector = k_matrix.i() * k_matrix.t() * X_vector;
        site.set_c_vector(c_vector[0], c_vector[1], c_vector[2]);
    }
}

void generate_half_spin_vector(std::vector<Square_Site> &test_lattice) {
    for (auto &site : test_lattice) {
        arma::vec X_vector(3);
        X_vector[0] = site.get_spin_x(); X_vector[1] = site.get_spin_y(); X_vector[2] = site.get_spin_z();
        arma::vec half_X_vector(3);
        half_X_vector[0] = site.get_c_x(); half_X_vector[1] = site.get_c_y(); half_X_vector[2] = site.get_c_z();
        arma::vec half_vector(3);
        half_vector = (X_vector + half_X_vector) * 0.5;
        site.set_i_spin_vector(half_vector[0], half_vector[1], half_vector[2]);
    }
}

void intermediate_output(std::vector<Square_Site> &test_lattice, arma::mat &intermediate_mat) {
    std::ofstream share_file("share_file.csv", std::ios::out);
    for (u_long i = 0; i < test_lattice.size(); ++i) {
        intermediate_mat.at(i, 6) = test_lattice[i].get_i_spin_x();
        intermediate_mat.at(i, 7) = test_lattice[i].get_i_spin_y();
        intermediate_mat.at(i, 8) = test_lattice[i].get_i_spin_z();
    }
    intermediate_mat.save(share_file, arma::csv_ascii);
    share_file.close();
}

void intermediate_passing(std::vector<Square_Site> &test_lattice, double h = 0.01, double alpha = 0.1) {
    arma::mat first_back = copy_from_ml(test_lattice);
    generate_t_vector(test_lattice, alpha);
    generate_c_vector(test_lattice, h);
    generate_half_spin_vector(test_lattice);
    intermediate_output(test_lattice, first_back);
}

arma::mat copy_from_ml_2(std::vector<Square_Site> &test_lattice) {
    arma::mat second_return_set_up(test_lattice.size(), 12);
    second_return_set_up.load("share_file.csv", arma::csv_ascii);
    for (u_long i = 0; i < test_lattice.size(); ++i) {
        test_lattice[i].set_spins(second_return_set_up.at(i, 0), second_return_set_up.at(i, 1), second_return_set_up.at(i, 2));
        test_lattice[i].set_force(second_return_set_up.at(i, 3), second_return_set_up.at(i, 4), second_return_set_up.at(i, 5));
        test_lattice[i].set_i_spin_vector(second_return_set_up.at(i, 6), second_return_set_up.at(i, 7), second_return_set_up.at(i, 8));
        test_lattice[i].set_i_force(second_return_set_up.at(i, 9), second_return_set_up.at(i, 10), second_return_set_up.at(i, 11));
    }
    return second_return_set_up;
}

void generate_i_t_vector(std::vector<Square_Site> &test_lattice, double alpha = 0.1) {
    for (auto &site: test_lattice) {
        arma::vec B_vector(3);
        arma::vec X_vector(3);
        B_vector[0] = site.get_i_force_x(); B_vector[1] = site.get_i_force_y(); B_vector[2] = site.get_i_force_z();
        X_vector[0] = site.get_i_spin_x(); X_vector[1] = site.get_i_spin_y(); X_vector[2] = site.get_i_spin_z();
        arma::vec T_vector = B_vector + alpha * arma::cross(X_vector, B_vector);
        site.set_i_t_vector(T_vector[0], T_vector[1], T_vector[2]);
    }
}

void generate_new_spin(std::vector<Square_Site> &test_lattice, double h = 0.01) {
    for (auto &site: test_lattice) {
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
}

void re_initialize(std::vector<Square_Site> &test_lattice) {
    std::ofstream share_file("share_file.csv", std::ios::out);
    arma::mat initial_set_up(test_lattice.size(), 12);
    initial_set_up.zeros();
    for (u_long i = 0; i < test_lattice.size(); ++i) {
        initial_set_up.at(i, 0) = test_lattice[i].get_spin_x();
        initial_set_up.at(i, 1) = test_lattice[i].get_spin_y();
        initial_set_up.at(i, 2) = test_lattice[i].get_spin_z();
    }
    initial_set_up.save(share_file, arma::csv_ascii);
    share_file.close();
}

void final_passing(std::vector<Square_Site> &test_lattice, double h = 0.01, double alpha = 0.1, int round = 0) {
    arma::mat second_back = copy_from_ml_2(test_lattice);
    calculate_mag(test_lattice, round);
    generate_i_t_vector(test_lattice, alpha);
    generate_new_spin(test_lattice, h);

    if (round % 10 == 0) {
        std::ofstream screen_shot("snapshot_save/screenshot_" + std::to_string(round) + ".csv", std::ios::out);
        second_back.save(screen_shot, arma::csv_ascii);
    }
    // re-initializing start
    re_initialize(test_lattice);

}

#endif //DE_C_PURE_ANALYSIS_H
