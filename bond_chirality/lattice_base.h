//
// Created by Puhan Zhang on 12/16/19.
//

#ifndef DE_C_TORCH_LATTICE_BASE_H
#define DE_C_TORCH_LATTICE_BASE_H
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <set>
#include "help_function.h"
#include <algorithm>
//#include <ATen/core/Tensor.h>

extern const int lattice_length = 30;

class Site {
protected:
    int x_coordinate, y_coordinate;
    double theta, phi;
    double spin[3] = {0};
    double force[3] = {0};
    double c_vector[3] = {0};
    double t_vector[3] = {0};
    double intermediate_spin[3] = {0};
    double intermediate_force[3] = {0};
    double intermediate_t_vector[3] = {0};

public:
    Site() {
        x_coordinate = y_coordinate = 0;
        theta = phi = 0.0;
    }
    virtual void set_x(int coordinate) {
        x_coordinate = coordinate;
    }
    virtual void set_y(int coordinate) {
        y_coordinate = coordinate;
    }

    virtual inline int get_x() { return x_coordinate; }
    virtual inline int get_y() { return y_coordinate; }

    virtual inline double get_theta() { return theta; }
    virtual inline double get_phi() { return phi; }

    virtual ~Site() = default;

};


class Square_Site : public Site {
private:
    int index;

public:
    Square_Site(const int lattice_index, const int latticeLength) {
        Site::set_x(lattice_index % latticeLength);
        Site::set_y(lattice_index / latticeLength);
        index = lattice_index;
        distance.reserve(lattice_length * lattice_length);
        true_x.reserve(lattice_length * lattice_length);
        true_y.reserve(lattice_length * lattice_length);
    }

    std::vector<double> distance;;
    std::vector<double> true_x;
    std::vector<double> true_y;
    std::vector< std::vector<int> > neighbors;

    inline int get_x() override { return Site::get_x(); };
    inline int get_y() override { return Site::get_y(); };
    inline double get_theta() override { return Site::get_theta(); }
    inline double get_phi() override { return Site::get_phi(); }
    inline int get_index() { return index; }

    inline void set_theta(double theta_parameter) { theta = theta_parameter; }

    inline void set_phi(double phi_parameter) { phi = phi_parameter; }

    void set_spins(double spin_x, double spin_y, double spin_z) {
        spin[0] = spin_x;
        spin[1] = spin_y;
        spin[2] = spin_z;
    }

    void set_force(double force_x, double force_y, double force_z) {
        force[0] = force_x;
        force[1] = force_y;
        force[2] = force_z;
    }

    void set_i_force(double force_x, double force_y, double force_z) {
        intermediate_force[0] = force_x;
        intermediate_force[1] = force_y;
        intermediate_force[2] = force_z;
    }

    void set_c_vector(double c_x, double c_y, double c_z) {
        c_vector[0] = c_x;
        c_vector[1] = c_y;
        c_vector[2] = c_z;
    }

    void set_t_vector(double t_x, double t_y, double t_z) {
        t_vector[0] = t_x;
        t_vector[1] = t_y;
        t_vector[2] = t_z;
    }

    void set_i_t_vector(double t_x, double t_y, double t_z) {
        intermediate_t_vector[0] = t_x;
        intermediate_t_vector[1] = t_y;
        intermediate_t_vector[2] = t_z;
    }

    void set_i_spin_vector(double s_x, double s_y, double s_z) {
        intermediate_spin[0] = s_x;
        intermediate_spin[1] = s_y;
        intermediate_spin[2] = s_z;
    }

    inline double get_spin_x() { return spin[0]; }
    inline double get_spin_y() { return spin[1]; }
    inline double get_spin_z() { return spin[2]; }

    inline double get_force_x() { return force[0]; }
    inline double get_force_y() { return force[1]; }
    inline double get_force_z() { return force[2]; }

    inline double get_c_x() { return c_vector[0]; }
    inline double get_c_y() { return c_vector[1]; }
    inline double get_c_z() { return c_vector[2]; }

    inline double get_t_x() { return t_vector[0]; }
    inline double get_t_y() { return t_vector[1]; }
    inline double get_t_z() { return t_vector[2]; }

    inline double get_i_t_x() { return intermediate_t_vector[0]; }
    inline double get_i_t_y() { return intermediate_t_vector[1]; }
    inline double get_i_t_z() { return intermediate_t_vector[2]; }

    inline double get_i_spin_x() { return intermediate_spin[0]; }
    inline double get_i_spin_y() { return intermediate_spin[1]; }
    inline double get_i_spin_z() { return intermediate_spin[2]; }

    inline double get_i_force_x() { return intermediate_force[0]; }
    inline double get_i_force_y() { return intermediate_force[1]; }
    inline double get_i_force_z() { return intermediate_force[2]; }

    std::vector<std::vector<int>> neighbors_1;
    std::vector<std::vector<int>> neighbors_2;
    std::vector<std::vector<int>> neighbors_3;
    std::vector<std::vector<int>> neighbors_4;
    std::vector<std::vector<int>> neighbors_5;
    std::vector<std::vector<int>> neighbors_6;
    std::vector<std::vector<int>> neighbors_7;
    std::vector<std::vector<int>> neighbors_8;

    std::vector<int> flatten_neighbors_0;
    std::vector<int> flatten_neighbors_1;
    std::vector<int> flatten_neighbors_2;
    std::vector<int> flatten_neighbors_3;
    std::vector<int> flatten_neighbors_4;
    std::vector<int> flatten_neighbors_5;
    std::vector<int> flatten_neighbors_6;
    std::vector<int> flatten_neighbors_7;
    std::vector<int> flatten_neighbors_8;

    std::vector<std::vector<std::vector<int>>> bond;

    ~Square_Site() override = default;
};

void flatten_neighbors(std::vector<Square_Site> &spin) {
    for (auto &site : spin) {
        for (const auto &group : site.neighbors) {
            for (const int &neighbor : group) {
                site.flatten_neighbors_0.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_1) {
            for (const int &neighbor : group) {
                site.flatten_neighbors_1.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_2) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_2.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_3) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_3.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_4) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_4.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_5) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_5.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_6) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_6.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_7) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_7.push_back(neighbor);
            }
        }
        for (const auto &group : site.neighbors_8) {
            for (const auto &neighbor :  group) {
                site.flatten_neighbors_8.push_back(neighbor);
            }
        }
    }

}

void generate_bonds(std::vector<Square_Site> &spin) {
    int n, r = 2;
    n = int(spin[0].flatten_neighbors_0.size());
    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + r, true);
    std::vector< std::vector<int> > pairs;
    do {
        std::vector<int> temp_pair;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                temp_pair.push_back(i);
            }
        }
        pairs.push_back(temp_pair);
    } while (std::prev_permutation(v.begin(), v.end()));

    for (auto &site : spin) {
        std::set< std::vector<int> > global_bond;
        for (auto &pair : pairs) {
            // std::cout << pair.at(0) << "," << pair.at(1) << std::endl;
            // std::cout << site.flatten_neighbors_0.at(pair.at(0)) << "," << site.flatten_neighbors_0.at(pair.at(1)) << std::endl;
            std::vector< std::vector<int> > bond_group;
            std::vector<int> bond[9];
            bond[0].push_back(site.flatten_neighbors_0.at(pair.at(0)));
            bond[0].push_back(site.flatten_neighbors_0.at(pair.at(1)));
            // std::sort(bond[0].begin(), bond[0].end());
            auto bond_0_dummy = bond[0];
            std::sort(bond_0_dummy.begin(), bond_0_dummy.end());
            if (global_bond.count(bond_0_dummy)) {
                continue;
            }

            if (global_bond.count(bond[0])) {
                continue;
            }

            bond[1].push_back(site.flatten_neighbors_1.at(pair.at(0)));
            bond[1].push_back(site.flatten_neighbors_1.at(pair.at(1)));
            // std::sort(bond[1].begin(), bond[1].end());

            bond[2].push_back(site.flatten_neighbors_2.at(pair.at(0)));
            bond[2].push_back(site.flatten_neighbors_2.at(pair.at(1)));
            // std::sort(bond[2].begin(), bond[2].end());

            bond[3].push_back(site.flatten_neighbors_3.at(pair.at(0)));
            bond[3].push_back(site.flatten_neighbors_3.at(pair.at(1)));
            // std::sort(bond[3].begin(), bond[3].end());

            bond[4].push_back(site.flatten_neighbors_4.at(pair.at(0)));
            bond[4].push_back(site.flatten_neighbors_4.at(pair.at(1)));
            // std::sort(bond[4].begin(), bond[4].end());

            bond[5].push_back(site.flatten_neighbors_5.at(pair.at(0)));
            bond[5].push_back(site.flatten_neighbors_5.at(pair.at(1)));
            // std::sort(bond[5].begin(), bond[5].end());

            bond[6].push_back(site.flatten_neighbors_6.at(pair.at(0)));
            bond[6].push_back(site.flatten_neighbors_6.at(pair.at(1)));
            // std::sort(bond[6].begin(), bond[6].end());

            bond[7].push_back(site.flatten_neighbors_7.at(pair.at(0)));
            bond[7].push_back(site.flatten_neighbors_7.at(pair.at(1)));
            // std::sort(bond[7].begin(), bond[7].end());

            bond[8].push_back(site.flatten_neighbors_8.at(pair.at(0)));
            bond[8].push_back(site.flatten_neighbors_8.at(pair.at(1)));
            // std::sort(bond[8].begin(), bond[8].end());

            std::vector<int> bond_dummy[9];
            for (int i = 0; i < 9; ++i) {
                auto k = bond[i];
                std::sort(k.begin(), k.end());
                bond_dummy[i] = k;
            }

            for (int i = 0; i < 9; ++i) {
                bool flag = true;
                for (unsigned int j = 0; j < bond_group.size(); ++j) {
                    auto pair_dummy = bond_group.at(j);
                    std::sort(pair_dummy.begin(), pair_dummy.end());
                    if (pair_dummy == bond_dummy[i]) {
                        flag = false;
                    }
                }

                if (flag) {
                    bond_group.push_back(bond[i]);
                    global_bond.insert(bond_dummy[i]);
                }
            }
            site.bond.push_back(bond_group);
        }
    }
}

void cut_off_bonds(std::vector<Square_Site> &spin, double cut_off_distance) {
    for (auto &site : spin) {
        std::vector< std::vector<std::vector<int>>>::iterator it;
        for (it = site.bond.begin(); it != site.bond.end(); ) {
            if (site.get_index() != spin[(*it)[0][0]].get_index()
                && (spin[(*it)[0][0]].distance.at(site.get_index()) > cut_off_distance
                || spin[(*it)[0][1]].distance.at(site.get_index()) > cut_off_distance)) {
                site.bond.erase(it);
            }
            else {
                it++;
            }
        }
    }
}

void cut_off_bonds_addition(std::vector<Square_Site> &spin, double cut_off_distance) {
    for (auto &site : spin) {
        std::vector< std::vector<std::vector<int>>>::iterator it;
        for (it = site.bond.begin(); it != site.bond.end(); ) {
            if (spin[(*it)[0][0]].distance.at((*it)[0][1]) > cut_off_distance
                && site.get_index() != spin[(*it)[0][0]].get_index()) {
                site.bond.erase(it);
            }
            else {
                it++;
            }
        }
    }
}

//remember the shortest distance between 2 lattice points and the vector between 2 points, legacy code, do not touch.
void true_distance(std::vector<Square_Site> &spin)
{
    for(int i=0;i<lattice_length*lattice_length;++i)
    {
        for(int j=0;j<lattice_length*lattice_length;++j)
        {
            std::vector<double> distance_temp;
            double L[9];
            L[0]=sqrt(pow(spin[i].get_x()-spin[j].get_x(),2.0)+pow(spin[i].get_y()-spin[j].get_y(),2.0));
            L[1]=sqrt(pow(spin[i].get_x()-spin[j].get_x(),2.0)+pow(spin[i].get_y()-(spin[j].get_y()+double(lattice_length)),2.0));
            L[2]=sqrt(pow(spin[i].get_x()-spin[j].get_x(),2.0)+pow(spin[i].get_y()-(spin[j].get_y()-double(lattice_length)),2.0));
            L[3]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()+double(lattice_length)),2.0)+pow(spin[i].get_y()-spin[j].get_y(),2.0));
            L[4]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()-double(lattice_length)),2.0)+pow(spin[i].get_y()-spin[j].get_y(),2.0));
            L[5]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()+double(lattice_length)),2.0)+pow(spin[i].get_y()-(spin[j].get_y()+double(lattice_length)),2.0));
            L[6]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()-double(lattice_length)),2.0)+pow(spin[i].get_y()-(spin[j].get_y()+double(lattice_length)),2.0));
            L[7]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()+double(lattice_length)),2.0)+pow(spin[i].get_y()-(spin[j].get_y()-double(lattice_length)),2.0));
            L[8]=sqrt(pow(spin[i].get_x()-(spin[j].get_x()-double(lattice_length)),2.0)+pow(spin[i].get_y()-(spin[j].get_y()-double(lattice_length)),2.0));
            for(double k : L)
            {
                distance_temp.push_back(k);
            }
            std::vector<double>::iterator it;
            it=min_element(distance_temp.begin(), distance_temp.end());
            spin[i].distance.emplace_back(distance_temp.at(it-distance_temp.begin()));
            switch (it-distance_temp.begin()) {
                case 0:
                    spin[i].true_x.emplace_back(spin[j].get_x()-spin[i].get_x());
                    spin[i].true_y.emplace_back(spin[j].get_y()-spin[i].get_y());
                    break;
                case 1:
                    spin[i].true_x.emplace_back(spin[j].get_x()-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()+double(lattice_length))-spin[i].get_y());
                    break;
                case 2:
                    spin[i].true_x.emplace_back(spin[j].get_x()-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()-double(lattice_length))-spin[i].get_y());
                    break;
                case 3:
                    spin[i].true_x.emplace_back((spin[j].get_x()+double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back(spin[j].get_y()-spin[i].get_y());
                    break;
                case 4:
                    spin[i].true_x.emplace_back((spin[j].get_x()-double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back(spin[j].get_y()-spin[i].get_y());
                    break;
                case 5:
                    spin[i].true_x.emplace_back((spin[j].get_x()+double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()+double(lattice_length))-spin[i].get_y());
                    break;
                case 6:
                    spin[i].true_x.emplace_back((spin[j].get_x()-double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()+double(lattice_length))-spin[i].get_y());
                    break;
                case 7:
                    spin[i].true_x.emplace_back((spin[j].get_x()+double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()-double(lattice_length))-spin[i].get_y());
                    break;
                case 8:
                    spin[i].true_x.emplace_back((spin[j].get_x()-double(lattice_length))-spin[i].get_x());
                    spin[i].true_y.emplace_back((spin[j].get_y()-double(lattice_length))-spin[i].get_y());
                    break;
                default:
                    break;
            }
        }
    }
}

// legacy code, do not touch
void arrange_for_neighbors(std::vector<Square_Site> &spin, const std::vector<double> &distance_range) {

//
//    std::vector<double> distance_range(arr, arr + sizeof(arr)/sizeof(arr[0]));

    for(int i=0;i<lattice_length*lattice_length;++i)
    {
        for(double k : distance_range)
        {
            std::vector<int> neighbor_temp;
            std::vector<int> neighbor_temp2;

            for(int j=0;j<lattice_length*lattice_length;++j)
            {
                if(spin[i].distance.at(j)==k)
                {
                    neighbor_temp.push_back(j);
                }
            }

            if(neighbor_temp.size()==4)
            {
                int temp_order[4];
                for(int n : neighbor_temp)
                {
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)>=0)
                    {
                        temp_order[0]=n;
                    }
                    if(spin[i].true_x.at(n)<=0 && spin[i].true_y.at(n)>0)
                    {
                        temp_order[1]=n;
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)<=0)
                    {
                        temp_order[2]=n;
                    }
                    if(spin[i].true_x.at(n)>=0 && spin[i].true_y.at(n)<0)
                    {
                        temp_order[3]=n;
                    }
                }
                for(unsigned long n=0; n != neighbor_temp.size(); ++n)
                {
                    neighbor_temp2.push_back(temp_order[n]);
                }
                spin[i].neighbors.push_back(neighbor_temp2);
            }

            if(neighbor_temp.size()==8)
            {
                int temp_order[8];
                for(int n : neighbor_temp)
                {
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)>0)
                    {
                        if(spin[i].true_x.at(n)>spin[i].true_y.at(n))
                        {
                            temp_order[0]=n;
                        }
                        else
                        {
                            temp_order[1]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)>0)
                    {
                        if(-spin[i].true_x.at(n)<spin[i].true_y.at(n))
                        {
                            temp_order[2]=n;
                        }
                        else
                        {
                            temp_order[3]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)<0)
                    {
                        if(-spin[i].true_x.at(n)>-spin[i].true_y.at(n))
                        {
                            temp_order[4]=n;
                        }
                        else
                        {
                            temp_order[5]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)<0)
                    {
                        if(spin[i].true_x.at(n)<-spin[i].true_y.at(n))
                        {
                            temp_order[6]=n;
                        }
                        else
                        {
                            temp_order[7]=n;
                        }
                    }
                }
                for(unsigned long n=0; n != neighbor_temp.size(); ++n)
                {
                    neighbor_temp2.push_back(temp_order[n]);
                }
                spin[i].neighbors.push_back(neighbor_temp2);
            }

            if(neighbor_temp.size()==12)
            {
                int temp_order[12];
                for(int n : neighbor_temp)
                {
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)==0)
                    {
                        temp_order[0]=n;
                    }
                    if(spin[i].true_x.at(n)==0 && spin[i].true_y.at(n)>0)
                    {
                        temp_order[3]=n;
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)==0)
                    {
                        temp_order[6]=n;
                    }
                    if(spin[i].true_x.at(n)==0 && spin[i].true_y.at(n)<0)
                    {
                        temp_order[9]=n;
                    }
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)>0)
                    {
                        if(spin[i].true_x.at(n)>spin[i].true_y.at(n))
                        {
                            temp_order[1]=n;
                        }
                        else
                        {
                            temp_order[2]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)>0)
                    {
                        if(-spin[i].true_x.at(n)<spin[i].true_y.at(n))
                        {
                            temp_order[4]=n;
                        }
                        else
                        {
                            temp_order[5]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)<0 && spin[i].true_y.at(n)<0)
                    {
                        if(-spin[i].true_x.at(n)>-spin[i].true_y.at(n))
                        {
                            temp_order[7]=n;
                        }
                        else
                        {
                            temp_order[8]=n;
                        }
                    }
                    if(spin[i].true_x.at(n)>0 && spin[i].true_y.at(n)<0)
                    {
                        if(spin[i].true_x.at(n)<-spin[i].true_y.at(n))
                        {
                            temp_order[10]=n;
                        }
                        else
                        {
                            temp_order[11]=n;
                        }
                    }
                }
                for(unsigned long n=0; n != neighbor_temp.size(); ++n)
                {
                    neighbor_temp2.push_back(temp_order[n]);
                }
                spin[i].neighbors.push_back(neighbor_temp2);
            }

            if(neighbor_temp.size()!=4 && neighbor_temp.size()!=8 && neighbor_temp.size()!=12)
            {
                spin[i].neighbors.push_back(neighbor_temp);
            }
        }
    }
}

void symmetry_generator(std::vector<Square_Site> &spin, int symmetry_type) {
    switch (symmetry_type) {
        case 1: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();
//                        std::cout << intermediate_x << ", " << intermediate_y << ", " << spin[neighbor].get_index() << std::endl;


                        int temp = intermediate_x;
                        intermediate_x = intermediate_y + site.get_x();
                        intermediate_y = temp + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

//                        std::cout << intermediate_x << ", " << intermediate_y << ", " << spin[neighbor].get_index() << std::endl;
//                        std::cout << "***************" << std::endl;

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_1 = temp_group;
            }
            break;
        }
        case 2: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_1) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        // int temp = intermediate_x;
                        intermediate_x = -intermediate_x + site.get_x();
                        intermediate_y = intermediate_y + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_2 = temp_group;
            }
            break;
        }
        case 3: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_2) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        int temp = intermediate_x;
                        intermediate_x = -intermediate_y + site.get_x();
                        intermediate_y = -temp + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_3 = temp_group;
            }
            break;
        }
        case 4: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_3) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        // int temp = intermediate_x;
                        intermediate_x = intermediate_x + site.get_x();
                        intermediate_y = -intermediate_y + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_4 = temp_group;
            }
            break;
        }
        case 5: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_4) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        int temp = intermediate_x;
                        intermediate_x = intermediate_y + site.get_x();
                        intermediate_y = temp + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_5 = temp_group;
            }
            break;
        }
        case 6: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_5) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        // int temp = intermediate_x;
                        intermediate_x = -intermediate_x + site.get_x();
                        intermediate_y = intermediate_y + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_6 = temp_group;
            }
            break;
        }
        case 7: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_6) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        int temp = intermediate_x;
                        intermediate_x = -intermediate_y + site.get_x();
                        intermediate_y = -temp + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_7 = temp_group;
            }
            break;
        }
        case 8: {
            for (auto &site : spin) {
                std::vector< std::vector<int> > temp_group;

                for (const auto &group : site.neighbors_7) {
                    std::vector<int> new_neighbor;
                    for (const int &neighbor : group) {
                        int intermediate_x, intermediate_y;

                        intermediate_x = spin[neighbor].get_x() - site.get_x();
                        intermediate_y = spin[neighbor].get_y() - site.get_y();

                        // int temp = intermediate_x;
                        intermediate_x = intermediate_x + site.get_x();
                        intermediate_y = -intermediate_y + site.get_y();

                        intermediate_x = mod(intermediate_x, lattice_length);
                        intermediate_y = mod(intermediate_y, lattice_length);

                        for (const int &neighbor_add : group) {
                            if (intermediate_x == spin[neighbor_add].get_x() &&
                                intermediate_y == spin[neighbor_add].get_y()) {
                                new_neighbor.push_back(spin[neighbor_add].get_index());
                            }
                        }
                    }
                    temp_group.push_back(new_neighbor);
                }
                site.neighbors_8 = temp_group;
            }
            break;
        }
        default: break;

    }
}

// exhibit neighbors, legacy code
void show_neighbors(std::vector<Square_Site> &spin) {
    std::ofstream new_file_1("neighbor.csv");
    for(int i=0;i<lattice_length*lattice_length;++i) {
//        std::cout << spin[i].true_x.size() << std::endl;
        for (const auto& group : spin[i].neighbors) {
            for (unsigned long k = 0; k < group.size(); ++k) {
                new_file_1 << group.at(k) ;
                // std::cout << group.at(k) ;
                if (k != group.size() - 1) {
                    // std::cout << ",";
                    new_file_1 << ",";
                }
            }
            // std::cout << std::endl;
            new_file_1 << std::endl;
        }
    }
    new_file_1.close();
}

#endif //DE_C_TORCH_LATTICE_BASE_H
