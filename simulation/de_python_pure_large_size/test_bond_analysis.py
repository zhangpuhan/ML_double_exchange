""" Analysis functions for generate irreducible representations bond """
import torch
import constant
import csv


def read_neighbor_list(file_name, neighbor_size):
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        neighbor_list = list(reader)

    convert_list = []
    temp_neighbor_list = []
    neighbor_counter = 0
    for i in range(len(neighbor_list)):
        temp_neighbor_list.append(list(map(int, neighbor_list[i])))
        neighbor_counter += 1
        if neighbor_counter == neighbor_size + 1:
            convert_list.append(temp_neighbor_list)
            neighbor_counter = 0
            temp_neighbor_list = []

    return convert_list


def generate_bond_features(spin_tensor, index, bond_list):
    bond_feature = []
    for i in range(1, len(bond_list[0])):
        test = torch.index_select(spin_tensor, 0, bond_list[index][i])
        bond_feature.append(torch.sum(test[::2] * test[1::2], 1))
        if bond_list[index][i][0] == index:
            bond_feature.append(torch.sum(test[::2] * torch.cross(test[1::2], torch.roll(test[1::2], 1, 0)), 1))

    return torch.cat(bond_feature).reshape((1, -1))


def new_create_bond_matrix(spin_tensor, bond_list):
    bond_feature_for_all = []
    for i in range(constant.lattice_n * constant.lattice_n):
        bond_feature_for_all.append(generate_bond_features(spin_tensor, i, bond_list))

    return torch.cat(bond_feature_for_all)


def new_process_bond_list(bond_list):
    all_site = []
    for i in range(len(bond_list)):
        each_site = []
        k = 0
        for j in range(len(bond_list[0])):
            each_site.append(torch.tensor(bond_list[i][j], device=constant.device))

        all_site.append(each_site)

    return all_site
