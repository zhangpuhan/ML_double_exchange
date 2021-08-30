import constant
import help
import torch
import math


def select_bond_data(data_tensor, lattice_index, neighbor_list):
    sub_neighbor_list = neighbor_list[lattice_index]
    selected_neighbors = []
    for i in range(len(sub_neighbor_list)):
        little_neighbor_list = []
        for j in range(len(sub_neighbor_list[i])):
            little_neighbor_list.append(torch.index_select(data_tensor, 0, torch.tensor(sub_neighbor_list[i][j])))
        selected_neighbors.append(little_neighbor_list)

    return selected_neighbors


def generate_bond_a1_feature(selected_neighbors):
    spin_sum = []
    for i in range(len(selected_neighbors)):
        if len(selected_neighbors[i]) != 2:
            if len(selected_neighbors[i]) == 4:
                spin_sum.append(torch.abs(torch.sum(torch.cat(tuple(selected_neighbors[i]), 1)[0] *
                                                    torch.cat(tuple(selected_neighbors[i]), 1)[1]).reshape(1, -1)))

            if len(selected_neighbors[i]) == 8:
                spin_sum.append((1.0 / (2.0 * math.sqrt(2.0))) *
                                torch.abs(torch.sum(torch.cat(tuple(selected_neighbors[i]), 1)[0] *
                                                    torch.cat(tuple(selected_neighbors[i]), 1)[1]).reshape(1, -1)))
    return torch.cat(tuple(spin_sum), 1)


def generate_bond_a2_feature(selected_neighbors):
    feature = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
    spin_sum = []
    for i in range(len(selected_neighbors)):
        if len(selected_neighbors[i]) != 2:
            if len(selected_neighbors[i]) == 8:
                spin_sum.append((1.0 / (2.0 * math.sqrt(2.0))) *
                                torch.abs(torch.sum(torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature,
                                                                         selected_neighbors[i])]), 1)[0] *
                                                    torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature,
                                                                         selected_neighbors[i])]),
                                                              1)[1]).reshape(1, -1)))

    return torch.cat(tuple(spin_sum), 1)


def generate_bond_b1_feature(selected_neighbors):
    feature_4 = [1.0, -1.0, 1.0, -1.0]
    feature_8 = [1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0]
    spin_sum = []
    for i in range(len(selected_neighbors)):
        if len(selected_neighbors[i]) != 2:
            if len(selected_neighbors[i]) == 4:
                spin_sum.append(0.5 *
                                torch.abs(torch.sum(torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature_4,
                                                                         selected_neighbors[i])]), 1)[0] *
                                                    torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature_4,
                                                                         selected_neighbors[i])]),
                                                              1)[1]).reshape(1, -1)))
            if len(selected_neighbors[i]) == 8:
                spin_sum.append((1.0 / (2.0 * math.sqrt(2.0))) *
                                torch.abs(torch.sum(torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature_8,
                                                                         selected_neighbors[i])]), 1)[0] *
                                                    torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature_8,
                                                                         selected_neighbors[i])]),
                                                              1)[1]).reshape(1, -1)))

    return torch.cat(tuple(spin_sum), 1)


def generate_bond_b2_feature(selected_neighbors):
    feature = [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0]
    spin_sum = []
    for i in range(len(selected_neighbors)):
        if len(selected_neighbors[i]) != 2:
            if len(selected_neighbors[i]) == 8:
                spin_sum.append((1.0 / (2.0 * math.sqrt(2.0))) *
                                torch.abs(torch.sum(torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature,
                                                                         selected_neighbors[i])]), 1)[0] *
                                                    torch.cat(tuple([a * b
                                                                     for a, b in
                                                                     zip(feature,
                                                                         selected_neighbors[i])]),
                                                              1)[1]).reshape(1, -1)))

    return torch.cat(tuple(spin_sum), 1)


def create_matrix(data_tensor, neighbor_list):
    a1_for_all = []
    b1_for_all = []
    a2_for_all = []
    b2_for_all = []
    for i in range(constant.lattice_n * constant.lattice_n):
        selected_neighbor = select_bond_data(data_tensor, i, neighbor_list)
        a1_for_all.append(generate_bond_a1_feature(selected_neighbor))
        a2_for_all.append(generate_bond_a2_feature(selected_neighbor))
        b1_for_all.append(generate_bond_b1_feature(selected_neighbor))
        b2_for_all.append(generate_bond_b2_feature(selected_neighbor))

    return torch.cat((torch.cat(a1_for_all, 0),
                      torch.cat(a2_for_all, 0),
                      torch.cat(b1_for_all, 0),
                      torch.cat(b2_for_all, 0)), 1)

