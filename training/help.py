import csv
from os import listdir
from os.path import isfile, join
import pandas as pd
import torch
import constant
from constant import device
from natsort import natsorted
import collections


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
