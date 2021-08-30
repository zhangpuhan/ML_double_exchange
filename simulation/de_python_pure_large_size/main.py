import torch
import torch.nn.functional as f
import torch.utils.data as data
import constant
import torchvision
import subprocess
import pandas as pd
from test_bond_analysis import new_create_bond_matrix
from test_bond_analysis import new_process_bond_list
from test_bond_analysis import read_neighbor_list

# preprocessing
torch.manual_seed(2293)
torch.set_default_dtype(torch.float64)

bond_list = read_neighbor_list("data_input/bond_sq26_3only_75_466_30.csv", 75)
bond_list = new_process_bond_list(bond_list)


class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.input = torch.nn.Linear(554, 1024)
        self.hidden_1 = torch.nn.Linear(1024, 512)
        self.hidden_2 = torch.nn.Linear(512, 256)
        self.hidden_3 = torch.nn.Linear(256, 128)
        self.hidden_4 = torch.nn.Linear(128, 64)
        self.output = torch.nn.Linear(64, 1)

    def forward(self, x):
        x = f.relu(self.input(x))
        x = f.relu(self.hidden_1(x))
        x = f.relu(self.hidden_2(x))
        x = f.relu(self.hidden_3(x))
        x = f.relu(self.hidden_4(x))
        x = self.output(x)
        return x


def generate_force(flag):
    left_bound, right_bound = 0, 3
    if flag == 0:
        left_bound, right_bound = 0, 3
    if flag == 1:
        left_bound, right_bound = 6, 9
    data = pd.read_csv("share_file.csv", header=None)
    spin_tensor = torch.tensor(data.iloc[:, left_bound:right_bound].values, device=constant.device).requires_grad_(True)
    x_temp = new_create_bond_matrix(spin_tensor, bond_list)
    energy_prediction_temp = net(x_temp)

    energy_prediction = torch.sum(energy_prediction_temp)

    force_prediction = -torch.autograd.grad(energy_prediction, spin_tensor, create_graph=True)[0]

    data.iloc[:, (left_bound + 3):(right_bound + 3)] = force_prediction.detach().numpy()
    data.to_csv("share_file.csv", header=False, index=False)
    return energy_prediction.detach().numpy()


print("Initializing net.")
net = Net()
net.load_state_dict(torch.load("data_input/final_model_save_at_2_task_1.pt"))
net.to(constant.device)
print("Initialized!")

# Initialize the first round
print("Generate first set up")
args = "./de_c_pure initialize".split()
args.append(str(constant.lattice_n))
print(args)
popen_status = subprocess.Popen(args)
exit_code = popen_status.wait()
print("exit_status: " + str(exit_code))
print("First set up done!")

# simulation
print("Simulation start")

for round_n in range(5000):
    energy = generate_force(0)
    with open("energy.csv", "a") as energy_file:
        energy_file.write(str(round_n) + "," + str(energy) + "\n")
    # print("Generate second set up for round " + str(round_n))
    args = "./de_c_pure intermediate".split()
    args.append(str(constant.lattice_n))
    args.append(str(constant.h))
    args.append(str(constant.alpha))
    args.append(str(round_n))
    print(args)
    popen_status = subprocess.Popen(args)
    exit_code = popen_status.wait()

    generate_force(1)
    print("Generate final pass set up for " + str(round_n))
    args = "./de_c_pure finalpass".split()
    args.append(str(constant.lattice_n))
    args.append(str(constant.h))
    args.append(str(constant.alpha))
    args.append(str(round_n))
    print(args)
    popen_status = subprocess.Popen(args)
    exit_code = popen_status.wait()
    print("exit_status: " + str(exit_code))
    print("Final pass set up done!")
    print("round " + str(round_n) + " $$$$$$$$$$$$$$$$$$$$$$$")

print("Simulation done!")
