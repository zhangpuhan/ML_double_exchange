from help import read_neighbor_list
from test_bond_analysis import new_create_bond_matrix
from test_bond_analysis import new_process_bond_list
from bond_analysis import create_matrix
import torch
import torch.nn.functional as f
import torch.utils.data as data
import constant

bench_mark_model = "model_save_2020032302.pt"
bond_list_input = "bond_sq26_3only_75_466_30.csv"
bond_ramp = 75
training_start = 1
start_model = "model_save_2020032302.pt"
from_bench_mark = True
task = 1
sample_period = 100

torch.manual_seed(0)
torch.set_default_dtype(torch.float64)


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


bond_list = read_neighbor_list("data_input/" + bond_list_input, bond_ramp)
bond_list = new_process_bond_list(bond_list)

net = Net()
net.to(constant.device)
if from_bench_mark:
    net.load_state_dict(torch.load("data_input/" + bench_mark_model))

optimizer = torch.optim.Adam(net.parameters(), lr=0.001)
loss_func = torch.nn.MSELoss()

while training_start <= 4:
    with open("data_output/file/training_" + str(training_start) + ".txt", "w") as outfile:
        pass
    spin_tensor_1 = torch.load("training_data/randomkpm_spin_" + str(training_start) + ".pt").to(constant.device)
    force_tensor_1 = torch.load("training_data/randomkpm_force_" + str(training_start) + ".pt").to(constant.device)
    spin_tensor_2 = torch.load("training_data/kpm_spin_" + str(training_start) + ".pt").to(constant.device)
    force_tensor_2 = torch.load("training_data/kpm_force_" + str(training_start) + ".pt").to(constant.device)
    spin_tensor = torch.cat((spin_tensor_1, spin_tensor_2))
    force_tensor = torch.cat((force_tensor_1, force_tensor_2))

    torch_data_set = data.TensorDataset(spin_tensor, force_tensor)
    loader = data.DataLoader(dataset=torch_data_set, batch_size=1, shuffle=True)

    for step, (spin, force) in enumerate(loader):
        temp_coordinate = spin[0].requires_grad_(True)

        optimizer.zero_grad()
        x_temp = new_create_bond_matrix(temp_coordinate, bond_list)
        energy_prediction_temp = net(x_temp)
        energy_prediction = torch.sum(energy_prediction_temp)
        force_prediction = -torch.autograd.grad(energy_prediction, temp_coordinate, create_graph=True)[0]

        loss = loss_func(force_prediction, force[0])

        output_string = "Task: " + str(task) + "| Step: " + str(step) + '| Loss: ' + str(loss.cpu().detach().numpy())
        with open("data_output/file/training_" + str(training_start) + ".txt", "a") as outfile:
            outfile.write(output_string)
            outfile.write("\n")
            outfile.write("***************\n")

        loss.backward()
        optimizer.step()
        if step != 0 and step % sample_period == 0:
            temp_model_save_string = "data_output/model/status_model_save_at_" + str(training_start) + "_task_" + str(
                task) + "_step_" + str(step) + ".pt"
            torch.save(net.state_dict(), temp_model_save_string)

    model_save_string = "data_output/model/final_model_save_at_" + str(training_start) + "_task_" + str(task) + ".pt"
    torch.save(net.state_dict(), model_save_string)

    with open("data_output/file/training_" + str(training_start) + ".txt", "a") as outfile:
        outfile.write("\n\n\n")
        outfile.write("Below is for test!")

    spin_tensor_test = torch.load("training_data/kpm_spin_5.pt").to(constant.device)
    force_tensor_test = torch.load("training_data/kpm_force_5.pt").to(constant.device)

    torch_data_set = data.TensorDataset(spin_tensor_test, force_tensor_test)
    loader = data.DataLoader(dataset=torch_data_set, batch_size=1, shuffle=False)

    for step, (spin, force) in enumerate(loader):
        temp_coordinate = spin[0].requires_grad_(True)

        x_temp = new_create_bond_matrix(temp_coordinate, bond_list)
        energy_prediction_temp = net(x_temp)

        energy_prediction = torch.sum(energy_prediction_temp)

        force_prediction = -torch.autograd.grad(energy_prediction, temp_coordinate, create_graph=True)[0]
        loss = loss_func(force_prediction, force[0])
        output_string = "Task: " + str(task) + "| Step: " + str(step) + '| Loss: ' + str(loss.cpu().detach().numpy())
        with open("data_output/file/training_" + str(training_start) + ".txt", "a") as outfile:
            outfile.write(output_string)
            outfile.write("\n")
            outfile.write("***************\n")

    outfile.close()

    training_start += 1

