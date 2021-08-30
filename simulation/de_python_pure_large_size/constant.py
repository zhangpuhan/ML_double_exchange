import math
import torch

lattice_n = 30
h = 0.01
alpha = 0.5
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

