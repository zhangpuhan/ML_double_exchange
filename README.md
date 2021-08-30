<img width="727" alt="截屏2021-08-30 下午4 50 52" src="https://user-images.githubusercontent.com/32048073/131403794-dda267ba-9528-4ebf-9688-37e7ba130f34.png">

# Machine Learning for Double Exchange Model

## Introduction
This repository includes codes, trained model samples and data samples to successfully run the machine learning spin dynamics. This will reproduce results in the paper https://arxiv.org/abs/2105.08221. Sub-directories in this repo are:
1. *training_data_sample*:

      training data samples, from 30x30 lattice spin double exchange simulations. The training data includes both from random spin condiguration calculation and quench simulation. The columns of each datafile corresponding to x, y, spin_x, spin_y, spin_x, torque_x, torque_y, torque_z. A matlab file to plot out the spin configuration and electron density-- "plot_spin.m" and its color function "rainbow.dat" are also included. 
  
2. *bond_chirality*: 
      
      generate neighboring information for each lattice site 
      
3. *training*:

      
4. *simulation*：

      enter the simulation and do following, a trained model is already included and this code will do machine learning dynamics directly in a 30x30 lattice. You can use *bond_chirality* code to generate features for larger lattice such as 100x100

      ```shell
      cd de_c_pure
      mkdir build
      cd build
      cmake ..
      make
      cp de_c_pure ../../de_python_pure_large_size
      ```
      the step above create a C++ exectable to help realizing simulation, then go to *de_python_pure_large_size* folder and do:
      ```python
      python main.py
      ```
      the calculated spin configurations will to saved to folder *snapshot_save*. The first three columns corresponds to spin_x, spin_y and spin_z. 
