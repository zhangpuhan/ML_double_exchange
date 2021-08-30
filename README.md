<img width="727" alt="截屏2021-08-30 下午4 50 52" src="https://user-images.githubusercontent.com/32048073/131403794-dda267ba-9528-4ebf-9688-37e7ba130f34.png">

# Machine Learning for Double Exchange Model

## Introduction
This repository includes codes, trained model samples and data samples to successfully run the machine learning spin dynamics. This will reproduce results in the paper https://arxiv.org/abs/2105.08221. Sub-directories in this repo are:
1. *training_data_sample*:

      training data samples, from 30x30 lattice spin double exchange simulations. The training data includes both from random spin condiguration calculation and quench simulation. The columns of each datafile corresponding to x, y, density, spin_x, spin_y, spin_x, torque_x, torque_y, torque_z. A matlab file to plot out the spin configuration and electron density-- "plot_spin.m" and its color function "rainbow.dat" are also included. 
  
2. *bond_chirality*: 
      
      generate neighboring information for each lattice site. Download the latest "libtorch" package into this folder and do following:
      ```shell
      mkdir build
      cd build
      cmake ..
      make
      ./de_c_torch
      ```
      to generate the neighboring information in "bond.csv". The size of the lattice "lattice_length" can be adjusted in the file "lattice_base.h".
      
3. *training*:

      this folder includes codes for training and testing model. The error decreasing can be illustrated with the training data included in this folder. Go to this folder and do:
      
      ```shell
      python main.py
      ```
      you can see the prediction error decrease as the training goes on. The trained model at some points will be saved for further use.

      
4. *simulation*：

      enter the *simulation* folder and do following, a trained model (this model is used to generate FIG. 2 in the paper) is already included and this code will do machine learning dynamics directly in a 30x30 lattice. You can use *bond_chirality* code to generate features for larger lattice such as 100x100.

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
      
 5. *simulation_result_sample*:
      
      this folder include simulation result from machine learning that is used to plot FIG. 3 and FIG. 4 in the paper. The model used to generate the spin condiguration data is included in bullet point 4.
      
When you reach this line, you have all you need to reproduce the results in https://arxiv.org/abs/2105.08221. If you have more questions, please contact pz4ee@virginia.edu for information.
