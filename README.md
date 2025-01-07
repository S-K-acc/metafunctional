# Metadensity functional theory for classical fluids: Extracting the pair potential

This repository contains code, data and a neural model for the methods presented in:

**Metadensity functional theory for classical fluids: Extracting the pair potential**  
*Stefanie M. Kampa, Florian Sammüller, Matthias Schmidt and Robert Evans; [arXiv:2411.06972](https://arxiv.org/abs/2411.06972).*


## Instructions

### Run Julia locally 

A recent version of [Julia](https://julialang.org/downloads/) needs to be installed on your system.
Launch the Julia interpreter within this directory and type `]` to enter the package manager.
Activate the environment and install the required packages as follows:

```julia
activate .
instantiate
```

Type backspace to exit the package manager.
Start a Jupyter server:

```julia
using IJulia
jupyterlab()
```

This should open JupyterLab in your browser where you can navigate to `Prediction.ipynb` or to `Training.ipynb`.

### Contents

`Prediction.ipynb` is a script in which the prediction of a density profile with the help of the trained neural metadensity functional (`model.bson`) is performed. This also containes a comparison with simulation data. 
The training procedure of such a neural network is shown in `Training.ipynb` where the training data can be downloaded. `GenerateVext.jl` contains functions with the purpose to generate a random external potential that leads to an inhomogeneous density profile which is used in `Prediction.ipynb`. The file `neural.jl` provides functions needed for the usage of the neural metadensity functional as well as for the DFT-minimization (Picard iteration) and preparation of training data. `simulation.jl` allows to carry out Grand Canonical Monte Carlo simulations which can serve as training data and verification of the accuracy of predictions made with the neural method. 

This code is build up on https://github.com/sfalmo/NeuralDFT-Tutorial. There, the underlying methods are explained in more detail. 