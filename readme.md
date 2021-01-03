# Interior Point Method

This is an implementation of IPM in CPP and PYTHON.


## Subfolders

- `code`: contains the source code of implementation.
    - `code/main.py`: source code of solver in PYTHON.
    - `code/test_extra.py`: full test process on extra data.
    - `code/cpp/*`: CPP version
- `data`: contains the 3 examples and extra data.
- `appendix`: the given *data1*, from whose folder structure I built other data. 

## Usage

### PYTHON 
recommended but slower

---
First, run `main.py` to see usage.
```
> python code/main.py -h
usage: main.py [-h] [--detail] [--solution] [-M] data_folder

positional arguments:
  data_folder  a folder containing `A.csv`, `b.csv` and `c.csv` [and
               `x_star.csv`]

optional arguments:
  -h, --help   show this help message and exit
  --detail     show detailed information each iteration
  --solution   show solution if solved
  -M           use the big-M form
```
---
To run *example 1* directly with Big-M, use the following command.
```c
> python code/main.py data/example1
total solving wall time = 0.029999732971191406 sec
status                  = DONE : successfully solved
primal-dual optimal solution:
optimal objective value = 3.0000
numbers of iterations   = 30
solver terminated successfully
```
---
To run *example 1* with Big-M in detail, use the following command, which will be slower.
```c
> python code/main.py data/example1 -M --detail 
k=  1,f(x)=273.231411,|dx|=1.144e+00,|dlam|=3.760e+01,|ds|=6.088e+01,|F(x)_d|=5.569e+01,|F(x)_p|=5.723e+00,|F(x)_0|=1.718e+02,alpha_p=5.905e-01,alpha_d=5.905e-01,
k=  2,f(x)=156.695938,|dx|=5.617e-01,|dlam|=5.663e+00,|ds|=2.540e+01,|F(x)_d|=2.281e+01,|F(x)_p|=2.343e+00,|F(x)_0|=7.777e+01,alpha_p=5.905e-01,alpha_d=5.905e-01,
k=  3,f(x)=82.995895,|dx|=2.964e-01,|dlam|=1.730e+00,|ds|=1.344e+01,|F(x)_d|=9.339e+00,|F(x)_p|=9.597e-01,|F(x)_0|=3.527e+01,alpha_p=5.905e-01,alpha_d=5.905e-01,
k=  4,f(x)=41.268022,|dx|=1.503e-01,|dlam|=1.531e+00,|ds|=6.088e+00,|F(x)_d|=3.824e+00,|F(x)_p|=3.930e-01,|F(x)_0|=1.558e+01,alpha_p=5.314e-01,alpha_d=4.812e-01,
k=  5,f(x)=22.637005,|dx|=2.205e-01,|dlam|=7.740e-01,|ds|=3.243e+00,|F(x)_d|=1.984e+00,|F(x)_p|=1.841e-01,|F(x)_0|=7.540e+00,alpha_p=5.314e-01,alpha_d=4.446e-01,
k=  6,f(x)=13.075237,|dx|=2.837e-01,|dlam|=2.870e-01,|ds|=1.825e+00,|F(x)_d|=1.102e+00,|F(x)_p|=8.628e-02,|F(x)_0|=3.601e+00,alpha_p=4.783e-01,alpha_d=3.855e-01,
k=  7,f(x)=8.685073,|dx|=2.792e-01,|dlam|=1.195e-01,|ds|=1.133e+00,|F(x)_d|=6.772e-01,|F(x)_p|=4.501e-02,|F(x)_0|=1.894e+00,alpha_p=4.783e-01,alpha_d=3.941e-01,
k=  8,f(x)=6.199979,|dx|=2.216e-01,|dlam|=9.355e-02,|ds|=6.933e-01,|F(x)_d|=4.103e-01,|F(x)_p|=2.348e-02,|F(x)_0|=1.033e+00,alpha_p=4.783e-01,alpha_d=4.153e-01,
k=  9,f(x)=4.789300,|dx|=2.491e-01,|dlam|=6.271e-02,|ds|=4.062e-01,|F(x)_d|=2.399e-01,|F(x)_p|=1.225e-02,|F(x)_0|=6.020e-01,alpha_p=4.783e-01,alpha_d=4.389e-01,
k= 10,f(x)=3.992972,|dx|=2.540e-01,|dlam|=3.623e-02,|ds|=2.242e-01,|F(x)_d|=1.346e-01,|F(x)_p|=6.392e-03,|F(x)_0|=3.346e-01,alpha_p=4.783e-01,alpha_d=4.573e-01,
k= 11,f(x)=3.545955,|dx|=2.216e-01,|dlam|=1.857e-02,|ds|=1.174e-01,|F(x)_d|=7.306e-02,|F(x)_p|=3.335e-03,|F(x)_0|=1.788e-01,alpha_p=4.783e-01,alpha_d=4.607e-01,
k= 12,f(x)=3.296468,|dx|=1.627e-01,|dlam|=9.660e-03,|ds|=6.070e-02,|F(x)_d|=3.940e-02,|F(x)_p|=1.740e-03,|F(x)_0|=9.501e-02,alpha_p=4.783e-01,alpha_d=4.649e-01,
k= 13,f(x)=3.158918,|dx|=1.043e-01,|dlam|=5.026e-03,|ds|=3.135e-02,|F(x)_d|=2.108e-02,|F(x)_p|=9.076e-04,|F(x)_0|=5.025e-02,alpha_p=4.783e-01,alpha_d=4.695e-01,
k= 14,f(x)=3.084294,|dx|=6.118e-02,|dlam|=2.587e-03,|ds|=1.623e-02,|F(x)_d|=1.119e-02,|F(x)_p|=4.735e-04,|F(x)_0|=2.645e-02,alpha_p=4.783e-01,alpha_d=4.731e-01,
k= 15,f(x)=3.044394,|dx|=3.403e-02,|dlam|=1.334e-03,|ds|=8.507e-03,|F(x)_d|=5.894e-03,|F(x)_p|=2.470e-04,|F(x)_0|=1.387e-02,alpha_p=4.783e-01,alpha_d=4.754e-01,
k= 16,f(x)=3.023281,|dx|=1.837e-02,|dlam|=6.903e-04,|ds|=4.506e-03,|F(x)_d|=3.092e-03,|F(x)_p|=1.289e-04,|F(x)_0|=7.259e-03,alpha_p=4.783e-01,alpha_d=4.767e-01,
k= 17,f(x)=3.012180,|dx|=9.758e-03,|dlam|=3.584e-04,|ds|=2.370e-03,|F(x)_d|=1.618e-03,|F(x)_p|=6.723e-05,|F(x)_0|=3.793e-03,alpha_p=4.783e-01,alpha_d=4.775e-01,
k= 18,f(x)=3.006363,|dx|=5.139e-03,|dlam|=1.865e-04,|ds|=1.242e-03,|F(x)_d|=8.455e-04,|F(x)_p|=3.508e-05,|F(x)_0|=1.981e-03,alpha_p=4.783e-01,alpha_d=4.779e-01,
k= 19,f(x)=3.003322,|dx|=2.695e-03,|dlam|=9.715e-05,|ds|=6.495e-04,|F(x)_d|=4.415e-04,|F(x)_p|=1.830e-05,|F(x)_0|=1.034e-03,alpha_p=4.783e-01,alpha_d=4.781e-01,
k= 20,f(x)=3.001734,|dx|=1.409e-03,|dlam|=5.065e-05,|ds|=3.393e-04,|F(x)_d|=2.304e-04,|F(x)_p|=9.547e-06,|F(x)_0|=5.394e-04,alpha_p=4.783e-01,alpha_d=4.782e-01,
k= 21,f(x)=3.000905,|dx|=7.363e-04,|dlam|=2.641e-05,|ds|=1.771e-04,|F(x)_d|=1.202e-04,|F(x)_p|=4.980e-06,|F(x)_0|=2.815e-04,alpha_p=4.783e-01,alpha_d=4.782e-01,
k= 22,f(x)=3.000472,|dx|=3.844e-04,|dlam|=1.378e-05,|ds|=9.243e-05,|F(x)_d|=6.274e-05,|F(x)_p|=2.598e-06,|F(x)_0|=1.468e-04,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 23,f(x)=3.000246,|dx|=2.006e-04,|dlam|=7.186e-06,|ds|=4.823e-05,|F(x)_d|=3.273e-05,|F(x)_p|=1.356e-06,|F(x)_0|=7.661e-05,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 24,f(x)=3.000129,|dx|=1.047e-04,|dlam|=3.749e-06,|ds|=2.516e-05,|F(x)_d|=1.708e-05,|F(x)_p|=7.072e-07,|F(x)_0|=3.997e-05,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 25,f(x)=3.000067,|dx|=5.462e-05,|dlam|=1.956e-06,|ds|=1.313e-05,|F(x)_d|=8.909e-06,|F(x)_p|=3.689e-07,|F(x)_0|=2.085e-05,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 26,f(x)=3.000035,|dx|=2.850e-05,|dlam|=1.020e-06,|ds|=6.850e-06,|F(x)_d|=4.648e-06,|F(x)_p|=1.925e-07,|F(x)_0|=1.088e-05,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 27,f(x)=3.000018,|dx|=1.487e-05,|dlam|=5.323e-07,|ds|=3.574e-06,|F(x)_d|=2.425e-06,|F(x)_p|=1.004e-07,|F(x)_0|=5.675e-06,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 28,f(x)=3.000010,|dx|=7.756e-06,|dlam|=2.777e-07,|ds|=1.864e-06,|F(x)_d|=1.265e-06,|F(x)_p|=5.239e-08,|F(x)_0|=2.961e-06,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 29,f(x)=3.000005,|dx|=4.046e-06,|dlam|=1.449e-07,|ds|=9.726e-07,|F(x)_d|=6.600e-07,|F(x)_p|=2.733e-08,|F(x)_0|=1.545e-06,alpha_p=4.783e-01,alpha_d=4.783e-01,
k= 30,f(x)=3.000003,|dx|=2.111e-06,|dlam|=7.558e-08,|ds|=5.074e-07,|F(x)_d|=3.443e-07,|F(x)_p|=1.426e-08,|F(x)_0|=8.059e-07,alpha_p=4.783e-01,alpha_d=4.783e-01,
total solving wall time = 0.11339664459228516 sec
status                  = DONE : successfully solved
primal-dual optimal solution:
optimal objective value = 3.0000
numbers of iterations   = 30
solver terminated successfully
```
---
To run *example 1* with Big-M and show the solution, use the following command.
```c
> python code/main.py data/example1 --solution
total solving wall time = 0.030968427658081055 sec
status                  = DONE : successfully solved
primal-dual optimal solution:
optimal objective value = 3.0000
numbers of iterations   = 30
solution x              = {
    1.0000,     1.6667,     1.3333,     0.0000,     0.0000,     0.0000
}
solver terminated successfully
```
---
To run an extra data piece *lpi_cplex2*, use the following command.
```c
> python code/main.py data/extra/lpi_cplex2 -M 
total solving wall time = 0.4212634563446045 sec
status                  = VIOLATED : infeasible
solver terminated successfully
```
---
To run an extra data piece *kb2*, use the following command.
```c
> python code/main.py data/extra/kb2 -M 
total solving wall time = 0.22266864776611328 sec
status                  = UNBOUNDED : feasible but unbounded
solver terminated successfully
```
---
To test all extra data, use the following command.
```c
> python code/test_extra.py
```
### CPP 
---
Super fast!

**May not support platforms other than LINUX**

First, run `make` and `./ipm` to build and run.
```
> cd code/cpp
> make
> ./ipm   
usage: ./ipm <data_folder> [--detail]
```
---
To run *example 1* directly, use the following command.
```c
> ./ipm ../../data/example1
total solving wall time = 0.0002217
status                  = DONE
primal-dual optimal solution:
optimal objective value = 3
numbers of iterations   = 32
solver terminated successfully
```
---
To run *example 1* in detail, use the following command, which will be slower.
```c
> ./ipm ../../data/example1 --detail
k=  1,f(x)=7.323404,|dx|=1.433e+00,|dlam|=1.241e+00,|ds|=1.896e+00,|F(x)_d|=1.795e+00,|F(x)_p|=7.112e+00,|F(x)_0|=5.066e+00,alpha_p=4.783e-01,alpha_d=4.783e-01,
k=  2,f(x)=5.960883,|dx|=8.590e-01,|dlam|=3.502e-01,|ds|=1.011e+00,|F(x)_d|=9.363e-01,|F(x)_p|=3.711e+00,|F(x)_0|=2.873e+00,alpha_p=4.783e-01,alpha_d=4.783e-01,
k=  3,f(x)=5.056838,|dx|=5.696e-01,|dlam|=8.293e-02,|ds|=7.453e-01,|F(x)_d|=4.885e-01,|F(x)_p|=1.936e+00,|F(x)_0|=1.646e+00,alpha_p=4.305e-01,alpha_d=4.305e-01,
k=  4,f(x)=4.449755,|dx|=4.072e-01,|dlam|=6.616e-02,|ds|=5.313e-01,|F(x)_d|=2.782e-01,|F(x)_p|=1.103e+00,|F(x)_0|=1.001e+00,alpha_p=4.305e-01,alpha_d=4.305e-01,
k=  5,f(x)=3.979717,|dx|=4.177e-01,|dlam|=4.801e-02,|ds|=3.476e-01,|F(x)_d|=1.584e-01,|F(x)_p|=6.279e-01,|F(x)_0|=6.011e-01,alpha_p=4.305e-01,alpha_d=4.269e-01,
k=  6,f(x)=3.635189,|dx|=4.005e-01,|dlam|=2.994e-02,|ds|=2.117e-01,|F(x)_d|=9.081e-02,|F(x)_p|=3.576e-01,|F(x)_0|=3.561e-01,alpha_p=4.305e-01,alpha_d=4.254e-01,
k=  7,f(x)=3.397363,|dx|=3.324e-01,|dlam|=1.637e-02,|ds|=1.221e-01,|F(x)_d|=5.218e-02,|F(x)_p|=2.037e-01,|F(x)_0|=2.081e-01,alpha_p=3.874e-01,alpha_d=3.779e-01,
k=  8,f(x)=3.257005,|dx|=2.519e-01,|dlam|=9.690e-03,|ds|=7.471e-02,|F(x)_d|=3.246e-02,|F(x)_p|=1.248e-01,|F(x)_0|=1.293e-01,alpha_p=3.874e-01,alpha_d=3.785e-01,
k=  9,f(x)=3.163411,|dx|=1.773e-01,|dlam|=6.709e-03,|ds|=4.739e-02,|F(x)_d|=2.017e-02,|F(x)_p|=7.643e-02,|F(x)_0|=8.000e-02,alpha_p=3.874e-01,alpha_d=3.805e-01,
k= 10,f(x)=3.102588,|dx|=1.186e-01,|dlam|=4.505e-03,|ds|=3.059e-02,|F(x)_d|=1.250e-02,|F(x)_p|=4.682e-02,|F(x)_0|=4.934e-02,alpha_p=3.874e-01,alpha_d=3.825e-01,
k= 11,f(x)=3.063838,|dx|=7.667e-02,|dlam|=2.940e-03,|ds|=1.944e-02,|F(x)_d|=7.717e-03,|F(x)_p|=2.868e-02,|F(x)_0|=3.036e-02,alpha_p=3.874e-01,alpha_d=3.842e-01,
k= 12,f(x)=3.039494,|dx|=4.855e-02,|dlam|=1.877e-03,|ds|=1.221e-02,|F(x)_d|=4.753e-03,|F(x)_p|=1.757e-02,|F(x)_0|=1.865e-02,alpha_p=3.874e-01,alpha_d=3.853e-01,
k= 13,f(x)=3.024343,|dx|=3.036e-02,|dlam|=1.181e-03,|ds|=7.596e-03,|F(x)_d|=2.921e-03,|F(x)_p|=1.076e-02,|F(x)_0|=1.145e-02,alpha_p=3.874e-01,alpha_d=3.861e-01,
k= 14,f(x)=3.014969,|dx|=1.883e-02,|dlam|=7.356e-04,|ds|=4.699e-03,|F(x)_d|=1.793e-03,|F(x)_p|=6.593e-03,|F(x)_0|=7.019e-03,alpha_p=3.874e-01,alpha_d=3.866e-01,
k= 15,f(x)=3.009191,|dx|=1.162e-02,|dlam|=4.553e-04,|ds|=2.897e-03,|F(x)_d|=1.100e-03,|F(x)_p|=4.039e-03,|F(x)_0|=4.303e-03,alpha_p=3.874e-01,alpha_d=3.869e-01,
k= 16,f(x)=3.005638,|dx|=7.153e-03,|dlam|=2.807e-04,|ds|=1.781e-03,|F(x)_d|=6.744e-04,|F(x)_p|=2.474e-03,|F(x)_0|=2.637e-03,alpha_p=3.874e-01,alpha_d=3.871e-01,
k= 17,f(x)=3.003457,|dx|=4.394e-03,|dlam|=1.727e-04,|ds|=1.094e-03,|F(x)_d|=4.133e-04,|F(x)_p|=1.516e-03,|F(x)_0|=1.616e-03,alpha_p=3.874e-01,alpha_d=3.872e-01,
k= 18,f(x)=3.002119,|dx|=2.697e-03,|dlam|=1.060e-04,|ds|=6.710e-04,|F(x)_d|=2.533e-04,|F(x)_p|=9.284e-04,|F(x)_0|=9.899e-04,alpha_p=3.874e-01,alpha_d=3.873e-01,
k= 19,f(x)=3.001298,|dx|=1.654e-03,|dlam|=6.504e-05,|ds|=4.114e-04,|F(x)_d|=1.552e-04,|F(x)_p|=5.687e-04,|F(x)_0|=6.065e-04,alpha_p=3.874e-01,alpha_d=3.873e-01,
k= 20,f(x)=3.000796,|dx|=1.014e-03,|dlam|=3.988e-05,|ds|=2.522e-04,|F(x)_d|=9.508e-05,|F(x)_p|=3.484e-04,|F(x)_0|=3.715e-04,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 21,f(x)=3.000487,|dx|=6.212e-04,|dlam|=2.444e-05,|ds|=1.545e-04,|F(x)_d|=5.825e-05,|F(x)_p|=2.134e-04,|F(x)_0|=2.276e-04,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 22,f(x)=3.000299,|dx|=3.806e-04,|dlam|=1.498e-05,|ds|=9.468e-05,|F(x)_d|=3.568e-05,|F(x)_p|=1.307e-04,|F(x)_0|=1.394e-04,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 23,f(x)=3.000183,|dx|=2.332e-04,|dlam|=9.178e-06,|ds|=5.800e-05,|F(x)_d|=2.186e-05,|F(x)_p|=8.008e-05,|F(x)_0|=8.541e-05,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 24,f(x)=3.000112,|dx|=1.429e-04,|dlam|=5.623e-06,|ds|=3.554e-05,|F(x)_d|=1.339e-05,|F(x)_p|=4.906e-05,|F(x)_0|=5.232e-05,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 25,f(x)=3.000069,|dx|=8.752e-05,|dlam|=3.445e-06,|ds|=2.177e-05,|F(x)_d|=8.203e-06,|F(x)_p|=3.005e-05,|F(x)_0|=3.205e-05,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 26,f(x)=3.000042,|dx|=5.361e-05,|dlam|=2.110e-06,|ds|=1.334e-05,|F(x)_d|=5.025e-06,|F(x)_p|=1.841e-05,|F(x)_0|=1.963e-05,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 27,f(x)=3.000026,|dx|=3.284e-05,|dlam|=1.293e-06,|ds|=8.169e-06,|F(x)_d|=3.078e-06,|F(x)_p|=1.128e-05,|F(x)_0|=1.203e-05,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 28,f(x)=3.000016,|dx|=2.012e-05,|dlam|=7.919e-07,|ds|=5.004e-06,|F(x)_d|=1.886e-06,|F(x)_p|=6.908e-06,|F(x)_0|=7.368e-06,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 29,f(x)=3.000010,|dx|=1.233e-05,|dlam|=4.851e-07,|ds|=3.066e-06,|F(x)_d|=1.155e-06,|F(x)_p|=4.232e-06,|F(x)_0|=4.513e-06,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 30,f(x)=3.000006,|dx|=7.550e-06,|dlam|=2.972e-07,|ds|=1.878e-06,|F(x)_d|=7.076e-07,|F(x)_p|=2.592e-06,|F(x)_0|=2.765e-06,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 31,f(x)=3.000004,|dx|=4.625e-06,|dlam|=1.820e-07,|ds|=1.150e-06,|F(x)_d|=4.335e-07,|F(x)_p|=1.588e-06,|F(x)_0|=1.694e-06,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 32,f(x)=3.000002,|dx|=2.833e-06,|dlam|=1.115e-07,|ds|=7.047e-07,|F(x)_d|=2.655e-07,|F(x)_p|=9.727e-07,|F(x)_0|=1.037e-06,alpha_p=3.874e-01,alpha_d=3.874e-01,
k= 33,f(x)=3.000001,|dx|=1.736e-06,|dlam|=6.831e-08,|ds|=4.317e-07,|F(x)_d|=1.627e-07,|F(x)_p|=5.959e-07,|F(x)_0|=6.355e-07,alpha_p=3.874e-01,alpha_d=3.874e-01,
total solving wall time = 0.0568735
status                  = DONE
primal-dual optimal solution:
optimal objective value = 3
numbers of iterations   = 32
solution x              = {
          1     1.66667     1.33333 1.73558e-06 3.63166e-07 3.21217e-07
}
solver terminated successfully
```
---
Since this is a previous version, it may not support some terminal conditions.