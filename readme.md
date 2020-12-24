# Interior Point Method (Linear Programming) Implementation

This is an implementation of IPM in CPP and PYTHON.

---
Subfolders are below:

- `code`: contains source code in python and cpp.

- `data`: contains several dataset of matrix A,b,c.
    
    - `data/LP_matlab`: contains raw matlab dataset from [COAP test  problems](http://users.clas.ufl.edu/hager/coap/format.html).

- `verification`: contains cplex script for verification and labeling.

- `report`: contains some notes.

- `appendix`: contains a special dataset.

---
## PYTHON demo
cd root folder (here), and run 

```bash
python code/main3.py
```

This will start benchmark test on all test cases. To run on a specific test case and show details, uncomment line 246, 247, 248
```python
    A, b, c, x_star = read_data(r"data/LP_MATLAB\lotfi.mat")
    x, tgt, res, tc, status = IPM(A, b, c, detail=True)
    exit()
```
and run again.

---
## CPP demo
cd to cpp folder and make
```bash
cd code/cpp
make
```
Then run
```bash
./ipm_v3  ../../data/data_matlab/stocfor1 --detail
```
to show the detailed solution.