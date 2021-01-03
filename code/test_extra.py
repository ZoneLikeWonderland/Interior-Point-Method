import os
import subprocess

f = open("code/label.txt").read().split("\n")

data_path = "data/extra"
for k, (i, j) in enumerate(zip(sorted(os.listdir(data_path)), f)):
    # if k < 26:
    #     continue
    print(k+1, i)
    continue
    j = j.split()[0]
    try:
        j = float(j)
    except:
        j = str(j)
    print("answer", j)
    output = subprocess.check_output("python code/main.py {} -M".format(os.path.join(data_path, i))).decode()
    print(output)

    yanzhong = True
    yanzhong = False

    if isinstance(j, float):
        if "DONE" not in output:
            print("WARNING not done")
            if yanzhong:
                exit(-1)
    else:
        if j not in output:
            print("WARNING not", j)
            if yanzhong:
                exit(-1)
