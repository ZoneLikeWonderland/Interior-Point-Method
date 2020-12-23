import time
import os
import glob
import subprocess

if __name__ == "__main__":
    d = "data/data_matlab/"
    for f in os.listdir(d):
        path = os.path.join(d, f)
        print("try", path)
        try:
            r = subprocess.check_output(["wsl", "-e", "code/cpp/ipm", path])
            print(r.decode())
            print()
        except Exception as e:
            print(e)
