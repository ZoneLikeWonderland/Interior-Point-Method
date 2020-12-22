import subprocess
import os 
import re
import time 
os.chdir(os.path.dirname(__file__))
path = os.listdir('data/')
fi = "cplex.run"
pattern = re.compile(r'(?<=/).*(\.dat)')
for p in path:

    file_data = ""
    with open(fi, "r", encoding="utf-8") as f:
      for line in f:
        line = re.sub(pattern, p, line)
        file_data += line
    with open(fi,"w",encoding="utf-8") as f:
      f.write(file_data)
    try  :
        start = time.time()
        output = subprocess.check_output("ampl ./cplex.run", shell = True)
        print('file name',p)    
        print('solution by cplex',output.decode().split('\n')[0].split(' ')[-1][:-1])
        print('time',time.time()-start,'s')
        print('-'*30)
    except :
        pass