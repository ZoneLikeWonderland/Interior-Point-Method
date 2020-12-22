import numpy as np
import os
import scipy.io as sio 

path = os.listdir('../data/LP_matlab/')
# path = ['adlittle.mat']
for p in path:
    if p[-3:] == 'mat':
        data = sio.loadmat(os.path.join('../data/LP_matlab/', p))
        try:
            A = data["Problem"]["A"][0][0].todense().tolist()
            b = data["Problem"]["b"][0][0].tolist()
            c = data["Problem"]["c"][0][0].tolist()
        except:
            A = data["A"].todense().tolist()
            b = data["b"].tolist()
            c = data["c"].tolist()
        with open('data/'+p[:-3]+'dat', 'w') as f :
            f.write('param m :=' + str(len(b)) +' ;\n')
            f.write('param n :=' + str(len(c)) +' ;\n\n\n')

            f.write('param: c := \n')
            for i in range(len(c)-1):
                f.write(str(i+1)  +  '  ' +str(c[i][0]) +'\n')
            f.write(str(i+2)  +  '  ' +str(c[i][0])+';\n\n\n')

            f.write('param: b := \n')
            for i in range(len(b)-1):
                f.write(str(i+1)  +  '  ' +str(b[i][0]) +'\n')
            f.write(str(i+2)  +  '  ' +str(b[i][0])+';\n\n\n')    
                    
            f.write('param A : \n')
            s = str([i+1 for i in range(len(c))]).replace(',',' ')[1:-1] + ':=\n'
            for i in range(len(b)-1):
                s+= str(i+1) +'  ' + str([j for j in A[i]]).replace(',',' ')[1:-1] +'\n'
            s+= str(i+2) +'  ' + str([j for j in A[i]]).replace(',',' ')[1:-1] +';'
            f.write(s)

