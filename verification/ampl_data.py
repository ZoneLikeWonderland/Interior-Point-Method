import numpy as np
import os
import scipy.io as sio 


def trans_csv(root):
    path = os.listdir(root)
    for p in path:
        if p[-3:] == 'mat':
            data = sio.loadmat(os.path.join(root, p))
            try:
                A = data["Problem"]["A"][0][0].todense()
                b = data["Problem"]["b"][0][0]
                c = data["Problem"]["c"][0][0]
            except:
                A = data["A"].todense().tolist()
                b = data["b"].tolist()
                c = data["c"].tolist()
            os.makedirs('../data/data_matlab/'+p[:-4], exist_ok=True)
            np.savetxt('../data/data_matlab/'+p[:-4]+'/A.csv',A, fmt='%.18e', delimiter=',', newline='\n')
            np.savetxt('../data/data_matlab/'+p[:-4]+'/b.csv',b, fmt='%.18e', delimiter=',', newline='\n')
            np.savetxt('../data/data_matlab/'+p[:-4]+'/c.csv',c, fmt='%.18e', delimiter=',', newline='\n')

def trans_dat(root):
    path = os.listdir(root)
    # root = ['adlittle.mat']
    for p in path:
        if p[-3:] == 'mat':
            data = sio.loadmat(os.path.join(root, p))
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
                for i in range(len(c)):
                    f.write(str(i+1)  +  '  ' +str(c[i][0]) +'\n')
                
                f.write(';\n\n\n')

                f.write('param: b := \n')
                for i in range(len(b)):
                    f.write(str(i+1)  +  '  ' +str(b[i][0]) +'\n')
                f.write(';\n\n\n')    
                        
                f.write('param A : \n')
                s = str([i+1 for i in range(len(c))]).replace(',',' ')[1:-1] + ':=\n'
                for i in range(len(b)):
                    s+= str(i+1) +'  ' + str([j for j in A[i]]).replace(',',' ')[1:-1] +'\n'
                s+= ';'
                f.write(s)


if __name__ == "__main__":
    root = '../data/LP_matlab/'
    trans_dat(root)
