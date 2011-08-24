from math import *             # any function could be used by set()
from numpy import *


fs = open('pitos.csv','rt')
datos=genfromtxt(fs, dtype=float, delimiter=' ', skip_header=2)
n1=datos.shape[0]
denom=zeros(n1)
results=zeros((n1, 5))
LHS=zeros(n1)
RHS=zeros(n1)
LHS[:]=datos[:, 2]+datos[:, 3]+datos[:, 4]+datos[:, 5]
RHS[:]=datos[:, 1]+datos[0, 2]+datos[0, 3]+datos[0, 4]+datos[0, 5]
for i in range(n1):
    denom[i]=max(1, (LHS[i]+RHS[i])/2.0 )
    results[i,4]=fabs(LHS[i]-RHS[i])/denom[i]
results[:, 0]=datos[:, 0]
results[:, 1]=datos[:, 3]
results[:, 2]=datos[:, 4]
results[:, 3]=datos[:, 5]
savetxt('myfile.csv', results)
