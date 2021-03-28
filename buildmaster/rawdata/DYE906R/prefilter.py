"""
script converting the data from Table 1 https://arxiv.org/pdf/2103.04024.pdf (reported in ./data_paper.dat)
differential in invarant mass M and Feynman-x xF
in distributions differential in invariant mass M and hadronic rapidity Y, 
using Eqs.(4.6),(4.7) of https://arxiv.org/pdf/1009.5691.pdf.
Output file data_E906.dat 
\sqrt(s)   Y   xt  M   ratio   stat    sys
"""
import numpy as np

xt, xb, M, pT, ratio, stat, sys = np.loadtxt("./data_paper.dat", usecols=(2,3,4,5,6,7,8), unpack=True)

# proton mass
m = 0.938
# beam energy
E = 160
# com energy
s = 2*m**2 + 2*E*m

xF = xb - xt
tau = M**2/s
pt = pT/M

# hadronic rapidity, Eq.(4.6) https://arxiv.org/pdf/1009.5691.pdf
Y = 0.5*np.log( (np.sqrt(xF**2 + 4*tau*(1+pt**2)) + xF)/(np.sqrt(xF**2 + 4*tau*(1+pt**2)) - xF) )

# jacobian, Eq.(4.7 https://arxiv.org/pdf/1009.5691.pdf)
j = np.sqrt(xF**2 +4*tau*(1+pt**2))
ratio *= j

output=open("data_E906.dat", "w+")
ndata = 6

for i in range (0,ndata):
    output.write(str(np.sqrt(s)) + "\t" + str(Y[i]) + "\t" + str(xt[i]) + "\t" + str(M[i]) + "\t" + str(ratio[i]) 
                   + "\t" + str(stat[i]) + "\t" + str(sys[i]) + "\t")
    output.write("\n")
output.close