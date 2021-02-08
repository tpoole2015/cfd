import matplotlib.pyplot as plt
import numpy as np

def areaMach(M, gamma):
    x=2*(1+(gamma-1)*0.5*M*M)/(gamma+1)
    return np.float_power(x, (gamma+1)*0.5/(gamma-1))/M

# given A/A* find M such that areaMach(M,gamma)=ratio, M0 is our initial guess
#def findMach(A_AstarRatio,gamma,M0):

gamma=1.4
x=np.arange(0.2,2,0.1)
plt.plot(x,[areaMach(M, gamma) for M in x])
plt.show()
