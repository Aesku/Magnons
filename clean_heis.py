import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

#program parameters
tol = 1e-5

#interaction parameters
J = [2.4, 0.6, -0.4] #Exchange interaction: 1st, 2nd, 3rd neighbors

#lattice parameters
a = 7.0
c = 2*a
nn = [6,12,18] #number of 1st, 2nd, 3rd... neighbors

#lattice vectors
    #Hexagonal
a1 = np.array([a,0.0,0.0])
a2 = np.array([-0.5*a,0.5*math.sqrt(3)*a,0.0])
a3 = np.array([0.0,0.0,c])
    #Square
#a1 = np.array([a,0.0,0.0])
#a2 = np.array([0.0,a,0.0])
#a3 = np.array([0.0,0.0,c])

lattice = [a1,a2,a3]
ng = [] #empty list that will contain the1st neighbors coordinates

#Creating reciprocal lattice
V = np.linalg.det(lattice) #volume
b1 = 2.0*math.pi*np.cross(a2,a3)/V
b2 = 2.0*math.pi*np.cross(a3,a1)/V
b3 = 2.0*math.pi*np.cross(a1,a2)/V

#Calculating critical points 1st Brillouin Zone
    #Hexagonal lattice
k0 = np.array([0,0,0])
km = 0.5*b1
kk = (2.0/3.0)*b1 + (1.0/3.0)*b2
    #Square lattice
#km = 0.5*b1 + 0.5*b2
#kx = 0.5*b2

#sampling neighbors
for m in np.arange(-5,5):
    for n in np.arange(-5,5):
        r = m*a1 + n*a2
        d = np.linalg.norm(r)
        if math.sqrt((a-d)*(a-d)) < tol: #1st
            ng.append(r)
n_1n = len(ng)
#print(ng)
plt.plot(0,0,'bo')
for i in range(n_1n):
    plt.plot(ng[i][0],ng[i][1],'ro')
plt.show()

#Calculating symmetry path
k = [] #Empty list that will store all the k's for the path
    #From 0 to Km
for i in range(50):
    k.append((i/50.0)*(km - k0))
    
    #From Km to Kk
for i in range(50):
    k.append(km + (i/50.0)*(kk - km))

    #From Kk to 0
for i in range(50):
    k.append(kk - (i/50.0)*kk)

#Defining funcion to calculate energy
def energy(k, ng): 
    s = 0.0
    for i in np.arange(nn[0]):    #Adding 1st neighbor contribution
        s += math.cos(np.dot(k,ng[i]))
    return 0.5*J[0]*(1.0*nn[0]/4.0 - nn[0] + s)

#Calculating dispersion
Xaxis = []
Yaxis = []
for count in range(len(k)):
    Yaxis.append(energy(k[count],ng))
    Xaxis.append(count)


#plotting dispersion
sym_hex = [r'$\Gamma$','M','K',r'$\Gamma$']
sym_sq = [r'$\Gamma$','M','X',r'$\Gamma$']
plt.scatter(Xaxis,Yaxis)
plt.xticks([])
ax = plt.gca()
Km = 50
Kk = 100
Kx = 100
# set the positions where we want ticks to appear
ax.xaxis.set_ticks([0,50,100,150])
# set what will actually be displayed at each tick.
ax.xaxis.set_ticklabels(sym_hex)
plt.axvline(x = 0, color = 'g', label = r'\Gamma')
plt.axvline(x = Km, color = 'g', label = 'M')
plt.axvline(x = Kx, color = 'g', label = 'X')
plt.axvline(x = 150, color = 'g', label = r'\Gamma')
plt.title("Hexagonal lattice")
plt.ylabel("meV")
plt.plot(Xaxis,Yaxis)
plt.savefig('Hex1.png', format='png')
plt.show()

print("program finished with 0 errors")
