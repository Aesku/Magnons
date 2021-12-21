import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

#program parameters
tol = 1e-5
deg = math.pi/180.0

#interaction parameters
J1 = 2.2 #Exchange interaction: 1st neighbors (meV)
S = 3.0/2.0 #Total spin
E0 = 10.12 #meV

#lattice parameters
#a = 18
a = 7.0
c = 2*a
nCr = 3 #number of Cr neighbors
nI = 6 #number of I neighbors

#lattice vectors
    #Hexagonal
a1 = np.array([a,0.0,0.0])
a2 = np.array([-a*math.cos(60*deg), a*math.sin(60*deg),0.0])
a3 = np.array([0.0,0.0,c])
    #Square
#a1 = np.array([a,0.0,0.0])
#a2 = np.array([0.0,a,0.0])
#a3 = np.array([0.0,0.0,c])

lattice = [a1,a2,a3]
print(lattice)
rCr = [] #empty list that will contain the 1st Cr neighbor coordinates
rI = [] #" 1st I neighbor coordinates

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
condCr1 = np.linalg.norm(1.0*a1/3.0 + 2.0*a2/3.0)
condCr2 = np.linalg.norm(5.0*a1/3.0 + 1.0*a2/3.0)
condCr3 = np.linalg.norm(4.0*a1/3.0 + 2.0*a2/3.0)

#Sampling Cr
for m in np.arange(-7,7):
    for n in np.arange(-7,7):
        r = 1.0*m*a1/3.0 + 1.0*n*a2/3.0
        d = np.linalg.norm(r)
        if math.sqrt((condCr1 - d)*(condCr1 - d)) < tol: #1st
            rCr.append(r)
        if math.sqrt((condCr2 - d)*(condCr2 - d)) < tol: #2nd
            rCr.append(r)
        if math.sqrt((condCr3 - d)*(condCr3 - d)) < tol: #3rd
            rCr.append(r)

#plotting structure
    #plotting lattice vectors
data = np.array([[a,0.0,0.0], [-a*math.cos(60*deg), a*math.sin(60*deg),0.0]])
plt.quiver([0,0,0],[0,0,0],data[:, 0], data[:, 1],units='xy', scale=1)
nCr = len(rCr)
for i in range(nCr):
    plt.plot(rCr[i][0], rCr[i][1], 'ro')
plt.title("CrI3 lattice (Only Cr)")
plt.savefig('CrI3.png', format='png')
plt.show()


#Calculating symmetry path
k = [] #Empty list that will store all the k's for the path
    #From 0 to KK
for i in range(50):
    k.append((1.0*i/50.0)*kk)

#Defining auxiliar functions
def complexmodule(z):
    return abs(z)

#Defining funcion to calculate energy
def energy1(k):
    z = 1 + complex(math.cos(np.dot(k,a1)),math.sin(np.dot(k,a1))) + complex(math.cos(np.dot(k,a2)),math.sin(np.dot(k,a2)))
    return E0 + J1*S*complexmodule(z)

def energy2(k):
    z = 1 + complex(math.cos(np.dot(k,a1)),math.sin(np.dot(k,a1))) + complex(math.cos(np.dot(k,a2)),math.sin(np.dot(k,a2)))
    return E0 - J1*S*complexmodule(z)    

#Calculating dispersion
Xaxis1 = []
Yaxis1 = []
Xaxis2 = []
Yaxis2 = []
for count in range(len(k)):
    #1st neighbors
    Yaxis1.append(energy1(k[count])/20.0)
    Xaxis1.append(count)
    Yaxis2.append(energy2(k[count])/20.0)
    Xaxis2.append(count)

#plotting dispersion
sym_hex = [r'$\Gamma$','K']
sym_sq = [r'$\Gamma$','M','X',r'$\Gamma$']
plt.scatter(Xaxis1,Yaxis1)
plt.scatter(Xaxis2,Yaxis2)
#plt.xticks([])
ax = plt.gca()
Kk = 50
# set the positions where we want ticks to appear
ax.xaxis.set_ticks([0,50])
# set what will actually be displayed at each tick.
ax.xaxis.set_ticklabels(sym_hex)
plt.axvline(x = 0, color = 'g', label = r'\Gamma')
plt.axvline(x = Kk, color = 'g', label = 'K')
plt.title("Lado2017")
plt.ylabel("meV")
plt.plot(Xaxis1,Yaxis1)
plt.plot(Xaxis2,Yaxis2)
plt.savefig('Hex1.png', format='png')
plt.show()

print("program finished with 0 errors")
