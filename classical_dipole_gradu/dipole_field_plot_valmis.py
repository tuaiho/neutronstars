#Plotting a dipole field

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Defining the functions
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**5     #Magnetic field at the equator on surface (Tesla)
R= 12000       # Stellar radius (m)
omega = 2.0   #Angular velocity of the star (rad/s)


# Cartesian to spherical, 2D slice
def sphericalize(y,z):
    r=np.sqrt(z**2+y**2)
    theta=np.arccos(z/np.sqrt(z**2+y**2))
    return (r,theta)  

#dipole field, cartesian
def B(z,y):
    r,theta=sphericalize(y,z)
    if r>R:
        B_r=B_0*R**3/r**3*2*np.cos(theta)
        B_theta=B_0*R**3/r**3*np.sin(theta)
    else:
        B_r=0
        B_theta=0
    B_z=B_r*np.cos(theta)-B_theta*np.sin(theta)
    B_y=B_r*np.sin(theta)+B_theta*np.cos(theta)
    return [B_z,B_y]

print(B(13000,0))


#corotation charge density

def rho(z,y):
    r,theta=sphericalize(y,z)
    k=-omega*B_0*R**3/(2*np.pi*c)
    if r>R:
        roo=k/r**3*(3*(np.cos(theta))**2-1)/(1.60217733*10**(-19))*10**(-6) # units should be e/cm^3
    else:
        roo=0
    return roo


#Creating the data
# This should be updated to 3D

z=np.linspace(-40000,40000,1000)
y=np.linspace(0,40000,500)

B_z=np.zeros((len(z),len(y))) #empty grid
B_y=np.zeros((len(z),len(y))) #empty grid
#B_abs=np.zeros((len(z),len(y))) #empty grid
for i in range(len(y)):
    for j in range(len(z)):
        B_z[j,i] = B(z[j],y[i])[0]
        B_y[j,i] = B(z[j],y[i])[1]
#        B_abs=np.sqrt((B_z[j,i])**2+(B_y[i,j])**2).item()



ourplot = plt.streamplot(y,z,B_y,B_z, broken_streamlines=False, color=B_z, linewidth=0.5,cmap='jet')
plt.xlabel('Y position (km)')
plt.ylabel('Z position (km)')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()