# Aligned dipole field + corotation charge density plotted


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Defining the functions
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**8     #Magnetic field at the equator on surface (Tesla)
R= 12000       # Stellar radius (m)
omega = 2*np.pi   #Angular velocity of the star (rad/s)


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


#corotation charge density

def rho(z,y):
    r,theta=sphericalize(y,z)
    e_0=8.8541878188*10**(-12)
    k=-2*e_0*omega*B_0*R**3
    if r>R:
        roo=k/r**3*(3*(np.cos(theta))**2-1)#/(1.60217733*10**(-19)*10**(6)) # units should be Q/m^3 / e/cm^3
    else:
        roo=0
    return roo*2998     #Converts to statC/cm^3

#print(rho(0,130000))
#Creating the data

z=np.linspace(-5*R,5*R,1000)
y=np.linspace(0,5*R,500)
#Z,Y=np.meshgrid(z,y)
#charge_density=rho(Z,Y)

#Slow method
charge_density=np.zeros((len(z),len(y))) #empty grid

for i in range(len(y)):
    for j in range(len(z)):
        charge_density[j,i] = rho(z[j],y[i])

B_z=np.zeros((len(z),len(y))) #empty grid
B_y=np.zeros((len(z),len(y))) #empty grid

for i in range(len(y)):
    for j in range(len(z)):
        B_z[j,i] = B(z[j],y[i])[0]
        B_y[j,i] = B(z[j],y[i])[1]



min=np.min(charge_density)
max=np.max(charge_density)
print(min)
print(max)
mylevels = np.linspace(min,-min,201)
fieldlines= plt.streamplot(y/R,z/R,B_y,B_z, broken_streamlines=False, linewidth=0.1, color='k')
ourplot = plt.contourf(y/R,z/R,charge_density,levels=mylevels,cmap='seismic')

plt.xlabel(r'x $(R_*)$')
plt.ylabel(r'z $(R_*)$')
#plt.title('Co-rotation charge density')
mycb = plt.colorbar(ourplot,label=r"$\rho_{\mathrm{co}} ~ (\mathrm{statC/cm^3})$")
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig("classical-co-rotationcahrge.png", dpi=300)
plt.show()
