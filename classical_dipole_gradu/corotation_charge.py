import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Defining the functions
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**5     #Magnetic field at the equator on surface (Tesla)
R= 12000       # Stellar radius (m)
omega = 2*np.pi   #Angular velocity of the star (rad/s)


# Cartesian to spherical, 2D slice
def sphericalize(y,z):
    r=np.sqrt(z**2+y**2)
    theta=np.arccos(z/np.sqrt(z**2+y**2))
    return (r,theta)  

#dipole field
def B(r,theta):
    if r>R:
        B_r=B_0*R**3/r**3*2*np.cos(theta)
        B_theta=B_0*R**3/r**3*np.sin(theta)
    else:
        B_r=0
        B_theta=0
    return [B_r,B_theta]

# vector components into cartesian coordinates
def cartesize(V_r, V_theta, theta):
    V_z=V_r*np.cos(theta)-V_theta*np.sin(theta)
    V_y=V_r*np.sin(theta)+V_theta*np.cos(theta)
    return [V_z,V_y]

def magnetic_field(z,y):
    r,theta=sphericalize(z,y)
    B_r,B_theta=B(r,theta)
    return [B_r,B_theta]


#corotation charge density

def rho(z,y):
    r,theta=sphericalize(y,z)
    k=-omega*B_0*R**3/(2*np.pi*c)
    if r>R:
        roo=k/r**3*(3*(np.cos(theta))**2-1)/(1.60217733*10**(-19))*10**(-6) # units should be e/cm^3
    else:
        roo=0
    return roo

print(rho(0,130000))
#Creating the data

z=np.linspace(-40000,40000,500)
y=np.linspace(0,40000,500)
#Z,Y=np.meshgrid(z,y)
#charge_density=rho(Z,Y)

#Slow method
charge_density=np.zeros((len(z),len(y))) #empty grid

for i in range(len(y)):
    for j in range(len(z)):
        charge_density[j,i] = rho(z[j],y[i])

#Testing
#print(charge_density)

#Plotting
#plt.contour(y,z,charge_density)
#plt.xlabel('Y position (m)')
#plt.ylabel('Z position (m)')
#plt.show()

min=np.min(charge_density)
max=np.max(charge_density)
print(min)
print(max)
mylevels = np.linspace(min,-min,201)

ourplot = plt.contourf(y/1000,z/1000,charge_density,levels=mylevels,cmap='seismic')
plt.xlabel('Y position (km)')
plt.ylabel('Z position (km)')
mycb = plt.colorbar(ourplot,label="rho")
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()