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
chi=0.3*np.pi
R_LC=c/omega
theta_pc=np.arcsin(np.sqrt(R/R_LC))
#print("Polar cap opening angle is ",theta_pc)
# Cartesian to spherical, 2D slice

def sphericalize(y,z):
    r=np.sqrt(z**2+y**2)
    theta=np.arccos(z/np.sqrt(z**2+y**2))
    return (r,theta)  
def cartezise(theta,r):
    z=r*np.sin(theta)
    y=r*np.cos(theta)
    return z,y

#theta_max=np.pi+np.arctan(-3+np.sqrt(9+8*(np.tan(chi))**2)/(4*np.tan(chi)))
#theta_max=2.17856
theta_max=1/2*(chi+np.pi+np.arcsin(1/3*np.sin(chi)))
print('thta_max ', theta_max)
s=np.sin(theta_max)*((np.sin(theta_max-chi))**2)
k=R_LC/s

theta_opening=np.arcsin(np.sqrt(s*R/R_LC))
print('Opening angle of the polar cap ', theta_opening)

#dipole field, cartesian
def radius(theta):
    r=k*(np.sin(theta-chi))**2
    return r

def aligned_r(theta):
    r=R_LC*(np.sin(theta))**2
    return r

theta=np.linspace(-np.pi,np.pi,2000)
r=radius(theta)

y,z=cartezise(theta,r)

r_aligned=aligned_r(theta)

y_a,z_a=cartezise(theta, r_aligned)
Light_cylinder=np.full(len(theta),R_LC)
#tilted=np.cos(np.pi/3)*z

plt.plot(y_a/R,z_a/R,'r',label='Aligned')
plt.plot(y/R,z/R,'b',label='Tilted')
#plt.plot(y/R,tilted/R, 'b', label='Tilted')
plt.plot(Light_cylinder/R,z/R, 'k', linewidth=0.6, label='Light cylinder')
plt.plot(-Light_cylinder/R,z/R,'k', linewidth=0.6)

plt.axis('scaled')
plt.xlabel(r'$x~(R_*)$')
plt.ylabel(r'$z~(R_*)$')
#plt.legend(shadow=True)
plt.legend(bbox_to_anchor=(1.04, 1.01))
#plt.title('Last closed field lines for tilted and aligned rotators')
plt.tight_layout()
plt.savefig("Tilted.png", dpi=300)
plt.show()

#f, ax = plt.subplots(1)
#ax.plot(y/R,z/R)
#ax.plot(Light_cylinder/R,z/R)
#ax.plot(-Light_cylinder/R,z/R)
#ax.plot(y/R,tilted/R)
#plt.show(f)

