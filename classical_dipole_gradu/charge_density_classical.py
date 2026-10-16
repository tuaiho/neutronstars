#just one field line

# Charge density along the last closed field line 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Defining the functions
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**5     #Magnetic field at the equator on surface (Tesla)
R= 12000       # Stellar radius (m)
n=2
omega = n*np.pi  #Angular velocity of the star (rad/s)
R_LC=c/omega

theta_pc_c=np.arcsin(np.sqrt(R/R_LC))
print("Polar cap opening angle is ",theta_pc_c)



def polararclength(r, a, b):
    bounds = (r >= a) & (r <= b)
    r=r[bounds]
    theta_b=theta[bounds]
    Dr_square=R_LC**2*4*np.sin(theta_b)**2*np.cos(theta_b)**2
    integrand= np.sqrt(r**2+Dr_square)
    line=np.trapezoid(integrand,theta_b)
    return line

def rho(r):     #Corotation charge density as a function of theta, along the last closed line
    e_0=8.8541878188*10**(-12)
    k=-2*e_0*omega*B_0*R**3
    theta=np.arcsin(np.sqrt(r/R_LC))
    roo=k/r**3*(3*(np.cos(theta))**2-1)#/(1.60217733*10**(-19))*10**(-6) # units should be e/cm^3
    return roo
print('surface charge', rho(R))

#
#theta=np.linspace(theta_pc_c,np.pi/2,1000) #Even spacing, logspace could be better
#theta=np.logspace(np.log10(theta_pc_c),np.log10(np.pi/2),1000)   #More resolution close to the star
r= np.logspace(np.log10(R),np.log10(R_LC),2000) 
theta=np.arcsin(np.sqrt(r/R_LC))

charge=rho(r)
l=np.zeros(len(r))
for i in range(1,len(r)):      #Start from >1 distance (index 1, not 0) to avoid negative due to numerics
    b=r[i]
    l[i] = polararclength(r, 0, b)     # The arclength to current r


#print("Length of the last open field line from the surface to light cylinder (km): ", polararclength(theta,r,theta_pc_c,np.pi/2)/1000)
print("Light cylinder radius (km): ", R_LC/1000)

#plt.semilogx(l,-charge)    #Only l logarithmic
plt.loglog(l, -charge)
plt.grid('minor')
plt.title(r"Corotation charge density along the last closed field line with $\Omega= %i\pi \, \frac{rad}{s}$" %n)
plt.xlabel("Distance (m) from the surface along the last closed field line")
plt.ylabel("charge density (Q/m^3) (absolute value)")
plt.show()