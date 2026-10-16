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
theta_pc=np.arcsin(np.sqrt(R/R_LC))
print("Polar cap opening angle is ",theta_pc)



def polararclength(theta, r, a, b):
    """
    Computes the arclength of the given curve
    defined 
    a: int
        The lower limit to integrate from
    b: int
        The upper limit to integrate to
    """
    if a>=b:
        return 0
    bounds = (theta >= a) & (theta <= b)

    return np.trapezoid(np.sqrt(r[bounds]**2+ (np.gradient(r[bounds], theta[bounds])) ** 2),theta[bounds])

def rho(theta):     #Corotation charge density as a function of theta, along the last closed line
    e_0=8.8541878188*10**(-12)
    k=-2*e_0*omega*B_0*R**3
    r=R_LC*(np.sin(theta))**2
    roo=k/r**3*(3*(np.cos(theta))**2-1)#/(1.60217733*10**(-19))*10**(-6) # units should be e/cm^3
    return roo
print('surface charge', rho(theta_pc))

#x=np.linspace(0,100,1000)
#y=x**2
#print(arclength(x,y,2,3))

#theta=np.linspace(theta_pc,np.pi/2,1000) #Even spacing, logspace could be better
theta=np.logspace(np.log10(theta_pc),np.log10(np.pi/2),1000)   #More resolution close to the star

#print(theta)

r=R_LC*(np.sin(theta))**2
charge=rho(theta)
l=np.zeros(len(theta))
for i in range(1,len(theta)):      #Start from >1 distance (index 1, not 0) to avoid negative due to numerics
    b=theta[i]
    l[i] = polararclength(theta,r,theta_pc,b)/R     # The arclength to current theta


#print("Length of the last open field line from the surface to light cylinder (km): ", polararclength(theta,r,theta_pc,np.pi/2)/1000)
print("Light cylinder radius (km): ", R_LC/1000)

#plt.semilogx(l,-charge)    #Only l logarithmic
plt.loglog(l, -charge)
plt.grid('minor')
plt.title(r"Corotation charge density along the last closed field line with $\Omega= %i\pi \, \frac{rad}{s}$" %n)
plt.xlabel("Distance (m) from the surface along the last closed field line")
plt.ylabel("charge density (Q/m^3) (absolute value)")
plt.show()