import numpy as np
from scipy.integrate import solve_ivp

### A program that solves the GR polar cap size and plots the charge density along last closed field line
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**8     #Magnetic field at the equator on surface (Tesla)
#B_0=10**12       # Magnetic field in Gauss
R_star= 12000       # Stellar radius (m)
epsilon=1/3     #typical pulsar compactness
n=2
omega = n*np.pi  # Angular velocity of the star (rad/s)
R_LC=c/omega # Light cylinder in m

r_s=1/3*R_star     #in units of stellar radius
# Limits for radius in units of stellar radius
R=R_star    
R_max=R_LC
#print('Light cylinder radius (m)', R_LC)

### GR ###

#intial conditions
G_0=B_0*R**3/R_max**3
dG_0=-3*B_0*R**3/R_max**4
A_0=0


def alpha(r):       #Lapse function
    return np.sqrt(1-r_s/r)

def fun(r, z):
    G, H, A = z
    F=-r*alpha(r)*H-G*(alpha(r)+r_s/(2*alpha(r)*r))
    dG=H
    dH= -(4/r+r_s/(alpha(r)**2*r**2))*H -(2/r**2+r_s/(alpha(r)**2*r**3)-r_s**2/(4*r**4*alpha(r)**4)-2/(alpha(r)**2*r**2))*G
    dA=G/(r*alpha(r)*F)
    return [dG, dH, dA]

sol = solve_ivp(fun, [R_max, R], [G_0, dG_0, A_0], dense_output=True)

#r = np.linspace(R, R_max, 20000)
r= np.logspace(np.log10(R),np.log10(R_max),2000) 

z = sol.sol(r)
G=z[0]
H=z[1]
A=z[2]
F=-r*alpha(r)*H-G*(alpha(r)+r_s/(2*alpha(r)*r))
theta=np.arcsin(np.exp(A))
theta_pc=np.arcsin(np.exp(A[0]))
print('GR polar cap opening angle: ', theta_pc)

def rterm(r):
    return (r**3-2/5*r_s*R**2)/(np.sqrt(1-r_s/r))
def drterm(r):
    return np.sqrt(1-r_s/r)*(30*r**4-35*r**3*r_s+2*r_s**2*R**2)/(10*(r-r_s)**2)

def divE(r,theta,G,H):
    result=alpha(r)*omega/r**2*np.sin(theta)**2*(drterm(r)*G+rterm(r)*H)-omega*(1-2/5*r_s*R**2/r**3)/alpha(r)*(3*np.cos(theta)**2-1)*F
    return result
e_0=8.8541878188*10**(-12)
dive=divE(r,theta,G, H)
rho=e_0*dive
print('GR surface charge', rho[0])
print('GR far away charge', rho[len(r)-1])

# l:S corresponding to r:s
def length(r,a,b):
    bounds = (r >= a) & (r <= b)
    r=r[bounds]
    theta_b=theta[bounds]
    G_b=G[bounds]
    F_b=F[bounds]
    H_b=H[bounds]
    A_b=A[bounds]
    DA=G_b/(r*alpha(r)*F_b)
    Dr_square=(np.exp(-2*A_b)-1)/DA**2
    integrand= np.sqrt(r**2+1/(alpha(r))**2*Dr_square)
    line=np.trapezoid(integrand,theta_b)
    return line


l=np.zeros(len(theta))
for i in range(1,len(theta)):      #Start from >1 distance (index 1, not 0) to avoid negative due to numerics
    b=r[i]
    l[i] = length(r,0,b)


### Classical: ###

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

def rho_c(r):     #Corotation charge density as a function of theta, along the last closed line
    e_0=8.8541878188*10**(-12)
    k=-2*e_0*omega*B_0*R**3
    theta=np.arcsin(np.sqrt(r/R_LC))
    roo=k/r**3*(3*(np.cos(theta))**2-1)#/(1.60217733*10**(-19))*10**(-6) # units should be e/cm^3
    return roo


rho_class=rho_c(r)
print('rho_class maximum ', np.max(rho_class))
print('Classical surface charge', rho_class[0])
print('Classical far away charge', rho_class[len(r)-1])
l_c=np.zeros(len(r))
for i in range(1,len(r)):      #Start from >1 distance (index 1, not 0) to avoid negative due to numerics
    b=r[i]
    l_c[i] = polararclength(r, 0, b)     # The arclength to current r




### Plot: ###
import matplotlib.pyplot as plt

plt.plot(figsize=(12,5))
plt.semilogx(l/R_star, rho)
plt.semilogx(l_c/R_star, rho_class)
#plt.xlabel(r'proper length from the surface of the star along the field line in $R_*$')
plt.xlabel(r'proper length ($R_*$)')
plt.ylabel(r'$\rho_{co}$  $\left(\frac{\mathrm{C}}{\mathrm{m}^3}\right)$')
plt.legend(['GR', 'Classical'], shadow=True)
plt.tight_layout()
#plt.savefig('charge_along_lcfl.eps', format='eps')
#plt.title(r'Corotation charge density $\rho_e$ along the last closed field line')
plt.show()


plt.semilogx(l_c/R_star, rho/rho_class)
plt.xlabel('proper length from the surface of the star along the field line (m)')
plt.ylabel(r'$\frac{\rho_{GR}}{\rho_{class}}$')
#plt.title(r'The ratio $\frac{\rho_{GR}}{\rho_{classical}}$ along the last closed field line')
#plt.savefig('charge_ratio.eps', format='eps')
plt.show()