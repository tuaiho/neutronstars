

import numpy as np
import matplotlib.pyplot as plt


# Defining the functions
c=299792458.0 #m/s
mu_0=1.256637061*10**(-6)   # N/A^2

B_0= 10**8     #Magnetic field at the equator on surface (Tesla)
R= 12000       # Stellar radius (m)
omega = 2*np.pi   #Angular velocity of the star (rad/s)
R_LC=c/omega
theta_pc_aligned=np.arcsin(np.sqrt(R/R_LC))



def opening_angle(chi):
    theta_max=1/2*(chi+np.pi+np.arcsin(1/3*np.sin(chi)))
    s=np.sin(theta_max)*((np.sin(theta_max-chi))**2)
    #k =R_LC/s
    theta_opening=np.arcsin(np.sqrt(s*R/R_LC))
    return theta_opening

chi=np.linspace(0,np.pi/2,1000)
openings= opening_angle(chi)
print(opening_angle(np.pi/2)/theta_pc_aligned)

plt.plot(chi,openings/theta_pc_aligned)
plt.xlabel(r'tilting angle $ \chi~\mathrm{(rad)}$')
plt.ylabel(r'$\theta_{\text{pc, tilted}} ~/~ \theta_{\text{pc, aligned}}$ ')
plt.tight_layout()
plt.savefig('pc_size_vs_tilting.png', dpi=300)
plt.show()