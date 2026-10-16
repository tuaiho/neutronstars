import numpy as np 
import matplotlib.pyplot as plt 

#Define the initial electromagnetic field
#Initial values

B_0 = 10**12        # equatorial surface magnetic field in gauss
R=1200000            # stellar radius in cm
X= np.pi/3          # tilt angle [0, pi/2]
Omega=2*np.pi       # angular velocity of the star in rad/s
c=29979245800       # speed of light cm/s

#Cartesian-spherical transformations from other files
#Is it smarter to define static B and update it with differenf function, or write time-dependent B?

#radial component:
def B_r(r,theta,phi,t):
    radial=B_0*R**3/r**3
    angular=2*(np.sin(X)*np.sin(theta)*np.cos(phi-Omega*t)+np.cos(X)*np.cos(theta))
    b_r=radial*angular
    return b_r
# poloidal component:
def B_theta(r,theta,phi,t):
    b_theta=B_0*R**3/r**3*(np.cos(X)*np.sin(theta)-np.sin(X)*np.cos(theta)*np.cos(phi-Omega*t))
    return b_theta
# toroidal component:
def B_phi(r,theta,phi,t):
    radial=B_0*R**3/r**3
    angular= np.sin(X)*np.sin(phi-Omega*t)
    return radial*angular




# co-rotation velocity
def v_co(r,theta):
    return Omega*r*np.sin(theta)

# co-rotation electric field (E_phi=0 always)
def E_r(r,theta,phi,t):
    e_r=v_co(r,theta)/c*B_theta(r,theta,phi,t)
    return e_r

def E_theta(r,theta,phi,t):
    e_theta=v_co(r,theta)/c*B_r(r,theta,phi,t)
    return e_theta


#testing everything in the same point
r=1.2*R
theta=np.pi/6
phi=0
t=1
print("In point r theta phi t  = ", r , theta , phi ,t)

print("B_r = ", B_r(r,theta,phi,t)/10**12, " 10^12 G")
print("B_theta = ", B_theta(r,theta,phi,t)/10**12, " 10^12 G")
print("B_phi = ", B_phi(r,theta,phi,t)/10**12, " 10^12 G")

print("v_co = ", v_co(r,theta), " cm/s")

print("E_r = ", E_r(r,theta,phi,t)/10**7, "10^7 G")
print("E_theta = ", E_theta(r,theta,phi,t)/10**7, "10^7 G")

# Create the cartesian grid points
# Express the points in spherical coords to calculate B, E and v_co in the point
# Express B, E and v_co in cartesian coords
# plot
