# import necessary modules/packages
import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename,type,number):
    '''
    this function returns the magnitude of the distance [kpc],
    magnitude of the velocity [km/s], and mass [solar masses]
    of any given particle 
    '''
    # call ReadFile function 
    time, total_parts, data = Read('MW_000.txt') 

    # conditions for identifying data points of selected type
    if type == 1:
        index = np.where(data['type']==1) # dark matter halo particles
    elif type == 2:
        index = np.where(data['type']==2) # disk particles
    elif type == 3:
        index = np.where(data['type']==3) # bulge particles

    # calculate distance magnitude of selected particle
    distance_mag = np.sqrt((data['x'][index][number-1])**2+(data['y'][index][number-1])**2+(data['z'][index][number-1])**2)
    dist = np.round(distance_mag,3)*u.kpc # round and identify kpc units
    dist_ly = np.round(dist.to(u.lightyear),3) # convert distance to light years
    print('Distance magnitude:', dist, '=', dist_ly) # print results
    
    # calculate velocity magnitude of selected particle
    velocity_mag = np.sqrt((data['vx'][index][number-1])**2+(data['vy'][index][number-1])**2+(data['vz'][index][number-1])**2)
    vel = np.round(velocity_mag,3) # round
    print('Velocity magnitude:', vel, 'km/s') # print

    # find mass of selected particle and print result
    mass = data['m'][index][number-1]
    print('Mass:', mass,'x 10^(10) solar masses')

    return dist, vel, mass
    
ParticleInfo('MW_000.txt',2,100)
