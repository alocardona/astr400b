# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
import matplotlib.pyplot as plt
from astropy import constants as const # import astropy constants
from CenterOfMass import CenterOfMass
from ReadFile import Read

class MassProfile:
# Class to define COM position and velocity properties 
# of a given galaxy and simulation snapshot

    def __init__(self, galaxy, snap):
        ''' 
        Class to calculate the mass enclosed by a given radius, Hernquist mass,
        and circular velocity.
         
        PARAMETERS
        ----------
        galaxy : `str`
        name of galaxy (i.e. MW,M31)
        snap : `int; 0 or 1`
        simulation timestamp
        '''
        # construct file name 
        ilbl = '000'+str(snap)
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy)+ilbl+'.txt'
        
        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                             
        
        # store the particle masses, positions, and galaxy name
        self.m = self.data['m']
        self.x = self.data['x']
        self.y = self.data['y']
        self.z = self.data['z']
        self.gname = galaxy

    def MassEnclosed(self, ptype, r):
        '''
        Function that calculates the mass enclosed within each
        radius in the r array.
        inputs: ptype= particle types, r= array of radii magnitudes
        outputs: masses= array of masses
        '''
        # assign units to the radius values
        r = r*u.kpc

        # call COM class
        COM = CenterOfMass(self.filename, 2)
        pos_COM = COM.COM_P(0.1)

        # deifing mask to select for particle type
        type = np.where(self.data['type'] == ptype)
        pdata = self.data[type]  # data of specified component
        
        # get positions from data and subtract COM components
        x = pdata['x']*u.kpc-pos_COM[0]
        y = pdata['y']*u.kpc-pos_COM[1]
        z = pdata['z']*u.kpc-pos_COM[2]

        # array to store the mass enclosed by each radius value
        masses = np.zeros(len(r))
        
        for i in range(len(r)):
            # looping through the radius array 
            # mask for selecting particles within the given radius
            radius = np.where(np.sqrt(x**2+y**2+z**2) <= r[i])
            
            pmass = self.data['m'][radius] # array of particle masses within given radius

            # sum enclosed particle masses
            m_sum = np.sum(pmass)
            masses[i] = m_sum # add summed mass to enclosed mass array

        # retrn array of masses as astropy quantities
        return masses*1e10*u.Msun

    def MassEnclosedTotal(self,r):
        '''
        Function that calculates the total mass eclosed by a radius.
        inputs: r = array of radii
        outputs: m_total = array of total mass at each radius value
        '''
        # condition to calculate total mass of M33 with no bulge
        if self.gname == 'M33':
            m_halo = self.MassEnclosed(1,r)
            m_disk = self.MassEnclosed(2,r)
 
            m_total = m_halo + m_disk

        # calculating total mass of MW and M31
        else:
            m_halo = self.MassEnclosed(1,r)
            m_disk = self.MassEnclosed(2,r)
            m_bulge = self.MassEnclosed(3,r)
            
            m_total = m_halo + m_disk + m_bulge

        return m_total
        
    def HernquistMass(self,R,a,mHalo):
        '''
        Function that calculates the Hernquist mass at a 
        given radius.
        inputs: R= radius, a= scale length, mHalo= mass of galaxy halo
        outputs: HernMass= Hernquist mass at R radius
        '''
        # breaking down equation into numerator and denominator
        num = mHalo*R**2
        den = (a+R)**2

        # Hernquist Mass
        HernMass = num/den

        return HernMass*u.Msun

    def CircularVelocity(self,ptype,r):
        '''
        Function that calculates the circular velocity of 
        any given galaxy component at specific radii.
        inputs: ptype= particle types, r= array of radii magnitudes
        outputs: v = array of circular velocities in km/s
        '''
        # define astropy gravitational constant
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # find mass of a given component
        m = self.MassEnclosed(ptype,r)
        r = r *u.kpc # astropy quantity

        # calculating circular velocity
        v = np.sqrt(G*m/r) 

        return v

    def CircularVelocityTotal(self,r):
        '''
        Function that calculates the circular velocity of each galaxy.
        inputs: r= array of radii magnitudes
        outputs: v_total= total circular velocity of all components
        '''
        # define circular velocity of each component
        v_halo = self.CircularVelocity(1,r)
        v_disk = self.CircularVelocity(2,r)
        v_bulge = self.CircularVelocity(3,r)

        # calculate total circular velocity (not just sum)
        v_total = np.sqrt(v_halo**2+v_disk**2+v_bulge**2)

        return v_total

    def HernquistVCirc(self,R,a,mHalo):
        '''
        Function that computes the circular speed using the Hernquist mass profile.
        inputs: R= radius, a= scale length, mHalo= mass of galaxy halo
        outputs: HernV= Hernquist circular velocity 
        '''
        # define astropy gravitational constant
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        # Hernquist mass 
        M = self.HernquistMass(R,a,mHalo)

        R = R*u.kpc # astropy quantity

        # calculate the Hernquist circular velocity
        HernV = np.sqrt(G*M/R)
        
        return HernV
