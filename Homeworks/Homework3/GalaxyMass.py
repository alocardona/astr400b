import numpy as np
import astropy.units as u
from ReadFile import Read

def ComponentMass(filename,type):
    '''
    this function returns the total mass of any given galaxy component
    inputs: filename = name of data file, type = particle type 
    designating a particular galaxy component
    '''
    # call Read function to get particle data 
    Tstart,particles,data = Read(filename)

    # mask to select for particles of specified type
    index = np.where(data['type'] == type) 
    mass = data['m'][index] # array of particle masses
    
    # loop through particle masses of a specific component and summing them
    comp_mass = 0 # initialize mass sum
    for val in mass:
        comp_mass += val

    galaxy = filename[0:3] # find name of galaxy for printing results
    M = np.round(comp_mass*0.01,3) # round and multiply by 0.01 to get units of 10^12 solmasses
    
    # print results based on galaxy component
    if type == 1:
        print(galaxy,'dark matter halo mass:',M,'x 10^12 solmasses')
    elif type == 2:
        print(galaxy,'disk mass:',M,'x 10^12 solmasses')
    elif type == 3:
        print(galaxy,'bulge mass:',M,'x 10^12 solmasses')
        
    return M

# call function to find masses of all galaxy components
# Milky Way
MW_bulge = ComponentMass('MW_000.txt',3) # bulge
MW_disk = ComponentMass('MW_000.txt',2) # disk 
MW_halo = ComponentMass('MW_000.txt',1) # dark matter 

# Andromeda (M31)
M31_bulge = ComponentMass('M31_000.txt',3) # bulge
M31_disk = ComponentMass('M31_000.txt',2) # disk
M31_halo = ComponentMass('M31_000.txt',1) # dark matter 

# Triangulum (M33)
M33_bulge = ComponentMass('M33_000.txt',3) # bulge 
M33_disk = ComponentMass('M33_000.txt',2) # disk 
M33_halo = ComponentMass('M33_000.txt',1) #dark matter




        
            


    
    
    