# import necessary modules/packages
import numpy as np
import astropy.units as u

def Read(filename):
    '''
    this function reads a file of values seperated by spaces, 
    finds values in first and second lines; returns the
    first and second values, and the rest of the data 
    '''
    file = open(filename, 'r') # open file inputted into the function
    line1 = file.readline() # reads first line
    label, value = line1.split() # splits values along space
    time = float(value)*u.Myr # convert the time value found to Myr

    line2 = file.readline() # reads second line
    label2, value2 = line2.split() # splits values along space
    
    file.close() # close the file

    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) # identify the data, excluding the header

    return time, value2, data