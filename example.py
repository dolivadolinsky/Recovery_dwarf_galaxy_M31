import numpy as np 
import random
import detection_dwarf_galaxy_M31 as detec #Importation of the module as detec

ra_refc,dec_refc=0.18648115837037746,0.7202828378876265 #right ascension (rad) and declination(rad) of M31

a=ra_refc+np.random.rand(10)/10 #Array of random right ascension (rad) around M31 
b=dec_refc+np.random.rand(10)/10 #Array of random declination (rad) around M31 
distance=600+np.random.rand(10)*300 #Array of random distance (kpc)
mv=-4.5-np.random.rand(10)*4 #Array of random magnitude 
rh=np.power(10,1.8+np.random.rand(10)*1.2) #Array of random half-light radius (pc)

print(detec.detection_dwarf_galaxy(a,b,distance,mv,rh))




	
