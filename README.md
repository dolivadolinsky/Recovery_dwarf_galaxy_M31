# Recovery_dwarf_galaxy_M31

This directory contains a script which from a list of properties of dwarf galaxies return their detection\non-detection and the probability that such dwarf galaxies would have been discovered in PAndAS. It using the recovery fraction model described in Doliva-Dolinsky and al, 2021.

# Content: 
The script is named "detection_dwarf_galaxy_M31.py". It contains a function "detection_dwarf_galaxy" which takes in input a list of right ascension (rad), a list of declination (rad), a list of distance (kpc), a list of absolute magnitude in the V-band and a list of half-light radius (pc). It returns a list of 1, 0 which correspond, respectively, to the galaxy is detected or to the galaxy is not detected. The second output is an array of the recovery rates. It also returns an array which contains an information on the position of the galaxies. Its value is 0 if the galaxy is not in the survey, 1 if it is in a field of the survey where the recovery fraction were calculated, 2 if the galaxy is in a masked field and 3 if the galaxy is in a hole of the survey.

The file "parameters2.txt" contains the number of the PAndAS field, the center of the field in the tangent plane rectangular coordinates (deg) and the corresponding recovery fractions model parameters if the dwarf galaxies are at the distance of M31. 

The file "field_corners.dat" and "field_middle_xkieta.dat" contain, respectively the coordinates of the corner of each PAndAS field in the tangent plane rectangular coordinates (deg) and the coordinates of the center of each PAndAS field in the tangent plane rectangular coordinates (deg). 

The file "example.py" presents a small example on how to use the module. 

# Principle:
The module first determines the fields in which the input dwarf galaxy would be located or if it is outside the PAndAS footprint. The parameters of the model described in Doliva-Dolinsky et al (2021), Mv_lim, α and σ are obtained for the field in which the galaxy is and for a distance of the galaxy equal to the one of M31.  Then, the relations described in Doliva-Dolinsky et al (2021), are used in order to determine the values of Mv_lim and α at the galaxy’s distance. Once the detection efficiencies are obtained, an acceptance/rejection method is used to decide if the galaxy is detected or not.

# Requirement:
The module was tested with Python 3.8.3 and the associate math and random module. It is also using numpy 1.19.5 and scipy 1.5.0. 

# Implementation:
To use the module, download the files, add import path_to_the_script/detection_dwarf_galaxy_M31 and call the function as detection_dwarf_galaxy_M31.detection_dwarf_galaxy(ra,dec,distance,Mv,rh). An example on how to use it is shown in the file example.py.
