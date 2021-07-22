# Recovery_dwarf_galaxy_M31

This directory contains a module which from a list of galaxies and a list of their properties return their detection\non-detection using the recovery fraction model described in Doliva-Dolinsky and al, 2021.

# Content: 
The module is named "detection_dwarf_galaxy_M31.py". It contains a function "detection_dwarf_galaxy" which takes in input a list of right ascension (rad), a list of declination (rad), a list of distance (kpc), a list of absolute magnitude in the V-band and a list of half-light radius (pc). It returns a list of 1, 0 which correspond, respectively, to the galaxy is detected or to the galaxy is not detected. It also returns an array which contains an information on the position of the galaxies. Its value is 0 if the galaxy is not in the survey, 1 if it is in a field of the survey where the recovery fraction were calculated, 2 if the galaxy in a masked field and 3 of the galaxy is in a hole of the survey.

The file "parameters2.txt" contains the number of the PAndAS field, the center of the field in the tangent plane rectangular coordinates (deg) and the corresponding recovery fractions model parameters if the dwarf galaxies are at the distance of M31. 

The file "field_corners.dat" and "field_middle_xkieta.dat" contain, respectively the coordinates of the corner of each PAndAS field in the tangent plane rectangular coordinates (deg) and the coordinates of the center of each PAndAS field in the tangent plane rectangular coordinates (deg). 

The file "example.py" presents a small example on how to use the module. 

# Principle:
The module works by first associating the center of the galaxy with its PAndAS field and its recovery fractions model parameters for a galaxy at a distance of M31. Then, the recovery fractions model parameters for a galaxy at the distance given in input is obtained thanks to the straight line model described in Doliva-Dolinsky and al, 2021. The efficiency of detection is then determined. Finally, thanks to an acceptance/rejection method the detection/non-detection of the galaxy is obtained. 

# Requirement:
The module was tested with Python 3.8.3 and the associate math and random module. It is also using numpy 1.19.5 and scipy 1.5.0. 

# Implementation:
To use the module, add import path_to_the_module/detection_dwarf_galaxy_M31 and call the function as detection_dwarf_galaxy_M31.detection_dwarf_galaxy(ra,dec,distance,Mv,rh). 
