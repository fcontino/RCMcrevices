# RCMcrevices
Script to evaluate if the size of the RCM crevices are well designed to absorb the boundary layer vortex.


## Launch the script: 

`python mStar.py csvFile.csv inputs.txt`

## Inputs to the script:

1. csv file (without header) and with the following columns
* time (s), 
* position of the piston (position(t0) = 0, position(end) = stroke) (m),
* pressure (Pa)
2. text file specifying the following (in this order)
* species name of the mixture (Name, in IUPAC form or common form or a synonym registered in PubChem, or CAS number)
* mass fractions for each species in the mixture
* time of beginning of compression (s), if None is specified, then it is detected automatically
* time of end of compression, if None is specied, this is by default the time of maximum pressure
* volume (m3) of the main chamber at TDC (without crevices but including other dead volumes)
* volume of the crevices (not including other dead voumes) (m3)
* radius of the combustion chamber (m)
* temperature at the beginning of compression (K)
* Stroke (m), if the stroke is obtained from the position vector, then indicate 0


## Output: 

The ratio between the volume of the crevices and the volume needed to absorb the boundary layer. The ratio should be bigger than 1 for the crevices to be effective.


## More information:

Details of the implementation are provided in the following article:
N. Bourgeois, H. Jeanmart, G. Winckelmans, O. Lamberts, F. Contino

[How to ensure the interpretability of experimental data in Rapid Compression Machines? A method to validate piston crevice designs
Combustion and Flame, 2018](https://www.sciencedirect.com/science/article/abs/pii/S0010218018304279)

## Contributors to this script: 

N. Bourgeois, M. Pochet, F. Contino