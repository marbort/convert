# convert
Scripts to convert from CP2K to dpdata and vice-versa

## CP2K_convert

Converts a CP2K files to dpdata format. In line arguments:

`--offset`  only extract every nth frame

`--box` set X, Y, Z box dimensions

`--feval` extract coordinates from a CP2K force evaluation calculation 


If used without arguments will look for each `.xyz` file  with ****pos**** in the name and extracts coordinate and energies from it. 
If an analogous `.xyz` file is found with the same name but ****frc**** instead of ****pos***** extracts forces from that file. 
If forces are not found creates a fake an empty array to be able to import the data as deepmd format. The box size needs to be specified explicitly with the `--box` argument.

If used with the `--feval` argument will look for each CP2K input file (`.cinp` extension) in the folder and extract coordinates and box size. Forces are extracted from the corresponding output which must have the same name as the input file with `.out` extension.

The output is written to _dpdata_ folder.




