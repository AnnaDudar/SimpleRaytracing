# SimpleRaytracing
SimpleRaytracing is a computer graphics program aiming to render a 3D scene on the screen.

## Description
The algorithm behind is raytracing. It is based on the algorithm provided by [Scratchapixel.com](https://www.scratchapixel.com/).\
Capabilities of the algorithm:
* rendering of several light sources;
* different material types, including diffuse, glossy, specular and transparent materials;
* setting camera options;
* reading object from a .geo file;
* writing output in a .ppm file.

The scene is at the moment predefined. It can be set either directly in source files or read from .geo files.\
For demonstration reasons the scene was set in the source files and one object is read from the .geo file located in the folder *scenes*.\
The rendered output is stored in the *outputSpheres.ppm* file (.ppm files can be read by a free image editor [GIMP](https://www.gimp.org/)).\
<img src="https://user-images.githubusercontent.com/76623102/110004203-73799600-7d17-11eb-92e8-ab1e782224f3.png" alt="" width="800"/>

## Installation
The files represented on GitHub are already archived source files. The project is intended for the CodeBlocks IDE. No furher development is supposed.\
If you would like to run the program and get the rendered output yourself, please download the *Raytracing_Release.zip* file.
The .zip file contains the following items:
* *log* folder for saving logs (one log is saved already as an example);
* *scenes* folder for keeping objects building up the scene;
* executable *Raytracing.exe*;
* output file *outputSpheres.ppm*.

## Contact
If you have any questions please do not hesitate to contact me via anna_dudar@hotmail.com
