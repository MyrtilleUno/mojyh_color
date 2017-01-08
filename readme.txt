Welcome to ** Mojyh_Color ** Colorimetry code
*********************************************
This software is under the GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007

This software was originaly created under the name "Yoshi_color" for 
the PhD thesis of M.O.J.Y. Hunault.
https://tel.archives-ouvertes.fr/tel-01137588/document

If you use this code, please do not forget to acknowledge the authors:
M.O.J.Y. Hunault, currently at Utrecht University, the Netherlands
V. Magron, currently at the Verimag Laboratory, Grenoble, France

Please use the Git web site to submit questions or comments.


General Purpose
***************
This program written in C and designed to calculate CIE-Lab values from an optical absorption spectrum.

This program calculates L*a*b* CIE coordinates for the D65 illuminant and CIE 1931 2˚ observer.

 
Install
*******
This program was designed for Unix platforms.
Use GIT to clone the repository or download the folder as it is.
Compilation: 
Open a Terminal window and go the directory mojyh_color.
Copy the "makefile.*" file that corresponds to the OS you are using, to a "makefile" file.
Then compile the code using the make file by typing "make".

Add the mojyh_color directory to your PATH in order to be able to run the code form anyother directory.

Running the software
********************
The Code run for two different kinds of input files :
1*/ Files formatted in columns :
-First column must be the wavelength from the spectrum and cover from 380 nm to 780 nm by steps of 5 nm.
-Next columns must be data in transmission from 0 to 1 (not in %).


To run the code, type :
mojyh_color filename spectro

2*/ Files output from Quanty software.
QUANTY is a software developped by Maurits Haverkort: https://www.thphys.uni-heidelberg.de/~haverkort/quanty/ 
This option was designed to deal directly with this type of files that contains spectral data in the imaginary part of the
Green's function calculated by QUANTY. 

Run :
mojyh_color filename quanty



License
*******




References
**********
References for the color matching functions used to analyse the spectrum
http://cvrl.ioo.ucl.ac.uk/database/text/cmfs/sbrgb2.htm
http://cvrl.ioo.ucl.ac.uk/database/text/cmfs/sbrgb10.htm

