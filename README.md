# Q2DTor


## Licensing and Distribution 

Q2DTor version 1.1

Copyright (c) 2019, Antonio Fernandez Ramos and David Ferro Costas 

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software 
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

## About Q2DTor

    Name of the Program: Q2DTor

    Program Version : 1.1

    Program Version Date: March 20, 2019

    Manual Version Date: March 20, 2018 

Q2DTor is a program designed to calculate partition functions and thermodynamic 
properties of molecular systems with two or more torsional modes. 
It allows calculating rotational-vibrational (or rovibrational) 
partition functions and thermodynamic functions by the multistructural harmonic 
oscillator (MS-HO) and the Extended Two-Dimensional Torsion (E2DT) methods. 
If the molecule has more than two torsions, the program treats only two of them as 
torsions, with the others harmonic. 
Q2DTor has to be executed in several steps, each one of them performing a different task. 
With this 'step-by-step' procedure the user can check the output file after a 
given task (and before proceeding to the next one).


## Version information

### version 1.1

List of known bugs and changes/additions in version 1.1:
- Important:
  - Corrected a bug in the function that reads the xyz file (problem with empty lines).
  - Corrected a bug related to the calculation of partition functions
    in systems with the 'a' symmetry condition.
  - Corrected a bug in the calculation of S and G for MSHO.
- Minor:
  - The version of Q2DTor is printed after the logo of the program.
  - The date is added at the beginning of each execution.
  - The Q2DTor.py executable and its modules are now separated.
    Modules are stored in the modsQ2DTor/ folder.
  - If tsigma1 * tsigma2 is not equal to 1, it is assumed that the system meets the
    'c'  type condition of symmetry even if that it is not specified in the input file.
  - For systems with 'b' type symmetry, odd symmetry terms of the Fourier series are
    automatically removed, even if they are specified by the user in the input file.
  - If tsigma_i is not 1, the corresponding Fourier terms of the phi_i dihedral
    are multiplies of tsigma_i. Other terms are automatically removed even if they are
    specified in the input file.
  - Now, by default, the scan of the PES starts at the closest point to the input
    geometry instead of at (phi1,phi2)=(0,0).
  - Option --units is now available.
    Energies in kJ/mol are printed with: --units joule
  - A new table for each thermodynamic function is added.


## How to cite

D. Ferro-Costas, M. N. D. S. Cordeiro, D. G. Truhlar, A. Fernández-Ramos, Comput.  Phys. Commun. 232, 190-205, 2018.


## Description of files

 Contents of the folders distributed in this version:
  - **src/**       : Q2DTor source files
  - **manual/**    : Manual of Q2DTor
  - **tests/**     : All the files related to the tests set

The Q2DTor source folder consists of:

   - **Q2DTor.py**:
     This is the main file of Q2DTor, i.e., the one to be executed

   - **mesc.txt**:
     It contains the path to the executables to the Electronic Structure Calculation software.
     THIS IS THE ONLY FILE THAT HAS TO BE MODIFIED BY THE USER.

   - Directory **mq2dtor/**, which contains the following files: 
     - **constants.py**:
        It contains different physical constants as well as atomic masses and covalent radii.
     - **classes.py**:
       It contains different Python classes that are used along the whole code. 
       For example, it contains different classes related to graph theory. 
       These allow to obtain the connectivity of a given geometry or a set of
       internal coordinates in just a few calls.
     - **gtsfile.py**:
       It contains functions that read and write a gts file (the format used by Q2Dtor). 
       These gts files contain information about stationary points: geometry, gradient, 
       Hessian matrix, etc.
     - **helpfns.py**:
       It contains different auxiliary functions (for example for example the function
       that converts Cartesian coordinates to mass-scaled Cartesian).
     - **mesc_gaussian.py**:
       It contains the interface to the Gaussian Package.  
     - **mesc_orca.py**:
       It contains the interface to the Orca Package.
     - **tesselation.py**:
       It contains functions related to the tesselation of the 2D potential energy surface.
     - **quotes.py**:
       It contains some quotes to print when the program ends.

The tests folder contains the output files of the tests set and an script to run them. 
The documents folder contains this manual.


## Installation

Q2Dtor is a program written in Python 2. Consequently, it does not need any kind 
of compilation, as it would be the case with C or Fortran programs.
The user should install Python 2 in order to use Q2DTor (version 2.7 is recommended), 
as well as the following Python libraries:
   - matplotlib
   - numpy
   - pylab
   - scipy

WARNING: do not use Python 3 to execute Q2DTor.


-----------------------
 Setting up the program
-----------------------

Before using Q2DTor, the user has to define the path to the executable(s) of the 
software for the electronic structure calculation.
This is done in file 'mesc.txt'.

For Gaussian users, this file has to contain these lines:

   mesc_gaussian  gauexe  "PATH_TO_GAUSSIANEXE"
   mesc_gaussian  fchk    "PATH_TO_FCHK"

where PATH_TO_GAUSSIANEXE and PATH_TO_FCHK are the paths to the Gaussian and formchk 
executables.
Example:

   mesc_gaussian  gauexe  "/home/programs/G09_64D/g09/g09"
   mesc_gaussian  fchk    "/home/programs/G09_64D/g09/formchk"

Notice that the path is between quotation marks(").

Similarly, for Orca users, this file has to contain the line:

   mesc_orca      orca  "PATH_TO_ORCAEXE"

where, again, the path to the Orca executable is between quotation marks. 
Example:

   mesc_orca      orca  "/home/programs/orca_4_0_1_2/orca"


## Execution

You can run Q2DTor by invoking the Python interpreter manually as follows:

   ```python2.7 Q2DTor.py```

If you prefer to avoid invoking the Python interpreter, you have to follow these
two simple steps:

(1) Add as first line in the Q2DTor.py the following:

      #!PATH_FOR_PYTHON python 

   where PATH_FOR_PYTHON indicates the location of the Python interpreter.

   Example:

     ```#!/usr/bin/env python```

   In this example Python is located in /usr/bin/env. 

(2) Make the main program Q2DTor.py executable:

      ```chmod u+x Q2DTor.py```

    This allows you to run Q2DTor just using:

      ```Q2DTor.py```

Before run Q2DTor, we recommend to read the help menu.
It can be displayed either by typing

   ```Q2DTor.py --help```

or 

   ```Q2DTor.py -h```

## Running the tests set

In the tests directory there is a script called Q2DTor_tests.py 
that facilitates running the tests set. 
                                                            
The script runs under Python version 2.x. It is executed by typing:

  ```python Q2DTor_tests.py```

This action opens the next interactive menu:                  

     | << Test creator for Q2DTor >>                        
     |                                                      
     |    (1) Create input files                            
     |    (2) Check results                                 
     |                                                      
     |    your choice:                                      
     |                                                      

where you should choose option 1. 

At this point, the script generates two folders:                                      
- GAUSSIAN                                              
- ORCA                                                  
and each of them contains 20 test folders (one per system). 
                                                            
For each system, you can find the xyz file with the         
reference geometry and the corresponding Q2DTor input file. 
                                                            
After executing Q2DTor.py for a given system (or for all of    
them), you can check your results by using the option 2     
of the Q2DTor_tests.py script. This checking process can    
be performed after every single step of Q2DTor or after     
executing it with all the options.                          
                                                            
IMPORTANT NOTICE: The searching algorithm for systems  
S10 and S19 is not able to find one transition state, and     
Q2DTor_tests.py complains if you execute it once      
Q2DTor ends the tasks associated to the --optsp argument.   
In these cases, the user should add the following line  
at the end of the SXX/IOfiles/SXX.splist file (XX=10,19). 

For XX=10:

   1 170.00 60.00 - NO S10_170_060

For XX=19

   1 130.00 70.00 - NO S19_130_070

Check the manual for more information about the search and
optimization procedures of the stationary points.                              
                                                            

