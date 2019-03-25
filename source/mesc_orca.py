'''
*----------------------------------*
 Q2DTor and TheRa Program Suites

 Copyright (c) 2018 Universidade de Santiago de Compostela

 This file is part of both Q2DTor and TheRa softwares.

 Q2DTor and TheRa arefree software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Q2DTor and TheRa are distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 inside both Q2DTor and TheRa manuals.  If not, see <http://www.gnu.org/licenses/>.
*----------------------------------*

This is a Module for Electronic Structure Calculations (mesc) using ORCA

The functions here are used to perform different
actions with the ORCA software

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
'''

import os, sys
angstrom = 1.0 / 1.8897261339e+00

#-------------------------------------------------#
# READING PATHS DEFINED BY THE USER               #
#-------------------------------------------------#
OEXEC = None
DIR_CODES = os.path.dirname(os.path.realpath(__file__))+"/"
txtfile   = DIR_CODES+"mesc.txt"
if os.path.isfile(txtfile):
   thefile = open(txtfile,'r')
   lines = thefile.readlines()
   for line in lines:
       line = line.split("#")[0].strip()
       if "mesc_orca " in line and " orca " in line: OEXEC = line.split('"')[1]
   thefile.close()
#-------------------------------------------------#
def check_path():
    abort = False
    if OEXEC is not None:
       if not os.path.exists(OEXEC):
          print "   ERROR: unable to find '%s'"%OEXEC
          abort = True
    else:
       abort = True
       if not os.path.exists(txtfile):
          print "   ERROR: unable to find '%s'"%txtfile
       else:
          print "   ERROR: unable to find 'orca' in '%s'"%txtfile
    if abort:
       sys.exit("   Aborting... Check output for more information.")
#-------------------------------------------------#




#------------------#
# READ BASIC FILES #
#------------------#

def read_output(ofile):

    # Read lines
    output = open(ofile,'r')
    lines  = [line.split('#')[0].strip() for line in output.readlines()]
    output.close()
    
    # Initialize
    ch, mtp, E = None, None, None
    xcc, symbols, atonums = [], [], []

    # Get ch, mtp, E
    pos = None
    for idx in range(len(lines)):
        line = lines[idx]
        if "Total Charge           Charge" in line: ch  = int(line.split()[-1])
        if "Multiplicity           Mult"   in line: mtp = int(line.split()[-1])
        if "FINAL SINGLE POINT ENERGY "    in line: E   = float(line.split()[-1])
        if "CARTESIAN COORDINATES (A.U.)"  in line: pos = idx

    # Get xcc, symbols, atonums
    if pos is None: sys.exit("Unable to find 'CARTESIAN COORDINATES (A.U.)' label in file!")
    pos = pos + 3
    for line in lines[pos:]:
        line = line.strip()
        if line == "": break
        idx, symbol, atonum, dummy, mass, x, y, z = line.split()
        xcc = xcc + [float(x),float(y),float(z)]
        symbols.append(symbol)
        atonums.append(int(atonum.split(".")[0]))

    # Return
    return xcc, symbols, atonums, ch, mtp, E

def read_engrad(engradfile):
    '''
    '''

    # Read lines
    engrad = open(engradfile,'r')
    lines  = [line.split('#')[0].strip() for line in engrad.readlines()]
    engrad.close()

    # Initialize
    natoms   = None
    atonums  = []   ; Etot = None
    x_cc     = []   ; g_cc = []

    # Read natoms, Etot, g_cc, atonums, xcc
    natoms = int(lines[3])
    Etot   = float(lines[7])
    g_cc   = [float(line)  for line in lines[11:11+3*natoms]]
    Zxyz   = [line.split() for line in lines[-natoms:]]
    for Z,x,y,z in Zxyz:
       x_cc = x_cc + [float(x),float(y),float(z)]
       atonums = atonums + [int(Z)]

    # Return data
    return g_cc, (x_cc, atonums, Etot)

def read_hess(hessfile):
    '''
    Input data:
       * hessfile
    Output data:
       * F_cc    : a list containing the lower-triangular part of the hessian matrix (None if not found)
    '''

    # Read lines
    hess = open(hessfile,'r')
    lines  = [line.split('#')[0].strip() for line in hess.readlines()]
    hess.close()

    # Initialize
    F_cc = []

    # Read hessian
    for idx in range(len(lines)):
        line = lines[idx]
        if "$hessian" in line:
           nrows = int(lines[idx+1])
           idx_hessian = idx

    F_cc = [ [0.0]*nrows for i in range(nrows)]

    fcol = int(lines[idx_hessian+2].split()[ 0]) # first column
    lcol = int(lines[idx_hessian+2].split()[-1]) # last  column
    row  = 0

    for line in lines[idx_hessian+3:]:
        if row == nrows-1:
           fcol = int(line.split()[ 0]) # first column
           lcol = int(line.split()[-1]) # last  column
           row  = 0
           continue
        # Split data in line
        line_data = line.split()
        # Save data
        row = int(line_data[0])
        F_cc[row][fcol:lcol+1] = [float(value) for value in line_data[1:]]
        # Finished?
        if row == nrows -1 and lcol == nrows-1: break

    # Promediate hessian
    for i in range(1,nrows,1):
        for j in range(i):
           Fij = F_cc[i][j]
           Fji = F_cc[j][i]
           average = 0.5*(Fij + Fji)
           F_cc[i][j] = average
           F_cc[j][i] = average

    # Get lower triangle
    F_lt = []
    for i in range(nrows):
        for j in range(0,i+1):
            F_lt.append(F_cc[i][j])
    F_cc = F_lt

    # Return data
    return F_cc

#------------------#
# SEND CALCULATION #
#------------------#
def sendcalc(ifile,ofile,err,folder=None):
    '''
    * queue: if True, calculation is sent with nohup
    '''
    check_path()

    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       ifile = folder + ifile
       ofile = folder + ofile
       err   = folder + err

    # Execute command
    command = "%s %s 1>%s 2>%s"%(OEXEC,ifile,ofile,err)
    status  = os.system(command)
    return status


#------------------#
# SOME BASIC NAMES #
#------------------#
def iofiles(name,folder=None):
    '''
    For a given name, it returns the name of all files
    for the Gaussian calculation
    '''

    if   folder is None          : folder = ""
    elif not folder.endswith("/"): folder = folder+"/"
    else                         : folder = folder
    ifile  = folder + name + ".inp"
    ofile  = folder + name + ".out"
    engrad = folder + name + ".engrad"
    hess   = folder + name + ".hess"
    gbw    = folder + name + ".gbw"
    err    = folder + name + ".err"
    return ifile, ofile, engrad, hess, gbw, err

#----------------------------------------#
#   ___   ____   ____  _____             #
#  / _ \ |___ \ |  _ \|_   _|___   _ __  #
# | | | |  __) || | | | | | / _ \ | '__| #
# | |_| | / __/ | |_| | | || (_) || |    #
#  \__\_\|_____||____/  |_| \___/ |_|    #
#                                        #
#----------------------------------------#
q2dtor_key1 = "[Q2DTor_geometry]"
q2dtor_key2 = "[Q2DTor_name]"
q2dtor_key3 = "[Q2DTor_MOread]"
q2dtor_key4 = "[Q2DTor_phi1]"
q2dtor_key5 = "[Q2DTor_phi2]"
#----------------------------------------#

def q2dtor_defaults(level,ch,mtp,torsion1,torsion2,ttype="min"):

    NPROC = 1
    the_string = ""

    #-----------------------------------#
    # Default string for PES point file #
    #-----------------------------------#
    the_string = the_string + "#-------------------------------------#\n"
    the_string = the_string + "# Geometry optimization at the 2D-PES #\n"
    the_string = the_string + "#-------------------------------------#\n"
    the_string = the_string + "start_scangeom orca\n"
    the_string = the_string + "%%pal nprocs %i end\n"%NPROC
    the_string = the_string + "! %s TightSCF\n"%level
    the_string = the_string + "! Opt\n"
    the_string = the_string + q2dtor_key3 + "\n"
    the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
    the_string = the_string + q2dtor_key1 + "\n"
    the_string = the_string + "*\n"
    the_string = the_string + "%geom\n"
    the_string = the_string + "   Constraints\n"
    the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion1)
    the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion2)
    the_string = the_string + "   end\n"
    the_string = the_string + "end\n"
    the_string = the_string + "\n"
    the_string = the_string + "end_scangeom\n"
    the_string = the_string + "#-------------------------------------#\n"
    the_string = the_string + "\n\n\n"

    #----------------------------#
    # Default string for minimum #
    #----------------------------#
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "start_sp0 orca\n"
    the_string = the_string + "#Optimization of a minimum in the 2D-PES\n"
    the_string = the_string + "%%pal nprocs %i end\n"%NPROC
    the_string = the_string + "! %s TightSCF\n"%(level)
    the_string = the_string + "! Opt\n"
    the_string = the_string + "! Freq\n"
    #the_string = the_string + "! NumFreq\n"
    the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
    the_string = the_string + q2dtor_key1 + "\n"
    the_string = the_string + "*\n"
   ##------
   #the_string = the_string + "%geom\n"
   #the_string = the_string + "   Constraints\n"
   #the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion1)+[q2dtor_key4])
   #the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion2)+[q2dtor_key5])
   #the_string = the_string + "   end\n"
   #the_string = the_string + "end\n"
   ##------
    the_string = the_string + "\n"
    the_string = the_string + "end_sp0\n"
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "\n\n\n"


    #----------------------------#
    # Default string for saddle  #
    #----------------------------#
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "start_sp1 orca\n"
    the_string = the_string + "#Optimization of a saddle point in the 2D-PES\n"
    the_string = the_string + "%%pal nprocs %i end\n"%NPROC
    the_string = the_string + "! %s TightSCF\n"%(level)
    the_string = the_string + "! OptTs\n"
    the_string = the_string + "! Freq\n"
   #the_string = the_string + "! NumFreq\n"
    the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
    the_string = the_string + q2dtor_key1 + "\n"
    the_string = the_string + "*\n"
    the_string = the_string + "%geom\n"
   #the_string = the_string + "%geom GDIISMaxEq 20\n"
   #the_string = the_string + "  InHess      Read       # Initial Hessian will be read\n"
   #the_string = the_string + '  InHessName  "%s.hess"  # File with Hessian matrix\n'%q2dtor_key2
    the_string = the_string + "  Calc_Hess   true       # Calculate Hessian in the beginning\n"
   #the_string = the_string + "  NumHess     true       # Request numerical Hessian (analytical not available)\n"
   #the_string = the_string + "  Recalc_Hess 5          # Recalculate the Hessian every 5 steps\n"
   #the_string = the_string + "  UseGDIIS    false      #\n"
    the_string = the_string + "end\n"
    the_string = the_string + "\n"
    the_string = the_string + "end_sp1\n"
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "\n\n\n"

    #----------------------------#
    # Default string for maximum #
    #----------------------------#
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "start_sp2 orca\n"
    the_string = the_string + "#Optimization of a maximum in the 2D-PES\n"
    the_string = the_string + "%%pal nprocs %i end\n"%NPROC
    the_string = the_string + "! %s TightSCF\n"%(level)
    the_string = the_string + "! Opt\n"
    the_string = the_string + "! Freq\n"
    #the_string = the_string + "! NumFreq\n"
    the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
    the_string = the_string + q2dtor_key1 + "\n"
    the_string = the_string + "*\n"
    #------
    the_string = the_string + "%geom\n"
    the_string = the_string + "   Constraints\n"
    the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion1)+[q2dtor_key4])
    the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion2)+[q2dtor_key5])
    the_string = the_string + "   end\n"
    the_string = the_string + "end\n"
    #------
    the_string = the_string + "\n"
    the_string = the_string + "end_sp2\n"
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "\n\n\n"

    return the_string

#def q2dtor_defaults(level,ch,mtp,torsion1,torsion2,ttype="min"):
#
#    NPROC = 1
#    the_string = ""
#
#    #-----------------------------------#
#    # Default string for PES point file #
#    #-----------------------------------#
#    the_string = the_string + "#-------------------------------------#\n"
#    the_string = the_string + "# Geometry optimization at the 2D-PES #\n"
#    the_string = the_string + "#-------------------------------------#\n"
#    the_string = the_string + "start_scangeom orca\n"
#    if ttype == "ts":
#       the_string = the_string + "# Firstly, hessian calculation\n"
#       the_string = the_string + "%%pal nprocs %i end\n"%NPROC
#       the_string = the_string + "! %s\n"%level
#       the_string = the_string + "! NumFreq\n"
#       the_string = the_string + q2dtor_key3 + "\n"
#       the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
#       the_string = the_string + q2dtor_key1 + "\n"
#       the_string = the_string + "*\n"
#       the_string = the_string + "# Now, optimization\n"
#       the_string = the_string + "$new_job\n"
#       the_string = the_string + "%%pal nprocs %i end\n"%NPROC
#       the_string = the_string + "! %s TightSCF\n"%level
#       the_string = the_string + "! OptTS\n"
#       the_string = the_string + q2dtor_key3 + "\n"
#       the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
#       the_string = the_string + q2dtor_key1 + "\n"
#       the_string = the_string + "*\n"
#       the_string = the_string + "%geom GDIISMaxEq 20\n"
#       the_string = the_string + "   UseGDIIS false\n"
#       the_string = the_string + "   InHess Read\n"
#       the_string = the_string + '   InHessName "%s.hess"\n'%q2dtor_key2
#       the_string = the_string + "   Constraints\n"
#       the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion1)
#       the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion2)
#       the_string = the_string + "   end\n"
#       the_string = the_string + "end\n"
#       the_string = the_string + "\n"
#    if ttype == "min":
#       the_string = the_string + "%%pal nprocs %i end\n"%NPROC
#       the_string = the_string + "! %s TightSCF\n"%level
#       the_string = the_string + "! Opt\n"
#       the_string = the_string + q2dtor_key3 + "\n"
#       the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
#       the_string = the_string + q2dtor_key1 + "\n"
#       the_string = the_string + "*\n"
#       the_string = the_string + "%geom\n"
#       the_string = the_string + "   Constraints\n"
#       the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion1)
#       the_string = the_string + "     {D %2i %2i %2i %2i C }\n"%tuple(torsion2)
#       the_string = the_string + "   end\n"
#       the_string = the_string + "end\n"
#       the_string = the_string + "\n"
#    the_string = the_string + "end_scangeom\n"
#    the_string = the_string + "#-------------------------------------#\n"
#
#    the_string = the_string + "\n\n\n"
#
#    #----------------------------#
#    # Default string for CP file #
#    #----------------------------#
#    for i in [0,1,2]:
#        the_string = the_string + "#-----------------------------#\n"
#        if i == 0: the_string = the_string + "start_sp0 orca\n#Optimization of a minimum in the 2D-PES\n"
#        if i == 1: the_string = the_string + "start_sp1 orca\n#Optimization of a saddle point in the 2D-PES\n"
#        if i == 2: the_string = the_string + "start_sp2 orca\n#Optimization of a maximum in the 2D-PES\n"
#
#        # If TS, first hessian calculation
#        if i == 1:
#          the_string = the_string + "# Firstly, hessian calculation\n"
#          the_string = the_string + "%%pal nprocs %i end\n"%NPROC
#          the_string = the_string + "! %s\n"%level
#          the_string = the_string + "! NumFreq\n"
#          the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
#          the_string = the_string + q2dtor_key1 + "\n"
#          the_string = the_string + "*\n"
#          the_string = the_string + "# Now, optimization\n"
#          the_string = the_string + "$new_job\n"
#
#        # Optimization + Frequency calculation
#        the_string = the_string + "%%pal nprocs %i end\n"%NPROC
#        the_string = the_string + "! %s TightSCF\n"%(level)
#        if i == 0 and ttype == "min": the_string = the_string + "! Opt\n"
#        if i == 1 and ttype == "min": the_string = the_string + "! OptTS\n"
#        if i == 2 and ttype == "min": the_string = the_string + "! Opt\n"
#        if i == 0 and ttype == "ts" : the_string = the_string + "! OptTS\n"
#        if i == 1 and ttype == "ts" : the_string = the_string + "! Opt\n"
#        if i == 2 and ttype == "ts" : the_string = the_string + "! Opt\n"
#        the_string = the_string + "! NumFreq\n"
#        the_string = the_string + "* xyz %i %i\n"%(ch,mtp)
#        the_string = the_string + q2dtor_key1 + "\n"
#        the_string = the_string + "*\n"
#
#        if (ttype == "min" and i == 1):
#          the_string = the_string + "%geom GDIISMaxEq 20\n"
#          the_string = the_string + "   UseGDIIS false\n"
#          the_string = the_string + "   InHess Read\n"
#          the_string = the_string + '   InHessName "%s.hess"\n'%q2dtor_key2
#          the_string = the_string + "end\n"
#
#        if (ttype == "min" and i == 2) or (ttype == "ts" and i > 0):
#           the_string = the_string + "%geom\n"
#           the_string = the_string + "   Constraints\n"
#           the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion1)+[q2dtor_key4])
#           the_string = the_string + "     {D %2i %2i %2i %2i %s C }\n"%tuple(list(torsion2)+[q2dtor_key5])
#           the_string = the_string + "   end\n"
#           the_string = the_string + "end\n"
#
#        the_string = the_string + "\n"
#        if i == 0: the_string = the_string + "end_sp0\n"
#        if i == 1: the_string = the_string + "end_sp1\n"
#        if i == 2: the_string = the_string + "end_sp2\n"
#        the_string = the_string + "#-----------------------------#\n"
#        the_string = the_string + "\n\n\n"
#
#    return the_string

def q2dtor_pespoint(ilines, xvec,symbols,name,previous=None):
    '''
    xcc = [(x1,y1,z1),(x2,y2,z2),...]
    mos  = main name for reference Molecular Orbitals
           if mos is "system", then "scangeom_system.gbw" will be used as
           MO guess; file mos.gbw is in folder
    '''


    #------------------------------#
    # Calculation folder and files #
    #------------------------------#
    mainname = "scangeom_"+name
    folder   = "TMP_CALCS/"
    if not os.path.exists(folder): os.mkdir(folder)
    ifile, ofile, engrad, hess, gbw, err = iofiles(mainname,folder)

    #----------#
    # MO guess #
    #----------#
    string_guess = ""
    if previous is not None:
       ref_ifile, ref_ofile, ref_engrad, ref_hess, ref_gbw, ref_err = iofiles("scangeom_"+previous,folder)
       if os.path.exists(ref_gbw): string_guess = '! MORead\n%%moinp "%s"\n'%ref_gbw

    #--------------#
    # Input string #
    #--------------#
    string_ifile = ""
    for line in ilines:
        if q2dtor_key1 in line:
           line = ""
           for idx in range(len(symbols)):
               symbol = symbols[idx]
               x,y,z  = [value*angstrom for value in xvec[3*idx:3*idx+3]]
               line   = line + "%2s   %+11.6f  %+11.6f  %+11.6f\n"%(symbol,x,y,z)
        if q2dtor_key2 in line:
           pos  = line.find(q2dtor_key2)
           line = line[0:pos] + "%s"%(folder+mainname) + line[pos+len(q2dtor_key2):]
        if q2dtor_key3 in line:
           line = string_guess
        string_ifile = string_ifile + line

    #-----------------------#
    # Send Orca calculation #
    #-----------------------#
    ff = open(ifile,"w")
    ff.write(string_ifile)
    ff.close()
    status = sendcalc(ifile,ofile,err)

    # Check ofile
    ff = open(ofile,'r')
    olines = ff.readlines()
    ff.close()
    if "ORCA TERMINATED NORMALLY" not in olines[-2]: return None

   ##--------------#
   ## Rename files #
   ##--------------#
   #j2_ifile, j2_ofile, j2_engrad, j2_hess, j2_gbw, j2_err = iofiles(mainname+"_job2",folder)
   #if os.path.exists(j2_engrad): os.rename(j2_engrad, engrad)
   #if os.path.exists(j2_gbw   ): os.rename(j2_gbw   , gbw   )

    #-----------#
    # Read data #
    #-----------#
    x_cc, symbols, atonums, ch, mtp, Etot = read_output(ofile)
    g_cc = None
    F_cc = None

    #-------------#
    # Remove data #
    #-------------#
    files = os.listdir(folder)
    files = [fff for fff in files if fff.startswith(mainname)]
    files = [fff for fff in files if not fff.endswith(".gbw")]
    files = [fff for fff in files if not fff.endswith(".inp")]
    files = [fff for fff in files if not fff.endswith(".out")]
   #for fff in files: os.remove(folder+fff)

    return x_cc, atonums, ch, mtp,  Etot, g_cc, F_cc


def q2dtor_optCP(ilines,xvec,symbols,mainname,phis=(None,None)):
    '''
    xcc = [(x1,y1,z1),(x2,y2,z2),...]
    phis in degrees
    '''

    #--------------------#
    # Calculation folder #
    #--------------------#
    folder = "TMP_CALCS/"
    if not os.path.exists(folder): os.mkdir(folder)


    #-------#
    # Files #
    #-------#
    phi1 = phis[0]
    phi2 = phis[1]

    ifile, ofile, engrad, hess, gbw, err = iofiles(mainname,folder)

    # The input string
    string_ifile = ""
    for line in ilines:
        if q2dtor_key1 in line:
           line = ""
           for idx in range(len(symbols)):
               symbol = symbols[idx]
               x,y,z  = [value*angstrom for value in xvec[3*idx:3*idx+3]]
               line   = line + "%2s   %+11.6f  %+11.6f  %+11.6f\n"%(symbol,x,y,z)
        if q2dtor_key2 in line:
           pos  = line.find(q2dtor_key2)
           line = line[0:pos] + "%s"%(folder+mainname) + line[pos+len(q2dtor_key2):]
        if q2dtor_key4 in line:
           pos  = line.find(q2dtor_key4)
           line = line[0:pos] + "%.2f"%phi1 + line[pos+len(q2dtor_key4):]
        if q2dtor_key5 in line:
           pos  = line.find(q2dtor_key5)
           line = line[0:pos] + "%.2f"%phi2 + line[pos+len(q2dtor_key5):]
        string_ifile = string_ifile + line

    #-----------------------#
    # Send Orca calculation #
    #-----------------------#
    ff = open(ifile,"w")
    ff.write(string_ifile)
    ff.close()
    print "        * Geometry optimization + Freq calculation with Orca..."
    status = sendcalc(ifile,ofile,err)

    # Check ofile
    ff = open(ofile,'r')
    olines = ff.readlines()
    ff.close()
    if "ORCA TERMINATED NORMALLY" not in olines[-2]: return None


   ##--------------#
   ## Rename files #
   ##--------------#
   #j2_ifile, j2_ofile, j2_engrad, j2_hess, j2_gbw, j2_err = iofiles(mainname+"_job2",folder)
   #if os.path.exists(j2_engrad): os.rename(j2_engrad, engrad)
   #if os.path.exists(j2_hess  ): os.rename(j2_hess  , hess  )
   #if os.path.exists(j2_gbw   ): os.rename(j2_gbw   , gbw   )

    #-----------#
    # Read data #
    #-----------#
    x_cc, symbols, atonums, ch, mtp, Etot = read_output(ofile)
    g_cc, dummy = read_engrad(engrad)
    F_cc = read_hess(hess)

    # Clean
   #print "        * Remove calculation files (inp, out, engrad, hess)..."
    files = os.listdir(folder)
    files = [fff for fff in files if fff.startswith(mainname)]
    files = [fff for fff in files if not fff.endswith(".inp")]
    files = [fff for fff in files if not fff.endswith(".out")]
   #for fff in files: os.remove(folder+fff)

    return x_cc, atonums, ch, mtp, Etot, g_cc, F_cc

#----------------------------------------#
#    _____  _            ____            #
#   |_   _|| |__    ___ |  _ \  __ _     #
#     | |  | '_ \  / _ \| |_) |/ _` |    #
#     | |  | | | ||  __/|  _ <| (_| |    #
#     |_|  |_| |_| \___||_| \_\\__,_|    #
#                                        #
#----------------------------------------#
thera_key1 = "[TheRa_geometry]"
thera_key2 = "[TheRa_name]"
thera_key3 = "[TheRa_gradhess]"
#----------------------------------------#

def thera_default(ch=0,mtp=1,nproc=4,mem=4):

    string = ""
    string = string + "%%pal nprocs %i end\n"%nproc
    string = string + "! hf sto-3g TightSCF\n"
    string = string + "%s\n"%thera_key3
    string = string + "* xyz %i %i\n"%(ch,mtp)
    string = string + "%s\n"%thera_key1
    string = string + "*\n"
    string = string + "\n"
    return string

def thera_userfolder(folder):
    return []
#
#   # List of fchk files
#   file_list = [filename[:-5] for filename in os.listdir(folder) if filename.lower().endswith("fchk")]
#   # Get data from fchk files
#   data = []
#   for mainname in file_list:
#       # read file
#       x_cc, atonums, ch, mtp,  E, g_cc, F_cc = readfchk(mainname,folder=folder)
#       # append data
#       data.append( [mainname, x_cc, atonums, ch, mtp, E, g_cc, F_cc] )
#   # Sort by energy
#   data.sort(key=lambda xx: xx[5])
#   # Return data
#   return data

def thera_spc(pointname,xvec,symbols,ref_lines,folder=None,hessian=False):
    pass



