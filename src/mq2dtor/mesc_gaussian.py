'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Q2DTor
Version     : 1.1
License     : MIT/x11

Copyright (c) 2019, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

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
---------------------------


*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*

 This is a Module for Electronic Structure Calculations (mesc)
 using GAUSSIAN

 The functions here are used to perform different
 actions with the GAUSSIAN software

'''
import os
import sys
import shutil
import time

angstrom = 1.0 / 1.8897261339e+00

#-------------------------------------------------#
# READING PATHS DEFINED BY THE USER               #
#-------------------------------------------------#
GFCHK,GEXEC = None, None
DIR_CODES = os.path.dirname(os.path.realpath(__file__))+"/../"
txtfile   = DIR_CODES+"mesc.txt"
if os.path.isfile(txtfile):
   thefile = open(txtfile,'r')
   lines = thefile.readlines()
   for line in lines:
       line = line.split("#")[0].strip()
       if "mesc_gaussian" in line and " fchk "    in line: GFCHK = line.split('"')[1]
       if "mesc_gaussian" in line and " gauexe "  in line: GEXEC = line.split('"')[1]
   thefile.close()
#-------------------------------------------------#
def check_path():
    abort = False
    if (GEXEC is not None) and (GFCHK is not None):
       if not os.path.exists(GEXEC):
          print "   ERROR: unable to find '%s'"%GEXEC
          abort = True
       if not os.path.exists(GFCHK):
          print "   ERROR: unable to find '%s'"%GFCHK
          abort = True
    else:
       abort = True
       if not os.path.exists(txtfile):
          print "   ERROR: unable to find '%s'"%txtfile
       elif GEXEC is None:
          print "   ERROR: unable to find 'gauexe' in '%s'"%txtfile
       elif GFCHK is None:
          print "   ERROR: unable to find 'fchk' in '%s'"%txtfile
    if abort:
       sys.exit("   Aborting... Check output for more information.")
#-------------------------------------------------#



#------------------#
# READ BASIC FILES #
#------------------#
def readfchk(fchkfile):
    '''
    Output data:
       * x_cc    : a list of atomic coordinates;
                      x_cc = [x1,y1,z1,x2,y2,z2,...]
                   where (xi,yi,zi) are the (x,y,z) coordinates of atom i
       * atonums : a list of atomic numbers
       * ch      : the change of the system
       * mtp     : the multiplicity of the system
       * E       : the energy of the system
       * g_cc    : similar to xyz_list, but it contains the gradient coordinates (None if not found)
       * F_cc    : a list containing the lower-triangular part of the hessian matrix (None if not found)
    '''

    #-------------------------#
    # Variable initialization #
    #-------------------------#
    natoms   = None ; ch       = None ; mtp  = None  ; E = None
    atonums  = []   ; masslist = []
    x_cc     = []   ; g_cc     = []   ; F_cc = []

    #------------------------#
    # The labels to look for #
    #------------------------#
    labels     = {}                              ; found_dict    = {}
    labels[0]  = "Number of atoms"               ; found_dict[0] = False
    labels[1]  = "Charge"                        ; found_dict[1] = False
    labels[2]  = "Multiplicity"                  ; found_dict[2] = False
    labels[3]  = "Total Energy"                  ; found_dict[3] = False
    labels[4]  = "Atomic numbers"                ; found_dict[4] = False
    labels[5]  = "Current cartesian coordinates" ; found_dict[5] = False
    labels[6]  = "Real atomic weights"           ; found_dict[6] = False
    labels[7]  = "Cartesian Gradient"            ; found_dict[7] = False
    labels[8]  = "Cartesian Force Constants"     ; found_dict[8] = False

    idx_basic_labels = [0,1,2,3,4,5,6]

    #-------------------------------------#
    # Introduce folder route in fchk name #
    #-------------------------------------#
    fchk = open(fchkfile,'r')

    #--------------------------#
    # Read info from fchk file #
    #--------------------------#
    for line in fchk:
        # Number of atoms
        if line.startswith(labels[0]):
            found_dict[0] = True
            natoms = int(line.split()[-1])
        # Charge
        elif line.startswith(labels[1]):
            found_dict[1] = True
            ch = int(line.split()[-1])
        # Spin multiplicity
        elif line.startswith(labels[2]):
            found_dict[2] = True
            mtp = int(line.split()[-1])
        # Total Energy
        elif line.startswith(labels[3]):
            found_dict[3] = True
            E = float(line.split()[-1])
        # Atomic Numbers
        elif line.startswith(labels[4]):
            found_dict[4] = True
            length = int(line.split()[-1])
            while len(atonums) != length:
                  nextline = next(fchk)
                  atonums += [int(i) for i in nextline.split()]
        # Cartesian Coordinates
        elif line.startswith(labels[5]):
            found_dict[5] = True
            length = int(line.split()[-1])
            while len(x_cc) != length:
                  nextline = next(fchk)
                  x_cc += [float(i) for i in nextline.split()]
        # List of atomic masses
        elif line.startswith(labels[6]):
            found_dict[6] = True
            length = int(line.split()[-1])
            while len(masslist) != length:
               nextline = next(fchk)
               masslist += [float(i) for i in nextline.split()]
        # Cartesian Gradient
        elif line.startswith(labels[7]) and natoms != 1:
            found_dict[7] = True
            length = int(line.split()[-1])
            while len(g_cc) != length:
                  nextline = next(fchk)
                  g_cc += [float(i) for i in nextline.split()]
        # Cartesian Force Constant Matrix
        elif line.startswith(labels[8]) and natoms != 1:
            found_dict[8] = True
            length = int(line.split()[-1])
            while len(F_cc) != length:
                  nextline = next(fchk)
                  F_cc += [float(i) for i in nextline.split()]
    fchk.close()

    #--------------------------------------------#
    # Some previous things before returning data #
    #--------------------------------------------#
    if g_cc == []: g_cc = None
    if F_cc == []: F_cc = None
    assert natoms == len(atonums),"Problems with number of atoms in %s!!"%fchkfile
    for idx in idx_basic_labels:
       if found_dict[idx] is False:
          key = labels[idx]
          exit("*--> Label '%s' not found <--*"%key)

    #-------------#
    # Return data #
    #-------------#
    return x_cc , atonums, ch , mtp ,  E , g_cc, F_cc

#------------------#
# SEND CALCULATION #
#------------------#
def sendcalc(ifile,ofile,err,folder=None):

    check_path()
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       ifile = folder + ifile
       ofile = folder + ofile
       err   = folder + err

    command = "%s <%s 1>%s 2>%s"%(GEXEC,ifile,ofile,err)
    status = os.system(command)
    # wait a while to be sure Gaussian writes everything
    time.sleep(0.5)
    return status

def genfchk(chk,fchk,err,folder=None):

    check_path()

    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       chk  = folder + chk
       fchk = folder + fchk
       err  = folder + err
    command = "%s %s %s 1>%s 2>&1"%(GFCHK,chk,fchk,err)
    status = os.system(command)
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
    ifile = folder + name + ".gjf"
    ofile = folder + name + ".out"
    chk   = folder + name + ".chk"
    fchk  = folder + name + ".fchk"
    err   = folder + name + ".err"
    return ifile, ofile, chk, fchk, err


#-----------------------#
# SOME USEFUL FUNCTIONS #
#-----------------------#
def fchk2gts(mainname):
    from helpfns import get_pgs
    from gtsfile import write_gtsfile
    import constants as cons
    x_cc , atonums, ch , mtp ,  E , g_cc, F_cc = readfchk(mainname+".fchk")
    masslist = [cons.dict_atomasses[atonum] for atonum in atonums]
    pgroup, sigma = get_pgs(atonums,masslist, x_cc)
    write_gtsfile(x_cc , atonums, ch , mtp ,  E , pgroup, sigma, g_cc, F_cc, mainname+".gts")


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
#----------------------------------------#

def q2dtor_defaults(level,ch,mtp,torsion1,torsion2,ttype="min"):

    NPROC, RAM = 1,1
    the_string = ""
    
    #-----------------------------------#
    # Default string for PES point file #
    #-----------------------------------#
    the_string = the_string + "#-----------------------------#\n"
    the_string = the_string + "start_scangeom gaussian\n"
    the_string = the_string + "%%nproc=%i\n"%NPROC
    the_string = the_string + "%%mem=%iGB\n"%RAM
    the_string = the_string + "%%chk=%s.chk\n"%q2dtor_key2
    the_string = the_string + "#p %s\n"%level
    the_string = the_string + "   scf=tight NoSymmetry %s\n"%q2dtor_key3
    if ttype=="min":
       the_string = the_string + "   opt=(tight,modredundant)\n"
    if ttype=="ts":
       the_string = the_string + "   opt=(tight,modredundant,ts,calcfc,noeigentest)\n"
    the_string = the_string + "\n"
    the_string = the_string + "Geometry optimization at the 2D-PES\n"
    the_string = the_string + "\n"
    the_string = the_string + "%i %i\n"%(ch,mtp)
    the_string = the_string + "%s\n"%q2dtor_key1
    the_string = the_string + "\n"
    the_string = the_string + "%2i %2i %2i %2i F\n"%tuple([ii+1 for ii in torsion1])
    the_string = the_string + "%2i %2i %2i %2i F\n"%tuple([ii+1 for ii in torsion2])
    the_string = the_string + "\n"
    the_string = the_string + "end_scangeom\n"
    the_string = the_string + "#-----------------------------#\n"

    the_string = the_string + "\n\n\n"

    #----------------------------#
    # Default string for CP file #
    #----------------------------#
    for i in [0,1,2]:
        the_string = the_string + "#-----------------------------#\n"
        if i == 0:
           the_string = the_string + "start_sp0 gaussian\n"
           if ttype=="min": opt = "opt=tight"
           if ttype=="ts":  opt = "opt=(tight,ts,calcfc,noeigentest)"
           comment = "Optimization of a minimum in the 2D-PES"
        if i == 1:
           the_string = the_string + "start_sp1 gaussian\n"
           if ttype=="min": opt = "opt=(tight,ts,calcfc,noeigentest)"
           if ttype=="ts":  opt = "opt=(tight,saddle=2,calcfc,noeigentest)"
           comment = "Optimization of a saddle point in the 2D-PES"
        if i == 2:
           the_string = the_string + "start_sp2 gaussian\n"
           if ttype=="min": opt = "opt=(tight,saddle=2,calcfc,noeigentest)"
           if ttype=="ts":  opt = "opt=(tight,saddle=3,calcfc,noeigentest)"
           comment = "Optimization of a maximum in the 2D-PES"
        the_string = the_string + "%%nproc=%i\n"%NPROC
        the_string = the_string + "%%mem=%iGB\n"%RAM
        the_string = the_string + "%%chk=%s.chk\n"%q2dtor_key2
        the_string = the_string + "#p %s\n"%level
        the_string = the_string + "   scf=verytight\n"
        the_string = the_string + "   %s\n"%opt
        the_string = the_string + "\n"
        the_string = the_string + "%s\n"%comment
        the_string = the_string + "\n"
        the_string = the_string + "%i %i\n"%(ch,mtp)
        the_string = the_string + "%s\n"%q2dtor_key1
        the_string = the_string + "\n"
        the_string = the_string + "--Link1--\n"
        the_string = the_string + "%%nproc=%i\n"%NPROC
        the_string = the_string + "%%mem=%iGB\n"%RAM
        the_string = the_string + "%%chk=%s.chk\n"%q2dtor_key2
        the_string = the_string + "#p %s\n"%level
        the_string = the_string + "   scf=verytight freq=noraman geom=allcheck\n"
        the_string = the_string + "\n"
        if i == 0:
           the_string = the_string + "end_sp0\n"
        if i == 1:
           the_string = the_string + "end_sp1\n"
        if i == 2:
           the_string = the_string + "end_sp2\n"
        the_string = the_string + "#-----------------------------#\n"
        the_string = the_string + "\n\n\n"

    return the_string


def q2dtor_pespoint(ilines, xvec,symbols,name,previous=None):
    '''
    xvec = [(x1,y1,z1),(x2,y2,z2),...]
    mos  = main name for reference Molecular Orbitals
           if mos is "system", then "system.chk" will be used as
           MO guess; file mos.chk is in folder
    '''


    #------------------------------#
    # Calculation folder and files #
    #------------------------------#
    mainname = "scangeom_"+name
    folder   = "TMP_CALCS/"
    if not os.path.exists(folder): os.mkdir(folder)
    ifile, ofile, chk, fchk, err = iofiles(mainname,folder)

    #----------#
    # MO guess #
    #----------#
    guessread = False
    if previous is not None:
       ref_ifile, ref_ofile, ref_chk, ref_fchk, ref_err = iofiles("scangeom_"+previous,folder)
       if (not os.path.exists(chk)) and (os.path.exists(ref_chk)):
          shutil.copyfile(ref_chk, chk)
       guessread = True

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
           pos = line.find(q2dtor_key2)
           line = line[0:pos] + folder + mainname + line[pos+len(q2dtor_key2):]
        if q2dtor_key3 in line:
           pos = line.find(q2dtor_key3)
           if guessread: line = line[0:pos] + "guess=read" + line[pos+len(q2dtor_key3):]
           else        : line = line[0:pos] + "          " + line[pos+len(q2dtor_key3):]
        string_ifile = string_ifile + line 

    #---------------------------#
    # Send Gaussian calculation #
    #---------------------------#
    ff = open(ifile,"w")
    ff.write(string_ifile)
    ff.close()
    status = sendcalc(ifile,ofile,err)

    # Check output
    ff = open(ofile, 'r')
    outlines = ff.readlines()
    ff.close()
    if "Normal termination" not in outlines[-1]: return None

    # Generate fchk
    status = genfchk(chk,fchk,err)

    #-----------#
    # Read data #
    #-----------#
    xcc, atonums, ch, mtp,  Etot, gcc, Fcc = readfchk(fchk)

    #-------------#
    # Remove data #
    #-------------#
    #os.remove(err)
    #os.remove(ifile)
    #os.remove(fchk)
    #os.remove(ofile)
    #os.remove(chk) # not removed, as it is useful for next calculation

    return xcc, atonums, ch, mtp,  Etot, gcc, Fcc


def q2dtor_optCP(ilines,xvec,symbols,mainname,phis=(None,None)):
    '''
    xcc = [(x1,y1,z1),(x2,y2,z2),...]
    * phis is included just to have the same arguments as in mesc_orca.py
    '''

    #------------------------------#
    # Calculation folder and files #
    #------------------------------#
    folder = "TMP_CALCS/"
    if not os.path.exists(folder): os.mkdir(folder)
    ifile, ofile, chk, fchk, err = iofiles(mainname,folder)

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
           pos = line.find(q2dtor_key2)
           line = line[0:pos] + folder + mainname + line[pos+len(q2dtor_key2):]
        string_ifile = string_ifile + line

    # Send calculation
    ff = open(ifile,"w")
    ff.write(string_ifile)
    ff.close()
    print "        * Geometry optimization + Freq calculation with Gaussian..."
    status = sendcalc(ifile,ofile,err)

    # Check output
    ff = open(ofile, 'r')
    outlines = ff.readlines()
    ff.close()
    if "Normal termination" not in outlines[-1]: return None

    # Get data
    print "        * Creating '%s' file"%fchk
    status = genfchk(chk,fchk,err)
    xcc, atonums, ch, mtp,  Etot, gcc, Fcc = readfchk(fchk)

    #print "        * Remove calculation files (gjf,out,chk,fchk)..."
    #os.remove(ifile)
    #os.remove(ofile)
    #os.remove(chk  )
    #os.remove(fchk )
    #os.remove(err  )

    return xcc, atonums, ch, mtp,  Etot, gcc, Fcc


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

def thera_default(ch=0,mtp=1,nproc=1,mem=1):
    string = ""
    string = string + "%%nproc=%i\n"%nproc
    string = string + "%%mem=%iGB\n"%mem
    string = string + "%%chk=%s.chk\n"%thera_key2
   #string = string + "#p hf sto-3g scf=tight NoSymm %s\n"%thera_key3
    string = string + "#p mpwb95/6-31+G(d,p) IOp(3/76=0560004400) scf=verytight int=ultrafine NoSymm %s\n"%thera_key3
    string = string + "\n"
    string = string + "Input file for MEP calculation\n"
    string = string + "\n"
    string = string + "%i %i\n"%(ch,mtp)
    string = string + "%s\n"%thera_key1
    string = string + "\n"
    return string

def thera_datainfolder(folder):

    if not folder.endswith("/"): folder = folder + "/"
    # List of fchk files
    file_list = [folder+filename for filename in os.listdir(folder) if filename.lower().endswith("fchk")]
    # et data from fchk files
    data = []
    for fchkfile in file_list:
        mainname = fchkfile[:-5]
        # read file
        x_cc, atonums, ch, mtp,  E, g_cc, F_cc = readfchk(fchkfile)
        # append data
        data.append( [mainname, x_cc, atonums, ch, mtp, E, g_cc, F_cc, []] )
    # Sort by energy
    data.sort(key=lambda xx: xx[5])
    # Return data
    return data

def thera_spc(ref_lines, xvec, symbols, pointname, hessian=False, folder=None):

    ifile, ofile, chk, fchk, err = iofiles(pointname,folder)

    #--------------#
    # Input string #
    #--------------#
    string_ifile = ""
    for line in ref_lines:
        if thera_key1 in line:
           line = ""
           for idx in range(len(symbols)):
               symbol = symbols[idx]
               x,y,z  = [value*angstrom for value in xvec[3*idx:3*idx+3]]
               line   = line + "%2s   %+11.6f  %+11.6f  %+11.6f\n"%(symbol,x,y,z)
        if thera_key2 in line:
           pos = line.find(thera_key2)
           line = line[0:pos] + chk[:-4] + line[pos+len(thera_key2):]
        if thera_key3 in line:
           pos = line.find(thera_key3)
           if hessian: line = line[0:pos] + "freq=noraman" + line[pos+len(thera_key3):]
           else      : line = line[0:pos] + "force       " + line[pos+len(thera_key3):]
        string_ifile = string_ifile + line

    #---------------------------#
    # Send Gaussian calculation #
    #---------------------------#
    ff = open(ifile,"w")
    ff.write(string_ifile)
    ff.close()
    status = sendcalc(ifile,ofile,err)

    # Check output
    ff = open(ofile, 'r')
    outlines = ff.readlines()
    ff.close()
    if "Normal termination" not in outlines[-1]: return None

    # Get data
    status = genfchk(chk,fchk,err)
    xcc, atonums, ch, mtp,  Etot, gcc, Fcc = readfchk(fchk)

    #-------------#
    # Remove data #
    #-------------#
    #os.remove(ifile)
    #os.remove(ofile)
    #os.remove(fchk)
    #os.remove(chk) # not removed, as it is useful for next calculation
    #os.remove(err)

    return xcc, atonums, ch, mtp,  Etot, gcc, Fcc
 

