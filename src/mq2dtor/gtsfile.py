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

 A module to read and write gts files

''' 

import os
import constants  as cons
import helpfns as hf


def read_gtsfile(gtsfile):

    gts_lines = hf.readfile(gtsfile)

    # Read cartesian coordinates
    info_cc = hf.select_lines(gts_lines,"start_cc","end_cc")
    xyz_list = []
    atonum_list = []
    for line in info_cc:
        atonum , x , y , z = line.split()
        xyz_list.append( (float(x),float(y),float(z)) )
        atonum_list.append(int(atonum))

    # Read basic information
    info_basic = hf.select_lines(gts_lines,"start_basic","end_basic")
    charge, multiplicity, energy, pgroup, rotsigma = None, None, None, None, None
    for line in info_basic:
        if line.startswith("charge"):       charge       = int(line.split()[1])
        if line.startswith("multiplicity"): multiplicity = int(line.split()[1])
        if line.startswith("energy"):       energy       = float(line.split()[1])
        if line.startswith("pointgroup"):   pgroup       = line.split()[1]
        if line.startswith("rotsigma"):     rotsigma     = int(line.split()[1])

    # Read cartesian gradient
    info_grad = hf.select_lines(gts_lines,"start_grad","end_grad")
    g_list = []
    for line in info_grad:
        gx , gy , gz = line.split()
        g_list.append( (float(gx),float(gy),float(gz)) )
    if g_list == []: g_list = None

    # Read Hessian matrix
    info_hess = hf.select_lines(gts_lines,"start_hess","end_hess")
    F_list = []
    for line in info_hess:
        F_list += [float(Fij) for Fij in line.split()]
    if F_list == []: F_list = None

    # Read frequencies (only special cases)
    info_freqs = hf.select_lines(gts_lines,"start_freqs","end_freqs")
    freq_list = []
    for line in info_freqs:
        freq_list += [float(freqs) for freqs in line.split()]
    if freq_list == []: freq_list = None

    # Return data
    return xyz_list , atonum_list , charge, multiplicity, energy, pgroup, rotsigma, g_list , F_list, freq_list

def write_gtsfile(xyz_list , atonum_list , charge, multiplicity, energy, \
                   pgroup, rotsigma , g_list , F_list, gtsfile, freqs=None):

    natoms   = len(atonum_list)
    xyz_list = hf.xvecformat(xyz_list,natoms,out="Nx3")
    if g_list is not None: g_list = hf.xvecformat(g_list,natoms,out="Nx3")


    gts = open(gtsfile,"w")

    # Write atomic numbers and cartesian coordinates
    gts.write("# Atomic number and non-scaled cartesian coordinates [bohr]\n")
    gts.write("start_cc\n")
    for atonum, xyz in zip(atonum_list,xyz_list):
        x,y,z = xyz
        gts.write("   %03i   %+15.8E  %+15.8E  %+15.8E\n"%(atonum,x,y,z))
    gts.write("end_cc\n\n")


    # Write basic data
    rotsigma = "%02i   "%rotsigma
    while len(pgroup) < 5 : pgroup = pgroup+" " 
    gts.write("# Charge, multiplicity, energy [hartree], point group, rotational symmetry number\n")
    gts.write("start_basic\n")
    gts.write("   charge        %i\n"%charge)
    gts.write("   multiplicity  %i\n"%multiplicity)
    gts.write("   energy       %+15.8E # Total energy in hartree\n"%energy)
    gts.write("   pointgroup    %5s          # Point group from automatic assignation\n"%pgroup)
    gts.write("   rotsigma      %5s          # Rotational sigma from automatic assignation\n"%rotsigma)
    gts.write("end_basic\n\n")


    # Write cartesian gradiend
    if g_list is not None and g_list != []:
       gts.write("# Non-scaled cartesian gradient [hartree/bohr]\n")
       gts.write("start_grad\n")
       for gx,gy,gz in g_list:
           gts.write("    %+15.8E  %+15.8E  %+15.8E\n"%(gx,gy,gz))
       gts.write("end_grad\n\n")


    # Write force constant matrix (i.e. hessian matrix)
    if F_list is not None and F_list != []:
       gts.write("# Low triangular part of symmetric force constant (hessian) matrix [hartree/bohr^2]\n")
       gts.write("# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...\n")
       gts.write("start_hess\n")
       for idx in range(0,len(F_list),5):
           line = "  ".join(["%+15.8E"%Fij for Fij in F_list[idx:idx+5]])
           gts.write("           %s\n"%line)
       gts.write("end_hess\n\n")

    # Write freqs
    if F_list is None and freqs is not None:
       gts.write("start_freqs\n")
       for idx in range(0,len(freqs),5):
           line = "  ".join(["%7.2f"%(freq) for freq in freqs[idx:idx+5]])
           gts.write("           %s\n"%line)
       gts.write("end_freqs\n\n")

    gts.close()


#----------------------------------------#
#    _____  _            ____            #
#   |_   _|| |__    ___ |  _ \  __ _     #
#     | |  | '_ \  / _ \| |_) |/ _` |    #
#     | |  | | | ||  __/|  _ <| (_| |    #
#     |_|  |_| |_| \___||_| \_\\__,_|    #
#                                        #
#----------------------------------------#

#def thera_default(ch=0,mtp=1,nproc=4,mem=4):
#    pass

def thera_datainfolder(folder):

    if not folder.endswith("/"): folder = folder + "/"
    # List of gts files
    file_list = [folder+gtsfile for gtsfile in os.listdir(folder) if gtsfile.lower().endswith("gts")]
    # Get data from gtsfile
    data = []
    for gtsfile in file_list:
        mainname = gtsfile[:-4]
        # read file
        x_cc, atonums, ch, mtp,  E, pgroup, rotsigma, g_cc, F_cc, freqs = read_gtsfile(gtsfile)
        # append data
        data.append( [mainname, x_cc, atonums, ch, mtp, E, g_cc, F_cc, freqs] )
    # Sort by energy
    data.sort(key=lambda xx: xx[5])
    # Return data
    return data

