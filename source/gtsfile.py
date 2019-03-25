'''
*----------------------------------*
 Q2DTor and TheRa Programs

 Copyright (c) 2018 Universidade de Santiago de Compostela

 This file is part of both Q2DTor and TheRa softwares.

 Q2DTor and TheRa are free softwares: you can redistribute it and/or modify
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

 A module to read and write gts files

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
''' 
'''
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

def getdata_folder(folder):

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

