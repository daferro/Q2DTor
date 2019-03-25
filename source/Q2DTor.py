#!/usr/bin/env python
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


#--------------------------------------------------------------------#
# Software name: Q2DTor                                              #
# Function     : Consideration of torsional anharmonicity due to two #
#                torsions in the calculation of partition functions. #
#--------------------------------------------------------------------#
# Main Author        :  David Ferro-Costas                           #
# Last Update        :  Mar-02nd-2018                                #
# Last Update (v1.1) :  Jan-09th-2019                                #
#--------------------------------------------------------------------#

Index to find things easily in the code to find things easily in the code:

  ( 0) Implemented Softwares
       * import_mes(software,function)

  ( 1) Main function and arguments
       * main()
       * get_arguments(head_message,hstring)

  ( 2) Q2DTor logo and help message
       * q2dtor_getlogo()
       * q2dtor_gethelpstring()

  ( 3) Q2DTor options
       * q2dtor_init(name,inpdata,files2read,files2write,software)
       * q2dtor_pes(name,inpdata,files2read,files2write,f_ics)
       * q2dtor_fourier(name,inpdata,files2read,files2write,TOTPESPOINTS)
         - string_fitting(fterms, popt, ignore, errors, fitting_time=None)
       * q2dtor_findsp(name,inpdata,files2read,files2write)
       * q2dtor_optsp(name,inpdata,files2read,files2write,mode="all")
       * q2dtor_dijfitting(name,tuple_fit,icoords,symmetry,tsigma1,tsigma2,f_pes,file_Dfit,masslist=None)
       * q2dtor_tor2dns(name,inpdata,files2read,files2write,Tlist)
       * q2dtor_rovibpf(name,inpdata,files2read,files2write,Tlist,thermo=False)

  ( 4) Q2DTor tools
       * tools(name,avail_args)
       * tool_pdf(name,ifile)
         - sp_colors(cpnames)
       * tool_gts(mainname,software)
       * tool_icoords(name,mode,f_ics)
         - gtsicoords(f_gts, f_ics)

  ( 5) Helper functions related to: symmetry
       * pessym_inregion(phi1, phi2, symmetry="", tsigma1=1, tsigma2=1)
       * pessym_apply(phi1, phi2, symmetry="", tsigma1=1, tsigma2=1)
       * pessym_points(dphi1,dphi2,symmetry="",tsigma1=1,tsigma2=1)
       * pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2,pprint=False)

  ( 6) Helper functions related to: stationary points
       * sp_optimization(cp_name,symbols,dict_CPs,dict_pes,functionOPT,ilines)
       * sp_checking(gtsfile,cp_name,cptype,torsion1,torsion2,tolerance,dict_warns,symmetry,tsigma1,tsigma2,ttype)
       * sp_sorting(dict_CPs)
       * sp_icdepuration(dict_CPs,torsion1,torsion2,xyzStruct,f_ics)
       * sp_allCPs(dict_CPs,symmetry,name,tsigma1,tsigma2)
       * sp_analysis(f_spinfo,dict_CPs,icoords,freqscal=1.0,masslist=None)

  ( 7) Helper functions related to: 2DNS
       * tor2DNS_hamiltonian(molname, V_terms, V_coefs, LkinD, dijtuple, coef_nums, kmax=100)
       * tor2DNS_diagonalization(kmax,fil,col,ham2d,LkinD,file_evalues,ENERMAX=10000.,EIGMAX=250,MAXCYCLE=100)

  ( 8) Helper functions related to: integration
       * integration2D(integrand,dphi=2.0,args=[])
       * intgrnd_ClasTor_pf(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0)
       * intgrnd_ClasTor_fd(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0)
       * intgrnd_ClasTor_sd(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0)
       * intgrnd_EHR_pf(phi1,phi2,T,Fourier_PES,rel_sqrtS,sqrtS0,delaunay,refEHR,bb_list)
       * intgrnd_EHR_fd(phi1,phi2,T,Fourier_PES,rel_sqrtS,sqrtS0,delaunay,refEHR,bb_list)
       * intgrnd_EHR_sd(phi1,phi2,T,Fourier_PES,rel_sqrtS,sqrtS0,delaunay,refEHR,bb_list)

  ( 9) Helper functions related to: strings
       * string_checkname(thename,exit=True)
       * string_getpname(phi1,phi2,molname="",units="rad")
       * string_getinfo(pname)
       * string_pfntable(Tlist,pf1WHO,pfMSHO,pfE2DT,pf2DNS,pfTorClas,pfEHR,ratio2DNS,Eref)
       * string_thermotable(Tlist,ther1WHO,therMSHO, therE2DT)

  (10) Helper functions related to: diverse things
       * ic_checktorsions(icoords,torsion1,torsion2)
       * clean_close(dict_cps,basic=[],diff_angle=4.0,diff_Ecm=10.0)
       * change_Eref(old_Eref,new_Eref,f_pes,f_vfit)
       * pfn_derivatives(structure,T,k="cc")

  (11) Functions for reading files
       * get_defaults():
       * deal1Dfourierstring(string1D):
       * deal2Dfourierstring(string2D):
       * readfile_inp(inputname):
       * readfiles(filenames,which):
       * readfile_ics(nricfile)
       * readfile_calcs(f_calcs,start,end)
       * readfile_pes(f_pes,exclude=True)
       * readfile_xfit(filename,which="V")
       * readfile_splist(spfile)
       * readfile_evals(f_evals)

  (12) Functions for writing files
       * writefile_inp(name)
       * writefile_ics(icoords,nricfile,bonds)
       * writefile_xfit(fterms,popt,ofile,mode="w")
       * writefile_spinfo(f_spinfo,dict_CPs,icoords,Tlist,freqscal=1.0,masslist=None)
'''

# Tesselation
TESSMODE     = 6
TESSMODEPLOT = 4

# Q2DTor infor
VERSION  = "1.1"
DATE     = "January, 2019"

# Q2DTor folders
DIRFILES = "IOfiles/"
DIRSP    = "IOfiles/SP/"

# Q2DTor error message
EXITMESS = "ERROR while executing Q2DTor!"

# updates
BOOL_v1_1 = True

# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #
#                  Python libraries              #
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #
import argparse
import datetime
import numpy as np
import os
import shutil
import random
import sys
import time
import warnings

from scipy.interpolate   import SmoothBivariateSpline
from scipy.interpolate   import interp2d
from scipy.optimize      import fmin
from scipy.optimize      import newton
from scipy.sparse        import csr_matrix
from scipy.sparse        import csc_matrix
from scipy.sparse.linalg import eigsh
from scipy.spatial       import Delaunay
#random.seed(3333)
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #


# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #
#                  Own libraries                 #
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #
import mq2dtor.constants   as     cons
import mq2dtor.helpfns     as     hf
from   mq2dtor.classes     import xyz2Struct
from   mq2dtor.classes     import basic2Struct
from   mq2dtor.classes     import gts2Struct
from   mq2dtor.classes     import Struct
from   mq2dtor.classes     import Fourier2D
from   mq2dtor.classes     import Logger
from   mq2dtor.gtsfile     import write_gtsfile
from   mq2dtor.quotes      import random_quote
import mq2dtor.tesselation as     tess
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #

# About energy units
HU1 = cons.kcalmol # human-units 1
HU2 = cons.calmol  # human-units 1
SU1 = "kcal/mol"   # string-units 1
SU2 = "cal/mol/K"  # string-units 1

# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 0) Implemented Softwares                       #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def import_mes(software,function):
    if software == "gaussian":
       if function == "defaults":
          global mes_defaults
          from mq2dtor.mesc_gaussian import q2dtor_defaults as mes_defaults
       if function == "pespoint":
          global mes_pespoint
          from mq2dtor.mesc_gaussian import q2dtor_pespoint as mes_pespoint
       if function == "optCP":
          global mes_optCP
          from mq2dtor.mesc_gaussian import q2dtor_optCP    as mes_optCP

    elif software == "orca":
       if function == "defaults":
          global mes_defaults
          from mq2dtor.mesc_orca     import q2dtor_defaults as mes_defaults
       if function == "pespoint":
          global mes_pespoint
          from mq2dtor.mesc_orca     import q2dtor_pespoint as mes_pespoint
       if function == "optCP":
          global mes_optCP
          from mq2dtor.mesc_orca     import q2dtor_optCP    as mes_optCP

    else:
       print "        - ERROR: Unknown software (%s)\n"%software
       sys.exit(EXITMESS)
#----------------------------------------------------------#
#----------------------------------------------------------#




# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 1) Main function and arguments                 #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def main():
    '''
    The main function, which is in charge
    of acting according any user option
    '''

    now = datetime.datetime.now()

    OPTSTRING1  = "   *==============================\n"
    OPTSTRING1 += "   | Executing Q2DTor with option: %s  [%s]\n"
    OPTSTRING1 += "   *==============================\n"
    OPTSTRING1 += "     Data: %s\n"%(now.strftime("%Y-%m-%d %H:%M"))

    OPTSTRING2  = "\n"
    OPTSTRING2 += "                                          | Elapsed time %.1f %s |\n"

    OPTSTRING3  = " "*3 + "\\"*23+" "+"/"*23 + "\n"
    OPTSTRING3 += " "*3 + "    End of execution of Q2DTor with %s  \n"
    OPTSTRING3 += " "*3 + "\\"*23+" "+"/"*23 + "\n"
    OPTSTRING3 += "\n\n"


    #---------------------#
    # Deal with arguments #
    #---------------------#
    head_message     = q2dtor_getlogo()
    help_string      = q2dtor_gethelpstring()
    name, avail_args = get_arguments(head_message,help_string)

    #-----------#
    # Any tool? #
    #-----------#
    if tools(name, avail_args): sys.exit()

    #-----------------#
    # Input file name #
    #-----------------#
    if name.endswith(".inp"):
       ifile = name
       name  = name[:-4]
    else:
       ifile = name + ".inp"
    string_checkname(name,exit=True)
    
    #---------------------#
    # Generate Input file #
    #---------------------#
    if avail_args.keys() == []:
      if not os.path.isfile(ifile):
         writefile_inp(name)
         sys.exit("Q2DTor input file '%s' has been generated!\nPlease, check it and modify it"%ifile)
      sys.exit("Q2DTor input file '%s' already exists!"%ifile)

    #-----------------#
    # Read Input file #
    #-----------------#
    if not os.path.isfile(ifile): sys.exit("Q2DTor input file '%s' does not exists!"%ifile)
    input_data   = readfile_inp(ifile)
    torsions     = input_data[0]
    calcs        = input_data[1]
    pes          = input_data[2]
    fourier      = input_data[3]
    statpoints   = input_data[4]
    tor2dns      = input_data[5]
    rovibpf      = input_data[6]
    Tlist        = input_data[7]
    iofiles      = input_data[8]
    # Expand each tuple
    (torsion1,torsion2,tsigma1,tsigma2) = torsions
    (ttype,level,charge,multiplicity) = calcs
    (t1step,t2step,symmetry) = pes
    (fterms,nums,weight,ignore) = fourier
    (tolerance,freqscal) = statpoints
    (dijvar,kmax,maxeigen) = tor2dns
    (interpolation,integrationstep) = rovibpf
    (f_xyz,f_calcs,f_ics,f_pes,f_vfit,f_dfit,f_splist,f_spinfo,f_evals,f_sqrt,f_tables,f_pdf,f_out) = iofiles
    # Add folder name to files
    f_ics     = DIRFILES + f_ics
    f_pes     = DIRFILES + f_pes
    f_vfit    = DIRFILES + f_vfit
    f_dfit    = DIRFILES + f_dfit
    f_splist  = DIRFILES + f_splist
    f_spinfo  = DIRFILES + f_spinfo
    f_evals   = DIRFILES + f_evals
    f_sqrt    = DIRFILES + f_sqrt
    f_tables  = DIRFILES + f_tables

    #-------------------------#
    # Print data in terminal? #
    #-------------------------#
    tprint = False
    if "--print" in avail_args.keys(): tprint = True
    sys.stdout = Logger(f_out,"a",tprint)

    #--------#
    # units? #
    #--------#
    if "--SI" in avail_args.keys():
       global HU1, HU2, SU1, SU2
       HU1 = cons.kjmol # human-units 1
       HU2 = cons.jmol  # human-units 1
       SU1 = "kJ/mol"   # string-units 1
       SU2 = "J/mol/K"  # string-units 1

    #-----------------------#
    # Now, option by option #
    #-----------------------#
    if "--init" in avail_args.keys():
      sys.stdout = Logger(f_out,"w",tprint)
      print q2dtor_getlogo()
      #---------------
      print OPTSTRING1%("--init", "Initialization")
      t1 = time.time()
      #---------------
      software = avail_args["--init"].lower()
      #---------------
      inpdata     = (level,charge,multiplicity,torsion1,torsion2,ttype)
      files2read  = (f_xyz,)
      files2write = (f_ics,f_calcs)
      #---------------
      q2dtor_init(name,inpdata,files2read,files2write,software)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--init")

    elif "--pes" in avail_args.keys():
      print OPTSTRING1%("--pes", "Calculation of the 2D-PES")
      t1 = time.time()
      #---------------
      inpdata     = (torsion1,torsion2,tsigma1,tsigma2,t1step,t2step,symmetry)
      files2read  = (f_xyz,f_calcs)
      files2write = (f_pes,)
      #---------------
      q2dtor_pes(name,inpdata,files2read,files2write,f_ics)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--pes")

    elif "--fourier" in avail_args.keys():
      print OPTSTRING1%("--fourier", "Fitting the numerical 2D-PES to Fourier series")
      t1 = time.time()
      #---------------
      dphi1 = t1step * cons.R2D
      dphi2 = t2step * cons.R2D
      TOTPESPOINTS = int(round(360.0/dphi1)+1) * int(round(360.0/dphi2)+1)
      #---------------
      inpdata     = (symmetry,tsigma1,tsigma2,fourier)
      files2read  = (f_pes,)
      files2write = (f_vfit,)
      #---------------
      q2dtor_fourier(name,inpdata,files2read,files2write,TOTPESPOINTS)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--fourier")

    elif "--findsp" in avail_args.keys():
      print OPTSTRING1%("--findsp", "Looking for stationary points in the fitted 2D-PES")
      t1 = time.time()
      #---------------
      inpdata     = (symmetry,tsigma1,tsigma2,tolerance)
      files2read  = (f_vfit,)
      files2write = (f_splist,)
      #---------------
      q2dtor_findsp(name,inpdata,files2read,files2write)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--findsp")

    elif "--optsp" in avail_args.keys():
      print OPTSTRING1%("--optsp" , "Optimizing the stationary points"   )
      t1 = time.time()
      #---------------
      mode = avail_args["--optsp"]
      #---------------
      inpdata     = (torsion1,torsion2,symmetry,tsigma1,tsigma2,tolerance,ttype,freqscal)
      files2read  = (f_calcs, f_xyz, f_pes, f_vfit, f_splist)
      files2write = (f_splist,f_ics,f_spinfo)
      #---------------
      q2dtor_optsp(name,inpdata,files2read,files2write,Tlist,mode=mode)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--optsp")

    elif '--tor2dns' in avail_args.keys():
      print OPTSTRING1%("--tor2dns" , "Calculating the 2D-NS torsional eigenvalues")
      t1 = time.time()
      #---------------
      inpdata     = (torsion1,torsion2,fourier,symmetry,tsigma1,tsigma2,dijvar,kmax,maxeigen)
      files2read  = (f_xyz,f_vfit,f_ics,f_pes,f_splist)
      files2write = (f_dfit,f_evals)
      #---------------
      q2dtor_tor2dns(name,inpdata,files2read,files2write,Tlist)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--tor2dns")

    elif '--rovibpf' in avail_args.keys():
      print OPTSTRING1%("--rovibpf" , "Calculating 1WHO, MSHO, and E2DT partition functions")
      t1 = time.time()
      #---------------
      thermo = False
      if avail_args["--rovibpf"].lower() == "thermo": thermo = True
      #---------------
      inpdata     = (torsion1,torsion2,ttype,symmetry,tsigma1,tsigma2,interpolation,integrationstep,freqscal)
      files2read  = (f_xyz,f_pes,f_vfit,f_splist,f_ics,f_evals)
      files2write = (f_sqrt,f_tables)
      #---------------
      q2dtor_rovibpf(name,inpdata,files2read,files2write,Tlist,thermo)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--rovibpf")
      print
      print random_quote("Q2DTor")

    #----------------------------------#
    # recalculate PES by interpolation #
    #----------------------------------#
    elif '--interpes' in avail_args.keys():
      print OPTSTRING1%("--interpes" , "Temporal option")
      t1 = time.time()
      #---------------
      inpdata = (symmetry,tsigma1,tsigma2)
      files2read  = (f_pes,f_vfit,f_splist)
      q2dtor_interpes(name,inpdata,files2read)
      #---------------
      t2 = time.time()
      print OPTSTRING2%(hf.time2human(t2-t1,"secs"))
      print OPTSTRING3%("--interpes")
      print
#----------------------------------------------------------#
def get_arguments(head_message,hstring):
    
    message = ""
    message = message+"basic usage: python Q2DTor.py system [argument]\n"
    message = message+"for help   : python Q2DTor.py -h"

    # Get args
    args = sys.argv[1:]

    # No arguments, help or version
    if len(args) == 0: sys.exit(message)
    if "-h" in args or "--help"    in args: sys.exit(head_message+hstring)
    if "-v" in args or "--version" in args: sys.exit("Q2DTor v%s (%s)"%(VERSION,DATE))
    
    # Get main file
    ifile = args.pop(0)
    if "--" in ifile: sys.exit(message)

    # Definition of arguments
    avail_args  = {}
    # Format: avail_args[key] = [boolean,k,default]
    # k is the number of options (k>0 mandatory, k<0 optional)
    avail_args["--init"    ] = (False ,-1 , ["gaussian"])
    avail_args["--pes"     ] = (False , 0 , [          ])
    avail_args["--fourier" ] = (False , 0 , [          ])
    avail_args["--findsp"  ] = (False , 0 , [          ])
    avail_args["--optsp"   ] = (False ,-1 , [ "all"    ])
    avail_args["--tor2dns" ] = (False , 0 , [          ])
    avail_args["--rovibpf" ] = (False ,-1 , ["pfs"     ])
    avail_args["--print"   ] = (False , 0 , [          ])
    avail_args["--interpes"] = (False , 0 , [          ])

    # Units 
    avail_args["--SI"      ] = (False , 0 , [          ])

    # Tools
    avail_args["--pdf"    ] = (False , 0 , [          ])
    avail_args["--gts"    ] = (False ,-1 , ["gaussian"])
    avail_args["--icoords"] = (False ,-2 , ["sp",None ])

    # Check arguments
    message_error1 = "Argument '%s':\n   * requires %i values\n   * user introduced: %s"
    message_error2 = "Argument '%s':\n   * accepts up to %i values\n   * user introduced: %s"

    flexible = False
    for idx in range(len(args)):
        coincidences = []
        for key,value in avail_args.items():
            boolean, k, default = value
            if flexible:
               if not key.startswith(args[idx]): continue
            else:
               if key != args[idx]: continue
            coincidences.append(key)
            # Get possible option
            options = []
            for arg in args[idx+1:]:
                if arg.startswith("--"): break
                options.append(arg)
            # No options
            if k==0:
               avail_args[key] = (True, k, None)
            # Mandatory options
            elif k>0:
               if len(options) != k: sys.exit(message_error1%(key,k," ".join(options)))
               if k == 1: avail_args[key] = (True, k, options[0])
               else     : avail_args[key] = (True, k, options)
            # Optional options
            elif k<0:
               if len(options) > abs(k): sys.exit(message_error2%(key,k," ".join(options)))
               for idx2 in range(len(options)): default[idx2] = options[idx2]
               if k == -1: avail_args[key] = (True, k, default[0])
               if k != -1: avail_args[key] = (True, k, default   )
        # Ambiguous?
        if len(coincidences) > 1:
           sys.exit("Ambiguous argument (%s): "%args[idx] + str(coincidences))

    # Check arguments and return
    toreturn = {}
    for key,(boolean, k, option) in avail_args.items():
        if boolean: toreturn[key] = option
    return ifile, toreturn
#----------------------------------------------------------#
#----------------------------------------------------------#




# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 2) Q2DTor logo and help message                #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def q2dtor_getlogo():
    LOGO = random.choice(["A","B","C","D","E","F","G"])

    head_message = "\n"
    head_message = head_message + "  #-----------------------------------------------------------------#\n"
    head_message = head_message + "  #                                                                 #\n"

    if LOGO=="A":
     head_message= head_message + "  #     .-')             _ .-') _   .-') _                _  .-')   #\n"
     head_message= head_message + "  #   .(  OO)           ( (  OO) ) (  OO) )              ( \( -O )  #\n"
     head_message= head_message + "  #  (_)---\_)   .-----. \     .'_ /     '._  .-'),-----. ,------.  #\n"
     head_message= head_message + "  #  '  .-.  '  / ,-.   \,`'--..._)|'--...__)( OO'  .-.  '|   /`. ' #\n"
     head_message= head_message + "  # ,|  | |  |  '-'  |  ||  |  \  ''--.  .--'/   |  | |  ||  /  | | #\n"
     head_message= head_message + "  #(_|  | |  |     .'  / |  |   ' |   |  |   \_) |  |\|  ||  |_.' | #\n"
     head_message= head_message + "  #  |  | |  |   .'  /__ |  |   / :   |  |     \ |  | |  ||  .  '.' #\n"
     head_message= head_message + "  #  '  '-'  '-.|       ||  '--'  /   |  |      `'  '-'  '|  |\  \  #\n"
     head_message= head_message + "  #   `-----'--'`-------'`-------'    `--'        `-----' `--' '--' #\n"

    if LOGO=="B":
     head_message= head_message + '  #     .d88888b.   .d8888b.  8888888b. 88888888888                 #\n'
     head_message= head_message + '  #    d88P" "Y88b d88P  Y88b 888  "Y88b    888                     #\n'
     head_message= head_message + '  #    888     888        888 888    888    888                     #\n'
     head_message= head_message + '  #    888     888      .d88P 888    888    888   .d88b.  888d888   #\n'
     head_message= head_message + '  #    888     888  .od888P"  888    888    888  d88""88b 888P"     #\n'
     head_message= head_message + '  #    888 Y8b 888 d88P"      888    888    888  888  888 888       #\n'
     head_message= head_message + '  #    Y88b.Y8b88P 888"       888  .d88P    888  Y88..88P 888       #\n'
     head_message= head_message + '  #     "Y888888"  888888888  8888888P"     888   "Y88P"  888       #\n'
     head_message= head_message + '  #           Y8b                                                   #\n'

    if LOGO=="C":
     head_message= head_message + "  #      $$$$$$\   $$$$$$\  $$$$$$$\ $$$$$$$$\                      #\n"
     head_message= head_message + "  #     $$  __$$\ $$  __$$\ $$  __$$\\\__$$  __|                     #\n"
     head_message= head_message + "  #     $$ /  $$ |\__/  $$ |$$ |  $$ |  $$ | $$$$$$\   $$$$$$\      #\n"
     head_message= head_message + "  #     $$ |  $$ | $$$$$$  |$$ |  $$ |  $$ |$$  __$$\ $$  __$$\     #\n"
     head_message= head_message + "  #     $$ |  $$ |$$  ____/ $$ |  $$ |  $$ |$$ /  $$ |$$ |  \__|    #\n"
     head_message= head_message + "  #     $$ $$\$$ |$$ |      $$ |  $$ |  $$ |$$ |  $$ |$$ |          #\n"
     head_message= head_message + "  #     \$$$$$$ / $$$$$$$$\ $$$$$$$  |  $$ |\$$$$$$  |$$ |          #\n"
     head_message= head_message + "  #      \___$$$\ \________|\_______/   \__| \______/ \__|          #\n"
     head_message= head_message + "  #          \___|                                                  #\n"

    if LOGO=="D":
     head_message=head_message +'''  #     .g8""8q.            `7MM"""Yb.  MMP""MM""YMM                #\n'''
     head_message=head_message +'''  #   .dP'    `YM.            MM    `Yb.P'   MM   `7                #\n'''
     head_message=head_message +'''  #   dM'      `MM  pd*"*b.   MM     `Mb     MM  ,pW"Wq.`7Mb,od8    #\n'''
     head_message=head_message +'''  #   MM        MM (O)   j8   MM      MM     MM 6W'   `Wb MM' "'    #\n'''
     head_message=head_message +'''  #   MM.      ,MP     ,;j9   MM     ,MP     MM 8M     M8 MM        #\n'''
     head_message=head_message +'''  #   `Mb.    ,dP'  ,-='      MM    ,dP'     MM YA.   ,A9 MM        #\n'''
     head_message=head_message +'''  #     `"bmmd"'   Ammmmmmm .JMMmmmdP'     .JMML.`Ybmd9'.JMML.      #\n'''
     head_message=head_message +'''  #         MMb                                                     #\n'''
     head_message=head_message +'''  #          `bood'                                                 #\n'''

    if LOGO=="E":
     head_message=head_message +'''  #        ___    ____     ____    _____   U  ___ u   ____          #\n'''
     head_message=head_message +'''  #       / " \  |___"\   |  _"\  |_ " _|   \/"_ \/U |  _"\ u       #\n'''
     head_message=head_message +'''  #      | |"| | U __) | /| | | |   | |     | | | | \| |_) |/       #\n'''
     head_message=head_message +'''  #     /| |_| |\\\/ __/ \U| |_| |\ /| |\.-,_| |_| |  |  _ <         #\n'''
     head_message=head_message +'''  #     U \__\_\u|_____|u |____/ uu |_|U \_)-\___/   |_| \_\        #\n'''
     head_message=head_message +'''  #        \\\//  <<  //    |||_   _// \\\_     \\\     //   \\\_       #\n'''
     head_message=head_message +'''  #       (_(__)(__)(__)  (__)_) (__) (__)   (__)   (__)  (__)      #\n'''

    if LOGO=="F":
     head_message=head_message +'''  #      ________   ________  ________ ___________                  #\n'''
     head_message=head_message +'''  #      \_____  \  \_____  \ \______ \\\__    ___/____ _______      #\n'''
     head_message=head_message +'''  #       /  / \  \  /  ____/  |    |  \ |    |  /  _ \\\_  __ \     #\n'''
     head_message=head_message +'''  #      /   \_/.  \/       \  |    `   \|    | (  <_> )|  | \/     #\n'''
     head_message=head_message +'''  #      \_____\ \_/\_______ \/_______  /|____|  \____/ |__|        #\n'''
     head_message=head_message +'''  #             \__>        \/        \/                            #\n'''

    if LOGO=="G":
     head_message=head_message +'''  #     _|_|        _|_|    _|_|_|    _|_|_|_|_|                    #\n'''
     head_message=head_message +'''  #   _|    _|    _|    _|  _|    _|      _|      _|_|    _|  _|_|  #\n'''
     head_message=head_message +'''  #   _|  _|_|        _|    _|    _|      _|    _|    _|  _|_|      #\n'''
     head_message=head_message +'''  #   _|    _|      _|      _|    _|      _|    _|    _|  _|        #\n'''
     head_message=head_message +'''  #     _|_|  _|  _|_|_|_|  _|_|_|        _|      _|_|    _|        #\n'''
                                                               
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #-----------------------------------------------------------------#\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #    Q2DTor v1.1 (2018-12)                                        #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #    A software for the calculation of torsional                  #\n"
    head_message = head_message + "  #    anharmonicity of two coupled rotors                          #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #    Main authors:                                                #\n"
    head_message = head_message + "  #         * Ferro-Costas, David       (1)                         #\n"
    head_message = head_message + "  #         * Fernandez-Ramos, Antonio  (1)                         #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #    In collaboration with:                                       #\n"
    head_message = head_message + "  #         * Cordeiro, M. Natalia D. S.(2)                         #\n"
    head_message = head_message + "  #         * Truhlar, Donald G.        (3)                         #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #    (1) Centro Singular de Investigacion en Quimica Bioloxica    #\n"
    head_message = head_message + "  #        e Materiais Moleculares (CIQUS),                         #\n"
    head_message = head_message + "  #        Universidade de Santiago de Compostela, Galicia, Spain   #\n"
    head_message = head_message + "  #    (2) Department of Chemistry and Biochemistry,                #\n"
    head_message = head_message + "  #        University of Porto, Portugal                            #\n"
    head_message = head_message + "  #    (3) Department of Chemistry and Supercomputer Institute,     #\n"
    head_message = head_message + "  #        University of Minnesota, Minneapolis, Minnesota          #\n"
    head_message = head_message + "  #                                                                 #\n"
    head_message = head_message + "  #-----------------------------------------------------------------#\n"
   #head_message = head_message + "  # Logo created with: http://patorjk.com/software/taag             #\n"
   #head_message = head_message + "  #-----------------------------------------------------------------#\n"
    head_message = head_message + "\n"
    return head_message
#----------------------------------------------------------#
def q2dtor_gethelpstring():
    hstring =           "\n"
    hstring = hstring + "  Program information:   \n"
    hstring = hstring + "     * 'python Q2DTor.py -h' (or --help   ) shows this help message and exits\n"
    hstring = hstring + "     * 'python Q2DTor.py -v' (or --version) shows the program's version and exits\n"
    hstring = hstring + "     \n"
    hstring = hstring + "  ---------------------------------------------------\n"
    hstring = hstring + "  | Basic usage: python Q2DTor.py system [argument] |\n"
    hstring = hstring + "  ---------------------------------------------------\n"
    hstring = hstring + "     \n"
    hstring = hstring + "     system                the input file name (without extension)\n"
    hstring = hstring + "     \n"
    hstring = hstring + "     Without argument      creates the input file if it does not exists\n"
    hstring = hstring + "                           example: 'python Q2DTor.py system' to create system.inp\n"
    hstring = hstring + "     \n"
    hstring = hstring + "     Available arguments:   \n"
    hstring = hstring + "      --init [software]    initialization step (requires input and xyz files)\n"
    hstring = hstring + "      --pes                calculates torsional 2D-PES\n"
    hstring = hstring + "      --fourier            fits the calculated torsional PES to Fourier series\n"
    hstring = hstring + "      --findsp             finds the stationary points in the Fourier series torsional PES   \n"
    hstring = hstring + "      --optsp  [spname]    optimizes the geometry of one/all stationary point \n"
    hstring = hstring + "      --tor2dns            calculates the 2D-NS eigenvalues \n"
    hstring = hstring + "      --rovibpf [thermo]   calculates the partition functions;  \n"
    hstring = hstring + "                           with 'thermo', also calculates the thermodynamic functions\n"
    hstring = hstring + "     \n"
    hstring = hstring + "     By default, no screen information is printed in the execution of Q2DTor.\n"
    hstring = hstring + "     To print the output information also in the screen \n"
    hstring = hstring + "     type --print as an additional argument\n"
    hstring = hstring + "     \n"
    hstring = hstring + "     Values for 'software' option:\n"
    hstring = hstring + "       * gaussian [DEFAULT]\n"
    hstring = hstring + "       * orca    \n"
    hstring = hstring + "     \n"
    hstring = hstring + "     ----------------- \n"
    hstring = hstring + "     | Diverse Tools | \n"
    hstring = hstring + "     ----------------- \n"
    hstring = hstring + "\n"
    hstring = hstring + "      (1) --pdf\n"
    hstring = hstring + "          This tool generates a pdf file with plots.\n"
    hstring = hstring + "          The information for the plots is collected from the Q2DTor output files.\n"
    hstring = hstring + "\n"
    hstring = hstring + "          Basic usage: python Q2DTor.py system --pdf\n"
    hstring = hstring + "\n"
    hstring = hstring + "      (2) --gts\n"
    hstring = hstring + "          This tool can be used to generate a gts file from the\n"
    hstring = hstring + "          output data of the following electronic structure softwares:\n"
    hstring = hstring + "              * gaussian\n"
    hstring = hstring + "              * orca    \n"
    hstring = hstring + "\n"
    hstring = hstring + "          Basic usage: python Q2DTor.py system --gts software\n"
    hstring = hstring + "\n"
    hstring = hstring + "          example for 'gaussian'\n"
    hstring = hstring + "          ----------------------\n"
    hstring = hstring + "          By executing:\n"
    hstring = hstring + "                python Q2DTor.py system --gts gaussian\n"
    hstring = hstring + "          the program extracts data from system.fchk\n"
    hstring = hstring + "\n"
    hstring = hstring + "          example for 'orca'\n"
    hstring = hstring + "          ----------------------\n"
    hstring = hstring + "          By executing\n"
    hstring = hstring + "                python Q2DTor.py system --gts orca\n"
    hstring = hstring + "          the program extracts data from:\n"
    hstring = hstring + "                system.out\n"
    hstring = hstring + "                system.engrad (if exists)\n"
    hstring = hstring + "                system.hess   (if exists)\n"
    hstring = hstring + "\n"
    hstring = hstring + "      (3) --icoords\n"
    hstring = hstring + "          This tool can be used to check if a given set of internal\n"
    hstring = hstring + "          coordinates is appropiate for the calculation of vibrational frequencies.\n"
    hstring = hstring + "\n"
    hstring = hstring + "          Basic usage: python Q2DTor.py system --icoords mode icsfile\n"
    hstring = hstring + "\n"
    hstring = hstring + "          mode = sp\n"
    hstring = hstring + "          ----------\n"
    hstring = hstring + "          Command: python Q2DTor.py system --icoords sp\n"
    hstring = hstring + "\n"
    hstring = hstring + "          Use this mode to check the internal coordinates of all the\n"
    hstring = hstring + "          stationary points listed in the Q2DTor_files_system/system.splist file\n"
    hstring = hstring + "\n"
    hstring = hstring + "          If icsfile is not defined, the internal coordinates will be\n"
    hstring = hstring + "          read from file Q2DTor_files_system/system.ics\n"
    hstring = hstring + "\n"
    hstring = hstring + "          To use the internal coordinates stored in file 'file.ics'\n"
    hstring = hstring + "          just execute:\n"
    hstring = hstring + "                 python Q2DTor.py system --icoords sp file.ics\n"
    hstring = hstring + "\n"
    hstring = hstring + "          mode = gts\n"
    hstring = hstring + "          ----------\n"
    hstring = hstring + "          Command: python Q2DTor.py system --icoords gts file.ics\n"
    hstring = hstring + "\n"
    hstring = hstring + "          Use this option to check the internal coordinates of the 'file.ics'\n"
    hstring = hstring + "          file for the stationary point defined by the 'system.gts' file\n"
    hstring = hstring + "\n"
    return hstring
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 3) Q2DTor options                              #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def q2dtor_init(name,inpdata,files2read,files2write,software):

    #---------------#
    # Expand tuples #
    #---------------#
    (level,ch,mtp,torsion1,torsion2,ttype) = inpdata
    (f_xyz,) = files2read
    (f_ics,f_calcs) = files2write

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information from input file:"
    print  "        * torsion1      %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print  "        * torsion2      %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print  "        * level         %s"%level
    print  "        * charge        %s"%ch
    print  "        * multiplicity  %s"%mtp
    print 

    #------------#
    # Read files #
    #------------#
    which = ["xyz"]
    dataifiles = readfiles(files2read,which)
    # Expand data
    structure,symbols,masslist = dataifiles["xyz"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print

    #======><======#

    #----------------------#
    # Internal coordinates #
    #----------------------#
    print "     Generation of (redundant) internal coordinates:"
    if os.path.exists(f_ics):
        print "        - File '%s' already exists"%f_ics
        ricoords, nics, bonds = readfile_ics(f_ics)
    else:
       structure.graph_autoconnect(lengthfactor=1.3)
       structure.graph_fragconnect()
       nbonds   = structure.graph_nbonds()
       poincare = 1 - structure.get("natoms") + nbonds
       print "          nbonds + 1 - natoms = %i"%poincare
       # Generate ricoords
       ricoords, nics = structure.gen_ricoords(torsions=[torsion1,torsion2],check=False)
       # Write ricoords
       bonds    = structure.graph_getbonds()
       writefile_ics(ricoords,f_ics,bonds)
       print "        - Internal coordinates stored in '%s' "%f_ics

    print "        - List of internal coordinates (%i)"%nics
    for idx in range(0,len(ricoords),5):
        line = ""
        for xx in range(5):
            pos = idx + xx
            if pos > len(ricoords)-1: continue
            tt, ic = ricoords[pos]
            if tt=="3": ic = "=".join("%i"%(a+1) for a in ic)
            else      : ic = "-".join("%i"%(a+1) for a in ic)
            line = line + "  %11s  "%ic
        print "     %s"%line
    print

    #-----------------------------#
    # Write reference input files #
    #-----------------------------#
    print "     File generation for electronic structure calculations:"
    software = software.lower()
    print "        - Selected software is '%s'"%software
    if os.path.exists(f_calcs):
       print "        - File '%s' already exists"%f_calcs
       ilines, software_in_file = readfile_calcs(f_calcs,"start_scangeom","end_scangeom")
       if software_in_file.lower() != software.lower():
          print "        - ERROR: software in %s is '%s'"%(f_calcs,software_in_file)
          print
          sys.exit(EXITMESS)
    else:
       import_mes(software,"defaults")
       string = mes_defaults(level,ch,mtp,torsion1,torsion2,ttype)
       ifile  = open(f_calcs,'w')
       ifile.write(string)
       ifile.close()
       print "        - File '%s' has been created"%f_calcs
    print
#----------------------------------------------------------#
def q2dtor_pes(name,inpdata,files2read,files2write,f_ics=None):
    global mes_pespoint

    #---------------#
    # Expand tuples #
    #---------------#
    (torsion1,torsion2,tsigma1,tsigma2,t1step,t2step,symmetry) = inpdata
    (f_xyz,f_calcs) = files2read
    (f_pes,) = files2write

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print "     Information from the input file:"
    print "        * torsion1      %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print "        * torsion2      %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print "        * tsigma1       %i"%tsigma1
    print "        * tsigma2       %i"%tsigma2
    print "        * steptorsion1  %.2f degrees"%(t1step*cons.R2D)
    print "        * steptorsion2  %.2f degrees"%(t2step*cons.R2D)
    print "        * symmetry      %s"%symmetry
    print
    torsional_bonds = ( (torsion1[1],torsion1[2]) , (torsion2[1],torsion2[2]) )

    #------------#
    # Read files #
    #------------#
    which = ["xyz","calcs_scangeom"]
    dataifiles = readfiles(files2read,which)
    # Expand data
    structure,symbols,masslist = dataifiles["xyz"]
    ilines, software = dataifiles["calcs_scangeom"]
    # Software
    import_mes(software,"pespoint")
    # Connectivity
    connect = True
    if f_ics is not None:
       if os.path.isfile(f_ics):
          print "      File '%s' exists: using its connectivity"%f_ics
          print
          ics, nics, all_bonds = readfile_ics(f_ics)
          connect = False
          for (idx1,idx2) in all_bonds: structure.graph_addbond(idx1,idx2)
    if connect:
       structure.graph_autoconnect(lengthfactor=1.3)
       structure.graph_fragconnect(maxfrags=1)
       all_bonds = structure.graph_getbonds()

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print

    #======><======#


    #--------------------------------#
    # File from previous calculation #
    #--------------------------------#
    if not os.path.exists(f_pes):
       print "      Creating a file with the PES energies: '%s'"%f_pes
       dict_pes   = {}
       calculated = []
       MINENERGY  = float("inf")
       pes_file   = open(f_pes,'w')
       rewrite    = False
    else:
       print "      File '%s' already exists:"%f_pes
       dict_pes, dummy, MINENERGY = readfile_pes(f_pes,exclude=False)
       print "        * %i PES points stored"%len(dict_pes.keys())
       print "        * Calculated points will be omitted"
       calculated = dict_pes.keys()
       pes_file   = open(f_pes,'a')
       rewrite    = True
    print

    #-------------------#
    # Go point by point #
    #-------------------#
    pespoints = pessym_points(t1step,t2step,symmetry,tsigma1,tsigma2)
    # v1.1 - get angles for input geometry and reorder points
    if BOOL_v1_1:
       iphi1 = hf.ranged_angle( structure.icvalue(torsion1) )
       iphi2 = hf.ranged_angle( structure.icvalue(torsion2) )
       minidx, mind2 = 0, float("inf")
       for idx,(phi1,phi2) in enumerate(pespoints):
           d2 = (phi1-iphi1)**2 + (phi2-iphi2)**2
           if d2 < mind2: minidx,mind2 = idx, d2
       pespoints = pespoints.tolist()
       pespoints = pespoints[minidx:]+pespoints[minidx::-1]
    # start calculations
    print "      Calculating PES (%i points)"%len(pespoints)
    MOs_previous = None
    pnames = []
    for phi1, phi2 in pespoints:
        point_name = string_getpname(phi1,phi2,molname=name,units="rad")
        pnames.append(point_name)
        if point_name in calculated:
           # v1.1 - read geometry in order to have a better guess
           if BOOL_v1_1:
              atonums = structure.get("atonums")
              Etot,xcc = dict_pes[point_name][2:4]
              structure = basic2Struct(point_name,xcc,atonums,0,1,Etot,None,None)
              for (idx1,idx2) in all_bonds: structure.graph_addbond(idx1,idx2)
           #print "         %s: omitting calculation..."%point_name
           continue
        print "         %s: calculating..."%point_name
        #-------------------------#
        # Generate geometry guess #
        #-------------------------#
        iphi1 = hf.ranged_angle( structure.icvalue(torsion1) )
        iphi2 = hf.ranged_angle( structure.icvalue(torsion2) )
        rotation1  = hf.angle_diff(iphi1,phi1)
        rotation2  = hf.angle_diff(iphi2,phi2)
        rot_angles = (rotation1,rotation2)
        xvec       = structure.graph_irotation( torsional_bonds , rot_angles )
        #-----------------#
        # Calculate point #
        #-----------------#
        idata = (ilines,list(xvec),symbols,point_name,MOs_previous)
        odata = mes_pespoint(*idata)
        if odata is None:
           pes_file.write(" 1\n")
           pes_file.write(" FAILED    0.0000000   %s\n"%point_name)
           pes_file.write(" XX       0.000000       0.000000       0.000000\n")
           dict_pes[point_name] = [phi1,phi2,None,None]
           continue
        xcc, atonums, ch, mtp,  Etot, gcc, Fcc = odata
        MOs_previous = point_name
        # v1.1 - add to calculated
        if BOOL_v1_1: calculated.append(point_name)
        #----------------------------------#
        # Generate structure and save data #
        #----------------------------------#
        structure = basic2Struct(point_name,xcc,atonums,ch,mtp,Etot,gcc,Fcc)
        for (idx1,idx2) in all_bonds: structure.graph_addbond(idx1,idx2)
        # Get data
        xvec   = structure.get("xcc")
        energy = structure.get("Etot")
        if energy < MINENERGY: MINENERGY = energy
        calc_phi1 = hf.ranged_angle( structure.icvalue(torsion1),v=1 )
        calc_phi2 = hf.ranged_angle( structure.icvalue(torsion2),v=1 )
        diff1 = abs(hf.angle_diff(calc_phi1,phi1))*cons.R2D
        diff2 = abs(hf.angle_diff(calc_phi2,phi2))*cons.R2D
        if diff1 > 0.2 or diff2 > 0.2:
           print "  --> Final dihedrals in %s are wrong..."%point_name
           print "      theta_1 is %5.1f and should be %5.1f"%(calc_phi1*cons.R2D,phi1*cons.R2D)
           print "      theta_2 is %5.1f and should be %5.1f"%(calc_phi2*cons.R2D,phi2*cons.R2D)
           sys.exit(EXITMESS)
        # Write in pes file
        pes_file.write(" %i\n"%len(symbols))
        pes_file.write(" Geometry  %+14.8f   %7.3f  %7.3f  %s  YES\n"%(energy, phi1*cons.R2D, phi2*cons.R2D, point_name))
        for idx in range(len(symbols)):
            symbol = symbols[idx]
            x,y,z = xvec[3*idx:3*idx+3] * cons.angstrom
            pes_file.write(" %-2s   %+13.6f  %+13.6f  %+13.6f\n"%(symbol,x,y,z))
        # Save in dict
        dict_pes[point_name] = [phi1,phi2,energy,xvec]
    print

    print "      Minimum energy (EMIN) in calculated points is %+14.8f hartree"%MINENERGY
    print
    #pes_file.write("MINENERGY         %+14.8f\n"%MINENERGY)
    pes_file.close()

    print "      phi1: torsion defined by atoms %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print "      phi2: torsion defined by atoms %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print
    print "      PES Summary Table:"
    print "                  ----------------------------"
    print "                    phi1  |  phi2  | E - EMIN "
    print "                  ----------------------------"
    for point_name in sorted(dict_pes.keys()):
        phi1,phi2,energy,xvec = dict_pes[point_name]
        if energy is None:
           linedata = (phi1*cons.R2D,phi2*cons.R2D,"   --   ")
           print "                   %6.2f | %6.2f | %8s "%linedata
        else:
           linedata = (phi1*cons.R2D,phi2*cons.R2D,(energy-MINENERGY)*cons.h2c)
           print "                   %6.2f | %6.2f | %8.2f "%linedata
    print "                  ----------------------------"
    print "                     phi1, phi2 in degrees    "
    print "                     E-EMIN in cm^-1          "
    print

    # Rewrite all data
    if rewrite:
        fails = 0
        print "      Rewriting '%s' file"%f_pes
        pes_file = open(f_pes,'w')
        for point_name in sorted(dict_pes.keys()):
            phi1,phi2,energy,xvec = dict_pes[point_name]
            if energy is None:
               fails += 1
               pes_file.write(" 1\n")
               pes_file.write(" FAILED    0.0000000   %s\n"%point_name)
               pes_file.write(" XX       0.000000       0.000000       0.000000\n")
            else:
               phi1,phi2,energy,xvec = dict_pes[point_name]
               pes_file.write(" %i\n"%len(symbols))
               if point_name in pnames:
                  pes_file.write(" Geometry  %+14.8f   %7.3f  %7.3f  %s  YES\n"%(energy, phi1*cons.R2D, phi2*cons.R2D, point_name))
               else:
                  pes_file.write(" Geometry  %+14.8f   %7.3f  %7.3f  %s  NO \n"%(energy, phi1*cons.R2D, phi2*cons.R2D, point_name))
               for idx in range(len(symbols)):
                   symbol = symbols[idx]
                   x,y,z = xvec[3*idx:3*idx+3] * cons.angstrom
                   pes_file.write(" %-2s   %+13.6f  %+13.6f  %+13.6f\n"%(symbol,x,y,z))
        pes_file.close()
        print "        * Number of 'FAIL' points: %i"%fails
        print
#----------------------------------------------------------#
def q2dtor_fourier(name,inpdata,files2read,files2write,TOTPESPOINTS):
    #------------------------#
    def string_fitting(fterms, popt, ignore, errors, fitting_time=None):
        nibs = 10
        rsquare, aae, aae_small = errors
        data_string = ""
        # Print data
        data_string = data_string + " "*nibs+"---------------------------------------------------------------\n"
        data_string = data_string + " "*nibs+"|                           Fitting                           |\n"
        data_string = data_string + " "*nibs+"---------------------------------------------------------------\n"
        for term, parameter in zip(fterms,popt):
            term_type, idx1, idx2 = term
            if   term_type == "const" :
               line = " "*(nibs+2) + "V(phi1,phi2) = %+11.4f                               +"%parameter
               data_string = data_string + line + "\n"
               continue
            # String of new line
            elif term_type == "cos"   :
                 if idx1 == "-"       : line = "%+11.4f * cos(%02i*phi2)                +"%(parameter,idx2)
                 if idx2 == "-"       : line = "%+11.4f * cos(%02i*phi1)                +"%(parameter,idx1)
            elif term_type == "sin"   :
                 if idx1 == "-"       : line = "%+11.4f * sin(%02i*phi2)                +"%(parameter,idx2)
                 if idx2 == "-"       : line = "%+11.4f * sin(%02i*phi1)                +"%(parameter,idx1)
            elif term_type == "coscos": line = "%+11.4f * cos(%02i*phi1) * cos(%02i*phi2) +"%(parameter,idx1,idx2)
            elif term_type == "cossin": line = "%+11.4f * cos(%02i*phi1) * sin(%02i*phi2) +"%(parameter,idx1,idx2)
            elif term_type == "sinsin": line = "%+11.4f * sin(%02i*phi1) * sin(%02i*phi2) +"%(parameter,idx1,idx2)
            elif term_type == "sincos": line = "%+11.4f * sin(%02i*phi1) * cos(%02i*phi2) +"%(parameter,idx1,idx2)
            else: print "Unknown term!"; sys.exit(EXITMESS)
            # Incorporate line
            data_string = data_string + " "*(nibs+17) + line
            if abs(parameter) < ignore: data_string = data_string + " [<--]"
            data_string = data_string + "\n"
        data_string = data_string[:-2]+"\n"
        data_string = data_string + " "*nibs+"---------------------------------------------------------------\n"
        data_string = data_string + " "*nibs+"  number of parameters is %i\n"%len(popt)
        data_string = data_string + " "*nibs+"  fitting correlation: (1.0 - r^2) = %.1e\n"%(1.0-rsquare)
        data_string = data_string + " "*nibs+"  average abs. errors:\n"
        data_string = data_string + " "*nibs+"       %.1e cm^-1\n"%(aae)
        data_string = data_string + " "*nibs+"       %.1e cm^-1 (for points below mean value)\n"%(aae_small)
        data_string = data_string + " "*nibs+"  elapsed time: %.1f seconds\n"%fitting_time
        data_string = data_string + " "*nibs+"---------------------------------------------------------------\n"
        return data_string
    #------------------------#

    #---------------#
    # Expand tuples #
    #---------------#
    (symmetry,tsigma1,tsigma2,tuple_fit) = inpdata
    fterms, coef_nums, weight, ignore    = tuple_fit
    (f_pes,) = files2read
    (f_vfit,) = files2write

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information read by the input file:"
    print  "        * symmetry      %s"%symmetry
    print  "        * tsigma1       %i"%tsigma1
    print  "        * tsigma2       %i"%tsigma2
    print  "        * weight        %.4f"%weight
    print  "        * ignore        %.4f"%ignore
    print

    #------------#
    # Read files #
    #------------#
    which = ["pes"]
    dataifiles = readfiles(files2read,which)
    # Expand data
    dict_pes, symbols, Eref = dataifiles["pes"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print

    #======><======#


    #---------------------------#
    # Expand PES using symmetry #
    #---------------------------#
    print  "     Applying symmetry conditions"
    dict_pes = pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2,pprint=True)
    print "                  ----------------------------"
    print "                    phi1  |  phi2  | E - EMIN "
    print "                  ----------------------------"
    for point_name in sorted(dict_pes.keys()):
        phi1,phi2,energy,xvec = dict_pes[point_name]
        linedata = (phi1*cons.R2D,phi2*cons.R2D,(energy-Eref)*cons.h2c)
        print "                   %6.2f | %6.2f | %8.2f "%linedata
    print "                  ----------------------------"
    print "                     phi1, phi2 in degrees    "
    print "                     E-EMIN in cm^-1          "
    print

    points   = dict_pes.keys()
    xdata    = np.array([tuple(dict_pes[point][0:2]) for point in points])
    ydata    = np.array([ (dict_pes[point][2]-Eref)*cons.h2c for point in points])

    # Limit values in ydata
    y_min  = min(ydata)
    y_mean = sum(ydata) / len(ydata)
    y_max  = max(ydata)
    print "        * Minimum energy in potential is: %+9.2f cm^-1"%y_min
    print "        * Average energy in potential is: %+9.2f cm^-1"%y_mean
    print "        * Maximum energy in potential is: %+9.2f cm^-1"%y_max
    if len(ydata) != TOTPESPOINTS:
       print "        * WARNING: the PES contains %i points instead of %i..."%(len(ydata),TOTPESPOINTS)
    print

    #------------------------#
    # Fourier series Fitting #
    #------------------------#
    # File with fitting already exists
    print 
    print  "     ===================================="
    print  "      Performing a fit to Fourier series "
    print  "     ===================================="
    print
    guess = [0.0]*len(fterms)
    if os.path.isfile(f_vfit):
       print "        * File '%s' exists. Using its data as guess..."%f_vfit
       print
       fterms2, parameters2, dummy = readfile_xfit(f_vfit)
       for idxA in range(len(fterms)):
           term = fterms[idxA]
           if term in fterms2:
              idxB = fterms2.index(term)
              guess[idxA] = parameters2[idxB]
    else:
      guess = None

    while True:
          # Print elements
          print "        * Elements in Fourier series (%i):"%len(fterms)
          count_const  = 0
          count_cos1   = 0
          count_cos2   = 0
          count_sin1   = 0
          count_sin2   = 0
          count_coscos = 0
          count_sincos = 0
          count_cossin = 0
          count_sinsin = 0
          for ft,idx1,idx2 in fterms:
              if ft == "const"            : count_const  += 1
              if ft == "cos" and idx2=="-": count_cos1   += 1
              if ft == "cos" and idx1=="-": count_cos2   += 1
              if ft == "sin" and idx2=="-": count_sin1   += 1
              if ft == "sin" and idx1=="-": count_sin2   += 1
              if ft == "coscos"           : count_coscos += 1
              if ft == "sincos"           : count_sincos += 1
              if ft == "cossin"           : count_cossin += 1
              if ft == "sinsin"           : count_sinsin += 1
          print "          number of constants               : %i"%count_const
          print "          number of cos(i*phi1)             : %i"%count_cos1
          print "          number of cos(j*phi2)             : %i"%count_cos2
          print "          number of sin(i*phi1)             : %i"%count_sin1
          print "          number of sin(j*phi2)             : %i"%count_sin2
          print "          number of cos(i*phi1).cos(j*phi2) : %i"%count_coscos
          print "          number of sin(i*phi1).cos(j*phi2) : %i"%count_sincos
          print "          number of cos(i*phi1).sin(j*phi2) : %i"%count_cossin
          print "          number of sin(i*phi1).sin(j*phi2) : %i"%count_sinsin
          # Create object and perform fitting
          FourierObject = Fourier2D(fterms)
          fterms2, fcoefs, fit_errors, fitting_time = FourierObject.fit(xdata,ydata,weight=weight,guess=guess)
          # Print result of fitting
          print string_fitting(fterms2, fcoefs, ignore, fit_errors, fitting_time)
          # Remove small values
          finish = True
          fterms, guess = [], []
          for ff, vv in zip(fterms2,fcoefs):
              if abs(vv) < ignore:
                finish = False
                continue
              fterms.append(ff)
              guess.append(vv)
          # Exit fitting
          if finish: break

    # write fitting parameters first
    vfit = open(f_vfit,'w')
    vfit.write("#----------------------------------------------#\n")
    vfit.write("# fitting correlation: (1.0 - r^2) = %7.1e   #\n"%(1.0-fit_errors[0]))
    vfit.write("# average abs. errors:                         #\n")
    vfit.write("#    %7.1e                                   #\n"%fit_errors[1])
    vfit.write("#    %7.1e (for points below mean value)     #\n"%fit_errors[2])
    vfit.write("# elapsed time: %6.1f seconds                 #\n"%fitting_time)
    vfit.write("#----------------------------------------------#\n")
    vfit.write("reference %.7f \n"%Eref)
    vfit.close()
    # Now, write fitting
    print "     Storing the fitting parameters at: %s"%f_vfit
    writefile_xfit(fterms2,fcoefs,f_vfit,mode='a')
    print
#----------------------------------------------------------#
def q2dtor_findsp(name,inpdata,files2read,files2write):

    #---------------#
    # Expand tuples #
    #---------------#
    (symmetry,tsigma1,tsigma2,tolerance) = inpdata
    (f_vfit,) = files2read
    (f_splist,) = files2write

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information read by the input file:"
    print  "        * symmetry      %s"%symmetry
    print  "        * tsigma1       %i"%tsigma1
    print  "        * tsigma2       %i"%tsigma2
    print  "        * tolerance     %.4f degrees"%(tolerance*cons.R2D)
    print

    #------------#
    # Read files #
    #------------#
    which = ["vfit"]
    dataifiles = readfiles(files2read,which)
    # Expand data
    fterms, parameters, Eref_vfit = dataifiles["vfit"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print

    #======><======#



    #-------------------------------#
    # Looking for stationary points #
    #-------------------------------#
    searchway, downhill = "A", True
   #searchway, downhill = "B", False
    print "     Searching for stationary points in the fitted 2D-PES"

    Fourier_PES = Fourier2D(fterms)
    Fourier_PES.set_coefs(parameters)

    # step in degrees
    phi1_list = np.arange(-3*tolerance , 2*np.pi+3*tolerance, tolerance)
    phi2_list = np.arange(-3*tolerance , 2*np.pi+3*tolerance, tolerance)

    npoints = len(phi1_list)

    if searchway == "A":
       print "        * Generating a grid with PES values (step of %.2f degrees)..."%(tolerance*cons.R2D)
       pes_matrix = [[ Fourier_PES.value(phi1,phi2) for phi2 in phi2_list] for phi1 in phi1_list]
       # get difference in x- and y- direction
       print "        * Generating gradient grid (x-component)..."
       grad_x = np.diff(pes_matrix,n=1,axis=0)
       print "        * Generating gradient grid (y-component)..."
       grad_y = np.diff(pes_matrix,n=1,axis=1)

    if searchway == "B":
       # get grid with gradient
       print "        * Generating gradient grid (x-component)..."
       grad_x = [[ Fourier_PES.derphi1(phi1, phi2) for phi2 in phi2_list] for phi1 in phi1_list]
       print "        * Generating gradient grid (y-component)..."
       grad_y = [[ Fourier_PES.derphi2(phi1, phi2) for phi2 in phi2_list] for phi1 in phi1_list]

    cpoints = []
    print "        * Looking for stationary points in the grid..."
    DX = 0.1*cons.D2R
    DY = 0.1*cons.D2R
    for i in range(2,len(phi1_list)-2):
        for j in range(2,len(phi2_list)-2):
            # check when the difference changes its sign
            gx_a, gy_a = grad_x[i-1][j], grad_y[i][j-1]
            gx_b, gy_b = grad_x[i  ][j], grad_y[i][j  ]
            cond1 = ((gx_a<0) != (gx_b<0)) or (gx_b==0.0)
            cond2 = ((gy_a<0) != (gy_b<0)) or (gy_b==0.0)
           #gx_c, gy_c = grad_x[i+1][j], grad_y[i][j+1]
           #cond1 = ((gx_a<0) != (gx_c<0)) or (gx_b==0.0)
           #cond2 = ((gy_a<0) != (gy_c<0)) or (gy_b==0.0)

            #---------------------------#
            # Classify stationary points#
            #---------------------------#
            if (cond1 and cond2):
               phi1 = phi1_list[i]
               phi2 = phi2_list[j]
               # Exclude those outside [0,2pi]
               if phi1*cons.R2D <   0.0: phi1 = 0.0
               if phi1*cons.R2D > 360.0: continue
               if phi2*cons.R2D <   0.0: phi2 = 0.0
               if phi2*cons.R2D > 360.0: continue
              #phi1 = hf.ranged_angle(phi1,v=1)
              #phi2 = hf.ranged_angle(phi2,v=1)
               # Calculate numeric hessian
               #
               # phi2
               # | 
               # | *  *  *  *  *  --> 40 41 42 43 44
               # | *  *  *  *  *  --> 30 31 32 33 34
               # | *  *  x  *  *  --> 20 21 22 23 24
               # | *  *  *  *  *  --> 10 11 12 13 14
               # | *  *  *  *  *  --> 00 01 02 03 04
               # ------------ phi1
               #
               if searchway == "A":
                  V22 = pes_matrix[i  ][j  ]
                  V20 = pes_matrix[i-2][j  ]
                  V24 = pes_matrix[i+2][j  ]
                  V02 = pes_matrix[i  ][j-2]
                  V42 = pes_matrix[i  ][j+2]
                  V33 = pes_matrix[i+1][j+1]
                  V11 = pes_matrix[i-1][j-1]
                  V13 = pes_matrix[i+1][j-1]
                  V31 = pes_matrix[i-1][j+1]
                  dxy = (V33+V11-V13-V31)/(2*tolerance*tolerance)
                  dxx = (V24+V20-2.0*V22)/(2*tolerance*tolerance)
                  dyy = (V42+V02-2.0*V22)/(2*tolerance*tolerance)

               if searchway == "B":
                  V22 = Fourier_PES.value(phi1     ,phi2     )
                  V20 = Fourier_PES.value(phi1-2*DX,phi2     )
                  V24 = Fourier_PES.value(phi1+2*DX,phi2     )
                  V02 = Fourier_PES.value(phi1     ,phi2-2*DY)
                  V42 = Fourier_PES.value(phi1     ,phi2+2*DY)
                  V33 = Fourier_PES.value(phi1+  DX,phi2+  DY)
                  V11 = Fourier_PES.value(phi1-  DX,phi2-  DY)
                  V13 = Fourier_PES.value(phi1+  DX,phi2-  DY)
                  V31 = Fourier_PES.value(phi1-  DX,phi2+  DY)
                  dxy = (V33+V11-V13-V31)/(2*DX*DY)
                  dxx = (V24+V20-2.0*V22)/(2*DX*DX)
                  dyy = (V42+V02-2.0*V22)/(2*DY*DY)

               hessian = np.matrix([[dxx,dxy],[dxy,dyy]])
               evalues, evectors = np.linalg.eigh(hessian)
              #print phi1*cons.R2D, phi2*cons.R2D, evalues
               if   evalues[0] > 0.0 and evalues[1] > 0.0: cpoints.append( (phi1,phi2,0,V22) )
               elif evalues[0] < 0.0 and evalues[1] < 0.0: cpoints.append( (phi1,phi2,2,V22) )
               else                                      : cpoints.append( (phi1,phi2,1,V22) )

    #------------------#
    # Print candidates #
    #------------------#
    print "        * A total of %i candidates were found"%len(cpoints)
    for phi1,phi2,cptype,V in cpoints:
        if cptype == 0: xxx = "minimum"
        if cptype == 1: xxx = "saddle "
        if cptype == 2: xxx = "maximum"
        print "            (%6.2f,%6.2f) ==> %s"%(phi1*cons.R2D,phi2*cons.R2D,xxx)
    print

    #------------------------#
    # Improve candidate list #
    #------------------------#
    print "     Improving list of candidates:"

    if downhill:
       print "        * Improving position of minima and maxima with downhill algorithm..."
       new_cpoints = []
       for phi1,phi2,cptype,V in cpoints:
           if cptype == 1:
              new_cpoints += [(phi1,phi2,cptype,V)]
              continue
           # (a) Minima
           if cptype == 0:
              nphi1, nphi2 = fmin(lambda x: Fourier_PES.value(x[0],x[1]),np.array([phi1,phi2]),disp=False)
           # (b) Maxima
           elif cptype == 2:
              nphi1, nphi2 = fmin(lambda x:-Fourier_PES.value(x[0],x[1]),np.array([phi1,phi2]),disp=False)
           # Save point
           V = Fourier_PES.value(nphi1,nphi2)
           new_cpoints += [(nphi1,nphi2,cptype,V)]
           phi1, nphi1 = phi1*cons.R2D, nphi1*cons.R2D
           phi2, nphi2 = phi2*cons.R2D, nphi2*cons.R2D
           print "            (%6.2f,%6.2f) --> (%6.2f,%6.2f)"%(phi1,phi2,nphi1,nphi2)
       cpoints = new_cpoints

    print "        * Applying symmetry conditions to each point..."
    new_cpoints = []
    for phi1,phi2,cptype,V in cpoints:
        equivalents = pessym_apply(phi1,phi2,symmetry,tsigma1,tsigma2)
        for p1,p2 in equivalents: new_cpoints += [ (p1,p2,cptype,V) ]
    cpoints = new_cpoints

    print "        * Removing stationary points outside the symmetry region..."
    dict_cps = {}
    for phi1,phi2,cptype,V in cpoints:
        cpname = string_getpname(phi1,phi2,"","rad")
        keep   = pessym_inregion(phi1, phi2, symmetry, tsigma1, tsigma2)
        if keep: dict_cps[cpname] = [phi1,phi2,cptype,V]

    print "        * Merging closed points..."
    dict_cps = clean_close(dict_cps,diff_angle=15.0,diff_Ecm=10.0)
    cpoints  = dict_cps.values()
    cpoints.sort(key=lambda x: x[3]) # sort by energy

    repeated = []
    if "b" in symmetry:
       print "        * Comparing energies of remaining points..."
       for idxA in range(len(cpoints)):
           phi1A, phi2A, cptypeA, VA = cpoints[idxA]
           for idxB in range(len(cpoints)):
            if idxA == idxB: continue
            if (idxA,idxB) in repeated: continue
            if (idxB,idxA) in repeated: continue
            phi1B, phi2B, cptypeB, VB = cpoints[idxB]
            if cptypeA != cptypeB: continue
            if abs(VA-VB) > 1.0  : continue
            tp1 = 2.0*np.pi/tsigma1
            tp2 = 2.0*np.pi/tsigma2
            cd1 = abs(phi1A + phi1B - tp1) < 2.0 * tolerance
            cd2 = abs(phi2A + phi2B - tp2) < 2.0 * tolerance
            if cd1 or cd2:
               distA = phi1A**2 + phi2A**2
               distB = phi1B**2 + phi2B**2
               if distA > distB: repeated.append( (idxB,idxA) )
               else            : repeated.append( (idxA,idxB) )
               phi1A,phi2A = phi1A*cons.R2D,phi2A*cons.R2D
               phi1B,phi2B = phi1B*cons.R2D,phi2B*cons.R2D
               print "            (%6.2f,%6.2f) <-> (%6.2f,%6.2f)"%(phi1A,phi2A,phi1B,phi2B)
       repeated = [idxB for idxA,idxB in repeated]
       print
          

    # Write list of points
    print "     Storing the location of each stationary point at: %s"%f_splist
    statpoints = open(f_splist,'w')
    statpoints.write("#                                                                   \n")
    statpoints.write("#  List of stationary points (SP) in the 2D-PES                     \n")
    statpoints.write("# ---------------------------------------------                     \n")
    statpoints.write("#         * Type = 0 for minima                                     \n")
    statpoints.write("#         * Type = 1 for saddle points                              \n")
    statpoints.write("#         * Type = 2 for maxima                                     \n")
    statpoints.write("#         * Phi1 and Phi2 in [degrees]                              \n")
    statpoints.write("#         * Energy in [1/cm] with regard to global the minimum      \n")
    statpoints.write("#                                                                   \n")
    statpoints.write("#  Type     Phi1     Phi2        Energy          opt OK?            SP name      \n")
    statpoints.write("# ------ | ------ | ------ | ---------------- | ---------- | --------------------\n")
    ncps = 0
    for idx in range(len(cpoints)):
        phi1, phi2, cptype, E = cpoints[idx]
        cp_name    = string_getpname(phi1,phi2,name,"rad")
        if idx in repeated: continue
        lineformat = "    %i      %6.2f   %6.2f     %+13.2f        NO        %18s  \n"
        statpoints.write( lineformat%(cptype,phi1*cons.R2D,phi2*cons.R2D,E,cp_name) )
        ncps += 1
    statpoints.write("\n")
    statpoints.write("# Number of stationary points: %i\n"%ncps)
    statpoints.close()
    print "                -------------------------------"
    print "                  phi1  |  phi2  |  E (cm^-1)  "
    print "                -------------------------------"
    num0 = 1
    num1 = 1
    num2 = 1
    for idx in range(len(cpoints)):
        if idx in repeated: continue
        phi1, phi2, cptype, E = cpoints[idx]
        if cptype == 0: xxx = "MIN%02i"%num0; num0 +=1
        if cptype == 1: xxx = "TS_%02i"%num1; num1 +=1
        if cptype == 2: xxx = "MAX%02i"%num2; num2 +=1
        line = " %7s : %6.2f | %6.2f | %+9.2f "%(xxx,phi1*cons.R2D,phi2*cons.R2D,E)
        print "      "+line
    print "                -------------------------------"
    print
#----------------------------------------------------------#
def q2dtor_optsp(name,inpdata,files2read,files2write,Tlist,mode="all"):
    global mes_optCP

    #---------------#
    # Expand tuples #
    #---------------#
    (torsion1,torsion2,symmetry,tsigma1,tsigma2,tolerance,ttype,freqscal) = inpdata
    (f_calcs, f_xyz, f_pes, f_vfit, f_splist) = files2read
    (f_splist,f_ics,f_spinfo) = files2write

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information read by the input file:"
    print  "        * ttype         %s"%ttype
    print  "        * torsion1      %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print  "        * torsion2      %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print  "        * symmetry      %s"%symmetry
    print  "        * tsigma1       %i"%tsigma1
    print  "        * tsigma2       %i"%tsigma2
    print  "        * tolerance     %.4f"%(tolerance*cons.R2D)
    print  "        * freqscal      %.3f"%freqscal
    print

    #------------#
    # Read files #
    #------------#
    which = ["calcs_sp","xyz","pes","vfit","splist"]
    dataifiles = readfiles(files2read,which)
    # Expand data
    [ilines0,software0,ilines1,software1,ilines2,software2] = dataifiles["calcs_sp"]
    [xyzStruct,symbols,masslist] = dataifiles["xyz"]
    [dict_pes, symbols, Eref] = dataifiles["pes"]
    [fterms, parameters, Eref_vfit] = dataifiles["vfit"]
    [dict_CPs, Eref_cp] = dataifiles["splist"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print
    if not os.path.exists(DIRSP):
       print  "     Creating folder: '%s'"%(DIRSP)
       os.mkdir(DIRSP)
       print


    #======><======#


    #----------------------------#
    # Deal with all or with one? #
    #----------------------------#
    list_CPs = sp_sorting(dict_CPs)
    if mode != "all":
       if mode not in list_CPs:
          print "     ERROR: '%s' is not defined in '%s' file"%(mode,f_splist)
          print
          sys.exit(EXITMESS)
       elif dict_CPs[mode][5] == "YES":
          message = "     '%s' has YES in the 'opt OK?' column of '%s' file"%(mode,f_splist)
          print message
          print
          sys.exit(message)
       elif dict_CPs[mode][5] != "NO":
          print "     ERROR: '%s' has %s in the 'opt OK?' column of '%s' file"%(mode,dict_CPs[mode][5],f_splist)
          print
          sys.exit(EXITMESS)
       print "     User has selected stationary point: %s"%(mode)
       list_CPs = [mode]
       # Remove gts file
       gtsfile = DIRSP + mode + ".gts"
       if os.path.exists(gtsfile):
          print "        * Removing file '%s'..."%gtsfile
          os.remove(gtsfile)
       print 

    #--------------#
    # Optimization #
    #--------------#
    dict_warns = {}
    min_sp  = float("inf")
    for cp_name in list_CPs:
        phi1, phi2, cptype, Ecm, gtsfile = dict_CPs[cp_name][0:5]
        dict_warns[cp_name] = []
        if cptype == 0: print  "     Optimization of %s (min)"%cp_name; ilines = ilines0; software = software0
        if cptype == 1: print  "     Optimization of %s (ts )"%cp_name; ilines = ilines1; software = software1
        if cptype == 2: print  "     Optimization of %s (max)"%cp_name; ilines = ilines2; software = software2
        # Software
        import_mes(software,"optCP")
        # Optimization
        if os.path.isfile(gtsfile): print "        * File '%s' exists; calculation omitted"%gtsfile
        else                      : sp_optimization(cp_name,symbols,dict_CPs,dict_pes,mes_optCP,ilines) 
        # Checking structure
        if not os.path.isfile(gtsfile):
           dict_warns[cp_name].append(1)
           continue
        idata_check = (gtsfile,cp_name,cptype,torsion1,torsion2,tolerance,dict_warns,symmetry,tsigma1,tsigma2,ttype)
        Etot, dict_warns = sp_checking(*idata_check)
        # Reference energy
        if Etot < min_sp : min_sp  = Etot
        print

    # Modify f_vfit and f_pes if minimum energy changes
    #change_Eref(Eref,min_sp,f_pes,f_vfit)

    #---------------------#
    # Rewrite splist file #
    #---------------------#
    sortedlist = sp_sorting(dict_CPs)
    print "     Summary table (relative energy in cm^-1):"
    print "                    -------------------------------"
    print "                      phi1  |  phi2  |   energy    "
    print "                    -------------------------------"

    dict_lines = {}
    lineformat = "    %i      %6.2f   %6.2f     %+13.2f       %3s        %18s  \n"
    lineformat1= "    %i      %6.2f   %6.2f     %+13.2f        NO        %18s  \n"
    lineformat2= "    %i      %6.2f   %6.2f     %+13.2f       YES        %18s  \n"
    num0, num1, num2 = 1, 1, 1
    for cp_name in sortedlist:
        phi1, phi2, cptype, Ecm, gtsfile, okcol = dict_CPs[cp_name]
        dict_lines[cp_name] = (Ecm,lineformat%(cptype,phi1*cons.R2D,phi2*cons.R2D,Ecm,okcol,cp_name))
        # gts file?
        if os.path.isfile(gtsfile):
           structure = gts2Struct(gtsfile,name=cp_name)
           phi1      = hf.ranged_angle(structure.icvalue(torsion1),v=1)
           phi2      = hf.ranged_angle(structure.icvalue(torsion2),v=1)
           Ecm       = (structure.get("Etot") - min_sp) * cons.h2c
           # Line for sp file
           if dict_warns.get(cp_name,[]) != []:
              dict_lines[cp_name] = (Ecm,lineformat1%(cptype,phi1*cons.R2D,phi2*cons.R2D,Ecm,cp_name))
           else:
              dict_lines[cp_name] = (Ecm,lineformat2%(cptype,phi1*cons.R2D,phi2*cons.R2D,Ecm,cp_name))
           # Print
           if cptype == 0: xxx = "MIN%02i"%num0; num0 += 1
           if cptype == 1: xxx = "TS_%02i"%num1; num1 += 1
           if cptype == 2: xxx = "MAX%02i"%num2; num2 += 1
           line = " %7s : %6.2f | %6.2f | %+9.2f "%(xxx,phi1*cons.R2D,phi2*cons.R2D,Ecm)
           print "          "+line
    print "                    -------------------------------"
    print

    print "     Storing updated stationary point positions at: %s"%f_splist
    string_sp  = ""
    string_sp += "#                                                                   \n"
    string_sp += "#  List of stationary points (SP) in the 2D-PES                     \n"
    string_sp += "# ---------------------------------------------                     \n"
    string_sp += "#         * Type = 0 for minima                                     \n"
    string_sp += "#         * Type = 1 for saddle points                              \n"
    string_sp += "#         * Type = 2 for maxima                                     \n"
    string_sp += "#         * Phi1 and Phi2 in [degrees]                              \n"
    string_sp += "#         * Energy in [1/cm] with regard to the global minimum      \n"
    string_sp += "#                                                                   \n"
    string_sp += "#  Type     Phi1     Phi2       Energy            opt OK?           SP name      \n"
    string_sp += "# ------ | ------ | ------ | ---------------- | ---------- | --------------------\n"
    for E,line in sorted(dict_lines.values()): string_sp += line
    string_sp += "# -------------------------------------------------------------------------------\n"
    string_sp += "eref %.7f # hartree\n"%min_sp
    string_sp += "# Number of stationary points: %i\n"%len(dict_lines.keys())

    finalsp = open(f_splist,"w")
    finalsp.write(string_sp)
    finalsp.close()
    print

    #----------------------#
    # Internal coordinates #
    #----------------------#
    icoords = sp_icdepuration(dict_CPs,torsion1,torsion2,xyzStruct,f_ics)

    #---------------------------------------#
    # Information for each stationary point #
    #---------------------------------------#
    writefile_spinfo(f_spinfo,dict_CPs,icoords,Tlist,freqscal,masslist)

    #-------------------#
    # Indicate warnings #
    #-------------------#
    for cpname,warns in dict_warns.items():
        if warns == []: continue
        print "     Warnings in stationary point '%s':"%cpname
        if 1 in warns: print "           - gts file not found"
        if 2 in warns: print "           - gts file does not contain hessian"
        if 3 in warns: print "           - dihedrals changed significantly"
        if 4 in warns: print "           - incorrect number of real/imag frequencies"
        print 

#----------------------------------------------------------#
def q2dtor_dijfitting(name,tuple_fit,icoords,symmetry,tsigma1,tsigma2,f_pes,file_Dfit,masslist=None):

    fterms, coef_nums, weight, ignore = tuple_fit

    #-----------#
    # Read data #
    #-----------#
    dict_pes, symbols, Eref = readfile_pes(f_pes) # incomplete data

    list_dij = []
    for pes_point in dict_pes.keys():
        xvec = dict_pes[pes_point][3]
        structure = Struct("",xvec,symbols,masslist,-1)
        structure.basic_setups([0,3])
        Dmatrix = structure.get_Dmatrix(icoords) * cons.amu * cons.angstrom**2
        # in [amu]*[angstrom]^2
        rI1  = + Dmatrix[0,0]
        rI2  = + Dmatrix[1,1]
        l12  = - Dmatrix[0,1]
        detD = rI1*rI2 - l12*l12
        d11  = 100*(rI2/detD)
        d22  = 100*(rI1/detD)
        d12  = 100*(l12/detD)
        dij_tuple = np.array( (d11,d22,d12) )
        dict_pes[pes_point] = dict_pes[pes_point] + [dij_tuple]
        
    dict_pes = pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2)

    list_pes = dict_pes.keys()
    xdata    = np.array([tuple(dict_pes[point][0:2]) for point in list_pes])

    ##########################################
    # (C) Fitting process                    #
    ##########################################

    output = open(file_Dfit,'w')
    output.close()

    for dij in ["d11","d22","d12"]:

       #print "      *===================================*"
       #print "      | Starting fitting for %3s elements |"%dij
       #print "      *===================================*"
        if dij == "d11": ydata = [dict_pes[pes_point][-1][0] for pes_point in list_pes]
        if dij == "d22": ydata = [dict_pes[pes_point][-1][1] for pes_point in list_pes]
        if dij == "d12": ydata = [dict_pes[pes_point][-1][2] for pes_point in list_pes]
        y_min  = min(ydata)
        y_mean = sum(ydata) / len(ydata)
        y_max  = max(ydata)
       #print "          - Min. 100*%s value is: %+9.2f"%(dij,y_min)
       #print "          - Avg. 100*%s value is: %+9.2f"%(dij,y_mean)
       #print "          - Max. 100*%s value is: %+9.2f"%(dij,y_max)
       #print

        FourierObject = Fourier2D(fterms, imag=False)
        fterms2, fcoefs, fit_errors, fitting_time = FourierObject.fit(xdata,ydata)
    
        # write fitting parameters first
        dfit = open(file_Dfit,'a')
        dfit.write("#----------------------------------------------#\n")
        if dij == "d11": dfit.write("#            Fitting data for d_11             #\n")
        if dij == "d12": dfit.write("#            Fitting data for d_12             #\n")
        if dij == "d22": dfit.write("#            Fitting data for d_22             #\n")
        dfit.write("#----------------------------------------------#\n")
        dfit.write("# Min. 100*%s value is: %+9.2f \n"%(dij,y_min))
        dfit.write("# Avg. 100*%s value is: %+9.2f \n"%(dij,y_mean))
        dfit.write("# Max. 100*%s value is: %+9.2f \n"%(dij,y_max))
        # Write data
        if dij == "d11": dfit.write("start_fitdata11\n")
        if dij == "d12": dfit.write("start_fitdata12\n")
        if dij == "d22": dfit.write("start_fitdata22\n")
        dfit.close()
        writefile_xfit(fterms2,fcoefs,file_Dfit,"a")
        if dij == "d11": dfit = open(file_Dfit,'a'); dfit.write("end_fitdata11\n")
        if dij == "d12": dfit = open(file_Dfit,'a'); dfit.write("end_fitdata12\n")
        if dij == "d22": dfit = open(file_Dfit,'a'); dfit.write("end_fitdata22\n")
        dfit.write("#----------------------------------------------#\n")
        dfit.write("# fitting correlation: (1.0 - r^2) = %7.1e   #\n"%(1.0-fit_errors[0]))
        dfit.write("# average abs. errors:                         #\n")
        dfit.write("#    %7.1e                                   #\n"%fit_errors[1])
        dfit.write("#    %7.1e (for points below mean value)     #\n"%fit_errors[2])
        dfit.write("# elapsed time: %6.1f seconds                 #\n"%fitting_time)
        dfit.write("#----------------------------------------------#\n\n\n")
        dfit.close()
#----------------------------------------------------------#
def q2dtor_tor2dns(name,inpdata,files2read,files2write,Tlist):

    #---------------#
    # Expand tuples #
    #---------------#
    (torsion1,torsion2,vfit,symmetry,tsigma1,tsigma2,dijvar,kmax,maxeigen) = inpdata
    (f_xyz,f_vfit,f_ics,f_pes,f_splist) = files2read
    (f_dfit,f_evals) = files2write
    (fterms,nums,weight,ignore) = vfit

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information from input file:"
    print  "        * torsion1      %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print  "        * torsion2      %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print  "        * symmetry      %s"%symmetry
    print  "        * tsigma1       %i"%tsigma1
    print  "        * tsigma2       %i"%tsigma2
    if dijvar: print  "        * dijvar        yes"
    else     : print  "        * dijvar        no"
    print  "        * kmax          %i"%kmax
    print  "        * maxeigen      %.2f cm^-1"%maxeigen
    print


    #------------#
    # Read files #
    #------------#
    files2read = (f_xyz,f_vfit,f_ics,f_pes,f_splist)
    which      = ["xyz","vfit","ics","pes","splist"]
    dataifiles = readfiles(files2read,which)

    [dict_pes, symbols, Eref] = dataifiles["pes"]
    [dict_CPs, Eref_cp] = dataifiles["splist"]
    [xyzStruct,symbols,masslist] = dataifiles["xyz"]
    [fterms, parameters, Eref_vfit] = dataifiles["vfit"]
    [icoords, nics] = dataifiles["ics"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print


    #======><======#


    #--------------------#
    # Deal with energies #
    #--------------------#
    print "     Analyzing minimum energy:"
    EMIN_pes = float("inf")
    for value in dict_pes.values():
        if value[2] < EMIN_pes : EMIN_pes = value[2]
    print "        - Mininum energy in numeric PES: %+13.7f hartree"%EMIN_pes
    EMIN_sp = float("inf")
    for value in dict_CPs.values():
        struct_instance = gts2Struct(value[4])
        if struct_instance.get("Etot") < EMIN_sp : EMIN_sp = struct_instance.get("Etot")
    print "        - Mininum energy in Fourier PES: %+13.7f hartree"%EMIN_sp
    if EMIN_sp < EMIN_pes:
       dE = (EMIN_pes-EMIN_sp) / cons.h / cons.c0 / cons.cm
       print "        - Reference energy in Fourier PES will be changed:"
       for idx in range(len(fterms)):
           if fterms[idx][0] == 'const':
              print "          * Constant term in Fourier potential was   : %+11.5f cm^-1"%(parameters[idx])
              parameters[idx] = parameters[idx] + dE
              print "          * Constant term in Fourier potential is now: %+11.5f cm^-1"%(parameters[idx])
    print

    #------------------------------------#
    # Check that torsions are in icoords #
    #------------------------------------#
    icoords = ic_checktorsions(icoords,torsion1,torsion2)

    #----------------#
    # Initialization #
    #----------------#
    CALCULATE = True
    if os.path.isfile(f_evals):
       print "     File '%s' already exists"%f_evals
       print "        - Data comparison (file,input):"
       datainfile = readfile_evals(f_evals)
       evalues    = datainfile[0]
       dif_kmax   = datainfile[1]
       dif_dijvar = datainfile[2]
       print "          kmax   ==> (%3i,%3i)"%(dif_kmax,kmax)
       if       dif_dijvar and     dijvar: print "          dijvar ==> (yes,yes)"
       elif     dif_dijvar and not dijvar: print "          dijvar ==> (yes, no)"
       elif not dif_dijvar and     dijvar: print "          dijvar ==> ( no,yes)"
       else                              : print "          dijvar ==> ( no, no)"
       if dif_dijvar == dijvar:
          if   dif_kmax == kmax: print "        - Calculation omitted!"; CALCULATE = False
          elif dif_kmax  > kmax: print "        - Calculation omitted! (kmax in '%s' is bigger)"%f_evals; CALCULATE = False
       print

    if CALCULATE:
       #-------------------------------#
       # (a) dij variation is required #
       #-------------------------------#
       if dijvar:
          print "     Dependence of dij terms with (phi1,phi2) is required"
          fitdij = True
          if os.path.exists(f_dfit):
             print "        * File '%s' already exists:"%(f_dfit)
             d11_terms, d11_coefs, dummy = readfile_xfit(f_dfit,"d11")
             d12_terms, d12_coefs, dummy = readfile_xfit(f_dfit,"d12")
             d22_terms, d22_coefs, dummy = readfile_xfit(f_dfit,"d22")
             if d11_coefs == [] or d12_coefs == [] or d22_coefs == []:
                print "          - Something wrong in file '%s'! File will overwritten...\n"%f_dfit
             else:
                fitdij = False
                print "          - dij fitting will be omitted"
          # Perform dij fitting to Fourier Series
          if fitdij:
             print "        * Fitting dij elements across 2D-PES..."
             q2dtor_dijfitting(name,vfit,icoords,symmetry,tsigma1,tsigma2,f_pes,f_dfit,masslist)
             print "        * File '%s' with fittings for dij elements was created"%f_dfit
          # read dfit file
          d11_terms, d11_coefs, dummy = readfile_xfit(f_dfit,"d11")
          d12_terms, d12_coefs, dummy = readfile_xfit(f_dfit,"d12")
          d22_terms, d22_coefs, dummy = readfile_xfit(f_dfit,"d22")
          if d11_coefs == [] or d12_coefs == [] or d22_coefs == []:
             print "        * Something wrong in file '%s'!\n"%f_dfit
             sys.exit(EXITMESS)
          dijtuple = [ (d11_terms,d11_coefs) , (d12_terms,d12_coefs) , (d22_terms,d22_coefs) ]
       #-----------------------------------#
       # (b) dij variation is not required #
       #-----------------------------------#
       else:
          print "     Dependence of dij terms with (phi1,phi2) is not required"
          E_min, cpname_min = sorted([ (dict_CPs[cp_name][3],cp_name) for cp_name in dict_CPs.keys() ])[0]
          print "        * dij values will be defined by '%s'"%(cpname_min)
          gtsfile = DIRSP+cpname_min+".gts"
          print "        * reading file: %s"%(gtsfile)
          structure = gts2Struct(gtsfile,"",masslist)
          structure.basic_setups([0,3])
          Dmatrix = structure.get_Dmatrix(icoords)
          # reduced moments of inertia in [amu]*[angstrom]^2
          rI1 = Dmatrix[0,0] * cons.amu * cons.angstrom**2
          rI2 = Dmatrix[1,1] * cons.amu * cons.angstrom**2
          dijtuple = [rI1,rI2]
       print

       #----------------------------#
       # 2DNS torsional Hamiltonian #
       #----------------------------#
       print "     Torsional Hamiltonian matrix:"
       print "        * Calculating..."
       fil,col,ham2d = tor2DNS_hamiltonian(name, fterms, parameters, dijvar, dijtuple, nums, kmax)
       print "        * Diagonalizing..."
       evalues  = tor2DNS_diagonalization(kmax,fil,col,ham2d,dijvar,f_evals,ENERMAX=maxeigen)

       #---------------#
       # Check evalues #
       #---------------#
       negative = False
       for evalue in evalues:
           if evalue<0: negative = True
       if negative:
         print "        * Negative evalues have been obtained... You should check Hamiltonian options..."
       print

    #---------------#
    # Print evalues #
    #---------------#
    print "     2DNS eigenvalues (cm^-1):"
    for idx in range(0,len(evalues),8):
        print "            " + "  ".join(["%8.2f"%evalue for evalue in evalues[idx:idx+8]])
    print

    #------------------------------#
    # Calculate partition function #
    #------------------------------#
    print "     2DNS partition function:"
    print "       ---------------------"
    print "         T (K)  | 2D-NS pfn "
    print "       ---------------------"
    evalues   = [ Ecm * cons.c2h for Ecm  in evalues]
    for T in sorted(Tlist):
        beta = 1.0 / cons.kB / T
        ptfn_2DNS = sum([np.exp(-E*beta) for E in evalues]) / float(tsigma1) / float(tsigma2)
        print "        %7.2f | %9.3E "%(T,ptfn_2DNS)
    print "       ---------------------"
    print
#----------------------------------------------------------#
def q2dtor_rovibpf(name,inpdata,files2read,files2write,Tlist,thermo=False):

    #---------------#
    # Expand tuples #
    #---------------#
    (torsion1,torsion2,ttype,symmetry,tsigma1,tsigma2,interpolation,integrationstep,freqscal) = inpdata
    (f_xyz,f_pes,f_vfit,f_splist,f_ics,f_evals) = files2read
    (f_sqrt, f_tables) = files2write
    sigma_tor  = float(tsigma1) * float(tsigma2)

    #-----------------------------------#
    # Print info needed from input file #
    #-----------------------------------#
    print  "     Information needed from input file:"
    print  "        * torsion1         %s"%("-".join(["%i"%(i+1) for i in torsion1]))
    print  "        * torsion2         %s"%("-".join(["%i"%(i+1) for i in torsion2]))
    print  "        * ttype            %s"%ttype
    print  "        * symmetry         %s"%symmetry
    print  "        * tsigma1          %i"%tsigma1
    print  "        * tsigma2          %i"%tsigma2
    print  "        * interpolation    %s"%interpolation
    print  "        * integrationstep  %.3f degrees"%(integrationstep*cons.R2D)
    print  "        * freqscal         %.3f"%freqscal
    print

    #------------#
    # Read files #
    #------------#
    which      = ["xyz","pes","vfit","splist","ics","evals"]
    dataifiles = readfiles(files2read,which)
    [xyzStruct,symbols,masslist] = dataifiles["xyz"]
    [dict_pes, symbols, EMIN_pes] = dataifiles["pes"]
    [fterms, parameters, EMIN_vfit] = dataifiles["vfit"]
    [dict_CPs, EMIN_sp] = dataifiles["splist"]
    [icoords, nics] = dataifiles["ics"]
    [evalues,zpe_2DNS] = dataifiles["evals"]

    #----------------------#
    # Dealing with folders #
    #----------------------#
    if not os.path.exists(DIRFILES):
       print  "     Creating folder: '%s'"%(DIRFILES)
       os.mkdir(DIRFILES)
       print

    #======><======#

    #--------------------#
    # Deal with energies #
    #--------------------#
    print "     Analyzing minimum energy:"
    if EMIN_vfit is None: EMIN_vfit = EMIN_pes
    if EMIN_sp is None:
       EMIN_sp = float("inf")
       for value in dict_CPs.values():
           struct_instance = gts2Struct(value[4])
           if struct_instance.get("Etot") < EMIN_sp : EMIN_sp = struct_instance.get("Etot") 

    print "        - Mininum energy in numeric PES         : %+13.7f hartree"%EMIN_pes
    print "        - Mininum energy in Fourier PES         : %+13.7f hartree"%EMIN_vfit
    print "        - Mininum energy among stationary points: %+13.7f hartree"%EMIN_sp

    if EMIN_sp < EMIN_vfit:
       dE = (EMIN_vfit-EMIN_sp) / cons.h / cons.c0 / cons.cm
       print "        - Correcting Fourier PES:"
       for idx in range(len(fterms)):
           if fterms[idx][0] == 'const':
              print "          * Constant term in Fourier potential was   : %+11.5f cm^-1"%(parameters[idx])
              parameters[idx] = parameters[idx] + dE
              print "          * Constant term in Fourier potential is now: %+11.5f cm^-1"%(parameters[idx])
    print


    #-------------------------------#
    # Fourier PES + Critical Points #
    #-------------------------------#
    Fourier_PES = Fourier2D(fterms, parameters)
    dict_CPs, list_CPs, initial_CPs = sp_allCPs(dict_CPs,symmetry,name,tsigma1,tsigma2)

    #------------------------------------#
    # Check that torsions are in icoords #
    #------------------------------------#
    icoords = ic_checktorsions(icoords,torsion1,torsion2)

    #---------------------------------------------#
    # Calculate D matrix, its determinant and the #
    # 3N-6-2 nontorsional frequencies for each CP #
    #---------------------------------------------#
    print "     Analyzing stationary point structures"

    # Get 'statpoint' string for tables
    ml = max([len(cpname) for cpname in list_CPs])
    spstring = "statpoint"
    while len(spstring) < ml: spstring = " "+spstring+" "
    if    len(spstring) > ml: spstring = spstring[:-1]

    # Initialize some variables
    ref_pes  = [float("inf")]
    ref_msho = [float("inf")]
    ref_e2dt  = [float("inf")]
    data4tables = []

    # Analyze each CP
    min_structs = {}
    for cpname in list_CPs:
        phi1     = dict_CPs[cpname][0]
        phi2     = dict_CPs[cpname][1]
        cptype   = dict_CPs[cpname][2]
        gtsfile  = dict_CPs[cpname][4]
        # Get data
        data_sp = sp_analysis(gtsfile,icoords,cptype,freqscal,masslist,name=cpname)
        (struct,energy,rij,dij,ccfreqs,zpe,ntfreqs,zpe_2b) = data_sp
        (rI1,rI2,l12) = rij
        (d11,d22,d12) = dij
        # Moments of inertia and sigma rot
        IA, IB, IC = struct.get("imoments")
        sigma_rot  = float(struct.get("rotsigma"))
        # D matrix (in a.u. <--> [mass]*[distance]^2)
        detD  = (rI1*rI2 - l12*l12)
        sqrtD = np.sqrt(detD)
        # S matrix
        detS  = detD * IA * IB * IC
        sqrtS = np.sqrt(detS)
        # Energies
        PES   = energy
        ZPE   = zpe
        ntZPE = zpe_2b
        V1    = PES + ZPE
        V2    = PES + ntZPE + zpe_2DNS
        # Vibrational partition functions and its derivatives (non-torsional)
        ntmatrix = np.zeros( (len(Tlist) , 3) ) # columns: part function, 1st derivative, 2nd derivative
        for idx1 in range(len(Tlist)):
            T = Tlist[idx1]
            ntmatrix[idx1][0] = 1.0
            for ff in ntfreqs:
                qvib_i, zpe_i = ff.get_qvib(T)
                fdln_i = ff.get_fdln(T)
                sdln_i = ff.get_sdln(T)
                ntmatrix[idx1][0] *= qvib_i
                ntmatrix[idx1][1] += fdln_i
                ntmatrix[idx1][2] += sdln_i
        # Save data
        dict_CPs[cpname] = dict_CPs[cpname] + [sqrtD,sqrtS,sigma_rot,ntmatrix,zpe_2b]
        if cptype == 0:
           min_structs[cpname] = struct
           if PES < ref_pes[0] : ref_pes  = [PES,PES,  ZPE,cpname]; sqrtD0,sqrtS0 = sqrtD, sqrtS
           if V1  < ref_msho[0]: ref_msho = [V1 ,PES,  ZPE,cpname]
           if V2  < ref_e2dt[0]: ref_e2dt = [V2 ,PES,ntZPE,cpname]
        # Data for table
        data4tables.append( (cpname,cptype,PES,ZPE,ntZPE,V1,V2) )

    # Apply reference
    ER2 = ref_pes[0]
    ref_pes  = [ref_pes[0] -ER2 , ref_pes[1] -ER2 , ref_pes[2]  , ref_pes[3]]
    ref_msho = [ref_msho[0]-ER2 , ref_msho[1]-ER2 , ref_msho[2] , ref_msho[3]]
    ref_e2dt = [ref_e2dt[0]-ER2 , ref_e2dt[1]-ER2 , ref_e2dt[2] , ref_e2dt[3]]

    # Print data
    print "        - Energy of PES global minimum: %.7f hartree"%ER2
    print "        - Summary table (energies in %s):"%SU1
    print "            -%s----------------------------------------------------------------"%("-"*ml)
    print "             %s |  type   |  Erel   |   ZPE   |  ntZPE  ||    V1    |    V2    "%spstring
    print "            -%s----------------------------------------------------------------"%("-"*ml)
    for cpname,cptype,PES,ZPE,ntZPE,V1,V2 in data4tables:
        PES   = HU1 * (PES-ER2)
        ZPE   = HU1 * ZPE
        ntZPE = HU1 * ntZPE
        V1    = HU1 * (V1 -ER2)
        V2    = HU1 * (V2 -ER2)
        if cptype == 0: dataline = (cpname,"minimum",PES,ZPE,ntZPE,V1,V2)
        if cptype == 1: dataline = (cpname,"saddle ",PES,ZPE,ntZPE,V1,V2)
        if cptype == 2: dataline = (cpname,"maximum",PES,ZPE,ntZPE,V1,V2)
        print "             %s | %s | %7.3f | %7.3f | %7.3f || %8.3f | %8.3f"%dataline
    print "            -%s----------------------------------------------------------------"%("-"*ml)
    print "          where:"
    print "            Erel is the potential energy (with regard to absolute minimum)"
    print "            ZPE is zero-point energy"
    print "            ntZPE is the non-torsional zero-point energy"
    print "            V1 is Erel + ZPE"
    print "            V2 is Erel + ntZPE + 2DNS zero-point energy"

    print "        - Structure with min(PES): %s (%7.3f %s)"%(ref_pes[3] ,ref_pes[0] *HU1,SU1)
    print "        - Structure with min(V1) : %s (%7.3f %s)"%(ref_msho[3],ref_msho[0]*HU1,SU1)
    print "        - Structure with min(V2) : %s (%7.3f %s)"%(ref_e2dt[3],ref_e2dt[0]*HU1,SU1)
    print
    
    #------------------------#
    # Teselation of 2D space #
    #------------------------#
    print "     Obtaining Delaunay tesselation using CPs"
    print 
    dict_CPs     = tess.replicate_points(name,dict_CPs)
    Dlist_CPs    = sp_sorting(dict_CPs)
    Dlist_points = [dict_CPs[cp_name][0:2] for cp_name in Dlist_CPs]
    delaunay     = Delaunay(Dlist_points)
    Vmidpoints   = tess.get_Vmidpoints(delaunay,Fourier_PES)

    #----------------------------#
    # |D| and |S| interpolations #
    #----------------------------#
    print "     Interpolating det(D) and det(S)"
    fitsqrtD = True
    fitsqrtS = True

    # Fourier interpolation and file exists
    if interpolation == "fourier" and os.path.isfile(f_sqrt):
       print "        * File '%s' already exists"%f_sqrt
       # Read reference values
       file_refD = None
       file_refS = None
       with open(f_sqrt,'r') as filewithsqrt:
           for line in filewithsqrt:
               if "ref_sqrtD" in line: file_refD = float(line.split()[1])
               if "ref_sqrtS" in line: file_refS = float(line.split()[1])
       # Read fittings
       ftD, fcD, dummy = readfile_xfit(f_sqrt,"sqrtD")
       ftS, fcS, dummy = readfile_xfit(f_sqrt,"sqrtS")
       # Is sqrt(D) fitting valid?
       if (file_refD is not None) and (ftD != []):
          diff = abs(100*(file_refD-sqrtD0)/sqrtD0)
          if diff < 0.1:
             print "          - contains information of sqrt(|D|) fitting"
             rel_sqrtD = Fourier2D(ftD)
             rel_sqrtD.set_coefs(fcD)
             fitsqrtD = False
          else:
             print "          - reference value for sqrt(D) differs (%.4E vs %.4E)"%(file_refD,sqrtD0)
       else:
          print "          - does not contain (valid) information of sqrt(|D|) fitting"
       # Is sqrt(S) fitting valid?
       if (file_refS is not None) and (ftS != []):
          diff = abs(100*(file_refS-sqrtS0)/sqrtS0)
          if diff < 0.1:
             print "          - contains information of sqrt(|S|) fitting"
             rel_sqrtS = Fourier2D(ftS)
             rel_sqrtS.set_coefs(fcS)
             fitsqrtS = False
          else:
             print "          - reference value for sqrt(S) differs (%.4E vs %.4E)"%(file_refS,sqrtS0)
       else:
          print "          - does not contain (valid) information of sqrt(|S|) fitting"

    if interpolation != "fourier" or fitsqrtD or fitsqrtS:
       print "        * Obtaining |D| and |S| in all PES points..."
       print "            reference value for sqrt(|D|) is %.2E au"%sqrtD0
       print "            reference value for sqrt(|S|) is %.2E au"%sqrtS0
       for pes_point in dict_pes.keys():
           xvec = dict_pes[pes_point][3]
           structure = Struct("",xvec,symbols,masslist,-1)
           structure.basic_setups([0,3,4])
           # Principal moments of inertia, D matrix, and S matrix
           IA, IB, IC = structure.get("imoments")
           Dmatrix = structure.get_Dmatrix(icoords)
           rI1  = +Dmatrix[0,0]
           rI2  = +Dmatrix[1,1]
           l12  = -Dmatrix[0,1]
           detD = rI1*rI2 - l12*l12
           detS = detD * IA*IB*IC
           # Square-root of |D| and |S| (with regards to reference)
           sqrtD = np.sqrt(detD) / sqrtD0
           sqrtS = np.sqrt(detS) / sqrtS0
           dict_pes[pes_point] = dict_pes[pes_point] + [sqrtD,sqrtS]
       dict_pes = pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2)

       # Data for interpolations
       list_pes   = dict_pes.keys()
       data_phi1  = [dict_pes[pname][ 0] for pname in list_pes]
       data_phi2  = [dict_pes[pname][ 1] for pname in list_pes]
       data_sqrtD = [dict_pes[pname][-2] for pname in list_pes]
       data_sqrtS = [dict_pes[pname][-1] for pname in list_pes]

       if interpolation == "fourier":
          data_XY = [dict_pes[pname][0:2] for pname in list_pes]
          if fitsqrtD:
             print "        * Fitting sqrt(|D|) to Fourier series..."
             rel_sqrtD = Fourier2D(fterms,imag=False)
             ftD, fcD, feD, ftimeD = rel_sqrtD.fit(data_XY,data_sqrtD)
          if fitsqrtS:
             print "        * Fitting sqrt(|S|) to Fourier series..."
             rel_sqrtS = Fourier2D(fterms,imag=False)
             ftS, fcS, feS, ftimeS = rel_sqrtS.fit(data_XY,data_sqrtS)
          # Write data
          print "        * Saving Fourier Series fittings to: %s"%f_sqrt
          filesqrt = open(f_sqrt,'w')
          filesqrt.write("ref_sqrtD %.4E\n"%sqrtD0)
          filesqrt.write("ref_sqrtS %.4E\n"%sqrtS0)
          filesqrt.write("\n")
          filesqrt.close()
          filesqrt = open(f_sqrt,'a'); filesqrt.write("start_sqrtD\n"); filesqrt.close()
          writefile_xfit(ftD,fcD,f_sqrt,"a")
          filesqrt = open(f_sqrt,'a'); filesqrt.write("end_sqrtD\n\n")  ; filesqrt.close()
          filesqrt = open(f_sqrt,'a'); filesqrt.write("start_sqrtS\n"); filesqrt.close()
          writefile_xfit(ftS,fcS,f_sqrt,"a")
          filesqrt = open(f_sqrt,'a'); filesqrt.write("end_sqrtS\n\n")  ; filesqrt.close()
       else:
          warnings.filterwarnings("ignore")
          kx, ky = interpolation, interpolation
          print "        * Creating (%i,%i)-spline for sqrt(|D|)..."%(kx,ky)
          rel_sqrtD = SmoothBivariateSpline(data_phi1, data_phi2, data_sqrtD, kx=kx, ky=ky, s=0)
          print "        * Creating (%i,%i)-spline for sqrt(|S|)..."%(kx,ky)
          rel_sqrtS = SmoothBivariateSpline(data_phi1, data_phi2, data_sqrtS, kx=kx, ky=ky, s=0)

    print "        * Checking interpolation: calculated  vs interpolated values for stationary points [sqrt(|D|)], [sqrt(|S|)]"
    wrong = False
    MAXERROR = 2.0
    for cpname in list_CPs:
        phi1     = dict_CPs[cpname][0]
        phi2     = dict_CPs[cpname][1]
        cp_sqrtD = dict_CPs[cpname][-5]/sqrtD0
        cp_sqrtS = dict_CPs[cpname][-4]/sqrtS0
        sp_sqrtD = rel_sqrtD(phi1,phi2)
        sp_sqrtS = rel_sqrtS(phi1,phi2)
        # relative error
        diffD = 100.0*abs(sp_sqrtD-cp_sqrtD)/cp_sqrtD
        diffS = 100.0*abs(sp_sqrtS-cp_sqrtS)/cp_sqrtS
        line1 = "            %s:  [%.3f vs %.3f (%3.1f%%)] , [%.3f vs %.3f (%3.1f%%)]"
        line2 = "            %s:  [%.3f vs %.3f (%3.1f%%)] , [%.3f vs %.3f (%3.1f%%)] <--"
        if diffD > MAXERROR or diffS > MAXERROR:
           wrong = True
           print line2%(cpname,cp_sqrtD,sp_sqrtD,diffD,cp_sqrtS,sp_sqrtS,diffS)
       #else: print line1%(cpname,cp_sqrtD,sp_sqrtD,diffD,cp_sqrtS,sp_sqrtS,diffS)
    if wrong: print "        * WARNING!! Problems with interpolation!"
    print

    # v1.1: get max sigma rot
    if BOOL_v1_1: rotsigma_max = max([dict_CPs[cpname][-3] for cpname in dict_CPs.keys()])

    #---------------------------------------#
    # Calculate (tilde) partition functions #
    #---------------------------------------#

    print "    ****************************************"
    print "     CALCULATION OF THE PARTITION FUNCTIONS "
    print "    ****************************************"
    print
   #print "    IMPORTANT: the partition functions are listed as:"
   #print "               (a) The reference energy is set to the lowest"
   #print "                   allowed energy (ZPE-exclusive)"
   #print "               (b) The reference energy is set to the potential"
   #print "                   energy of the global minimum (PES)"
   #print
    # Get firstly some constants (check manual)
    A1 = 1.0 / (2.0 * np.pi * cons.hbar * cons.hbar * sigma_tor)
    A2 = 8.0 * np.pi**2 / sigma_tor / (2.0*np.pi*cons.hbar**2)**2.5 # sigmarot included in integrand
    DPHI = integrationstep*cons.R2D

    # Initialize lists
    pf1WHO, pfMSHO = [], []
    pf2DNS, pfTorClas, pfEHR  = [], [], []
    pfTRA, pfELE = [], []
    ratio2DNS = []

    # Begin calculation
    Tlist.sort()
    for idx in range(len(Tlist)):
        t1   = time.time()
        T    = Tlist[idx]
        beta = 1.0 / cons.kB / T

        #-----------------------------------------------#
        # Volume at 1E5 Pa = 1 bar; i.e. standard state #
        #-----------------------------------------------#
        p0      = 1E5 * (cons.meter**3 / cons.joule)
        volume  = cons.kB * T / p0 # volume of a single molecule

        #-------------------------------#
        # Calculate Q_1W-HO and Q_MS-HO #
        #-------------------------------#
        ptfn_MSHO   = 0.0 # rovibrational partition function (MSHO)
        for cpname,structure in min_structs.items():
            # Partition functions
            phi_tra, qrot, qvib, qele, E00 = structure.get_pfns(T,k="cc")
            # Save translational and electronic 
            ptfn_tr  = phi_tra * volume
            ptfn_ele = qele
            # 1WHO rovibrational part
            ptfn_1WHO_j = qrot * qvib
            if cpname == ref_msho[3]: ptfn_1WHO = ptfn_1WHO_j
            # v1.1 - remove contribution of rotsigma to j-th structure
            if BOOL_v1_1:
               rotsigma = structure.get("rotsigma")
               ptfn_1WHO_j *= rotsigma
            # add to MSHO rovibrational part
            dj  = (E00-ER2) - ref_msho[0]
            exp = np.exp(-dj*beta)
            ptfn_MSHO += ptfn_1WHO_j * exp
        ptfn_MSHO = ptfn_MSHO / sigma_tor

        #----------------------#
        # calculate Q_tor^2DNS #
        #----------------------#
        ptfn_2DNS = + sum([np.exp(-(E-zpe_2DNS)*beta) for E in evalues]) / sigma_tor

        #------------------------#
        # Calculate Q_cl,tor^(C) #
        #------------------------#
        args       = [T,Fourier_PES,rel_sqrtD,sqrtD0]
        pf_torClas = + (A1/beta) * integration2D(intgrnd_ClasTor_pf,args=args,dphi=DPHI)

        #-----------------#
        # Calculate Q^EHR #
        #-----------------#
        data4EHR = np.zeros( (5,len(Dlist_CPs)) )
        for cp_idx in range(len(Dlist_CPs)):
            cpname   = Dlist_CPs[cp_idx]
            sigmarot = dict_CPs[cpname][-3]
            ntmatrix = dict_CPs[cpname][-2]
            ntzpe    = dict_CPs[cpname][-1]
            ntpf     = float(ntmatrix[idx][0])
            ntfdln   = float(ntmatrix[idx][1])
            ntsdln   = float(ntmatrix[idx][2])
            # Add data
            data4EHR[0][cp_idx] = sigmarot
            data4EHR[1][cp_idx] = ntzpe-ref_e2dt[2]
            data4EHR[2][cp_idx] = ntpf
            data4EHR[3][cp_idx] = ntfdln
            data4EHR[4][cp_idx] = ntsdln
        tupleS  = (sqrtS0,rel_sqrtS)
        tupleNT = (delaunay,Vmidpoints,data4EHR[:3,:])
        args    = [T,Fourier_PES,tupleS,tupleNT,ref_e2dt[1]]
        rv_EHR  = + (A2 * beta**(-2.5)) * integration2D(intgrnd_EHR_pf,args=args,dphi=DPHI)
        
        # v1.1: divide by rotsigma_max
        if BOOL_v1_1:
           ptfn_MSHO  /= rotsigma_max
           rv_EHR     /= rotsigma_max

        #-------------#
        # Append data #
        #-------------#
        pf1WHO.append(float(ptfn_1WHO))
        pfMSHO.append(float(ptfn_MSHO))
        pf2DNS.append(float(ptfn_2DNS))
        pfEHR.append( float(rv_EHR   ))
        pfTRA.append(ptfn_tr)
        pfELE.append(ptfn_ele)
        pfTorClas.append(float(pf_torClas))
        expZPE = np.exp(-beta*zpe_2DNS)
        ratio2DNS.append(float(expZPE*ptfn_2DNS/pf_torClas))

        t2 = time.time()
        print "        * T = %7.2f K done! (%.1f seconds)"%(T,t2-t1)

    #------------------------------------#
    # Correct ratio in case it decreases #
    #------------------------------------#
    for idx in range(1,len(ratio2DNS)):
        value_a = ratio2DNS[idx-1]
        value_b = ratio2DNS[idx]
        if value_b == 1.0   : ratio2DNS[idx] = "unit"
        if value_b < value_a: ratio2DNS[idx] = "unit"

    #----------------#
    # Calculate E2DT #
    #----------------#
    pfE2DT = []
    for idx in range(len(Tlist)):
        T      = Tlist[idx]
        beta   = 1.0 / cons.kB / T
        expZPE = np.exp(-beta*zpe_2DNS)
        pfn_TorClas = pfTorClas[idx]
        Fq2DNS      = ratio2DNS[idx]
        if Fq2DNS  == "unit": Fq2DNS = 1.0
        pfn_EHR     = pfEHR[idx]
        pfn_2DNS    = Fq2DNS*pfn_TorClas/expZPE
        # Correct 2DNS
        pf2DNS[idx] = pfn_2DNS
        # Get E2DT
        pfn_E2DT = (pfn_2DNS/pfn_TorClas)*pfn_EHR
        pfE2DT.append(pfn_E2DT)

    #--------#
    # Tables #
    #--------#
    Eref_1WHO = ref_msho[0]
    Eref_MSHO = ref_msho[0]
    Eref_E2DT = ref_e2dt[0]
    tuple_Eref = (Eref_1WHO,Eref_MSHO,zpe_2DNS,Eref_E2DT-zpe_2DNS,Eref_E2DT)
    tables     = string_pfntable(Tlist,pf1WHO,pfMSHO,pfE2DT,pf2DNS,pfTorClas,pfEHR,ratio2DNS,tuple_Eref,pfTRA,pfELE)
    print
    for line in  tables.split("\n")[4:-2]: print line
    output = open(f_tables,'w')
    output.write(tables)
    output.close()
    print "     Tables were stored at '%s' file"%(f_tables)
    print

    if thermo is False: return

    #----------------#
    # Thermodynamics #
    #----------------#
    print "    ********************************************"
    print "     CALCULATION OF THE THERMODYNAMIC FUNCTIONS "
    print "    ********************************************"
    print
    ther1WHO, therMSHO, therE2DT = [], [], []
    for idx in range(len(Tlist)):
        t1   = time.time()
        T    = Tlist[idx]
        beta = 1.0 / cons.kB / T

        #-----------------------------#
        # Partition functions (tilde) #
        #-----------------------------#
        ptfn_tr      = pfTRA[idx]
        ptfn_torClas = pfTorClas[idx]
        ptfn_2DNS    = pf2DNS[idx]
        Fq2DNS       = ratio2DNS[idx]
        ptfn_EHR     = pfEHR[idx]
        ptfn_E2DT    = pfE2DT[idx]

        #-------------#
        # 1WHO & MSHO #
        #-------------#
        ptfn_MSHO   = 0.0
        fdln_MSHO   = 0.0 # (dlnQ/dbeta)     for rovibration (MSHO)
        sdln_MSHO   = 0.0 # (d^2lnQ/dbeta^2) for rovibration (MSHO)
        fd_MSHO     = 0.0 # (dQ/dbeta)       for rovibration (MSHO)
        sd_MSHO     = 0.0 # (d^2Q/dbeta^2)   for rovibration (MSHO)
        for cpname,structure in min_structs.items():
            # Partition functions
            phi_tra, qrot, qvib, qele, E00 = structure.get_pfns(T,k="cc")
            # Derivatives
            tuple_traV, tuple_traP, tuple_rot, tuple_vib, tuple_ele = pfn_derivatives(structure,qele,T,k="cc")
            # Translational part
            fdln_trV = tuple_traV[0]
            sdln_trV = tuple_traV[1]
            fdln_trP = tuple_traP[0]
            sdln_trP = tuple_traP[1]
            # Rotational part
            ptfn_rot = qrot
            fdln_rot = tuple_rot[0]
            sdln_rot = tuple_rot[1]
            # Vibrational part
            ptfn_vib = qvib
            fdln_vib = tuple_vib[0]
            sdln_vib = tuple_vib[1]
            # Electronic part 
            ptfn_ele  = qele
            fdln_ele = tuple_ele[0]
            sdln_ele = tuple_ele[1]
            # 1WHO rovibrational part
            ptfn_1WHO_j = ptfn_rot * ptfn_vib
            # v1.1 - remove rotsigma contribution
            if BOOL_v1_1:
               rotsigma     = structure.get("rotsigma")
               ptfn_1WHO_j *= rotsigma
            # derivative of log
            fdln_1WHO_j = fdln_rot + fdln_vib
            sdln_1WHO_j = sdln_rot + sdln_vib
            fd_1WHO_j   = fdln_1WHO_j*ptfn_1WHO_j
            sd_1WHO_j   = ptfn_1WHO_j*(sdln_1WHO_j+(fdln_1WHO_j)**2)
            if cpname == ref_msho[3]:
               ptfn_1WHO = ptfn_1WHO_j
               fdln_1WHO = fdln_1WHO_j
               sdln_1WHO = sdln_1WHO_j
               # v1.1 - apply rotsigma to 1WHO
               if BOOL_v1_1: ptfn_1WHO /= rotsigma
            # add to MSHO rovibrational part
            dj  = (E00-ER2) - ref_msho[0]
            exp = np.exp(-dj*beta)
            ptfn_MSHO += ptfn_1WHO_j * exp
            fd_MSHO   += (fd_1WHO_j-dj*ptfn_1WHO_j) * exp
            sd_MSHO   += (sd_1WHO_j+dj*dj*ptfn_1WHO_j-2*dj*fd_1WHO_j) * exp
        # v1.1 - divide by sigma_tor
        if BOOL_v1_1:
           ptfn_MSHO /= sigma_tor
           fd_MSHO   /= sigma_tor
           sd_MSHO   /= sigma_tor
        # v1.1 - correct with rotsigma_max
        if BOOL_v1_1:
           ptfn_MSHO /= rotsigma_max
           fd_MSHO   /= rotsigma_max
           sd_MSHO   /= rotsigma_max
        # MSHO rovibrational
        fdln_MSHO = fd_MSHO/ptfn_MSHO
        sdln_MSHO = sd_MSHO/ptfn_MSHO - (fd_MSHO/ptfn_MSHO)**2

        #----------------------------#
        # 2DNS, TorClas and Fq(2DNS) #
        #----------------------------#
        if Fq2DNS == "unit":
           fdln_Fq2DNS = 0.0
           sdln_Fq2DNS = 0.0
        else:
           args         = [T,Fourier_PES,rel_sqrtD,sqrtD0]
           fd_torClas   = - (A1/beta) * integration2D(intgrnd_ClasTor_fd,args=args,dphi=DPHI)
           sd_torClas   = + (A1/beta) * integration2D(intgrnd_ClasTor_sd,args=args,dphi=DPHI)
           # get derivative of log
           fdln_torClas = fd_torClas / ptfn_torClas
           sdln_torClas = sd_torClas / ptfn_torClas - (fd_torClas / ptfn_torClas)**2

           fd_2DNS   = - sum([  (E-zpe_2DNS)    *np.exp(-(E-zpe_2DNS)*beta) for E in evalues]) / sigma_tor
           sd_2DNS   = + sum([ ((E-zpe_2DNS)**2)*np.exp(-(E-zpe_2DNS)*beta) for E in evalues]) / sigma_tor
           # get derivative of log
           fdln_2DNS = fd_2DNS/ptfn_2DNS   
           sdln_2DNS = sd_2DNS/ptfn_2DNS - (fd_2DNS/ptfn_2DNS)**2

           fdln_Fq2DNS = fdln_2DNS - fdln_torClas - zpe_2DNS
           sdln_Fq2DNS = sdln_2DNS - sdln_torClas


        #-----------------#
        # Calculate Q^EHR #
        #-----------------#
        data4EHR = np.zeros( (5,len(Dlist_CPs)) )
        for cp_idx in range(len(Dlist_CPs)):
            cpname   = Dlist_CPs[cp_idx]
            sigmarot = dict_CPs[cpname][-3]
            ntmatrix = dict_CPs[cpname][-2]
            ntzpe    = dict_CPs[cpname][-1]
            ntpf     = float(ntmatrix[idx][0])
            ntfdln   = float(ntmatrix[idx][1])
            ntsdln   = float(ntmatrix[idx][2])
            # Add data
            data4EHR[0][cp_idx] = sigmarot
            data4EHR[1][cp_idx] = ntzpe-ref_e2dt[2]
            data4EHR[2][cp_idx] = ntpf
            data4EHR[3][cp_idx] = ntfdln
            data4EHR[4][cp_idx] = ntsdln
        tupleS  = (sqrtS0,rel_sqrtS)
        tupleNT = (delaunay,Vmidpoints,data4EHR)
        args    = [T,Fourier_PES,tupleS,tupleNT,ref_e2dt[1]]
        rv_EHR  = pfEHR[idx]
        fd_EHR  = + (A2 * beta**(-2.5)) * integration2D(intgrnd_EHR_fd,args=args,dphi=DPHI)
        sd_EHR  = + (A2 * beta**(-2.5)) * integration2D(intgrnd_EHR_sd,args=args,dphi=DPHI)
        # v1.1 - correct with rotsigma_max
        if BOOL_v1_1:
           fd_EHR /= rotsigma_max
           sd_EHR /= rotsigma_max
        # get derivative of log
        fdln_EHR = fd_EHR / rv_EHR    
        sdln_EHR = sd_EHR / rv_EHR - (fd_EHR / rv_EHR)**2

        #------#
        # E2DT #
        #------#
        fdln_E2DT = fdln_Fq2DNS + fdln_EHR + zpe_2DNS
        sdln_E2DT = sdln_Fq2DNS + sdln_EHR

        #-----------------------------------#
        # Calculate thermodynamic functions #
        #-----------------------------------#
        log_1WHO = np.log(ptfn_tr) + np.log(ptfn_1WHO) + np.log(ptfn_ele)
        U_1WHO   = -fdln_trV - fdln_1WHO - fdln_ele
        H_1WHO   = -fdln_trP - fdln_1WHO - fdln_ele
        cp_1WHO  = cons.kB * beta * beta * ( sdln_trP + sdln_1WHO + sdln_ele )
        S_1WHO   = cons.kB * log_1WHO - cons.kB * beta * (fdln_trP + fdln_1WHO + fdln_ele)
        G_1WHO   =  - log_1WHO / beta
        ther1WHO.append( (U_1WHO,H_1WHO,cp_1WHO,S_1WHO,G_1WHO) )

        log_MSHO = np.log(ptfn_tr) + np.log(ptfn_MSHO) + np.log(ptfn_ele)
        U_MSHO   = -fdln_trV - fdln_MSHO - fdln_ele
        H_MSHO   = -fdln_trP - fdln_MSHO - fdln_ele
        cp_MSHO  = cons.kB * beta * beta * ( sdln_trP + sdln_MSHO + sdln_ele )
        S_MSHO   = cons.kB * log_MSHO - cons.kB * beta * (fdln_trP + fdln_MSHO + fdln_ele)
        G_MSHO   =  - log_MSHO / beta
        therMSHO.append( (U_MSHO,H_MSHO,cp_MSHO,S_MSHO,G_MSHO) )

        log_E2DT = np.log(ptfn_tr) + np.log(ptfn_E2DT) + np.log(ptfn_ele)
        U_E2DT   = -fdln_trV - fdln_E2DT - fdln_ele
        H_E2DT   = -fdln_trP - fdln_E2DT - fdln_ele
        cp_E2DT  = cons.kB * beta * beta * ( sdln_trP + sdln_E2DT + sdln_ele )
        S_E2DT   = cons.kB * log_E2DT - cons.kB * beta * (fdln_trP + fdln_E2DT + fdln_ele)
        G_E2DT   =  - log_E2DT / beta
        therE2DT.append( (U_E2DT,H_E2DT,cp_E2DT,S_E2DT,G_E2DT) )

        t2 = time.time()
        print "        * T = %7.2f K done! (%.2f seconds)"%(T,t2-t1)

    #--------#
    # Tables #
    #--------#
    tables = string_thermotable(Tlist,ther1WHO, therMSHO, therE2DT)
    print
    for line in  tables.split("\n")[4:-2]: print line
    output = open(f_tables,'a')
    output.write(tables)
    output.close()
    print "     Tables were stored at '%s' file"%(f_tables)
    print
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 4) Q2DTor tools                                #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def tools(name,avail_args):
    TOOLMESS = ""
    TOOLMESS = TOOLMESS + " ====================\n"
    TOOLMESS = TOOLMESS + " | Using Q2DTor tool: %s\n"
    TOOLMESS = TOOLMESS + " ====================\n"

    #----------------------#
    # Any tool to be used? #
    #----------------------#
    if '--pdf' in avail_args.keys():
       print TOOLMESS%"--pdf"
       if name is None: sys.exit("     Main name needed")
       if name.endswith(".inp"): ifile = name; name  = name[:-4]
       else: ifile = name + ".inp"
       if not os.path.isfile(ifile): sys.exit("     Q2DTor input file '%s' does not exists!\n"%ifile)
       tool_pdf(name,ifile)
       return True

    if "--gts" in avail_args.keys():
       print TOOLMESS%"--gts"
       software = avail_args["--gts"]
       tool_gts(name,software)
       return True

    if '--icoords' in avail_args.keys():
       print TOOLMESS%"--icoords"
       mode, ics = avail_args["--icoords"]
       tool_icoords(name,mode,ics)
       return True

    return False
#----------------------------------------------------------#
def tool_pdf(name,ifile):
    #--------------------------------
    def sp_colors(cpnames):
        # Some starting colors
        mcolors = []
        mcolors.append( np.array( [0.25,0.25,0.25] ) )
        mcolors.append( np.array( [0.50,0.50,0.50] ) )
        mcolors.append( np.array( [0.75,0.75,0.75] ) )
        mcolors.append( np.array( [0.25,0.25,0.75] ) )
        mcolors.append( np.array( [0.25,0.75,0.25] ) )
        mcolors.append( np.array( [0.75,0.25,0.25] ) )
        mcolors.append( np.array( [0.25,0.75,0.75] ) )
        mcolors.append( np.array( [0.75,0.25,0.75] ) )
        mcolors.append( np.array( [0.75,0.75,0.25] ) )
        # Get one color per stationary point
        assignation = {}
        color = random.choice(mcolors)
        for cpname in cpnames:
            while True:
               # direction of variation in color
               dcx = random.choice([-1,+1]) * random.random()
               dcy = random.choice([-1,+1]) * random.random()
               dcz = random.choice([-1,+1]) * random.random()
               dcolor = np.array([dcx,dcy,dcz])
               dcolor = dcolor / np.linalg.norm(dcolor)
               # amount in variation
               magnitude = random.choice([0.25,0.50,0.75])
               # new color
               color = color + magnitude*dcolor
               # Valid color?
               if color[0] < 0.0 or color[0] > 1.0: color = random.choice(mcolors); continue
               if color[1] < 0.0 or color[1] > 1.0: color = random.choice(mcolors); continue
               if color[2] < 0.0 or color[2] > 1.0: color = random.choice(mcolors); continue
               break
            # Associate color to stationary point
            assignation[cpname] = list(color)
        return assignation
    #--------------------------------

    #-----------#
    # Libraries #
    #-----------#
    import pylab
    from   matplotlib        import cm
    from   matplotlib.backends.backend_pdf import PdfPages
    #from   scipy.spatial import Voronoi

    #-------------------------#
    # Variables for the plots #
    #-------------------------#
    dpi       = 250
    ncl       = 40 # number of contour lines
    CMAP      = cm.bwr # options: (cm.plasma, cm.bwr)
    cpmarker  = {0:"co",1:"g^",2:"ms"}
    xlabel    = r"$\phi_1$ (degrees)"
    ylabel    = r"$\phi_2$ (degrees)"
    xsize     = 16
    ysize     = 16
    axis      = [0.0 , 360., 0.0 , 360.]
    ticks     = [float(i) for i in range(0,360+1,60)]
    plots     = { i:False for i in range(6)}
    cbarlabel = r"cm$^{-1}$"
    title     = {}
    title[0]  = "Numerical PES"
    title[1]  = "Numerical PES"
    title[2]  = "Numerical PES"
    title[3]  = "Fourier PES"
    title[4]  = "Delaunay Teselation"
    title[5]  = "Delaunay Teselation"
    titlesize = 22

    #-----------------#
    # Read Input file #
    #-----------------#
    input_data   = readfile_inp(ifile)
    torsions     = input_data[0]
    calcs        = input_data[1]
    pes          = input_data[2]
    fourier      = input_data[3]
    statpoints      = input_data[4]
    tor2dns      = input_data[5]
    rovibpf      = input_data[6]
    Tlist        = input_data[7]
    iofiles      = input_data[8]
    # Expand each tuple
    (torsion1,torsion2,tsigma1,tsigma2) = torsions
    (ttype,level,charge,multiplicity) = calcs
    (t1step,t2step,symmetry) = pes
    (fterms,nums,weight,ignore) = fourier
    (tolerance,freqscal) = statpoints
    (dijvar,kmax,maxeigen) = tor2dns
    (interpolation,integrationstep) = rovibpf
    (f_xyz,f_calcs,f_ics,f_pes,f_vfit,f_dfit,f_splist,f_spinfo,f_evals,f_sqrt,f_tables,f_pdf,f_out) = iofiles
    # Add folder name to files
    f_ics     = DIRFILES + f_ics
    f_pes     = DIRFILES + f_pes
    f_vfit    = DIRFILES + f_vfit
    f_dfit    = DIRFILES + f_dfit
    f_splist = DIRFILES + f_splist
    f_spinfo = DIRFILES + f_spinfo
    f_evals   = DIRFILES + f_evals
    f_sqrt    = DIRFILES + f_sqrt
    f_tables  = DIRFILES + f_tables

    #-------------------#
    # Read Q2DTor files #
    #-------------------#
    # Read PES
    if os.path.exists(f_pes):
       dict_pes, symbols, Eref = readfile_pes(f_pes)
       pes_calc = dict_pes.keys()
       dict_pes = pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2)
       pes_all  = dict_pes.keys()
       plots[0] = True
       plots[1] = True
       plots[2] = True
    # Read Fourier fitting
    if os.path.exists(f_vfit):
       fterms, parameters, Eref_fourier = readfile_xfit(f_vfit)
       Fourier_PES = Fourier2D(fterms)
       Fourier_PES.set_coefs(parameters)
       plots[3] = True
    # Read critical points
    if os.path.exists(f_splist):
       dict_CPs, Eref_cps  = readfile_splist(f_splist)
       plots[4] = True
      #plots[5] = True
    else:
       dict_CPs = {}
       
    #------------------------------#
    # Critical points: replication #
    #------------------------------#
    # Set color for the region of each CP
    colors = sp_colors(dict_CPs.keys())
    for cp in colors.keys(): dict_CPs[cp] = dict_CPs[cp] + [colors[cp]]
    # Symmetry replication
    dict_CPs, list_CPs, initial_CPs = sp_allCPs(dict_CPs,symmetry,name,tsigma1,tsigma2)
    # Delaunay tesselation
    if dict_CPs != {}:
       dict_CPs     = tess.replicate_points(name,dict_CPs)
       Dlist_CPs    = sp_sorting(dict_CPs)
       Dlist_points = [dict_CPs[cp_name][0:2] for cp_name in Dlist_CPs]
       delaunay     = Delaunay(Dlist_points)
       Vmidpoints   = tess.get_Vmidpoints(delaunay,Fourier_PES)

    #-----------------#
    # Create pdf file #
    #-----------------#
    print "     Storing plots in file: %s"%f_pdf
    pp       = PdfPages(f_pdf)

    if plots[0]:
       #---------------#
       # Plot number 0 #
       #---------------#
       plotnum = 0
       fig     = pylab.figure(plotnum)
       # numerical PES
       xx = np.array([tuple(dict_pes[point][0:2]) for point in pes_all])*cons.R2D
       zz = np.array([(dict_pes[point][2]-Eref)*cons.h2c for point in pes_all])
       pylab.tricontourf(xx[:,0],xx[:,1],zz,ncl,alpha=0.75,cmap=CMAP,antialiased=True,zorder=1)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)

    if plots[1]:
       #---------------#
       # Plot number 1 #
       #---------------#
       plotnum = 1
       fig     = pylab.figure(plotnum)
       # numerical PES
       xx = np.array([tuple(dict_pes[point][0:2]) for point in pes_all])*cons.R2D
       zz = np.array([(dict_pes[point][2]-Eref)*cons.h2c for point in pes_all])
       pylab.tricontourf(xx[:,0],xx[:,1],zz,ncl,alpha=0.75,cmap=CMAP,antialiased=True,zorder=1)
       # calculated points
       for point in pes_calc:
           phi1, phi2 = dict_pes[point][0:2]
           pylab.plot([phi1*cons.R2D],[phi2*cons.R2D],'kx',alpha=0.50,zorder=2)
       # basic stationary points
       for cp_name in initial_CPs:
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=10, markeredgewidth=2,zorder=3)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)

    if plots[2]:
       #---------------#
       # Plot number 2 #
       #---------------#
       plotnum = 2
       fig     = pylab.figure(plotnum)
       # numerical PES
       xx = np.array([tuple(dict_pes[point][0:2]) for point in pes_all])*cons.R2D
       zz = np.array([(dict_pes[point][2]-Eref)*cons.h2c for point in pes_all])
       pylab.tricontourf(xx[:,0],xx[:,1],zz,ncl,alpha=0.75,cmap=CMAP,antialiased=True,zorder=1)
       # calculated points
       for point in pes_calc:
           phi1, phi2 = dict_pes[point][0:2]
           pylab.plot([phi1*cons.R2D],[phi2*cons.R2D],'kx',alpha=0.50,zorder=2)
       # all stationary points
       for cp_name in list_CPs:
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=8, markeredgewidth=1,alpha=0.50,zorder=3)
       # basic stationary points
       for cp_name in initial_CPs:
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=10, markeredgewidth=2,zorder=4)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)

    if plots[3]:
       #---------------#
       # Plot number 3 #
       #---------------#
       plotnum = 3
       fig     = pylab.figure(plotnum)
       # Plot Fourier PES
       Y, X = np.mgrid[2*np.pi:0.0:100j, 2*np.pi:0.0:100j]
       Z = Fourier_PES.value(X,Y)
       pylab.contourf(X*cons.R2D,Y*cons.R2D,Z,ncl, alpha=0.75, cmap=CMAP,zorder=1)
       # all stationary points
       for cp_name in list_CPs:
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=10, markeredgewidth=2,zorder=2)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)

    if plots[4]:
       #---------------#
       # Plot number 4 #
       #---------------#
       plotnum = 4
       fig     = pylab.figure(plotnum)
       # Plot Fourier PES
       Y, X = np.mgrid[2*np.pi:0.0:100j, 2*np.pi:0.0:100j]
       Z = Fourier_PES.value(X,Y)
       pylab.contourf(X*cons.R2D,Y*cons.R2D,Z,ncl, alpha=0.75, cmap=CMAP,zorder=1)
       # Plot Teselation
       plotmode = 1
       if plotmode == 1:
          for triangle in delaunay.simplices:
              p1, p2, p3 = [delaunay.points[ii]*cons.R2D for ii in triangle]
              pylab.plot( [p1[0],p2[0],p3[0],p1[0]] , [p1[1],p2[1],p3[1],p1[1]] , 'k-', zorder=2)
       if plotmode == 2:
           pylab.triplot(points[:,0]*cons.R2D, points[:,1]*cons.R2D, delaunay.simplices.copy(),'-k',lw=1.5,zorder=2)
       # all stationary points (including those generated for Delaunay)
       for cp_name in dict_CPs.keys():
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=10, markeredgewidth=2,zorder=3)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)

    if plots[5]:
       #---------------#
       # Plot number 5 #
       #---------------#
       plotnum = 5
       fig     = pylab.figure(plotnum)
       # Plot Fourier PES
       Y, X = np.mgrid[2*np.pi:0.0:100j, 2*np.pi:0.0:100j]
       Z = Fourier_PES.value(X,Y)
       pylab.contour(X*cons.R2D,Y*cons.R2D,Z,ncl, alpha=0.75, cmap='Greys',zorder=2)
       # Plot Teselation
       plotmode = 1
       if plotmode == 1:
          for triangle in delaunay.simplices:
              p1, p2, p3 = [delaunay.points[ii]*cons.R2D for ii in triangle]
              pylab.plot( [p1[0],p2[0],p3[0],p1[0]] , [p1[1],p2[1],p3[1],p1[1]] , 'k-', zorder=2)
       if plotmode == 2:
           pylab.triplot(points[:,0]*cons.R2D, points[:,1]*cons.R2D, delaunay.simplices.copy(),'-k',lw=1.5,zorder=2)
       # Get regions
       for xp in range(X.shape[0]):
           for yp in range(X.shape[1]):
               phi1, phi2 = X[xp,yp], Y[xp,yp]
               p0     = np.array([phi1,phi2])
               idx    = tess.which_statpoint(p0,delaunay,method=TESSMODEPLOT,args=[Fourier_PES,Vmidpoints])
               cpname = Dlist_CPs[idx]
               color  = dict_CPs[cpname][-1]
               pylab.plot( [phi1*cons.R2D],[phi2*cons.R2D],'o',color=color,mec=color,alpha=1.00,zorder=1)
       # all stationary points (including those generated for Delaunay)
       for cp_name in dict_CPs.keys():
           phi1   = dict_CPs[cp_name][0]*cons.R2D
           phi2   = dict_CPs[cp_name][1]*cons.R2D
           cptype = dict_CPs[cp_name][2]
           pylab.plot([phi1],[phi2],cpmarker[cptype], markersize=10, markeredgewidth=2,zorder=4)
       # Format
       pylab.title(title[plotnum],fontsize=titlesize)
       pylab.xlabel(xlabel,fontsize=xsize)
       pylab.ylabel(ylabel,fontsize=ysize)
       pylab.axis(axis)
       pylab.xticks(ticks)
       pylab.yticks(ticks)
       pylab.colorbar().ax.set_title(cbarlabel)
       # Save figure and close
       pp.savefig(fig,dpi=dpi) ; pylab.close(fig)


    pylab.close("all")
    pp.close()
    print
#----------------------------------------------------------#
def tool_gts(mainname,software):
    '''
    Function to create a gts file from files coming from other sofware
    '''
    software = software.lower()
    print "     Software: '%s'"%software

    #---------#
    if   software == "gaussian":
         ff = [mainname+"."+ext for ext in ["fchk"]]
         print "     Looking for files:"
         found = []
         for f in ff:
             if os.path.isfile(f):  print "         %s  (yes)"%f; found.append(True )
             else                :  print "         %s  (no)"%f ; found.append(False)
         # Read fchk file
         if found[0]:
            from mq2dtor.mesc_gaussian import readfchk
            xcc , atonums, ch , mtp ,  E , gcc, Fcc = readfchk(ff[0])
         else:
            sys.exit("     Problems reading file: %s\n"%ff[0])
    #---------#
    elif software == "orca":
         ff = [mainname+"."+ext for ext in ["out","engrad","hess"]]
         print "     Looking for files:"
         found = []
         for f in ff:
             if os.path.isfile(f):  print "         (yes) %s"%f; found.append(True )
             else                :  print "         (no)  %s"%f ; found.append(False)
         # Read output
         if found[0]:
            from mq2dtor.mesc_orca import read_output
            xcc, symbols, atonums, ch, mtp, E = read_output(ff[0])
         else:
            print
            sys.exit("     Problems reading file: %s\n"%ff[0])
         # Read engrad file
         if found[1]:
            from mq2dtor.mesc_orca import read_engrad
            gcc = read_engrad(ff[1])[0]
         else:
            gcc = None
         # Read hess file
         if found[2]:
            from mq2dtor.mesc_orca import read_hess
            Fcc = read_hess(ff[2])
         else:
            Fcc = None
    #---------#
    else: sys.exit("     Error: software not implemented\n")
    #---------#

    masslist = [cons.dict_atomasses[atonum] for atonum in atonums]
    print "     Calculatin point group:"
    pgroup, rotsigma = hf.get_pgs(atonums,masslist,xcc)
    print "        * %s"%pgroup
    print "     Creating file: %s"%(mainname+".gts")
    write_gtsfile(xcc , atonums , ch, mtp, E, pgroup, rotsigma , gcc , Fcc, mainname+".gts")
    # Write molden
    structure = gts2Struct(mainname+".gts")
    if structure.get("Fcc") is not None:
       structure.basic_setups([0,3])
       structure.calc_ccfreqs()
       structure.write_molden(mainname+".molden")
       print "     Creating file: %s"%(mainname+".molden")
    print
#----------------------------------------------------------#
def tool_icoords(name,mode,f_ics):
    #-------------------------------------------------------------
    def gtsicoords(f_gts, f_ics):
        # files?
        if not os.path.isfile(f_gts): sys.exit("  Unable to find file: %s\n"%f_gts)
        if not os.path.isfile(f_ics): sys.exit("  Unable to find file: %s\n"%f_ics)
        # read files
        print "  Reading file: %s"%f_gts
        struct = gts2Struct(f_gts,"")
        nvib = struct.get("nvib")
        print "    - number of atoms: %i"%(struct.get("natoms"))
        print "    - vibrational dof: %i"%nvib
        print "  Reading file: %s"%f_ics
        icoords, nics, bonds = readfile_ics(f_ics)
        print "    - number of ics  : %i"%nics
        # Calculate cc-frequencies
        print "  Calculating cc-frequencies"
        if struct.get("Fcc") is None: sys.exit("  Unable to find hessian matrix")
        struct.basic_setups([0,3])
        struct.calc_ccfreqs()
        for freq in struct.get("ccfreqs"): print "    %s cm^-1"%str(freq)
        # Calculate ic-frequencies
        print "  Calculating ic-frequencies"
        same_values = struct.check_icoords(icoords)
        for freq in struct.get("icfreqs"): print "    %s cm^-1"%str(freq)
        
        if same_values: print "  The set of internal coordinates is: VALID"
        else          : print "  The set of internal coordinates is: NOT VALID"

        return struct, icoords, nvib, nics, same_values
    #-------------------------------------------------------------

    if mode == "gts":
       if f_ics is None: sys.exit("  ERROR: define the file with internal coordinates!\n")
       f_gts = name+".gts"
       struct, icoords, nvib, nics, valid = gtsicoords(f_gts,f_ics)
       print
       if valid and nics > nvib:
         print "  Trying to debug internal coordinates"
         icoords, nics = struct.purify_ricoords(icoords,show=True,rbonds=True)
         print "  This set (%i) would be also valid:"%nics
         for ictype, ic in icoords:
             if ictype == "3": print "       "+"=".join([str(ii+1) for ii in ic])
             else            : print "       "+"-".join([str(ii+1) for ii in ic])

    elif mode == "sp":
       iofiles   = readfile_inp(name+".inp")[8]
       f_splist = DIRFILES + iofiles[6]
       if f_ics is None: f_ics = DIRFILES + iofiles[2]
       # read splist
       dict_CPs, Eref_cp = readfile_splist(f_splist)
       for sp in dict_CPs.keys():
           print " ----------------------------------" + "-"*len(sp)
           print "  Testing internal coordinates in: %s"%sp
           print " ----------------------------------" + "-"*len(sp)
           f_gts = DIRSP + sp + ".gts"
           gtsicoords(f_gts, f_ics)
           print " ----------------------------------" + "-"*len(sp)
           print
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 5) Helper functions related to: symmetry       #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def pessym_inregion(phi1, phi2, symmetry="", tsigma1=1, tsigma2=1):
    '''
    Given a point, it indicates if it
    belongs to the non-redundant region
    (that what is used to replicate using
    symmetry)
    '''

    # Define 2pi limits
    if "c" not in symmetry:
       mphi1 = cons.TWOPI
       mphi2 = cons.TWOPI
    if "c"     in symmetry:
       mphi1 = cons.TWOPI/tsigma1
       mphi2 = cons.TWOPI/tsigma2
    m12 = float(tsigma1)/float(tsigma2)

    # Check 2pi and (c) conditions
    condC1 = mphi1 - phi1 > cons.ZERO_RAD
    condC2 = mphi2 - phi2 > cons.ZERO_RAD
    condC  = condC1 and condC2
    if not condC: return False

    # Check condition A
    if "a" in symmetry:
       condA = phi1-phi2 > -cons.ZERO_RAD
       if not condA: return False

    # Check condition B
    if "b" in symmetry:
       condB1 = phi1 - mphi1/2.0 > -cons.ZERO_RAD
       if condB1:
          condB = phi2 - m12*(mphi1-phi1) < +cons.ZERO_RAD
       else:
          condB = phi2 - m12*(mphi1-phi1) < -cons.ZERO_RAD
       if not condB: return False

    # Point is inside
    return True
#----------------------------------------------------------#
def pessym_apply(phi1, phi2, symmetry="", tsigma1=1, tsigma2=1):
    '''
    Replicate a point using the PES symmetry
    inside the [0,2pi) range
    The returned list does not included the given point
    '''

    #----------------------------#
    # (0) Point in [0,2pi/sigma) #
    #----------------------------#
    if "c" in symmetry:
       maxphi1 = cons.TWOPI/float(tsigma1)
       maxphi2 = cons.TWOPI/float(tsigma2)
    else:
       maxphi1 = cons.TWOPI
       maxphi2 = cons.TWOPI

    while phi1 <      0.0: phi1 += maxphi1
    while phi2 <      0.0: phi2 += maxphi2
    while phi1 >= maxphi1: phi1 -= maxphi1
    while phi2 >= maxphi2: phi2 -= maxphi2

    points = [ (phi1,phi2) ]

    #-------------------------#
    # (1) Apply condition (a) #
    #-------------------------#
    if "a" in symmetry:
       toadd = []
       for phi1,phi2 in points:
           if abs(phi1-phi2) < cons.ZERO_RAD: continue
           newphi1 = phi2
           newphi2 = phi1
           toadd.append( (newphi1,newphi2) )
       points += toadd

    #-------------------------#
    # (2) Apply condition (b) #
    #-------------------------#
    if "b" in symmetry:
       toadd = []
       for phi1,phi2 in points:
           newphi1   = maxphi1 - phi1
           newphi2   = maxphi2 - phi2
           toadd.append( (newphi1,newphi2) )
       points += toadd

    #-------------------------#
    # (3) Apply condition (c) #
    #-------------------------#
    if "c" in symmetry:
       toadd = []
       for phi1,phi2 in points:
           for s1 in range(tsigma1):
               for s2 in range(tsigma2):
                   if s1 == 0 and s2 == 0: continue
                   newphi1 = phi1 + s1*maxphi1
                   newphi2 = phi2 + s2*maxphi2
                   toadd.append( (newphi1,newphi2) )
       points += toadd

    #-------------------------#
    # (4) Remove 2pi points   #
    #-------------------------#
    final_list = []
    for phi1,phi2 in points:
        if abs(phi1-cons.TWOPI) < cons.ZERO_RAD: continue
        if abs(phi2-cons.TWOPI) < cons.ZERO_RAD: continue
        final_list.append( (phi1,phi2) )

    return final_list
#----------------------------------------------------------#
def pessym_points(dphi1,dphi2,symmetry="",tsigma1=1,tsigma2=1):

    if "c" in symmetry:
       mphi1 = cons.TWOPI/tsigma1
       mphi2 = cons.TWOPI/tsigma2
    else:
       mphi1 = cons.TWOPI
       mphi2 = cons.TWOPI

    # Basic 1D list of points
    list_phi1 = hf.frange(0.0,mphi1,dphi1,False)
    list_phi2 = hf.frange(0.0,mphi2,dphi2,False)

    # Generate list of all points
    points = []
    for phi1 in list_phi1:
        # Select points
        for phi2 in list_phi2:
            add = pessym_inregion(phi1,phi2,symmetry,tsigma1,tsigma2)
            if add: points.append( (phi1,phi2) )
        # Reverse order
        list_phi2 = list_phi2[::-1]

    # Convert into array and return
    points = np.array(points)
    return points
#----------------------------------------------------------#
def pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2,pprint=False):

    dict_whole = {}

    np1 = len(dict_pes.keys())
    # Replicate point due to PES symmetry
    for pname1 in dict_pes.keys():
        phi1,phi2 = dict_pes[pname1][0:2]
        toadd     = pessym_apply(phi1, phi2, symmetry, tsigma1, tsigma2)
        for p1,p2 in toadd:
            pname2 = string_getpname(p1,p2,name,"rad")
            # Add point to new dict
            dict_whole[pname2] = [p1,p2]+dict_pes[pname1][2:]
        # Empty original dict
        dict_pes[pname1]   = None
    np2 = len(dict_whole.keys())
    # Replicate point due to 2pi symmetry
    for pname1 in dict_whole.keys():
        phi1,phi2 = dict_whole[pname1][0:2]
        if phi1 == 0.0:
           p1, p2 = cons.TWOPI, phi2
           pname2 = string_getpname(p1,p2,name,"rad")
           dict_whole[pname2] = [p1,p2]+dict_whole[pname1][2:]
        if phi2 == 0.0:
           p1, p2 = phi1, cons.TWOPI
           pname2 = string_getpname(p1,p2,name,"rad")
           dict_whole[pname2] = [p1,p2]+dict_whole[pname1][2:]
        if phi1 == 0.0 and phi2 == 0.0:
           p1, p2 = cons.TWOPI, cons.TWOPI
           pname2 = string_getpname(p1,p2,name,"rad")
           dict_whole[pname2] = [p1,p2]+dict_whole[pname1][2:]
    np3 = len(dict_whole.keys())

    if pprint:
       print "        * PES contains %4i points (calculated ones)"%np1
       print "        * PES contains %4i points after applying symmetry conditions"%np2
       print "        * PES contains %4i points after applying 2pi-periodicity\n"%np3

    return dict_whole
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>>> ### <<<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                             #
# SECTION ( 6) Helper functions related to: stationary points #
#                                                             #
# >>>>>>>>>>>>>>>>>>>>>>>>>>> ### <<<<<<<<<<<<<<<<<<<<<<<<<<< #
def sp_optimization(cp_name,symbols,dict_CPs,dict_pes,functionOPT,ilines):
        
   #phi1    = dict_CPs[cp_name][0]
   #phi2    = dict_CPs[cp_name][1]
    phi1, phi2 = [float(angle)*cons.D2R for angle in cp_name.split("_")[-2:]]
    cptype  = dict_CPs[cp_name][2]

    # Perform calculation
    try:
       # Select closest geometry from scan calculations
       target    = None
       min_norm2 = +float("inf")
       for pes_point in dict_pes.keys():
           ref_phi1, ref_phi2 = [float(angle)*cons.D2R for angle in pes_point.split("_")[-2:]]
           dphi1 = hf.angle_diff(phi1,ref_phi1)
           dphi2 = hf.angle_diff(phi2,ref_phi2)
           norm2 = dphi1**2 + dphi2**2
           if norm2 < min_norm2:
              min_norm2 = norm2
              target = pes_point
       xvec = dict_pes[target][3]

       # Adjust dihedral angles


       # Perfom calculation
       phis  = (phi1*cons.R2D,phi2*cons.R2D)
       odata = functionOPT(ilines,xvec,symbols,cp_name,phis)
       if odata is None: return None

       # Create gtsfile
       gtsfile = DIRSP+cp_name+".gts"
       xcc, atonums, ch, mtp,  Etot, gcc, Fcc = odata
       cp_energy = Etot
       masslist = [cons.dict_atomasses[atonum] for atonum in atonums]
       pgroup, rotsigma = hf.get_pgs(atonums,masslist,xcc)
       gts_tuple = (xcc, atonums, ch, mtp, Etot, pgroup, rotsigma, gcc,Fcc,gtsfile)
       print "        * Creating file '%s' with SP information"%(gtsfile)
       write_gtsfile(*gts_tuple)

    except: pass
#----------------------------------------------------------#
def sp_checking(gtsfile,cp_name,cptype,torsion1,torsion2,tolerance,dict_warns,symmetry,tsigma1,tsigma2,ttype):
    # Create structure
    structure = gts2Struct(gtsfile,name=cp_name)
    structure.basic_setups()
    if structure.get("Fcc") is None:
       dict_warns[cp_name].append(2)
       return structure.get("Etot"), dict_warns
    structure.calc_ccfreqs()

    # Check dihedral angles
    iphi1, iphi2 = cp_name.split("_")[-2:]
    iphi1 = float(iphi1)*cons.D2R
    iphi2 = float(iphi2)*cons.D2R
    phi1  = hf.ranged_angle(structure.icvalue(torsion1),v=1)
    phi2  = hf.ranged_angle(structure.icvalue(torsion2),v=1)
    dphi1 = hf.angle_diff(iphi1,phi1)
    dphi2 = hf.angle_diff(iphi2,phi2)
    diff  = np.sqrt(dphi1**2 + dphi2**2)
    if diff > 4.0*tolerance:
       dict_warns[cp_name].append(3)

    # Check number of imaginary frequencies
    if ttype == "min": expected = cptype
    if ttype == "ts" : expected = cptype + 1
    nimag = 0
    for ccfreq in structure.get("ccfreqs"):
        if ccfreq.isItImag(): nimag += 1
    if nimag != expected: dict_warns[cp_name].append(4)

    # Generating molden file with freqs
    molden  = DIRSP+cp_name+".molden"
    print "        * Creating Molden file: '%s'"%(molden)
    structure.write_molden(molden)

    return structure.get("Etot"), dict_warns
#----------------------------------------------------------#
def sp_sorting(dict_CPs):
    E_names = [ (dict_CPs[cpname][3],cpname) for cpname in dict_CPs.keys()]
    E_names.sort()
    list_CPs = [cpname for E,cpname in E_names]
    return list_CPs
#----------------------------------------------------------#
def sp_icdepuration(dict_CPs,torsion1,torsion2,xyzStruct,f_ics):

    list_CPs = sp_sorting(dict_CPs)

    #-----------------------------------------------#
    # Check that there exists at least one gts file #
    #-----------------------------------------------#
    list_CPs = [cpname for cpname in list_CPs if os.path.exists(dict_CPs[cpname][4])]
    if len(list_CPs) == 0:
        print "     ERROR: Not a single gts file... Cannot continue...\n"; sys.exit(EXITMESS)

    #-----------------------#
    # Is there an ics file? #
    #-----------------------#
    print "     Generation of non-redundant internal coordinates"
    if not os.path.isfile(f_ics):
       # Generate ricoords
       xyzStruct.graph_autoconnect(lengthfactor=1.3)
       xyzStruct.graph_fragconnect()
       bonds = xyzStruct.graph_getbonds()
       icoords, nics = xyzStruct.gen_ricoords(torsions=[torsion1,torsion2],check=False)
       # Write ricoords
       print "        - Internal coordinates generated and stored at '%s' "%f_ics
       writefile_ics(icoords,f_ics,bonds)
    else:
       print "        - File '%s' already exists"%f_ics
       icoords, nics, bonds = readfile_ics(f_ics)
       print "          number of internal coordinates: %i"%(nics)
    # Copy original
    f_ics_original = str(f_ics)+"_original"
    if not os.path.exists(f_ics_original):
       print "        - Creating copy of original set: %s"%f_ics_original
       shutil.copyfile(f_ics, f_ics_original)
    print

    #------------------------------------#
    # Check that torsions are in icoords #
    #------------------------------------#
    icoords = ic_checktorsions(icoords,torsion1,torsion2)

    #--------------------------------------#
    # Depurating icoords and checking them #
    #--------------------------------------#
    torsions = [torsion1,torsion2]
    if xyzStruct.get("nvib") < nics: targets = list(list_CPs)
    else                           : targets = []
    depured      = False
    icoords_copy = list(icoords)
    EMIN = float("inf")
    masslist = xyzStruct.get("masslist")

    while True:
      #---------------------------------------#
      # Depure using a given stationary point #
      #---------------------------------------#
      if targets != []:
         # Select an stationary point
         idx     = random.choice( range(len(targets)) )
         cpname  = targets.pop(idx)
         cptype  = dict_CPs[cpname][2]
         gtsfile = dict_CPs[cpname][4]

         # Depure icoords
         print "        - Depurating using '%s' critical point... "%cpname
         struct = gts2Struct(gtsfile,cpname,masslist,stype=cptype)
         if struct.get("Etot") < EMIN: EMIN = struct.get("Etot")
         struct.basic_setups([0,3])
         struct.calc_ccfreqs()
         icoords, nics = struct.purify_ricoords(icoords,torsions=torsions,show=True)
      else:
         cpname = None
         print "        - Skipping step: number of icoords = number of vib d.o.f."
      print

      #------------------------------#
      # Check in the rest of systems #
      #------------------------------#
      print "        - Checking set of non-redundant internal coordinates"
      num0, num1, num2 = 1, 1, 1
      for cpname2 in list_CPs:
          if cpname2 == cpname: continue
          cptype2  = dict_CPs[cpname2][2]
          gtsfile2 = dict_CPs[cpname2][4]
          struct2 = gts2Struct(gtsfile2,cpname2,masslist,stype=cptype2)
          if struct2.get("Etot") < EMIN: EMIN = struct2.get("Etot")
          struct2.basic_setups([0,3])
          struct2.calc_ccfreqs()
          try:
             depured = struct2.check_icoords(icoords)
          except:
             depured = False
          if cptype2 == 0: xxx = "MIN%02i"%num0; num0 += 1
          if cptype2 == 1: xxx = "TS_%02i"%num1; num1 += 1
          if cptype2 == 2: xxx = "MAX%02i"%num2; num2 += 1
          if not depured:
             print "             * %s --> %s --> fails!"%(xxx,cpname2)
             break
          else:
             print "             * %s --> %s --> OK"%(xxx,cpname2)
      print 

      # Finished?
      if len(targets) == 0: break
      if depured          : break
      icoords = list(icoords_copy)

    #-----------#
    # Success?? #
    #-----------#
    if not depured:
      print "        - Failure!!"
      print "          Removing file '%s'"%f_ics
      os.remove(f_ics)
      print
      sys.exit(EXITMESS)

    print "        - Success!!"
    print
    print "        - Storing non-redundant set at: '%s'"%f_ics
    writefile_ics(icoords,f_ics,bonds)
    print
    return icoords
#----------------------------------------------------------#
def sp_allCPs(dict_CPs,symmetry,name,tsigma1,tsigma2):

    MAX_dE   = 1.0 # in cm^-1
    MAX_dist = 2.5 # in degrees
    initial = dict_CPs.keys()

    #---------------------------------#
    # Check redundancy in initial set #
    #---------------------------------#
    initial = []
    for cpname in dict_CPs.keys():
        phi1, phi2, cptype, E = dict_CPs[cpname][0:4]
        keep = True
        for basic in initial:
            bphi1, bphi2, bcptype, bE = dict_CPs[basic][0:4]
            # Compare type
            if cptype != bcptype: continue
            # Compare energy
            if abs(E-bE) > MAX_dE: continue
            # Compare angles
            diff1 = hf.angle_diff(phi1,bphi1)
            diff2 = hf.angle_diff(phi2,bphi2)
            dist  = np.sqrt(diff1**2 + diff2**2)
            if dist*cons.R2D > MAX_dist: continue
            # The point is repeated
            keep = False
            break
        if keep: initial.append(cpname)
        else   : dict_CPs.pop(cpname)

    #--------------------------#
    # Apply symmetry condition #
    #--------------------------#
    for cpname in initial:
        phi1, phi2  = dict_CPs[cpname][0:2]
        equivalents = pessym_apply(phi1,phi2,symmetry,tsigma1,tsigma2)
        for p1, p2 in equivalents:
            # skip same point
            diff = np.array([phi1,phi2])-np.array([p1,p2])
            if np.linalg.norm(diff) < cons.ZERO_RAD: continue
            # Save point
            cpname2 = string_getpname(p1,p2,name,"rad")
            dict_CPs[cpname2] = [p1,p2] + dict_CPs[cpname][2:]

    #-------------#
    # Check again #
    #-------------#
    close = set([])
    for cpnameA in dict_CPs.keys():
        for cpnameB in dict_CPs.keys():
            if cpnameA == cpnameB: continue
            Aphi1, Aphi2, Acptype, AE = dict_CPs[cpnameA][0:4]
            Bphi1, Bphi2, Bcptype, BE = dict_CPs[cpnameB][0:4]
            if Acptype != Bcptype: continue
            if abs(AE-BE) > MAX_dE: continue
            diff1 = hf.angle_diff(Aphi1,Bphi1)
            diff2 = hf.angle_diff(Aphi2,Bphi2)
            dist  = np.sqrt(diff1**2 + diff2**2)
            if dist*cons.R2D > MAX_dist: continue
            close.add( tuple(sorted([cpnameA,cpnameB])) )
    for cpnameA , cpnameB in close:
        if   cpnameA in initial: dict_CPs.pop(cpnameB)
        elif cpnameB in initial: dict_CPs.pop(cpnameA)
        else                   : dict_CPs.pop( random.choice([cpnameA,cpnameB]) )

    # Return data
    list_CPs = sp_sorting(dict_CPs)
    return dict_CPs, list_CPs, initial
#----------------------------------------------------------#
def sp_analysis(gtsfile,icoords,cptype,freqscal=1.0,masslist=None,name=""):
    # Generate structure
    struct  = gts2Struct(gtsfile,name,masslist,stype=cptype)
    struct.set("freqscal",freqscal)
    if masslist is None: struct.basic_setups([0,  2,3,4])
    else               : struct.basic_setups([0,1,2,3,4])
    struct.calc_ccfreqs()
    # Energy
    energy = struct.get("Etot")
    # D matrix
    Dmatrix = struct.get_Dmatrix(icoords)
    rI1 =  Dmatrix[0,0]
    rI2 =  Dmatrix[1,1]
    l12 = -Dmatrix[0,1]
    detD = rI1*rI2 - l12*l12
    d11  = 100*(rI2/detD)
    d22  = 100*(rI1/detD)
    d12  = 100*(l12/detD)
    # All vibrational freqs
    struct.calc_ccfreqs()
    ccfreqs = struct.get("ccfreqs")
    zpe     = sum([freq.get("zpe") for freq in ccfreqs])
    # Non-torsional freqs
    struct.calc_icfreqs(icoords[:-2])
    ntfreqs = struct.get("icfreqs")
    zpe_2b     = sum([freq.get("zpe") for freq in ntfreqs])
    # Return data
    rij = (rI1,rI2,l12)
    dij = (d11,d22,d12)
    data_sp = (struct,energy,rij,dij,ccfreqs,zpe,ntfreqs,zpe_2b)
    return data_sp
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 7) Helper functions related to: 2DNS           #
#              These functions were written by A.F.R.      #
#              and re-adapted by D.F.C.                    #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def tor2DNS_hamiltonian(molname, V_terms, V_coefs, LkinD, dijtuple, coef_nums, kmax=100):

    '''
    This function calculates matrix
    elements of the type <m,j|H|n,k>

    In general, the matrix is really large.
    Therefore, data is returned as a sparse matrix
    * Elements for phi1 are <m|n>
    * Elements for phi2 are <j|k>
    MJ is the row of a given element
    NK is the column of a given element

    Input:
      * hdc
      * hd1
      * hd2
      * LkinD
      * nv1d=201 # odd to include zero -100-->0-->100

    Returns:
      * fil
      * col
      * ham2d
    '''

    nv1d = 2*kmax + 1
    npar1, npar2, ncoup, npar1c, npar2c, ncoupc = coef_nums

    #--------------------------------------#
    # Format change: from D.F.C. to A.F.R. #
    #--------------------------------------#
    # Initialize arrays for Real coefs
    ap  = np.zeros(npar1 )            # 1D coef. associated with cos(x)
    bp  = np.zeros(npar2 )            # 1D coef. associated with cos(y)
    cp  = np.zeros( (ncoup ,ncoup ) ) # 2D coef. associated with cos(x)*cos(y)
    dp  = np.zeros( (ncoup ,ncoup ) ) # 2D coef. associated with sin(x)*sin(y)
    # Initialize arrays for Complex coefs
    apc = np.zeros(npar1c)            # 1D coef. associated with sin(x)
    bpc = np.zeros(npar2c)            # 1D coef. associated with sin(y)
    cpc = np.zeros( (ncoupc,ncoupc) ) # 2D coef. associated with cos(x)*sin(y)
    dpc = np.zeros( (ncoupc,ncoupc) ) # 2D coef. associated with sin(x)*cos(y)

    for V_term, V_coef in zip(V_terms,V_coefs):
        term_type,idx_phi1,idx_phi2 = V_term

        if term_type == "const": ap0 = V_coef

        if term_type == "cos" and idx_phi2 == "-": ap[idx_phi1-1] = V_coef
        if term_type == "cos" and idx_phi1 == "-": bp[idx_phi2-1] = V_coef
        if term_type == "coscos":                  cp[idx_phi1-1][idx_phi2-1] = V_coef
        if term_type == "sinsin":                  dp[idx_phi1-1][idx_phi2-1] = V_coef

        if term_type == "sin" and idx_phi2 == "-": apc[idx_phi1-1] = V_coef
        if term_type == "sin" and idx_phi1 == "-": bpc[idx_phi2-1] = V_coef
        if term_type == "sincos":                  cpc[idx_phi1-1][idx_phi2-1] = V_coef
        if term_type == "cossin":                  dpc[idx_phi1-1][idx_phi2-1] = V_coef

    #-----------------------#
    # Define some constants #
    #-----------------------#
    CDIS=0.02966
    UMAKG=1.661E-27
    CONSTI=UMAKG*1.E-20
    CONSB=1.438833
    hdc=cons.hbar_SI/(cons.c0_SI*100)/(cons.TWOPI)/2./100./CONSTI

    #---------------------------------#
    # Read coefficients of K operator #
    #---------------------------------#
    if LkinD:
       [ (d11_terms,d11_coefs) , (d12_terms,d12_coefs) , (d22_terms,d22_coefs) ] = dijtuple
       for dd in ["d11","d12","d22"]:
         if dd == "d11": dd_terms, dd_coefs = d11_terms, d11_coefs
         if dd == "d12": dd_terms, dd_coefs = d12_terms, d12_coefs
         if dd == "d22": dd_terms, dd_coefs = d22_terms, d22_coefs
         #--------------------------------------#
         # Format change: from D.F.C. to A.F.R. #
         #--------------------------------------#
         # Initialize arrays for Real coefs
         ak  = np.zeros(npar1 )            # 1D coef. associated with cos(x)
         bk  = np.zeros(npar2 )            # 1D coef. associated with cos(y)
         ck  = np.zeros( (ncoup ,ncoup ) ) # 2D coef. associated with cos(x)*cos(y)
         dk  = np.zeros( (ncoup ,ncoup ) ) # 2D coef. associated with sin(x)*sin(y)
         # Initialize arrays for Complex coefs
         akc = np.zeros(npar1c)            # 1D coef. associated with sin(x)
         bkc = np.zeros(npar2c)            # 1D coef. associated with sin(y)
         ckc = np.zeros( (ncoupc,ncoupc) ) # 2D coef. associated with cos(x)*sin(y)
         dkc = np.zeros( (ncoupc,ncoupc) ) # 2D coef. associated with sin(x)*cos(y)

         for dd_term, dd_coef in zip(dd_terms,dd_coefs):
             term_type,idx_phi1,idx_phi2 = dd_term

             if term_type == "const": ak0 = dd_coef

             if term_type == "cos" and idx_phi2 == "-": ak[idx_phi1-1] = dd_coef
             if term_type == "cos" and idx_phi1 == "-": bk[idx_phi2-1] = dd_coef
             if term_type == "coscos":                  ck[idx_phi1-1][idx_phi2-1] = dd_coef
             if term_type == "sinsin":                  dk[idx_phi1-1][idx_phi2-1] = dd_coef

             if term_type == "sin" and idx_phi2 == "-": akc[idx_phi1-1] = dd_coef
             if term_type == "sin" and idx_phi1 == "-": bkc[idx_phi2-1] = dd_coef
             if term_type == "sincos":                  ckc[idx_phi1-1][idx_phi2-1] = dd_coef
             if term_type == "cossin":                  dkc[idx_phi1-1][idx_phi2-1] = dd_coef
         if dd == "d11": ak011,ak11,bk11,ck11,dk11 = ak0, ak, bk, ck, dk
         if dd == "d12": ak012,ak12,bk12,ck12,dk12 = ak0, ak, bk, ck, dk
         if dd == "d22": ak022,ak22,bk22,ck22,dk22 = ak0, ak, bk, ck, dk
    #------------------------------------------------#
    # Reduced moments of inertia for most stable min #
    #------------------------------------------------#
    else:
       [rI1,rI2] = dijtuple
       # Apply constants
       rI1=CONSTI*rI1
       rI2=CONSTI*rI2
       hd1=cons.hbar_SI/(cons.c0_SI*100)/(cons.TWOPI)/rI1/2.
       hd2=cons.hbar_SI/(cons.c0_SI*100)/(cons.TWOPI)/rI2/2.


    nfil=nv1d**2
    nrow=nfil
    imax=(nv1d-1)/2
    istart=-imax
    tolx=1.e-20
    fil=[]
    col=[]
    ham2d=[]
    Lcompx = False
    if npar1c > 0 or npar2c > 0 or ncoupc > 0: Lcompx = True

    #----------------------------#
    # Start generation of matrix #
    #----------------------------#
    for n in range(nv1d):
      for k in range(nv1d):
         for m in range(nv1d):
            for j in range(nv1d):
               nh=istart+n
               kh=istart+k
               mh=istart+m
               jh=istart+j
               mj=m*nv1d+j
               nk=n*nv1d+k
               hkvd=0.
               hdn=hdc*nh**2
               hdkn=2.*hdc*nh*kh
               hdk=hdc*kh**2
               # ------------------
               # Diagonal terms (kinetic energy + independent term)
               # ------------------
               if mj==nk:
                 if LkinD==True:
                    hkin=hdn*ak011 + hdkn*ak012 + hdk*ak022
                    hkvd=(hkin+ap0)*hf.delta_ij(mh,nh)*hf.delta_ij(jh,kh)
                 else:
                    hkvd=(hd1*nh**2+hd2*kh**2+ap0)*hf.delta_ij(mh,nh)*hf.delta_ij(jh,kh)
                 fil.append(mj)
                 col.append(nk)
                 if Lcompx==True:
                   ham2d.append(complex(hkvd,0.))
                 else:
                   ham2d.append(hkvd)
                 #print mj+1,nk+1,jh,mh,kh,nh,hkvd
               else:
                 # ------------------
                 # Off diagonal terms
                 # ------------------
                 hap=0.
                 hbp=0.
                 hcp=0.
                 hapc=0.
                 hbpc=0.
                 hcpc=0.
                 sumh=0.
                 sumhc=0.
                 tmp1=0
                 tmp2=0
                 tmp1=abs(mh-nh)
                 iap=tmp1-1
                 tmp2=abs(jh-kh)
                 ibp=tmp2-1
                 ic1=tmp1-1
                 ic2=tmp2-1
                 # ////////////\\\\\\\\\\\\\\\\\
                 # --Real part of the potential--
                 # \\\\\\\\\\\\/////////////////
                 # Potential 1 + kinetic energy terms (if LkinD)
                 if tmp1 <= npar1:
                   if LkinD==True:
                     hakin=ak11[iap]*hdn+ak12[iap]*hdkn+ak22[iap]*hdk
                     hakd=hdc*tmp1*hf.sign(mh-nh)*(ak11[iap]*nh+ak12[iap]*kh)
                   else:
                     hakin=0.
                     hakd=0.
                   hap=(hakin+hakd+ap[iap])*hf.delta_ij(jh,kh)/2.
                 # Potencial 2 + kinetic energy terms (if LkinD)
                 if tmp2 <= npar2:
                   if LkinD==True:
                     hbkin=bk11[ibp]*hdn+bk12[ibp]*hdkn+bk22[ibp]*hdk
                     hbkd=hdc*tmp2*hf.sign(jh-kh)*(bk12[ibp]*nh+bk22[ibp]*kh)
                   else:
                     hbkin=0.
                     hbkd=0.
                   hbp=(hbkin+hbkd+bp[ibp])*hf.delta_ij(mh,nh)/2.
                 # Coupling + kinetic energy terms (if LkinD)
                 if tmp1<=ncoup and tmp1>0 and tmp2<=ncoup and tmp2>0:
                   if LkinD==True:
                      hckin=ck11[ic1][ic2]*hdn+ck12[ic1][ic2]*hdkn+ck22[ic1][ic2]*hdk
                      #
                      hdkin=-hf.sign(mh-nh)*hf.sign(jh-kh)*\
                      (dk11[ic1][ic2]*hdn+dk12[ic1][ic2]*hdkn+dk22[ic1][ic2]*hdk)
                      # 
                      hcdkd=hdc*nh*tmp1*\
                      (ck11[ic1][ic2]*hf.sign(mh-nh)-dk11[ic1][ic2]*hf.sign(jh-kh))\
                      +hdc*kh*tmp1*\
                      (ck12[ic1][ic2]*hf.sign(mh-nh)-dk12[ic1][ic2]*hf.sign(jh-kh))\
                      +hdc*nh*tmp2*\
                      (ck12[ic1][ic2]*hf.sign(jh-kh)-dk12[ic1][ic2]*hf.sign(mh-nh))\
                      +hdc*kh*tmp2*\
                      (ck22[ic1][ic2]*hf.sign(jh-kh)-dk22[ic1][ic2]*hf.sign(mh-nh))
                   else:
                      hckin=0.
                      hdkin=0.
                      hcdkd=0.
                   hcp1=cp[ic1][ic2]
                   hcp2=-hf.sign(mh-nh)*hf.sign(jh-kh)*dp[ic1][ic2]
                   hcp=(hckin+hdkin+hcdkd+hcp1+hcp2)/4.
                 sumh=hap+hbp+hcp
                 # /////////////\\\\\\\\\\\\\\\\\
                 # Imaginary part of the potential
                 # \\\\\\\\\\\\\/////////////////
                 # Potential 1
                 if tmp1 <= npar1c:
                   if npar1c > 0:
                     hapc=-hf.sign(mh-nh)*apc[iap]*hf.delta_ij(jh,kh)/2.
                 # Potencial 2
                 if tmp2 <= npar2c:
                   if npar2c > 0:
                     hbpc=-hf.sign(jh-kh)*bpc[ibp]*hf.delta_ij(mh,nh)/2.
                 # Coupling
                 if ncoupc > 0:
                   if tmp1<=ncoupc and tmp1>0 and tmp2<=ncoupc and tmp2>0:
                     hcp1c=-hf.sign(jh-kh)*cpc[ic1][ic2]
                     hcp2c=-hf.sign(mh-nh)*dpc[ic1][ic2]
                     hcpc=(hcp1c+hcp2c)/4.
                 sumhc=hapc+hbpc+hcpc
                 if abs(sumh) > tolx or abs(sumhc) > tolx:
                   fil.append(mj)
                   col.append(nk)
                   if Lcompx==True:
                     ham2d.append(complex(sumh,sumhc))
                   else:
                     ham2d.append(sumh)

    return fil,col,ham2d
#----------------------------------------------------------#
def tor2DNS_diagonalization(kmax,fil,col,ham2d,LkinD,file_evalues,ENERMAX=10000.,EIGMAX=250,MAXCYCLE=100):
   # -----------------------------
   # This function diagonalizes the hamiltonian of _hamvar2d_
   # Subroutine _eigsh_ can diagonalize real or complex matrices
   # The eigenvalues are sorted and stored in eigenvalues.dat 
   # Defaults:
   #   ENERMAX : Maximum energy at which eigenvalues are calculated
   #   EIGMAX  : Number of eigenvalues required at each cycle
   #   The program stops if _eigsh_ has been executed more than MAXCYCLE times 
   #             The counter is icount.
   # match: True if the last eigen. is different from the before last eigen.
   #        This is important to avoid cutting degenerate eigenvalues.
   # oldeig: It is the last eigenvalue of the previous run
   #         In the new run of _eigsh_ oldeig has to be larger than some new eigen.
   #         The search stops when that condicion is not met; the position
   #         is given by list_pos and the eigenvalues are stored in list_eig. 

   nv1d = 2*kmax + 1

   list_eig=[]
   nfil=nv1d**2
   nrow=nfil
   fila=np.array(fil)
   cola=np.array(col)

   exit_calc = False #### new
   negative_evals = False

   # Start calculation
   icount=0
   numsig=0.
   while numsig <= ENERMAX:
      if icount == MAXCYCLE:
         print 'Number of cycles exceeded'
         break
      end__=-1
      init_=-2
      match=False
      while match==False:
        if icount > 0:
           if list_eig[init_] != list_eig[end__]:
              oldeig=list_eig[end__]+1.e-3
              match=True
           else:
              end__ -= 1
              init_ -= 1
        else:
           oldeig=0.
           match=True
      # Diagonalizing Compressed Sparse Row matrix (eigsh may require CS Column matrix instead)
      warnings.filterwarnings('error')
      try:
          sparsematrix = csr_matrix((ham2d,(fila,cola)), shape=(nfil,nrow))
          res2 = eigsh(sparsematrix,k=EIGMAX,which='LM',sigma=numsig,return_eigenvectors=False)
      except:
          sparsematrix = csc_matrix((ham2d,(fila,cola)), shape=(nfil,nrow))
          res2 = eigsh(sparsematrix,k=EIGMAX,which='LM',sigma=numsig,return_eigenvectors=False)
      #res2 = eigsh(csr_matrix((ham2d,(fila,cola)), shape=(nfil,nrow)),k=EIGMAX,which='LM',\
      #sigma=numsig,return_eigenvectors=False)
      res2.sort()
      list_pos=0
      if oldeig != 0.0: print "              last eigenvalue is %8.2f cm^-1"%oldeig
      if icount > 0:
         for i in range(EIGMAX):
            neweig=res2[i]
            if oldeig>=neweig:
               list_pos += 1
               if i==EIGMAX-1:
                  print '              Unable to find new eigenvalues above',oldeig
                  exit_calc = True #### new
                  break            #### new
                 #sys.exit()           #### old
            else:
               break
      newlasteig=res2[EIGMAX-1]
      numsig=newlasteig
      for i in range(list_pos,EIGMAX):
         list_eig.append(res2[i])
      icount += 1
      if exit_calc: break #### new

   # Write data in file
   print "        * Storing eigenvalues in '%s'..."%file_evalues
   neigen=len(list_eig)
   ofile = open(file_evalues,'w')
   ofile.write("kmax   %i\n"%kmax)
   if LkinD: ofile.write("dijvar yes\n")
   else    : ofile.write("dijvar no \n")
   for i in range(neigen):
       ofile.write('%5i %10.5f\n'%(i,list_eig[i]))
   ofile.close()

   return list_eig
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 8) Helper functions related to: integration    #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def integration2D(integrand,dphi=2.0,args=[]):
    '''
    * dphi in degrees
    '''
    nsteps   = int(round(360.0/dphi))
    dphi_int = 360.0 / nsteps * cons.D2R
    dA       = dphi_int * dphi_int

    integral = 0.0
    for i in range(nsteps):
        for j in range(nsteps):
            phi1 = i * dphi_int
            phi2 = j * dphi_int
            integrand_input = tuple( [phi1,phi2]+args )
            integral += integrand( *integrand_input )
    return integral * dA
#----------------------------------------------------------#
def intgrnd_ClasTor_pf(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0):
    # beta
    beta = 1.0 / cons.kB / T
    # Energy in hartree
    V = Fourier_PES.value(phi1,phi2) * cons.c2h
    # |D| in a.u.
    sqrtD = rel_sqrtD(phi1,phi2)*sqrtD0
    # F
    F = sqrtD * np.exp(-beta*V)
    # Return integrand
    return F
#----------------------------------------------------------#
def intgrnd_ClasTor_fd(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0):
    # beta
    beta = 1.0 / cons.kB / T
    # Energy in hartree
    V = Fourier_PES.value(phi1,phi2) * cons.c2h
    # |D| in a.u.
    sqrtD = rel_sqrtD(phi1,phi2)*sqrtD0
    # F
    F = sqrtD * np.exp(-beta*V)
    # integrand
    integrand = F * (V+1.0/beta)
    # Return integrand
    return integrand
#----------------------------------------------------------#
def intgrnd_ClasTor_sd(phi1,phi2,T,Fourier_PES,rel_sqrtD,sqrtD0):
    # beta
    beta = 1.0 / cons.kB / T
    # Energy in hartree
    V = Fourier_PES.value(phi1,phi2) * cons.c2h
    # |D| in a.u.
    sqrtD = rel_sqrtD(phi1,phi2)*sqrtD0
    # F
    F = sqrtD * np.exp(-beta*V)
    # integrand
    integrand = F * (1.0/beta/beta + (V+1.0/beta)**2)
    # Return integrand
    return integrand
#----------------------------------------------------------#
def intgrnd_EHR_pf(phi1,phi2,T,FPES,tupleS,tupleNT,ref_Vehr):
    '''
    * phi1, phi2: dihedral angles for torsions (in radians)
    * T: temperature (in Kelvin)
    * FPES: Fourier2D instance for the PES
    * tupleS = (sqrtS0,rel_sqrtS)
      sqrtS0   : reference value for sqrt(|D|)
      rel_sqrtS: 2D-interpolation (Fourier2D or spline) for the relative value
    * tupleNT = (delaunay,Vmidpoints,spdata)
      delaunay: a Delaunay instance
      Vmidpoints: a dictionary containing the mid point (in PES) for each triangle edge in the tesselation
      spdata: a list with the data to be interpolated using the tesselation
               spdata = [(sigmarot,ntzpe,ntpf,ntfdln,ntsdln)], length equal to the number of stationary points
               ntpf: non-torsional partition function
               ntfdln: non-torsional 1st derivative of ln(pf)
               ntsdln: non-torsional 2nd derivative of ln(pf)
    * ref_ehf = the PES part of the reference energy for EHR
    '''
    beta = 1.0/cons.kB/T
    p0   = np.array( (phi1,phi2) )

    # (a) Interpolation of potential
    U0 = FPES(phi1,phi2) * cons.c2h - ref_Vehr

    # (b) Interpolation of sqrt(|S|)
    sqrtS  = tupleS[0] * tupleS[1](phi1,phi2)

    # (c) Interpolations for non-torsional pfs (and sigma_rot)
    delaunay, Vmidpoints, spdata = tupleNT
    args = [FPES,Vmidpoints]
    interpolated_data, idx = tess.q2dtor_interpolation(p0,delaunay,spdata[1:3],method=TESSMODE,args=args)
    # (c.1) sigma_rot
    sigmarot = spdata[0][idx]
    # (c.2) diff(zpe)
    ntzpe   = interpolated_data[0]
    # (c.3) tilde partition function
    ntpf    = interpolated_data[1]

    # (d) Get integrand
    integrand = sqrtS*np.exp(-U0*beta)*np.exp(-ntzpe*beta)*ntpf / sigmarot
    # v1.1: ignore sigmarot (it'll be considered later)
    if BOOL_v1_1: integrand *= sigmarot

    return integrand

#----------------------------------------------------------#
def intgrnd_EHR_fd(phi1,phi2,T,FPES,tupleS,tupleNT,ref_Vehr):
    beta = 1.0/cons.kB/T
    p0   = np.array( (phi1,phi2) )

    # (a) Interpolation of potential
    U0 = FPES(phi1,phi2) * cons.c2h - ref_Vehr

    # (b) Interpolation of sqrt(|S|)
    sqrtS  = tupleS[0] * tupleS[1](phi1,phi2)

    # (c) Interpolations for non-torsional pfs (and sigma_rot)
    delaunay, Vmidpoints, spdata = tupleNT
    args = [FPES,Vmidpoints]
    interpolated_data, idx = tess.q2dtor_interpolation(p0,delaunay,spdata[1:4],method=TESSMODE,args=args)
    # (c.1) sigma_rot
    sigmarot = spdata[0][idx]
    # (c.2) diff(zpe)
    ntzpe   = interpolated_data[0]
    # (c.3) tilde partition function
    ntpf    = interpolated_data[1]
    # (c.4) 1st derivative of ln(tilde partition function)
    ntfdln  = interpolated_data[2]

    # (d) Get integrand
    integrand  = sqrtS*np.exp(-U0*beta)*np.exp(-ntzpe*beta)*ntpf/sigmarot
    integrand *= (ntfdln-2.5/beta-(U0+ntzpe))
    # v1.1: ignore sigmarot (it'll be considered later)
    if BOOL_v1_1: integrand *= sigmarot

    return integrand
#----------------------------------------------------------#
def intgrnd_EHR_sd(phi1,phi2,T,FPES,tupleS,tupleNT,ref_Vehr):
    beta = 1.0/cons.kB/T
    p0   = np.array( (phi1,phi2) )

    # (a) Interpolation of potential
    U0 = FPES(phi1,phi2) * cons.c2h - ref_Vehr

    # (b) Interpolation of sqrt(|S|)
    sqrtS  = tupleS[0] * tupleS[1](phi1,phi2)

    # (c) Interpolations for non-torsional pfs (and sigma_rot)
    delaunay, Vmidpoints, spdata = tupleNT
    args = [FPES,Vmidpoints]
    interpolated_data, idx = tess.q2dtor_interpolation(p0,delaunay,spdata[1:5],method=TESSMODE,args=args)
    # (c.1) sigma_rot
    sigmarot = spdata[0][idx]
    # (c.2) diff(zpe)
    ntzpe   = interpolated_data[0]
    # (c.3) tilde partition function
    ntpf    = interpolated_data[1]
    # (c.4) 1st derivative of ln(tilde partition function)
    ntfdln  = interpolated_data[2]
    # (c.5) 2nd derivative of ln(tilde partition function)
    ntsdln  = interpolated_data[3]


    # (d) Get integrand
    integrand  = sqrtS*np.exp(-U0*beta)*np.exp(-ntzpe*beta)*ntpf/sigmarot
    integrand *= (ntsdln+2.5/beta/beta+(2.5/beta+U0+ntzpe-ntfdln)**2)
    # v1.1: ignore sigmarot (it'll be considered later)
    if BOOL_v1_1: integrand *= sigmarot

    return integrand
#----------------------------------------------------------#
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION ( 9) Helper functions related to: strings        #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def string_checkname(thename,exit=True):
    large_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    small_letters = "abcdefghijklmnopqrstuvwxyz"
    numbers       = "0123456789"
    symbols       = ".-"
    allowed       = large_letters+small_letters+numbers+symbols
    check_string = ""
    for character in thename:
        if character not in allowed:
           check_string = check_string + "'%s',"%(character)
    if check_string != "":
       message =           "Name for input file : '%s'\n"%thename
       message = message + "Error in name due to: %s\n"%check_string[:-1]
       message = message + "Allowed characters  : '%s'"%allowed
       if exit: sys.exit(message)
       else   :   return message
    else: return None
#----------------------------------------------------------#
def string_getpname(phi1,phi2,molname="",units="rad"):
    '''
    This function generates a name for each point in the 2D-torsional PES
    * units: indicates the units of phi1 and phi2 ('rad' or 'deg')
    PS: points that are less than a degree apart cannot be differentiated by the generated name
    '''
    if units == "rad": int_phi1, int_phi2 = int(round(phi1*cons.R2D)), int(round(phi2*cons.R2D))
    if units == "deg": int_phi1, int_phi2 = int(round(phi1    )), int(round(phi2    ))
    return molname+"_%03i_%03i"%(int_phi1,int_phi2)
#----------------------------------------------------------#
def string_getinfo(pname):
    name, phi1, ph2 = pname.split("_")
    phi1 = float(phi1) * cons.D2R
    phi2 = float(phi2) * cons.D2R
    return phi1, phi2
#----------------------------------------------------------#
def string_pfntable(Tlist,pf1WHO,pfMSHO,pfE2DT,pf2DNS,pfTorClas,pfEHR,ratio2DNS,Eref,qtra,qele):

    tables = ""
    tables = tables + " #---------------------#\n"
    tables = tables + " # Partition Functions #\n"
    tables = tables + " #---------------------#\n"
    tables = tables + "\n"
    # Rovibrational partition functions (tilde)
    table1 = "     (a) Rovibrational partition functions using as zero of energy\n"       +\
             "         the lowest zero point level of the torsional PES         \n" +\
             "                                                              \n" +\
             "            T (K)  |  rv(1WHO)  |  rv(MSHO)  |  rv(E2DT)  | E2DT/MSHO \n" +\
             "          ------------------------------------------------------------\n"

    table2 = "         Components of E2DT\n"       +\
             "                                                       \n" +\
             "            T (K)  |    2DNS    |  TorsClas  |     EHR    \n" +\
             "          ------------------------------------------------\n"


    table3 = "     (b) Rovibrational partition functions obtained using as zero of energy\n"       +\
             "         the bottom of the torsional PES                   \n" +\
             "                                                       \n" +\
             "            T (K)  |  rv(1WHO)  |  rv(MSHO)  |  rv(E2DT)  | E2DT/MSHO \n" +\
             "          ------------------------------------------------------------\n"

    table4 = "         Components of E2DT\n"       +\
             "                                                        \n" +\
             "            T (K)  |    2DNS    |  TorsClas  |     EHR    | Fq(2DNS) \n" +\
             "          -----------------------------------------------------------\n"

    table5 = "     Translational and electronic partition functions \n"       +\
             "                                                        \n" +\
             "            T (K)  |   Qtrans   |   Qelec    \n" +\
             "          -----------------------------------\n"

    table6 = "     Total partition functions, from rovibrational partition functions in (b)\n"       +\
             "                                                       \n" +\
             "            T (K)  |    1WHO    |    MSHO    |    E2DT    \n" +\
             "          ------------------------------------------------\n"


    for idx in range(len(Tlist)):
        T = Tlist[idx]
        beta = 1.0 / cons.kB / T
        # Tilde version
        pfn_1WHO = pf1WHO[idx]
        pfn_MSHO = pfMSHO[idx]
        pfn_E2DT = pfE2DT[idx]
        ratio    = pfn_E2DT/pfn_MSHO
        line     = (T,pfn_1WHO,pfn_MSHO,pfn_E2DT,ratio)
        table1   = table1 + "           %7.2f | %10.3E | %10.3E | %10.3E |  %7.5f  \n"%line
        # Components E2DT
        pfn_2DNS = pf2DNS[idx]
        torclas = pfTorClas[idx]
        pfn_EHR = pfEHR[idx]
        line    = (T,pfn_2DNS,torclas,pfn_EHR)
        table2 = table2 + "           %7.2f | %10.3E | %10.3E | %10.3E \n"%line
        # Non-tilde
        pfn_1WHO = pf1WHO[idx] * np.exp(-beta*Eref[0])
        pfn_MSHO = pfMSHO[idx] * np.exp(-beta*Eref[1])
        pfn_E2DT = pfE2DT[idx] * np.exp(-beta*Eref[4])
        ratio    = pfn_E2DT/pfn_MSHO
        line     = (T,pfn_1WHO,pfn_MSHO,pfn_E2DT,ratio)
        table3 = table3 + "           %7.2f | %10.3E | %10.3E | %10.3E |  %7.5f  \n"%line
        # Components E2DT
        Fq2DNS   = ratio2DNS[idx]
        if Fq2DNS == "unit": Fq2DNS = 1.0
        torclas  = pfTorClas[idx]
        pfn_2DNS = Fq2DNS*torclas
        pfn_EHR  = pfEHR[idx] * np.exp(-beta*Eref[3])
        line     = (T,pfn_2DNS,torclas,pfn_EHR,Fq2DNS)
        table4 = table4 + "           %7.2f | %10.3E | %10.3E | %10.3E | %7.5f \n"%line
        # translation and electronic
        pfn_trans = qtra[idx]
        pfn_elect = qele[idx]
        line     = (T,pfn_trans,pfn_elect)
        table5 = table5 + "           %7.2f | %10.3E | %10.3E \n"%line
        # Total partition functions
        pfn_1WHO = pf1WHO[idx] * np.exp(-beta*Eref[0]) * pfn_trans * pfn_elect
        pfn_MSHO = pfMSHO[idx] * np.exp(-beta*Eref[1]) * pfn_trans * pfn_elect
        pfn_E2DT = pfE2DT[idx] * np.exp(-beta*Eref[4]) * pfn_trans * pfn_elect
        ratio    = pfn_E2DT/pfn_MSHO
        line     = (T,pfn_1WHO,pfn_MSHO,pfn_E2DT)
        table6 = table6 + "           %7.2f | %10.3E | %10.3E | %10.3E \n"%line

    tables = tables + "     Energy of the lowest zero point level of the torsional PES:\n"
    tables = tables + "       * 1WHO => %7.3f %s\n"%(Eref[0]*HU1,SU1)
    tables = tables + "       * MSHO => %7.3f %s\n"%(Eref[1]*HU1,SU1)
    tables = tables + "       * 2DNS => %7.3f %s\n"%(Eref[2]*HU1,SU1)
    tables = tables + "       * EHR  => %7.3f %s\n"%(Eref[3]*HU1,SU1)
    tables = tables + "       * E2DT => %7.3f %s\n"%(Eref[4]*HU1,SU1)
    tables = tables + "\n"
    tables = tables + table1 + "\n"
    tables = tables + table2 + "\n"
    tables = tables + table3 + "\n"
    tables = tables + table4 + "\n"
    tables = tables + table5 + "\n"
    tables = tables + table6 + "\n"
    tables = tables + "\n"

    # Add translational and electronic


    return tables
#----------------------------------------------------------#
def string_thermotable(Tlist,ther1WHO,therMSHO, therE2DT):
    tables  = ""
    tables += " #-------------------------#\n"
    tables += " # Thermodynamic Functions #\n"
    tables += " #-------------------------#\n"
    tables += "\n"
    tables += "     Note: All thermodynamic functions are calculated using as\n"
    tables += "     zero of energy the lowest zero point level of the torsional PES\n"
    tables += "\n"
    tables += "     Units:\n"
    tables += "        * U, H, G => %s\n"%SU1
    tables += "        * S, Cp   =>  %s\n"%SU2
    tables += "\n"

    # Table for each method
    for pfn,thermo in [ ("1WHO",ther1WHO) , ("MSHO",therMSHO), ("E2DT",therE2DT)]:
        table = "     Thermodynamic functions for %s partition function\n"%pfn +\
                "                                                                              \n" +\
                "           T (K)  |    U^o    |    H^o    |    S^o    |    G^o    |    Cp     \n" +\
                "         ---------------------------------------------------------------------\n"
        for idx in range(len(Tlist)):
            T = Tlist[idx]
            U, H, Cp, S, G = thermo[idx]
            idata  = (T, U*HU1, H*HU1, S*HU2, G*HU1, Cp*HU2)
            table += "          %7.2f | %9.3f | %9.3f | %9.3f | %9.3f | %9.3f \n"%idata
        tables += table + "\n"
    tables += "\n"

    # Table for each thermodynamic function
    tableU   = "     Internal Energy (U^o)                               \n"
    tableU  += "                                                         \n"
    tableU  += "           T (K)  |    1WHO    |    MSHO    |    E2DT    \n"
    tableU  += "         ------------------------------------------------\n"

    tableH   = "     Enthalpy (H^o)                                      \n"
    tableH  += "                                                         \n"
    tableH  += "           T (K)  |    1WHO    |    MSHO    |    E2DT    \n"
    tableH  += "         ------------------------------------------------\n"

    tableS   = "     Entropy (S^o)                                       \n"
    tableS  += "                                                         \n"
    tableS  += "           T (K)  |    1WHO    |    MSHO    |    E2DT    \n"
    tableS  += "         ------------------------------------------------\n"

    tableG   = "     Gibbs Free Energy (G^o)                             \n"
    tableG  += "                                                         \n"
    tableG  += "           T (K)  |    1WHO    |    MSHO    |    E2DT    \n"
    tableG  += "         ------------------------------------------------\n"

    tableCp  = "     Heat Capacity at constant pressure (Cp)             \n"
    tableCp += "                                                         \n"
    tableCp += "           T (K)  |    1WHO    |    MSHO    |    E2DT    \n"
    tableCp += "         ------------------------------------------------\n"

    for idx,T in enumerate(Tlist):
        list_U  = (T,ther1WHO[idx][0]*HU1,therMSHO[idx][0]*HU1,therE2DT[idx][0]*HU1)
        list_H  = (T,ther1WHO[idx][1]*HU1,therMSHO[idx][1]*HU1,therE2DT[idx][1]*HU1)
        list_Cp = (T,ther1WHO[idx][2]*HU2,therMSHO[idx][2]*HU2,therE2DT[idx][2]*HU2)
        list_S  = (T,ther1WHO[idx][3]*HU2,therMSHO[idx][3]*HU2,therE2DT[idx][3]*HU2)
        list_G  = (T,ther1WHO[idx][4]*HU1,therMSHO[idx][4]*HU1,therE2DT[idx][4]*HU1)
        tableU  += "          %7.2f | %10.3f | %10.3f | %10.3f \n"%list_U
        tableH  += "          %7.2f | %10.3f | %10.3f | %10.3f \n"%list_H
        tableS  += "          %7.2f | %10.3f | %10.3f | %10.3f \n"%list_S
        tableG  += "          %7.2f | %10.3f | %10.3f | %10.3f \n"%list_G
        tableCp += "          %7.2f | %10.3f | %10.3f | %10.3f \n"%list_Cp
    tables += tableU +"\n"
    tables += tableH +"\n"
    tables += tableS +"\n"
    tables += tableG +"\n"
    tables += tableCp+"\n"
    tables += "\n"
    return tables
#----------------------------------------------------------#
#----------------------------------------------------------#
#----------------------------------------------------------#




# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION (10) Helper functions related to: diverse things #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def ic_checktorsions(icoords,torsion1,torsion2):
    ics = []
    t1 = list(torsion1)
    t2 = list(torsion2)
    t1in = False
    t2in = False
    if t1[0] > t1[3]: t1 = t1[::-1]
    if t2[0] > t2[3]: t2 = t2[::-1]
    for ic_type, ic in icoords:
        if ic_type == "4":
           t0 = list(ic)
           if t0[0] > t0[3]: t0 = t0[::-1]
           if t0 == t1: t1in = True; continue
           if t0 == t2: t2in = True; continue
        ics.append( (ic_type,ic) )
    ics.append( ("4",torsion1) )
    ics.append( ("4",torsion2) )
    if not t1in:
       print "     ERROR: internal coordinates do not contain torsion 1"
       print
       sys.exit(EXITMESS)
    if not t2in:
       print "     ERROR: internal coordinates do not contain torsion 2"
       print
       sys.exit(EXITMESS)
    return ics
#----------------------------------------------------------#
def clean_close(dict_cps,basic=[],diff_angle=4.0,diff_Ecm=10.0):
    '''
    '''
    dict_cps2 = {}

    # The basic ones will remain
    for ppname in basic:
        dict_cps2[ppname] = dict_cps.pop(ppname)

    # Now compare the rest
    for ppname in dict_cps.keys():
        phi1, phi2, cptype, E = dict_cps[ppname][0:4]
        keep = True
        for rname in dict_cps2.keys():
            rp1, rp2, rt, rE = dict_cps2[rname][0:4]
            # Same point?
            diff1 = hf.angle_diff(phi1,rp1)
            diff2 = hf.angle_diff(phi2,rp2)
            dist = np.sqrt(diff1**2 + diff2**2)
            cond1 = (dist*cons.R2D < diff_angle)
            cond2 = abs(rE-E) < diff_Ecm
            cond3 = (rt == cptype)
            if cond1 and cond2 and cond3:
               keep = False
              #if cptype == 0 and E < rE:
              #   keep = True
              #   dict_cps2.pop(rname)
              #   if rname in basic:
              #      basic.remove(rname)
              #      basic.append(ppname)
              #if cptype == 2 and E > rE:
              #   keep = True
              #   dict_cps2.pop(rname)
              #   if rname in basic:
              #      basic.remove(rname)
              #      basic.append(ppname)
              #     #print rname, ppname
               break
        if keep: dict_cps2[ppname] = dict_cps[ppname]

    return dict_cps2
#----------------------------------------------------------#
def change_Eref(old_Eref,new_Eref,f_pes,f_vfit):

    #-------------------------#
    # Change reference Energy #
    #-------------------------#
    if old_Eref > new_Eref:
       print "     Reference energy has to be updated:"
       print "        - old reference: %+14.7f hartree"%old_Eref
       print "        - new reference: %+14.7f hartree"%new_Eref
       print "        - Modifying file: %s"%f_pes
       lines = []
       # Read file
       ff = open(f_pes,'r')
       for line in ff: lines.append(line)
       ff.close()
       # Modify file
       ff = open(f_pes,'w')
       for line in lines:
           if not line.startswith("MINENERGY"): ff.write(line)
       ff.write("MINENERGY         %+14.8f\n"%new_Eref)
       ff.close()

       print "        - Modifying file: %s"%f_vfit
       lines = []
       # Read file
       ff = open(f_vfit,'r')
       for line in ff:
           lines.append(line)
           if line.startswith("const"): A0 = float(line.split()[3])
       ff.close()
       # Modify file
       ff = open(f_vfit,'w')
       dE = (new_Eref - old_Eref) / cons.h / cons.c0 / cons.cm
       ff.write("const      %02s    %02s   %+11.4f\n"%("-","-",A0-dE))
       for line in lines:
           if not line.startswith("const"): ff.write(line)
       ff.close()
       print
#----------------------------------------------------------#
def pfn_derivatives(structure,qele,T,k="cc"):

    # Number of vibrational freqs to consider
    nimag = 0
    if structure._type >= 0: nimag = structure._type

    # Get beta and partition functions
    beta = 1.0 / cons.kB / T

    #---------------#
    # Translational #
    #---------------#
    # at constant V
    fdln_tra_V = -1.5 / beta
    sdln_tra_V = +1.5 / beta / beta
    # at constant p
    fdln_tra_P = -2.5 / beta
    sdln_tra_P = +2.5 / beta / beta

    #---------------#
    # Rotational    #
    #---------------#
    fdln_rot = -1.5 / beta
    sdln_rot = +1.5 / beta / beta

    #---------------#
    # Vibrational   #
    #---------------#
    if k=="cc":
       fdln_vib = sum([ff.get_fdln(T) for ff in structure._ccfreqs[nimag:]])
       sdln_vib = sum([ff.get_sdln(T) for ff in structure._ccfreqs[nimag:]])
    if k=="ic":
       fdln_vib = sum([ff.get_fdln(T) for ff in structure._icfreqs[nimag:]])
       sdln_vib = sum([ff.get_sdln(T) for ff in structure._icfreqs[nimag:]])

    #---------------#
    # Electronic    #
    #---------------#
    fd_ele = - sum( [(relE   )*mtp*np.exp(-beta*relE) for relE,mtp in structure._elstates]  )
    sd_ele = + sum( [(relE**2)*mtp*np.exp(-beta*relE) for relE,mtp in structure._elstates]  )
    fdln_ele = fd_ele / qele
    sdln_ele = sd_ele / qele - (fd_ele/qele)**2

    tuple_tra1 = (fdln_tra_V,sdln_tra_V) # at constant V
    tuple_tra2 = (fdln_tra_P,sdln_tra_P) # at constant p
    tuple_rot  = (fdln_rot,sdln_rot)
    tuple_vib  = (fdln_vib,sdln_vib)
    tuple_ele  = (fdln_ele,sdln_ele)

    return tuple_tra1, tuple_tra2, tuple_rot, tuple_vib, tuple_ele
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION (11) Functions for reading files                 #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def get_defaults():
    '''
    Default values as would be read from input file
    This means that angles are in degrees
    '''
    #---------------#
    torsion1        = None
    torsion2        = None
    tsigma1         = None
    tsigma2         = None
    #---------------#
    level           = "hf sto-3g"
    charge          = "0"
    multiplicity    = "1"
    ttype           = "min"
    #---------------#
    t1step          = "10.0"
    t2step          = "10.0"
    symmetry        = "none"
    #---------------#
    weight          = "0.9"
    ignore          = "0.0"
    cos1            = "1-6"
    cos2            = "1-6"
    sin1            = "none"
    sin2            = "none"
    cos1cos2        = "1-6 , 1-6"
    sin1sin2        = "1-6 , 1-6"
    cos1sin2        = "none , none"
    sin1cos2        = "none , none"
    #---------------#
    tolerance       = "2.0"
    freqscal        = "1.0"
    #---------------#
    dijvar          = "yes"
    kmax            = "100"
    maxeigen        = "1e4"
    #---------------#
    interpolation   = "fourier"
    integrationstep = "1.0"
    #---------------#

    return torsion1     , torsion2       , tsigma1     , tsigma2,  \
           level        , charge         , multiplicity, ttype  ,  \
           t1step       , t2step         , symmetry    ,           \
           weight       , ignore         ,                         \
           cos1         , cos2           , cos1cos2    , sin1sin2, \
           sin1         , sin2           , cos1sin2    , sin1cos2, \
           tolerance    , freqscal       ,                         \
           dijvar       , kmax           , maxeigen    ,           \
           interpolation, integrationstep
#----------------------------------------------------------#
def deal1Dfourierstring(string1D):
    if "none" in string1D: return [], 0
    indices = []
    if "-" in string1D:
       start, end = string1D.split("-")
       indices = range(int(start),int(end)+1,1)
    else:
       indices = list(set([int(ii) for ii in string1D.split()]))
    indices.sort()
    if len(indices) == 0: last_idx = 0
    else                : last_idx = indices[-1]
    return indices, last_idx
#----------------------------------------------------------#
def deal2Dfourierstring(string2D):
    if "none" in string2D: return [], [], 0

    if   "," in string2D: string1, string2 = string2D.split(",")
    elif ";" in string2D: string1, string2 = string2D.split(";")
    elif ":" in string2D: string1, string2 = string2D.split(":")
    else: print " ERROR: Problems reading 2D Fourier terms...\n"; sys.exit(EXITMESS)

    indices1, n1 = deal1Dfourierstring(string1)
    indices2, n2 = deal1Dfourierstring(string2)

    return indices1,indices2,max(n1,n2)
#----------------------------------------------------------#
def readfile_inp(inputname):

    mainname  = inputname.split(".inp")[0]
    f_xyz     = mainname+".xyz"
    f_calcs   = mainname+".calcs"
    f_ics     = mainname+".ics"
    f_pes     = mainname+".pes"
    f_vfit    = mainname+".vfit"
    f_dfit    = mainname+".dfit"
    f_splist  = mainname+".splist"
    f_spinfo  = mainname+".spinfo"
    f_evals   = mainname+".evals"
    f_sqrt    = mainname+".sqrt"
    f_tables  = mainname+".tables"
    f_pdf     = mainname+".pdf"
    f_out     = mainname+".out"

    #--------------------#
    # Get default values #
    #--------------------#
    torsion1     , torsion2       , tsigma1     , tsigma2,  \
    level        , charge         , multiplicity, ttype  ,  \
    t1step       , t2step         , symmetry    ,           \
    weight       , ignore         ,                         \
    cos1         , cos2           , cos1cos2    , sin1sin2, \
    sin1         , sin2           , cos1sin2    , sin1cos2, \
    tolerance    , freqscal       ,                         \
    dijvar       , kmax           , maxeigen    ,           \
    interpolation, integrationstep                            = get_defaults()

    #----------------------#
    # Read from input file #
    #----------------------#
    lines = hf.readfile(inputname,hashtag=True,strip=True,skipblank=True)
    for line in lines:
        if "{" in line and "}" in line:
          print "The input file has to be modified. Error reading line:\n%s"%line
          sys.exit(EXITMESS)

    for line in hf.select_lines(lines,"start_torsions" ,"end_torsions" ,True):
        if line.startswith("torsion1 "       ): torsion1        = " ".join(line.split()[1:])
        if line.startswith("torsion2 "       ): torsion2        = " ".join(line.split()[1:])
        if line.startswith("tsigma1 "        ): tsigma1         = " ".join(line.split()[1:])
        if line.startswith("tsigma2 "        ): tsigma2         = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_calcs"    ,"end_calcs"    ,True):
        if line.startswith("level "          ): level           = " ".join(line.split()[1:])
        if line.startswith("charge "         ): charge          = " ".join(line.split()[1:])
        if line.startswith("multiplicity "   ): multiplicity    = " ".join(line.split()[1:])
        if line.startswith("ttype "          ): ttype           = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_pes"      ,"end_pes"      ,True):
        if line.startswith("t1step "         ): t1step          = " ".join(line.split()[1:])
        if line.startswith("t2step "         ): t2step          = " ".join(line.split()[1:])
        if line.startswith("symmetry "       ): symmetry        = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_fourier"  ,"end_fourier"  ,True):
        if line.startswith("weight "         ): weight          = " ".join(line.split()[1:])
        if line.startswith("ignore "         ): ignore          = " ".join(line.split()[1:])
        if line.startswith("cos1 "           ): cos1            = " ".join(line.split()[1:])
        if line.startswith("cos2 "           ): cos2            = " ".join(line.split()[1:])
        if line.startswith("sin1 "           ): sin1            = " ".join(line.split()[1:])
        if line.startswith("sin2 "           ): sin2            = " ".join(line.split()[1:])
        if line.startswith("cos1cos2 "       ): cos1cos2        = " ".join(line.split()[1:])
        if line.startswith("sin1sin2 "       ): sin1sin2        = " ".join(line.split()[1:])
        if line.startswith("cos1sin2 "       ): cos1sin2        = " ".join(line.split()[1:])
        if line.startswith("sin1cos2 "       ): sin1cos2        = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_statpoint","end_statpoint",True):
        if line.startswith("tolerance "      ): tolerance       = " ".join(line.split()[1:])
        if line.startswith("freqscal "       ): freqscal        = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_tor2dns"  ,"end_tor2dns"  ,True):
        if line.startswith("dijvar "         ): dijvar          = " ".join(line.split()[1:])
        if line.startswith("kmax "           ): kmax            = " ".join(line.split()[1:])
        if line.startswith("maxeigen "       ): maxeigen        = " ".join(line.split()[1:])
    for line in hf.select_lines(lines,"start_rovibpf"  ,"end_rovibpf"  ,True):
        if line.startswith("interpolation "  ): interpolation   = " ".join(line.split()[1:])
        if line.startswith("integrationstep "): integrationstep = " ".join(line.split()[1:])

    #--------------#
    # Temperatures #
    #--------------#
    temperatures = []
    for line in hf.select_lines(lines,"start_temperatures","end_temperatures",True):
        if "range" in line:
           Ti, Tf, dT = [float(value) for value in line.split()[1:]]
           nsteps = int(round((Tf-Ti)/dT))+1
           temperatures += [Ti + j*dT for j in range(nsteps)]
        else:
           temperatures += [float(T) for T in line.split()]
    if temperatures == []:
       temperatures = [273.15, 298.15]
    temperatures.sort()

    #---------------------------------#
    # Check that basic data was given #
    #---------------------------------#
    if torsion1 is None: print "keyword 'torsion1' has no value!\n"; sys.exit(EXITMESS)
    if torsion2 is None: print "keyword 'torsion2' has no value!\n"; sys.exit(EXITMESS)
    if tsigma1  is None: print "keyword 'tsigma1'  has no value!\n"; sys.exit(EXITMESS)
    if tsigma2  is None: print "keyword 'tsigma2'  has no value!\n"; sys.exit(EXITMESS)

    #--------------#
    # Convert data #
    #--------------#
    torsion1 = tuple([int(ii)-1 for ii in torsion1.split("-")])
    torsion2 = tuple([int(ii)-1 for ii in torsion2.split("-")])
    if torsion1[0] > torsion1[3]: torsion1 = torsion1[::-1]
    if torsion2[0] > torsion2[3]: torsion2 = torsion2[::-1]
    tsigma1  = int(tsigma1)
    tsigma2  = int(tsigma2)
    #--------------#
    charge   = int(charge)
    multiplicity = int(multiplicity)
    #--------------#
    t1step   = float(t1step)*cons.D2R
    t2step   = float(t2step)*cons.D2R
    sym      = ""
    if "a" in symmetry: sym += "a"
    if "b" in symmetry: sym += "b"
    if "c" in symmetry: sym += "c"
    if sym == "": symmetry = "none"
    else        : symmetry = sym
    #--------------#
    tolerance = max(float(tolerance),1.0)*cons.D2R
    freqscal  = float(freqscal)
    #--------------#
    if dijvar.lower() in ["y","yes","true"]: dijvar = True
    else                                   : dijvar = False
    kmax = int(kmax)
    maxeigen = float(maxeigen)
    #--------------#
    try   : interpolation = int(interpolation)
    except: interpolation = "fourier"
    integrationstep = float(integrationstep)*cons.D2R
    #--------------#
    weight = float(weight)
    ignore = float(ignore)
    #--------------#
    idx_c1, num_cos1 = deal1Dfourierstring(cos1)
    idx_c2, num_cos2 = deal1Dfourierstring(cos2)
    idx_s1, num_sin1 = deal1Dfourierstring(sin1)
    idx_s2, num_sin2 = deal1Dfourierstring(sin2)
    #--------------#
    idx_cc1, idx_cc2, num_coscos = deal2Dfourierstring(cos1cos2)
    idx_ss1, idx_ss2, num_sinsin = deal2Dfourierstring(sin1sin2)
    idx_cs1, idx_cs2, num_cossin = deal2Dfourierstring(cos1sin2)
    idx_sc1, idx_sc2, num_sincos = deal2Dfourierstring(sin1cos2)
    #--------------#
    nums = (num_cos1, num_cos2, max(num_coscos,num_sinsin), num_sin1, num_sin2, max(num_sincos,num_cossin))
    #--------------#


    #--------------------------------#
    # v1.1 - symmetry considerations #
    #--------------------------------#
    if BOOL_v1_1:
       # add "c" symmetry if totaltsigma != 1
       if (tsigma1*tsigma2 != 1) and ("c" not in symmetry): symmetry += "c"
       # remove Fourier terms considering tsigma
       if tsigma1 != 1:
           idx_c1   = [idx for idx in idx_c1  if idx%tsigma1 == 0]
           idx_s1   = [idx for idx in idx_s1  if idx%tsigma1 == 0]
           idx_cc1  = [idx for idx in idx_cc1 if idx%tsigma1 == 0]
           idx_ss1  = [idx for idx in idx_ss1 if idx%tsigma1 == 0]
           idx_sc1  = [idx for idx in idx_sc1 if idx%tsigma1 == 0]
           idx_cs1  = [idx for idx in idx_cs1 if idx%tsigma1 == 0]
       if tsigma2 != 1:
           idx_c2   = [idx for idx in idx_c2  if idx%tsigma2 == 0]
           idx_s2   = [idx for idx in idx_s2  if idx%tsigma2 == 0]
           idx_cc2  = [idx for idx in idx_cc2 if idx%tsigma2 == 0]
           idx_ss2  = [idx for idx in idx_ss2 if idx%tsigma2 == 0]
           idx_sc2  = [idx for idx in idx_sc2 if idx%tsigma2 == 0]
           idx_cs2  = [idx for idx in idx_cs2 if idx%tsigma2 == 0]
       # remove odd Fourier terms if 'b' in symmetry
       if "b" in symmetry:
           idx_s1  = []
           idx_s2  = []
           idx_sc1 = []
           idx_cs1 = []
           idx_sc2 = []
           idx_cs2 = []

       # updata num variables
       if len(idx_c1) != 0: num_cos1 = max(idx_c1)
       else               : num_cos1 = 0
       if len(idx_s1) != 0: num_sin1 = max(idx_s1)
       else               : num_sin1 = 0
       if len(idx_c2) != 0: num_cos2 = max(idx_c2)
       else               : num_cos2 = 0
       if len(idx_s2) != 0: num_sin2 = max(idx_s2)
       else               : num_sin2 = 0
       if len(idx_cc1+idx_cc2) != 0: num_coscos = max(idx_cc1+idx_cc2)
       else                        : num_coscos = 0
       if len(idx_sc1+idx_sc2) != 0: num_sincos = max(idx_sc1+idx_sc2)
       else                        : num_sincos = 0
       if len(idx_cs1+idx_cs2) != 0: num_cossin = max(idx_cs1+idx_cs2)
       else                        : num_cossin = 0
       if len(idx_ss1+idx_ss2) != 0: num_sinsin = max(idx_ss1+idx_ss2)
       else                        : num_sinsin = 0
    #--------------------------------#

    #--------------#
    # Get fterms   #
    #--------------#
    fterms  = [ ("const" ,'-' ,'-' )]
    fterms += [ ("cos"   ,idx ,'-' ) for idx in idx_c1]
    fterms += [ ("cos"   ,"-" ,idx ) for idx in idx_c2]
    fterms += [ ("sin"   ,idx ,'-' ) for idx in idx_s1]
    fterms += [ ("sin"   ,"-" ,idx ) for idx in idx_s2]
    fterms += [ ("coscos",idx1,idx2) for idx2 in idx_cc2 for idx1 in idx_cc1]
    fterms += [ ("sinsin",idx1,idx2) for idx2 in idx_ss2 for idx1 in idx_ss1]
    fterms += [ ("cossin",idx1,idx2) for idx2 in idx_cs2 for idx1 in idx_cs1]
    fterms += [ ("sincos",idx1,idx2) for idx2 in idx_sc2 for idx1 in idx_sc1]
    #--------------#


    #---------------#
    # Create tuples #
    #---------------#
    torsions   = (torsion1,torsion2,tsigma1,tsigma2)
    calcs      = (ttype,level,charge,multiplicity)
    pes        = (t1step,t2step,symmetry)
    fourier    = (fterms,nums,weight,ignore)
    statpoints = (tolerance,freqscal)
    tor2dns    = (dijvar,kmax,maxeigen)
    rovibpf    = (interpolation,integrationstep)
    iofiles    = (f_xyz   ,f_calcs ,f_ics  ,f_pes ,f_vfit  ,f_dfit,\
                  f_splist,f_spinfo,f_evals,f_sqrt,f_tables,f_pdf ,f_out)

    return torsions, calcs, pes, fourier, statpoints, tor2dns, rovibpf, temperatures, iofiles
#----------------------------------------------------------#
def readfiles(filenames,which):

    data = {}

    print "     Reading files:"
    for ff, string in zip(filenames,which):

        # Does file exists?
        if not os.path.isfile(ff):
           print"         * %s  | ERROR: unable to find it\n"%ff; sys.exit(EXITMESS)
        else:
           print "        - file name: %s"%ff

        #--------------#
        # Reading file #
        #--------------#
        info = ""
        if string == "xyz":
           structure, masslist = xyz2Struct(ff)
           symbols   = structure.get("symbols")
           print "          Number of atoms      %i"%structure.get("natoms")
           print "          Vibrational d.o.f    %i"%structure.get("nvib")
           data[string] = [structure,symbols,masslist]
           print "          Geometry (Cartesian coordinates in Angstrom):"
           xvec = structure.get("xcc")
           for idx in range(structure.get("natoms")):
               x,y,z = xvec[3*idx:3*idx+3] * cons.angstrom
               symbol = symbols[idx]
               mass   = structure.get("masslist")[idx] * cons.amu
               print "            %2s   %+10.5f  %+10.5f  %+10.5f   (%7.3f amu)"%(symbol,x,y,z,mass)

        if string == "calcs_scangeom":
           ilines, software = readfile_calcs(ff,"start_scangeom","end_scangeom")
           software = software.lower()
           print "          Software: %s"%software
           print "          Reference lines for the input file:"
           print "          -----------------------%"
           for line in ilines: print "          |"+line[:-1]
           print "          -----------------------%"
           data[string] = [ilines,software]

        if string == "calcs_sp":
           ilines0, software0 = readfile_calcs(ff,"start_sp0","end_sp0")
           ilines1, software1 = readfile_calcs(ff,"start_sp1","end_sp1")
           ilines2, software2 = readfile_calcs(ff,"start_sp2","end_sp2")
           software0 = software0.lower()
           software1 = software1.lower()
           software2 = software2.lower()
           print "          Reference '%s' input for sp0:"%software0
           print "          -----------------------%"
           for line in ilines0: print "          |"+line[:-1]
           print "          -----------------------%"
           print "          Reference '%s' input for sp1:"%software1
           print "          -----------------------%"
           for line in ilines1: print "          |"+line[:-1]
           print "          -----------------------%"
           print "          Reference '%s' input for sp2:"%software2
           print "          -----------------------%"
           for line in ilines2: print "          |"+line[:-1]
           print "          -----------------------%"
           data[string] = [ilines0,software0,ilines1,software1,ilines2,software2]

        if string == "pes":
           dict_pes, symbols, Eref = readfile_pes(ff)
           print "          number of points    %i"%len(dict_pes.keys())
           print "          reference energy    %.6f hartree"%Eref
           data[string] = [dict_pes, symbols, Eref]

        if string == "vfit":
           fterms, parameters, Eref_vfit = readfile_xfit(ff)
           print "          number of terms in the Fourier series  %i"%len(fterms)
           data[string] = [fterms, parameters, Eref_vfit]

        if string == "splist":
           dict_CPs, Eref_cp = readfile_splist(ff)
           list_CPs = dict_CPs.keys()
           list0    = [(dict_CPs[cpname][3],cpname) for cpname in list_CPs if dict_CPs[cpname][2] == 0]
           list1    = [(dict_CPs[cpname][3],cpname) for cpname in list_CPs if dict_CPs[cpname][2] == 1]
           list2    = [(dict_CPs[cpname][3],cpname) for cpname in list_CPs if dict_CPs[cpname][2] == 2]
           print "          number of minima : %2i"%len(list0)
           print "          number of saddles: %2i"%len(list1)
           print "          number of maxima : %2i"%len(list2)
           data[string] = [dict_CPs,Eref_cp]

        if string == "ics":
           icoords, nics, bonds = readfile_ics(ff)
           data[string] = [icoords, nics]
           for idx in range(0,len(icoords),5):
               line = ""
               for xx in range(5):
                   pos = idx + xx
                   if pos > len(icoords)-1: continue
                   tt, ic = icoords[pos]
                   if tt=="3": ic = "=".join("%i"%(a+1) for a in ic)
                   else      : ic = "-".join("%i"%(a+1) for a in ic)
                   line = line + "  %11s  "%ic
               print "     %s"%line

        if string == "evals":
           dataevals = readfile_evals(ff)
           evalues   = [ Ecm * cons.c2h for Ecm  in dataevals[0]]
           zpe_2DNS  = min(evalues)
           print "          zpe = %5.3f %s"%(zpe_2DNS*HU1,SU1)
           data[string] = [evalues,zpe_2DNS]

        print

    return data
#----------------------------------------------------------#
def readfile_ics(nricfile):
    '''
    Read file of non-redundant internal coordinates
    '''
    lines  = hf.readfile(nricfile)
    lines1 = hf.select_lines(lines,"start_connectivity","end_connectivity")
    lines2 = hf.select_lines(lines,"start_ics"         ,"end_ics")

    # Read connectivity
    connectivity = []
    for line in lines1:
        line = line.split("#")[0].strip()
        if line == "": continue
        bonds = line.split()
        for bond in bonds:
            bond = tuple([int(at)-1 for at in bond.split("-")])
            connectivity.append(bond)

    # Read internal coordinates
    ics_stretch, ics_abend, ics_lbend, ics_tors = [], [], [], []
    for line in lines2:
        line = line.split("#")[0].strip()
        if line == "": continue
        ics = line.split()
        for ic in ics:
           if "=" in ic:
              at1, at2, at3 = ic.split("=")
              ics_lbend.append( (int(at1)-1,int(at2)-1,int(at3)-1) )
           else:
              icoord = tuple([int(at)-1 for at in ic.split("-")])
              if len(icoord) == 2: ics_stretch.append(icoord)
              if len(icoord) == 3: ics_abend.append(icoord)
              if len(icoord) == 4: ics_tors.append(icoord)

    icoords = []
    nics = len(ics_stretch)+len(ics_abend)+2*len(ics_lbend)+len(ics_tors)
    for stretching in ics_stretch: icoords.append( ("1",stretching) )
    for abending   in ics_abend  : icoords.append( ("2",abending  ) )
    for lbending   in ics_lbend  : icoords.append( ("3",lbending  ) )
    for torsion    in ics_tors   : icoords.append( ("4",torsion   ) )

    return icoords, nics, connectivity
#----------------------------------------------------------#
def readfile_calcs(f_calcs,start,end):

    lines_ifile = []

    lines = open(f_calcs,'r')
    record = False
    for line in lines:
        # Check line and save
        if line.startswith(end): break
        if record: lines_ifile.append(line)
        if line.startswith(start):
           software = line.split()[1].lower()
           record = True
    lines.close()
    return lines_ifile, software
#----------------------------------------------------------#
def readfile_pes(f_pes,exclude=True):

    dict_pes  = {}
    symbols   = None
    list_Etot = []

    # Open file
    pesfile = open(f_pes,'r')

    # Read data
    for line in pesfile:
        data_in_line = line.split()
        # Read geometries
        if   len(data_in_line) == 1:
             natoms = int(data_in_line[0])
             symbols = []
             xvec    = []
        elif "Geometry" in data_in_line:
             Etot   = float(data_in_line[1])
             phi1   = float(data_in_line[2])*cons.D2R
             phi2   = float(data_in_line[3])*cons.D2R
             ppname = data_in_line[4]
             yesno  = data_in_line[5].lower()
             list_Etot.append(Etot)
        elif "FAILED" in data_in_line:
             continue
        else:
             symbols.append( data_in_line[0] )
             xvec += [float(x)/cons.angstrom for x in data_in_line[1:]]
             if len(xvec) == 3*natoms:
                if natoms == 1: continue
                if exclude and yesno == "no": continue
                dict_pes[ppname] = [phi1,phi2,Etot,np.array(xvec)]

    # Minimum energy?
    minE = min(list_Etot+[float("inf")])

    # Close file and return data
    pesfile.close()
    
    return dict_pes, symbols, minE
#----------------------------------------------------------#
def readfile_xfit(filename,which="V"):
    '''
    Read parameters of Fourier fitting
    '''
    if filename is None or not os.path.exists(filename): return None, None

    lines = hf.readfile(filename)
    if which == "d11"  : lines = hf.select_lines(lines,"start_fitdata11","end_fitdata11")
    if which == "d12"  : lines = hf.select_lines(lines,"start_fitdata12","end_fitdata12")
    if which == "d22"  : lines = hf.select_lines(lines,"start_fitdata22","end_fitdata22")
    if which == "sqrtD": lines = hf.select_lines(lines,"start_sqrtD"    ,"end_sqrtD"    )
    if which == "sqrtS": lines = hf.select_lines(lines,"start_sqrtS"    ,"end_sqrtS"    )

    fterms    = []
    parameters = []

    ref = None
    for line in lines:

        # Read line
        line = line.split("#")[0].strip()
        if line == "": continue
        if line.lower().startswith("reference"):
           ref = float(line.split()[1])
           continue
        term_type, idx_phi1, idx_phi2, value = line.split()
        value = float(value)

        # Assure the constant value is saved
        if term_type == "const":
           fterms.append( (term_type,idx_phi1, idx_phi2) )
           parameters.append(value)
           continue

        # If value is zero, ignore it
        if value == 0.0: continue
        if idx_phi1 != "-": idx_phi1 = int(idx_phi1)
        if idx_phi2 != "-": idx_phi2 = int(idx_phi2)

        # Save value
        fterms.append( (term_type,idx_phi1, idx_phi2) )
        parameters.append(value)

    return fterms, parameters, ref
#----------------------------------------------------------#
def readfile_splist(spfile):

    dict_CPs = {}

    Eref = None
    sp = open(spfile,'r')
    for line in sp:
        line = line.split("#")[0].strip()
        if line == "": continue
        if line.lower().startswith("eref "):
           Eref = float(line.split()[1])
           continue
        cp_type, phi1, phi2, E, okcol, cp_name = line.split()
        try:    E = float(E)
        except: E = float("inf")
        gtsfile = DIRSP + cp_name + ".gts"
        dict_CPs[cp_name] = [float(phi1)*cons.D2R,float(phi2)*cons.D2R,int(cp_type),E,gtsfile,okcol.upper()]
    sp.close()
    return dict_CPs, Eref
#----------------------------------------------------------#
def readfile_evals(f_evals):
    # Read lines
    iofile = open(f_evals,'r')
    lines = iofile.readlines()
    iofile.close()
    # Get kmax and dij
    kmax   = int(lines[0].split()[1])
    dijvar = lines[1].split()[1].lower()
    if dijvar == "yes": dijvar = True
    if dijvar == "no" : dijvar = False
    # Get evalues
    evals = []
    for line in lines[2:]:
        line = line.split("#")[0].strip()
        if line == "": continue
        evals.append( float(line.split()[1]) )
    # return
    return evals, kmax, dijvar
#----------------------------------------------------------#
#----------------------------------------------------------#



# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
#                                                          #
# SECTION (12) Functions for writing files                 #
#                                                          #
# >>>>>>>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<<<<<<< #
def string_inp(torsion1     , torsion2       , tsigma1     , tsigma2,  \
               level        , charge         , multiplicity, ttype  ,  \
               t1step       , t2step         , symmetry    ,           \
               weight       , ignore         ,                         \
               cos1         , cos2           , cos1cos2    , sin1sin2, \
               sin1         , sin2           , cos1sin2    , sin1cos2, \
               tolerance    , freqscal       ,                         \
               dijvar       , kmax           , maxeigen    ,           \
               interpolation, integrationstep,                         \
               Tlist                                                   ):

    string  = "#----------------------------------#\n"
    string += "# Torsional information            #\n"
    string += "#----------------------------------#\n"
    string += "start_torsions\n"
    string += "  torsion1         %-11s     # atoms involved in torsion 1\n"%torsion1
    string += "  torsion2         %-11s     # atoms involved in torsion 2\n"%torsion2
    string += "  tsigma1          %-3s             # torsional symmetry number of hindered rotor 1\n"%tsigma1
    string += "  tsigma2          %-3s             # torsional symmetry number of hindered rotor 2\n"%tsigma2
    string += "end_torsions\n"
    string += "#----------------------------------#\n"
    string += "# Calculations                     #\n"
    string += "#----------------------------------#\n"
    string += "start_calcs\n"
   #string += "  ttype            %-3s             # type of structure to analyze torsional anharmonicity (min,ts)\n"%ttype
    string += "  level            %-11s     # the calculation level\n"%level
    string += "  charge           %-2s              # charge\n"%charge
    string += "  multiplicity     %-2s              # spin multiplicity\n"%multiplicity
    string += "end_calcs\n"
    string += "#----------------------------------#\n"
    string += "# Torsional PES                    #\n"
    string += "#----------------------------------#\n"
    string += "start_pes\n"
    string += "  t1step           %-4s            # step in phi1 for scan calculation [degrees]\n"%t1step
    string += "  t2step           %-4s            # step in phi2 for scan calculation [degrees]\n"%t2step
    string += "  symmetry         %-4s            # Symmetry condition for PES: [a,b,c,ab,ac,bc,abc] or none\n"%symmetry
    string += "end_pes\n"
    string += "#----------------------------------#\n"
    string += "# Fitting details                  #\n"
    string += "#----------------------------------#\n"
    string += "start_fourier\n"
    string += "  weight           %-5s           #\n"%weight
    string += "  ignore           %-5s           # Set to zero coefficients with smaller absolute value (in cm^-1)\n"%ignore
    string += "  # Fourier Terms (Even)           #\n"
    string += "  cos1             %-4s            # i values in cos(i*Phi_1)\n"%cos1
    string += "  cos2             %-4s            # j values in cos(j*Phi_2)\n"%cos2
    string += "  cos1cos2         %-11s     # i,j values in cos(i*Phi_1) * cos(j*Phi_2)\n"%cos1cos2
    string += "  sin1sin2         %-11s     # i,j values in sin(i*Phi_1) * sin(j*Phi_2)\n"%sin1sin2
    string += "  # Fourier Terms (Odd)            #\n"
    string += "  sin1             %-4s            # i values in sin(i*Phi_1)\n"%sin1
    string += "  sin2             %-4s            # j values in sin(j*Phi_2)\n"%sin2
    string += "  cos1sin2         %-11s     # i,j values in cos(i*Phi_1) * sin(j*Phi_2)\n"%cos1sin2
    string += "  sin1cos2         %-11s     # i,j values in sin(i*Phi_1) * cos(j*Phi_2)\n"%sin1cos2
    string += "end_fourier\n"
    string += "#----------------------------------#\n"
    string += "# Search and Opt stationary points #\n"
    string += "#----------------------------------#\n"
    string += "start_statpoint\n"
    string += "  tolerance        %-5s           # step (in degrees) to explore torsional PES when looking for CPs\n"%tolerance
    string += "  freqscal         %-5s           # scaling factor for frequencies\n"%freqscal
    string += "end_statpoint\n"
    string += "#----------------------------------#\n"
    string += "# 2D-NS Hamiltonian                #\n"
    string += "#----------------------------------#\n"
    string += "start_tor2dns\n"
    string += "  dijvar           %-3s             # yes (dij not constant) or no (dij constant)\n"%dijvar
    string += "  kmax             %-3s             # check 2013-JChemPhys_138_134112, eq (14)\n"%kmax
    string += "  maxeigen         %-3s             # threshold for H eigenvalues (in cm^-1)\n"%maxeigen
    string += "end_tor2dns\n"
    string += "#----------------------------------#\n"
    string += "# Partition functions              #\n"
    string += "#----------------------------------#\n"
    string += "start_rovibpf\n"
    string += "  interpolation    %-7s         # fourier or spline order (1,3,5)\n"%interpolation
    string += "  integrationstep  %-7s         # integration dphi\n"%integrationstep
    string += "end_rovibpf\n"
    string += "#----------------------------------#\n"
    string += "# Working temperatures             #\n"
    string += "#----------------------------------#\n"
    string += "start_temperatures                 #\n"
    for idx in range(0,len(Tlist),3):
        line = "  ".join(["%6.1f"%T for T in Tlist[idx:idx+3]])
        line = "  " + line + " "*(33-len(line))+"#\n"
        string += line
    string += "end_temperatures                   #\n"
    string += "#----------------------------------#\n"
    string += "\n"
    return string

def writefile_inp(name):
    # Defaults
    torsion1     , torsion2       , tsigma1     , tsigma2,  \
    level        , charge         , multiplicity, ttype  ,  \
    t1step       , t2step         , symmetry    ,           \
    weight       , ignore         ,                         \
    cos1         , cos2           , cos1cos2    , sin1sin2, \
    sin1         , sin2           , cos1sin2    , sin1cos2, \
    tolerance    , freqscal       ,                         \
    dijvar       , kmax           , maxeigen    ,           \
    interpolation, integrationstep                            = get_defaults()
    # 
    torsion1 = "{1-2-3-4}"
    torsion2 = "{2-3-4-5}"
    tsigma1  = "{1}"
    tsigma2  = "{1}"
    level    = "{"+level   +"}"
    Tlist    = range(100,2001,100)
    #----------------------#
    string = string_inp(torsion1     , torsion2       , tsigma1     , tsigma2,  \
                        level        , charge         , multiplicity, ttype  ,  \
                        t1step       , t2step         , symmetry    ,           \
                        weight       , ignore         ,                         \
                        cos1         , cos2           , cos1cos2    , sin1sin2, \
                        sin1         , sin2           , cos1sin2    , sin1cos2, \
                        tolerance    , freqscal       ,                         \
                        dijvar       , kmax           , maxeigen    ,           \
                        interpolation, integrationstep,                         \
                        Tlist                                                   )
    #----------------------#
    inputname = name+".inp"
    thefile = open(inputname,'w')
    thefile.write(string)
    thefile.close()
#----------------------------------------------------------#
def writefile_ics(icoords,nricfile,bonds):

    nepl  = 5
    bonds = sorted(list(bonds))
    ic1   = sorted([ic for kind,ic in icoords if kind=="1"])
    ic2   = sorted([ic for kind,ic in icoords if kind=="2"])
    ic3   = sorted([ic for kind,ic in icoords if kind=="3"])
    ic4   = sorted([ic for kind,ic in icoords if kind=="4"])

    nric = open(nricfile,'w')

    nric.write("start_connectivity\n")
    for idx in range(0,len(bonds),nepl):
        line = ""
        for bond in bonds[idx:idx+nepl]:
            bond = "-".join(["%i"%(ii+1) for ii in bond])
            line = line + " %-7s "%bond
        nric.write("    %s\n"%line)
    nric.write("end_connectivity\n")
    nric.write("\n")

    nric.write("start_ics\n")

    nric.write("    #--------\n")
    nric.write("    # Stretches (%i):\n"%len(ic1))
    nric.write("    #--------\n")
    for idx in range(0,len(ic1),nepl):
        line = ""
        for ic in ic1[idx:idx+nepl]:
            ic   = "-".join(["%i"%(ii+1) for ii in ic])
            line = line + "%-11s "%ic
        nric.write("    %s\n"%line)

    nric.write("    #--------\n")
    nric.write("    # Bent bond angles (%i):\n"%len(ic2))
    nric.write("    #--------\n")
    for idx in range(0,len(ic2),nepl):
        line = ""
        for ic in ic2[idx:idx+nepl]:
            ic   = "-".join(["%i"%(ii+1) for ii in ic])
            line = line + "%-11s "%ic
        nric.write("    %s\n"%line)

    nric.write("    #--------\n")
    nric.write("    # Linear bond angles (2x%i):\n"%len(ic3))
    nric.write("    #--------\n")
    for idx in range(0,len(ic3),nepl):
        line = ""
        for ic in ic3[idx:idx+nepl]:
            ic   = "=".join(["%i"%(ii+1) for ii in ic])
            line = line + "%-11s "%ic
        nric.write("    %s\n"%line)

    nric.write("    #--------\n")
    nric.write("    # Improper torsions & Torsions (%i):\n"%len(ic4))
    nric.write("    #--------\n")
    for idx in range(0,len(ic4),nepl):
        line = ""
        for ic in ic4[idx:idx+nepl]:
            ic   = "-".join(["%i"%(ii+1) for ii in ic])
            line = line + "%-11s "%ic
        nric.write("    %s\n"%line)

    nric.write("end_ics\n")

    nric.close()
#----------------------------------------------------------#
def writefile_xfit(fterms,popt,ofile,mode="w"):

    output = open(ofile,mode)

    for term,parameter in zip(fterms,popt):
        term_type, idx1, idx2 = term
        if   term_type == "const" : output.write("const      %02s    %02s   %+12.5f\n"%("-","-",parameter))
        elif term_type == "cos"   :
             if idx1 == "-"       : output.write("cos        %02s    %02i   %+12.5f\n"%("-",idx2,parameter))
             if idx2 == "-"       : output.write("cos        %02i    %02s   %+12.5f\n"%(idx1,"-",parameter))
        elif term_type == "sin"   :
             if idx1 == "-"       : output.write("sin        %02s    %02i   %+12.5f\n"%("-",idx2,parameter))
             if idx2 == "-"       : output.write("sin        %02i    %02s   %+12.5f\n"%(idx1,"-",parameter))
        elif term_type == "coscos": output.write("coscos     %02i    %02i   %+12.5f\n"%(idx1,idx2,parameter))
        elif term_type == "cossin": output.write("cossin     %02i    %02i   %+12.5f\n"%(idx1,idx2,parameter))
        elif term_type == "sinsin": output.write("sinsin     %02i    %02i   %+12.5f\n"%(idx1,idx2,parameter))
        elif term_type == "sincos": output.write("sincos     %02i    %02i   %+12.5f\n"%(idx1,idx2,parameter))
        else: print "Unknown term!"; sys.exit(EXITMESS)
    output.close()
#----------------------------------------------------------#
def writefile_spinfo(f_spinfo,dict_CPs,icoords,Tlist,freqscal=1.0,masslist=None):

    output = open(f_spinfo,'w')
    output.write(" #---------------------------------------------------#\n")
    output.write(" # Information associated with each stationary point #\n")
    output.write(" #---------------------------------------------------#\n")
    output.write("\n")

    # Get min energy
    list_CPs = sp_sorting(dict_CPs)
    list_CPs = [cpname for cpname in list_CPs if os.path.exists(dict_CPs[cpname][4])]

    EMIN = float("inf")
    for cpname in list_CPs:
        cptype  = dict_CPs[cpname][2]
        gtsfile = dict_CPs[cpname][4]
        struct = gts2Struct(gtsfile,cpname,masslist,stype=cptype)
        if struct.get("Etot") < EMIN: EMIN = struct.get("Etot") 

    # Now, calculate some things at CPs
    print "     -------------------------------------------------- "
    print "     INFORMATION ASSOCIATED WITH EACH STATIONARY POINT: "
    print "     -------------------------------------------------- "
    print
    print "     Qrot: rotational partition function"
    print "     Qvib: harmonic-oscillator partition function from ZPE"
    print "     Erel: relative potential energy with respect to the global minimum"
    print "     Qrv = Qrot * Qvib * exp[-(Erel+ZPE)/(kT)]"
    print
    num0, num1, num2 = 1, 1, 1
    for cpname in list_CPs:
        cptype   = dict_CPs[cpname][2]
        gtsfile  = dict_CPs[cpname][4]
        # Get data
        data_sp = sp_analysis(gtsfile,icoords,cptype,freqscal,masslist,name=cpname)
        # Expand data
        (struct,energy,rij,dij,ccfreqs,zpe,ntfreqs,zpe_2b) = data_sp
        (rI1,rI2,l12) = rij
        (d11,d22,d12) = dij
        # To amu*angstrom^2
        rI1 = rI1 * cons.amu * cons.angstrom**2
        rI2 = rI2 * cons.amu * cons.angstrom**2
        l12 = l12 * cons.amu * cons.angstrom**2
        # To (amu*angstrom^2)^-1
        d11 = d11 / cons.amu / cons.angstrom**2
        d22 = d22 / cons.amu / cons.angstrom**2
        d12 = d12 / cons.amu / cons.angstrom**2
        # Generate string
        if cptype == 0: spinfo = " ==> MIN%02i --> '%s'\n\n"%(num0,cpname); num0 += 1
        if cptype == 1: spinfo = " ==> TS_%02i --> '%s'\n\n"%(num1,cpname); num1 += 1
        if cptype == 2: spinfo = " ==> MAX%02i --> '%s'\n\n"%(num2,cpname); num2 += 1
        # Energy
        energy  = ((energy-EMIN)*cons.h2c,(energy-EMIN)*HU1,SU1)
        spinfo = spinfo + "     * Erel: %.2f cm^-1 = %.2f %s\n"%energy
        spinfo = spinfo + "\n"
        # Moments of inertia
        imoments = tuple(struct.get("imoments"))
        spinfo = spinfo + "     * Principal moments of inertia (amu * angstrom^2):\n"
        spinfo = spinfo + "              %.3E  %.3E  %.3E\n"%imoments
        spinfo = spinfo + "\n"
        # Normal freqs
        spinfo = spinfo + "     * Vibrational frequencies (cm^-1)\n"
        for idx in range(0,len(ccfreqs),5):
            line = "  ".join([str(freq) for freq in ccfreqs[idx:idx+5]])
            spinfo = spinfo + "              %s\n"%line
        spinfo = spinfo + "\n"
        spinfo = spinfo + "     * ZPE = %.3f cm^-1 = %.3f %s\n"%(zpe*cons.h2c,zpe*HU1,SU1)
        spinfo = spinfo + "\n"
        # Rovibrational partition functions
        ptra,qrot,qvib,qele,Vadi = struct.get_pfns(Tlist)
        spinfo = spinfo + "     * Rotational and vibrational partition functions (tilde):\n"
        spinfo = spinfo + "            -----------------------------------------------------------------\n"
        spinfo = spinfo + "              T (K)  |    Qrot    |    Qvib    | (PES+ZPE)/(kT) |    Qrv     \n"
        spinfo = spinfo + "            -----------------------------------------------------------------\n"
        for idx in range(len(Tlist)):
            Qrv = qrot[idx]*qvib[idx]*np.exp(-(Vadi - EMIN)/cons.kB/Tlist[idx])
            linedata = (Tlist[idx], qrot[idx], qvib[idx],(Vadi - EMIN)/cons.kB/Tlist[idx],Qrv)
            spinfo = spinfo + "             %7.2f | %10.3E | %10.3E |    %8.2E    | %10.3E \n"%linedata
        spinfo = spinfo + "            -----------------------------------------------------------------\n"
        spinfo = spinfo + ""
        spinfo = spinfo + "\n"
        # D matrix
        spinfo = spinfo + "     * D matrix for the two torsions (amu*angstrom^2):\n"
        spinfo = spinfo + "           %+9.4f   %+9.4f\n"%(rI1,-l12)
        spinfo = spinfo + "           %+9.4f   %+9.4f\n"%(-l12,rI2)
        spinfo = spinfo + "\n"
        spinfo = spinfo + "     * det(D) = %+8.3f\n"%(rI1*rI2-l12*l12)
        spinfo = spinfo + "\n"
        spinfo = spinfo + "     * reduced I1 = %+9.4f\n"%rI1
        spinfo = spinfo + "     * reduced I2 = %+9.4f\n"%rI2
        spinfo = spinfo + "     * lambda_12  = %+9.4f\n"%l12
        spinfo = spinfo + "\n"
        spinfo = spinfo + "     * inv(D)*100 matrix:\n"
        spinfo = spinfo + "           %+9.4f   %+9.4f\n"%(d11,d12)
        spinfo = spinfo + "           %+9.4f   %+9.4f\n"%(d12,d22)
        spinfo = spinfo + "\n"
        # Non-torsional freqs
        spinfo = spinfo + "     * non-torsional frequencies (cm^-1)\n"
        for idx in range(0,len(ntfreqs),5):
            line = "  ".join([str(freq) for freq in ntfreqs[idx:idx+5]])
            spinfo = spinfo + "              %s\n"%line
        spinfo = spinfo + "\n"
        linetuple = (zpe_2b*cons.h2c,zpe_2b*HU1,SU1)
        spinfo = spinfo + "     * non-torsional ZPE = %.3f cm^-1 = %.3f %s\n"%linetuple
        spinfo = spinfo + "\n\n"
        # Print data
        for line in spinfo.split("\n")[:-1]: print "         "+line
        # Write data
        output.write(spinfo)

    output.close()
    print "     The information about each stationary point is stored at: %s"%f_spinfo
    print 
#----------------------------------------------------------#
#----------------------------------------------------------#



#-.-.-.-.-.-.-.-.-.-.-.-.#
# some provisional stuff #
#-.-.-.-.-.-.-.-.-.-.-.-.#
def q2dtor_interpes(name,inpdata,files2read):


    # Read files
    print "Read files"
    which      = [ "pes",  "vfit", "splist"]
    (symmetry,tsigma1,tsigma2) = inpdata
    dataifiles = readfiles(files2read,which)

    # PES and Fourier series at low level
    [dict_pes, symbols, Eref_pes] = dataifiles["pes"]
    [fterms, parameters, Eref_vfit] = dataifiles["vfit"]
    pnames = dict_pes.keys()
    if Eref_vfit is None: Eref_vfit = Eref_pes

    Fourier_PES = Fourier2D(fterms, parameters)

    # Expand PES, so interpolation works properly
    dict_pes = pessym_whole(name,dict_pes,symmetry,tsigma1,tsigma2,pprint=True)

    # SPs at high level
    [dict_CPs, Eref_cp] = dataifiles["splist"]

    # Expand SPs
    dict_CPs, list_CPs, initial_CPs = sp_allCPs(dict_CPs,symmetry,name,tsigma1,tsigma2)

    # Min Energy at low level
    minE_low = Eref_pes
    for sp in dict_CPs.keys():
        phi1, phi2, sptype, E = dict_CPs[sp][0:4]
        Vcm = Fourier_PES(phi1,phi2)*cons.c2h + Eref_vfit
        if Vcm < minE_low: minE_low = Vcm

    # Tesselation
    dict_CPs     = tess.replicate_points(name,dict_CPs)
    Dlist_CPs    = sp_sorting(dict_CPs)
    Dlist_points = [dict_CPs[cp_name][0:2] for cp_name in Dlist_CPs]
    delaunay     = Delaunay(Dlist_points)

    # Recalculate each point in PES
    new_dict_pes = {}
    for ppoint in dict_pes.keys():
        phi1, phi2, E_0, xcc = dict_pes[ppoint]
        p0 = np.array([phi1,phi2])

        # Look in delaunay
        simplex  = delaunay.find_simplex( (p0[0],p0[1]) )
        triangle = list(delaunay.simplices[simplex])
        idxA, idxB, idxC = triangle

        # Really close to a vertex
        spA = Dlist_CPs[idxA]
        spB = Dlist_CPs[idxB]
        spC = Dlist_CPs[idxC]
        pA  = np.array(delaunay.points[idxA])
        pB  = np.array(delaunay.points[idxB])
        pC  = np.array(delaunay.points[idxC])

        #---------#
        # Lambdas #
        #---------#
        area_A = hf.shoelace_triangle(pB,pC,p0)
        area_B = hf.shoelace_triangle(pA,pC,p0)
        area_C = hf.shoelace_triangle(pA,pB,p0)
        tot_area = area_A + area_B + area_C

        lambda_A = area_A / tot_area
        lambda_B = area_B / tot_area
        lambda_C = area_C / tot_area

        #----------------------------------#
        # Calculate factor using low-level #
        #----------------------------------#
        E_0 =  (E_0-minE_low) * cons.h2c
        E_A =  (Fourier_PES(pA[0],pA[1]) * cons.c2h + Eref_vfit - minE_low) * cons.h2c
        E_B =  (Fourier_PES(pB[0],pB[1]) * cons.c2h + Eref_vfit - minE_low) * cons.h2c
        E_C =  (Fourier_PES(pC[0],pC[1]) * cons.c2h + Eref_vfit - minE_low) * cons.h2c
        E_i = E_A*lambda_A + E_B*lambda_B + E_C*lambda_C

        # interpolated value can be zero
        if E_i != 0.0: factor = E_0 / E_i
        else         : factor = 1.0

        #-----------------------------#
        # Interpolation to high-level #
        #-----------------------------#
        E_A = dict_CPs[spA][3]
        E_B = dict_CPs[spB][3]
        E_C = dict_CPs[spC][3]
        E_i = (E_A*lambda_A + E_B*lambda_B + E_C*lambda_C)
        E_0 = E_i * factor
       #E_0 = E_i * 1.0

        # In hartree
        E_0 = Eref_cp + E_0*cons.c2h

        new_dict_pes[ppoint] = phi1, phi2, E_0, xcc
    # Expand data
    f_copy = files2read[0]+"_copy"
    if not os.path.exists(f_copy): shutil.copy(files2read[0],f_copy)
    writefile_pes(files2read[0],new_dict_pes,symbols,pnames)


def writefile_pes(f_pes,dict_pes,symbols,pnames):
    fails  = 0
    string = ""
    for point_name in sorted(dict_pes.keys()):
        phi1,phi2,energy,xvec = dict_pes[point_name]
        if energy is None:
           fails += 1
           string += " 1\n"
           string += " FAILED    0.0000000   %s\n"%point_name
           string += " XX       0.000000       0.000000       0.000000\n"
        else:
           phi1,phi2,energy,xvec = dict_pes[point_name]
           string += " %i\n"%len(symbols)
           if point_name in pnames:
              string += " Geometry  %+14.8f   %7.3f  %7.3f  %s  YES\n"%(energy, phi1*cons.R2D, phi2*cons.R2D, point_name)
           else:
              string += " Geometry  %+14.8f   %7.3f  %7.3f  %s  NO \n"%(energy, phi1*cons.R2D, phi2*cons.R2D, point_name)
           for idx in range(len(symbols)):
               symbol = symbols[idx]
               x,y,z = xvec[3*idx:3*idx+3] * cons.angstrom
               string += " %-2s   %+13.6f  %+13.6f  %+13.6f\n"%(symbol,x,y,z)
    pes_file = open(f_pes,'w')
    pes_file.write(string)
    pes_file.close()
    return fails




if __name__ == '__main__': main()

