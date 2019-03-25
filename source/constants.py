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

 This module contains different constants and dictionaries.
 These are the dictionaries:
      * dict_z2symbol    ;        Z      --> symbol
      * dict_symbol2z    ;        symbol --> Z
      * dict_covradii    ;        Z      --> cov_radious
      * dict_atomasses   ;        Z      --> ato_mass
      * dict_isomasses   ;        string --> mass
      * dict_UFURC       ;  molecularity --> conversion factor

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
'''

import math

#------------------------------------#
# Radians --> Degrees and vice versa #
#------------------------------------#
R2D      = 180.0 / math.pi           #
D2R      = math.pi / 180.0           #
TWOPI    = 2.0*math.pi               #
LINEAR   = 165.0 * D2R               #
#------------------------------------#

#----------------------------#
# Physical constants in S.I. #
#----------------------------#------------------------------------#
h_SI    = 6.62606896E-34     # Planck constant    [J*s]           #
NA_SI   = 6.022140857E+023   # Avogradro's number [mol^-1]        #
c0_SI   = 2.99792458E8       # speed of light     [m/s]           #
kB_SI   = 1.3806504E-23      # Boltzmann const    [J/K]           #
me_SI   = 9.10938356E-031    # mass of electron   [kg]            #
mC12_SI = 1.66053904E-027    # mass of C12        [kg]            #
cal_SI  = 4.184              # 1 cal = 4.184 J                    #
#----------------------------#------------------------------------#
R_SI    = kB_SI *  NA_SI     # Ideal gas constant [J*K^-1*mol^-1] #
hbar_SI =  h_SI / TWOPI      # Planck constant divided by 2pi     #
#----------------------------#------------------------------------#

#-----------------------------------------------#
# Conversion from a.u. to the given unit        #
# eg: a distance of 5 a.u. in meters is 5*meter #
#-----------------------------------------------#
joule    = 4.35974465E-018                      #
meter    = 0.52917721067E-10                    #
second   = 1.0 / 4.1341373337e+16               #
#-----------------------------------------------#
angstrom = meter * 1E10                         #
amu      = me_SI/mC12_SI                        #
kg       = me_SI                                #
cal      = joule / cal_SI                       #
kcal     = cal /1000.0                          #
kjmol    = joule * NA_SI / 1000.0               #
kcalmol  = kjmol / cal_SI                       #
calmol   = kcalmol * 1000.0                     #
meter3   = meter**3                             #
cm       = 100 * meter                          #
mL       = cm**3                                #
#-----------------------------------------------#

#-----------------------------#
# Physical constants in a.u.  #
#-----------------------------#
h    = 2*math.pi              #
hbar = 1.0                    #
NA   = NA_SI                  #
c0   = c0_SI / (meter/second) #
kB   = kB_SI / (joule)        #
me   = me_SI / (kg)           #
R    = kB * NA                #
h2c  = 1.0 / (h * c0) / cm    # hartree to 1/cm
c2h  = 1.0/h2c                # 1/cm to hartree
#-----------------------------#

#-----------------------------#
# Some ZERO definitions       #
#-----------------------------#
ZERO_hartree = 5E-6           # value of zero for energies
ZERO_geom    = 5E-6           # value to compare structures
ZERO_wavenum = 1.0 / (1.0/cm) # 1 cm^-1
ZERO         = 5e-07          # a less strict zero value
ZERO2        = 1e-10          # a more strict zero value
ZERO_DEG     = 0.05           #
ZERO_RAD     = 0.05 * D2R     #
#-----------------------------#

#--------------------------#
# Atomic number --> symbol #
#--------------------------#
# alkali metals, alkaline earth metals, post-transition metals, metalloids, non-metals, halogens, and noble gases
dict_z2symbol = { 1:"H ",                                                        2:"He",
                  3:"Li",  4:"Be",  5:"B ",  6:"C ",  7:"N ",  8:"O ",  9:"F ", 10:"Ne",
                 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P ", 16:"S ", 17:"Cl", 18:"Ar",
                 19:"K ", 20:"Ca", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr",
                 37:"Rb", 38:"Sr", 49:"In", 50:"Sn", 51:"Sb", 52:"Te", 53:"I ", 54:"Xe"}
# Transition metals
dict_z2symbol.update({21:"Sc",22:"Ti",23:"V ",24:"Cr",25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",
                      39:"Y ",40:"Zr",41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",47:"Ag",48:"Cd"})
# The 6th row
dict_z2symbol.update({55:"Cs",56:"Ba",57:"La",72:"Hf",73:"Ta",74:"W ",75:"Re",76:"Os",77:"Ir",78:"Pt",
                      79:"Au",80:"Hg",81:"Tl",82:"Pb",83:"Bi",84:"Po",85:"At",86:"Rn"})

# Lanthanides
dict_z2symbol.update({58:"Ce",59:"Pr",60:"Nd",61:"Pm",62:"Sm",63:"Eu",64:"Gd",65:"Tb",66:"Dy",67:"Ho",
                      68:"Er",69:"Tm",70:"Yb",71:"Lu"})

# The 7th row
dict_z2symbol.update({87:"Fr",88:"Ra",89:"Ac",104:"Ku",105:"Ha"})

# Actinides
dict_z2symbol.update({ 90:"Th", 91:"Pa", 92:"U ", 93:"Np", 94:"Pu", 95:"Am", 96:"Cm",
                       97:"Bk", 98:"Cf", 99:"Es",100:"Fm",101:"Md",102:"No",103:"Lr"})
#--------------------------#

#--------------------------#
# Symbol --> atomic number #
#--------------------------#
dict_symbol2z = dict((v.strip(),k) for k,v in dict_z2symbol.iteritems())

#-----------------#
# Reference radii #
#-----------------#
covraddi = []
covraddi = covraddi + [0.31,0.28,1.28,0.96,0.84,0.73,0.71,0.66,0.57,0.58] # atoms Z =  1 to 10
covraddi = covraddi + [1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76] # atoms Z = 11 to 20
covraddi = covraddi + [1.70,1.60,1.53,1.39,1.50,1.42,1.38,1.24,1.32,1.22] # atoms Z = 21 to 30
covraddi = covraddi + [1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75] # atoms Z = 31 to 40
covraddi = covraddi + [1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,1.42,1.39] # atoms Z = 41 to 50
covraddi = covraddi + [1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01] # atoms Z = 51 to 60
covraddi = covraddi + [1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87] # atoms Z = 61 to 70
covraddi = covraddi + [1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32] # atoms Z = 71 to 80
covraddi = covraddi + [1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,2.06] # atoms Z = 81 to 90
covraddi = covraddi + [2.00,1.96,1.90,1.87,1.80,1.69                    ] # atoms Z = 91 to 96
# to bohr
covraddi = [radius/angstrom for radius in covraddi]
# as a dict
dict_covradii = {}
for Z in range(96): dict_covradii[Z+1] = covraddi[Z]

#---------------#
# Atomic masses #
#---------------#
atomasses = []
atomasses = atomasses + [1.007825E+00,4.002600E+00,7.016000E+00,9.012180E+00,1.100931E+01] # atoms Z =   1 to   5
atomasses = atomasses + [1.200000E+01,1.400307E+01,1.599491E+01,1.899840E+01,1.999244E+01] # atoms Z =   6 to  10
atomasses = atomasses + [2.298980E+01,2.398504E+01,2.698153E+01,2.797693E+01,3.097376E+01] # atoms Z =  11 to  15
atomasses = atomasses + [3.197207E+01,3.496885E+01,3.994800E+01,3.896371E+01,3.996259E+01] # atoms Z =  16 to  20
atomasses = atomasses + [4.495592E+01,4.790000E+01,5.094400E+01,5.194050E+01,5.493810E+01] # atoms Z =  21 to  25
atomasses = atomasses + [5.593490E+01,5.893320E+01,5.793530E+01,6.292980E+01,6.392910E+01] # atoms Z =  26 to  30
atomasses = atomasses + [6.892570E+01,7.392190E+01,7.492160E+01,7.991650E+01,7.891830E+01] # atoms Z =  31 to  35
atomasses = atomasses + [8.391150E+01,8.491170E+01,8.790560E+01,8.990540E+01,8.990430E+01] # atoms Z =  36 to  40
atomasses = atomasses + [9.290600E+01,9.790550E+01,9.700000E+01,1.019037E+02,1.029048E+02] # atoms Z =  41 to  45
atomasses = atomasses + [1.059032E+02,1.069041E+02,1.139036E+02,1.149041E+02,1.199022E+02] # atoms Z =  46 to  50
atomasses = atomasses + [1.209038E+02,1.299067E+02,1.269044E+02,1.319042E+02,1.329054E+02] # atoms Z =  51 to  55
atomasses = atomasses + [1.379052E+02,1.389063E+02,1.399054E+02,1.409076E+02,1.419077E+02] # atoms Z =  56 to  60
atomasses = atomasses + [1.449127E+02,1.519197E+02,1.529212E+02,1.579241E+02,1.589253E+02] # atoms Z =  61 to  65
atomasses = atomasses + [1.639292E+02,1.649303E+02,1.659303E+02,1.689342E+02,1.739389E+02] # atoms Z =  66 to  70
atomasses = atomasses + [1.749408E+02,1.799465E+02,1.809480E+02,1.839509E+02,1.869557E+02] # atoms Z =  71 to  75
atomasses = atomasses + [1.919615E+02,1.929629E+02,1.949648E+02,1.969665E+02,2.019706E+02] # atoms Z =  76 to  80
atomasses = atomasses + [2.049744E+02,2.079766E+02,2.089804E+02,2.089824E+02,2.099871E+02] # atoms Z =  81 to  85
atomasses = atomasses + [2.220176E+02,2.230197E+02,2.260254E+02,2.270278E+02,2.320381E+02] # atoms Z =  86 to  90
atomasses = atomasses + [2.310359E+02,2.380508E+02,2.370482E+02,2.440642E+02,2.430614E+02] # atoms Z =  91 to  95
atomasses = atomasses + [2.470703E+02,2.470703E+02,2.510796E+02,2.520829E+02,2.570751E+02] # atoms Z =  96 to 100
atomasses = atomasses + [2.580986E+02,2.591009E+02,2.601053E+02                          ] # atoms Z = 101 to 103
# in atomic units
atomasses = [atomass/amu for atomass in atomasses]
# as a dict
dict_atomasses = {}
for Z in range(103): dict_atomasses[Z+1] = atomasses[Z]

#-----------------#
# Isotopic masses #
#-----------------#
dict_isomasses = {}
dict_isomasses[ "2H"] =  2.01410178 / amu # mass of deuterium [a.m.u.]
dict_isomasses[ "3H"] =  3.01604928 / amu # mass of tritium   [a.m.u.]
dict_isomasses["13C"] = 13.003355   / amu # mass of C13       [a.m.u.]

#----------------------------------------#
# User friendly units for rate constants #
#----------------------------------------#
dict_UFURC    = {}                                         
dict_UFURC[1] = 1.0/second # s^-1                 
dict_UFURC[2] =  mL/second # cm^3 / (molecule * s)

