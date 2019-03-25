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

 Module with diverse helper functions

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
''' 

#------------------------------------------------#
# >>           Importation  section           << #
#------------------------------------------------#
import sys
import numpy as np
from   numpy.linalg      import norm
from   scipy.interpolate import UnivariateSpline
from   scipy.optimize    import fmin
from   scipy.integrate   import quad
#------------------------------------------------#
import constants         as cons
#------------------------------------------------#
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #


#-----------------------#
# Cartesian Coordinates #
#-----------------------#
def cc2ms_x(x_cc,masslist,mu):
    x_ms = []
    for idx in range(len(masslist)):
        f = np.sqrt(1.0 * masslist[idx] / mu) # mass-scaled factor
        x,y,z = x_cc[3*idx:3*idx+3]
        x_ms += [ x*f , y*f , z*f ] 
    return np.array(x_ms)

def ms2cc_x(x_ms,masslist,mu):
    x_cc = []
    for idx in range(len(masslist)):
        f = np.sqrt(1.0 * masslist[idx] / mu) # mass-scaled factor
        x,y,z = x_ms[3*idx:3*idx+3]
        x_cc += [ x/f , y/f , z/f ] 
    return np.array(x_cc)

#-----------------------#
# Cartesian Gradient    #
#-----------------------#
def cc2ms_g(g_cc,masslist,mu):
    if g_cc is None: return None
    g_ms = []
    for idx in range(len(masslist)):
        f = np.sqrt(1.0 * masslist[idx] / mu) # mass-scaled factor
        gx,gy,gz = g_cc[3*idx:3*idx+3]
        g_ms += [ gx/f , gy/f , gz/f ]
    return np.array(g_ms)

def ms2cc_g(g_ms,masslist,mu):
    if g_ms is None: return None
    g_cc = []
    for idx in range(len(masslist)):
        f = np.sqrt(1.0 * masslist[idx] / mu) # mass-scaled factor
        gx,gy,gz = g_ms[3*idx:3*idx+3]
        g_cc += [ gx*f , gy*f , gz*f ]
    return np.array(g_cc)

#-----------------------#
# Force Constant Matrix #
#-----------------------#
def cc2ms_F(F_cc,masslist,mu):
    if F_cc is None: return None
    natoms = len(masslist)
    shape  = (3*natoms,3*natoms)
    F_ms = np.matrix( np.zeros( shape ) )
    for i in range(3*natoms):
        mi = masslist[int(i/3)]
        for j in range(3*natoms):
            mj = masslist[int(j/3)]
            f = mu / np.sqrt(mi*mj)
            F_ms[i,j] = F_cc[i,j] * f
    return F_ms

def ms2cc_F(F_ms,masslist,mu):
    if F_ms is None: return None
    natoms = len(masslist)
    shape  = (3*natoms,3*natoms)
    F_cc = np.matrix( np.zeros( shape ) )
    for i in range(3*natoms):
        mi = masslist[int(i/3)]
        for j in range(3*natoms):
            mj = masslist[int(j/3)]
            f = mu / np.sqrt(mi*mj)
            F_cc[i,j] = F_ms[i,j] / f
    return F_cc

def frange(start,end,dx,include_end=True):
    '''
    returns a list of floats
    '''
    start = float(start)
    end   = float(end  )
    dx    = float(dx   )
    nsteps = int( round( (end-start)/dx ) )
    if include_end: nsteps = nsteps + 1
    return [ start+i*dx for i in range(nsteps)]

def xvecformat(xvec,natoms,out="Nx3"):
    '''
    out: "Nx3" or "3Nx1"
    '''
    # Convert to list in case of numpy array
    try:
       xvec = xvec.tolist()
    except:
       pass

    # xvec in Nx3
    if len(xvec) == natoms:
       if out == "Nx3":
          return np.array(xvec,copy=True)
       if out == "3Nx1":
          final_xvec = []
          for x,y,z in xvec: final_xvec += [x,y,z]
          return np.array(final_xvec,copy=True)
              
    # xvec in 3Nx1
    elif len(xvec) == 3*natoms:
       if out == "Nx3":
          final_xvec = []
          for idx in range(natoms):
              x,y,z = xvec[3*idx:3*idx+3]
              final_xvec.append( [x,y,z] )
          return np.array(final_xvec,copy=True)
       if out == "3Nx1":
          return np.array(xvec,copy=True)

    sys.exit("Problems in xvecformat function")

def hessianformat(hessian,natoms):
    if hessian is None : return None
    if hessian is False: return None

    try   : hessian = hessian.tolist()
    except: pass

    lowtriangle  = (len(hessian) == 3*natoms*(3*natoms+1)/2)
    whole_matrix = (len(hessian) == 3*natoms)
    if lowtriangle : return ltriangle2matrix(hessian)
    if whole_matrix: return np.matrix(hessian)

    sys.exit("Problems in hessianformat function")


def shift2cm(x_cc,masslist):
    '''
    Function to shift to center of mass
    '''
  
    natoms = len(masslist)
    xvec   = xvecformat(x_cc,natoms,"3Nx1")
    # Get center of mass
    center_of_mass = np.array([0.0,0.0,0.0])
    for i in range(len(masslist)):
       center_of_mass = center_of_mass + masslist[i] * xvec[3*i:3*i+3]
    center_of_mass = center_of_mass / sum(masslist)
    # Shift the system so the center of mass is the origin
    shifted_xvec = np.copy(xvec)
    for i in range(len(masslist)):
        shifted_xvec[i*3+0] = xvec[i*3+0] - center_of_mass[0]
        shifted_xvec[i*3+1] = xvec[i*3+1] - center_of_mass[1]
        shifted_xvec[i*3+2] = xvec[i*3+2] - center_of_mass[2]
    return shifted_xvec

def ltriangle2matrix(ltriangle_list):
    '''
    * ltriangle_list: a list with all the lower triangular elements of a symmetrix matrix
    '''
    l = len(ltriangle_list)
    N = int( (-1 + np.sqrt(1+8*l)) / 2)
    index = 0
    matrix = np.zeros( (N,N) )
    for i in range(1,N+1):
        for j in range(1,i+1):
            matrix[i-1,j-1] = ltriangle_list[index]
            matrix[j-1,i-1] = ltriangle_list[index]
            index += 1
    return np.matrix(matrix)

def sign(number):
    if number >= 0.0: return +1.0
    else:             return -1.0

def delta_ij(i,j):
   if i==j: return 1.0
   else:    return 0.0

def isitlinear(x_cc,eps=1e-8):
    linear = True
    natoms = len(x_cc) / 3
    if natoms > 2:
       # 2nd atom from reference atom (1st one)
       ref_vec = x_cc[3:6] - x_cc[0:3]
       ref_vec = ref_vec / np.linalg.norm(ref_vec)
       # Rest of atoms from reference atom
       for atom in range(2,natoms):
           idx = 3*atom
           nvec = x_cc[idx:idx+3] - x_cc[0:3]
           nvec = nvec / np.linalg.norm(nvec)
           absdot = abs(np.dot(nvec,ref_vec))
           if absdot < (1.0 - eps): linear = False; break
    return linear


#--------------------------------------#
#   Functions for the calculation of   #
#   distances, angles and dihedrals    #
#--------------------------------------#
def ranged_angle(angle,v=1):
    '''
    angle: in rads
    if v == 1: returns angle in [  0,2pi]
    if v == 2: returns angle in [-pi, pi]
    '''

    # First, in [0,2pi]
    N = int( abs(angle / (2.0*np.pi)) )           # These three lines are useful

    if angle < 0.0: angle = angle + N*(2.0*np.pi) # if a big value in inserted (e.g. 1e8)
    else:           angle = angle - N*(2.0*np.pi) # as the conversion is faster

    while angle <  0.0                    : angle += 2.0*np.pi
    while angle >= 2.0*np.pi-cons.ZERO_RAD: angle -= 2.0*np.pi

    if abs(angle) <= cons.ZERO_RAD: angle = 0.0

    if v == 1: return angle

    # Now, in [-pi,pi]
    if angle > np.pi: angle -= 2.0*np.pi

    if v == 2: return angle

def angle_diff(angle1,angle2):
    angle1 = ranged_angle(angle1)
    angle2 = ranged_angle(angle2)

    d1 = angle2-2.0*np.pi-angle1
    d2 = angle2          -angle1
    d3 = angle2+2.0*np.pi-angle1

    if abs(d1) <= abs(d2) and abs(d1) <= abs(d3): return d1
    if abs(d2) <= abs(d1) and abs(d2) <= abs(d3): return d2
    if abs(d3) <= abs(d1) and abs(d3) <= abs(d2): return d3

def calc_distance(x1,x2):
    '''
    Calculates the distance 1-2
    Input:
      * x1: a numpy array (x,y,z) for point 1
      * x2: a numpy array (x,y,z) for point 2
    Returns:
      * dist: the distance (units of x1 and x2)
    '''
    # Convert to numpy array
    x1 = np.array(x1)
    x2 = np.array(x2)
    # Calculate distance
    vector = x1 - x2
    dist   = np.linalg.norm(vector)
    return dist

def calc_angle(x1,x2,x3):
    '''
    Calculates the angle 1-2-3
    Input:
      * x1: a numpy array (x,y,z) for point 1
      * x2: a numpy array (x,y,z) for point 2
      * x3: a numpy array (x,y,z) for point 3
    Returns:
      * theta: the angle (in radians)
    '''
    # Convert to numpy array
    x1 = np.array(x1)
    x2 = np.array(x2)
    x3 = np.array(x3)
    # Vectors from middle atom
    u = x1 - x2; norm_u = np.linalg.norm(u)
    v = x3 - x2; norm_v = np.linalg.norm(v)
    # Cosine and angle (from dot product)
    cosine = np.dot(u,v) / (norm_u * norm_v)
    theta  = np.arccos(cosine)
    return theta

def calc_dihedral(x1,x2,x3,x4):
    '''
    Calculates the dihedral angle 1-2-3-4
    Input:
      * x1: a numpy array (x,y,z) for point 1
      * x2: a numpy array (x,y,z) for point 2
      * x3: a numpy array (x,y,z) for point 3
      * x4: a numpy array (x,y,z) for point 4
    Returns:
      * theta: the dihedral angle (in radians) [-pi,pi] <-- Same criterium as gview and molden
    '''
    zero_rads = 0.001 * (np.pi/180.0)

    # Convert to numpy array
    x1 = np.array(x1)
    x2 = np.array(x2)
    x3 = np.array(x3)
    x4 = np.array(x4)

    # Vector 1-->2
    vec_12 = x2 - x1 ; vec_12 /= np.linalg.norm(vec_12)
    # Vector 2-->3
    vec_23 = x3 - x2 ; vec_23 /= np.linalg.norm(vec_23)
    # Vector 3-->4
    vec_34 = x4 - x3 ; vec_34 /= np.linalg.norm(vec_34)

    # Compute plane vectors
    n1 = np.cross(vec_12,vec_23); n1 /= np.linalg.norm(n1)
    n2 = np.cross(vec_23,vec_34); n2 /= np.linalg.norm(n2)

    # Vector perpendicular to (n1,vec23)
    m1 = np.cross(n1,vec_23)

    # Coordinates of n2 in this frame (n1,vec23,m1)
    x = np.dot(n1,n2)
    y = np.dot(m1,n2)

    # Angle
    theta = - np.arctan2(y,x)

    if abs(theta) < zero_rads: theta = 0.0

    return theta


def gen_rotmatrix(axis,theta):
    '''
    Generates the rotation matrix around axis.
    Input:
      * axis: the (x,y,z) coordinates of the
              direction vector of the axis
      * theta: the rotation angle (radians)
    Returns:
      * rot_matrix: a 3x3 numpy matrix

    Note: The rotation considers that the axis is
    situated at the origin

    Direction of rotation given by right-hand rule
    '''

    # Theta in [0,2*pi]
    while theta < 0.0:
          theta = theta + 2.0*np.pi
    while theta > 2*np.pi:
          theta = theta - 2.0*np.pi
    # Theta in [-pi,pi]
    if theta > np.pi:
       theta = -(2.0*np.pi - theta)

    ux, uy, uz = axis
    st = np.sin(theta)
    ct = np.cos(theta)

    R11 = ux*ux*(1.- ct) +    ct
    R12 = ux*uy*(1.- ct) - uz*st
    R13 = ux*uz*(1.- ct) + uy*st

    R21 = uy*ux*(1.- ct) + uz*st
    R22 = uy*uy*(1.- ct) +    ct
    R23 = uy*uz*(1.- ct) - ux*st

    R31 = uz*ux*(1.- ct) - uy*st
    R32 = uz*uy*(1.- ct) + ux*st
    R33 = uz*uz*(1.- ct) +    ct

    rot_matrix = np.matrix( [ [R11,R12,R13],[R21,R22,R23],[R31,R32,R33] ] )
    return rot_matrix, theta
#--------------------------------------#




#----------------------------------#
#   Functions related to strings   #
#----------------------------------#
def get_molformula(symbols):
    '''
    Returns the molecular formula of a given molecule
    Input:
      * symbols: the list of atomic symbols of the molecule
    Returns:
      * molformula: the string with the molecular formula
    '''
    # See symbols and count repetitions
    formula_dict = {symbol:symbols.count(symbol) for symbol in symbols}
    molformula = ""
    for key,value in sorted(formula_dict.items()):
        if value != 1: molformula = molformula + "%s(%i)"%(key,value)
        if value == 1: molformula = molformula + "%s"%(key)
    return molformula

def increase_string(the_string,length):
    if len(the_string) >= length: return the_string
    while len(the_string) < length:
          the_string = " "+the_string+" "
    if len(the_string) > length: the_string = the_string[:-1]
    return the_string

def remove_duplicates(a_list_of_floats,error=1e-8):
    # Copy list and sort it
    the_list = list(a_list_of_floats)
    the_list.sort()
    # Get position of duplicates
    repeated = [] 
    if len(the_list) > 1: 
       for idx in range(1,len(the_list)):
           previous_float = the_list[idx-1]
           current_float  = the_list[idx]
           if abs(current_float - previous_float) < error: repeated.append(idx)
    # Change values of positions to Nones
    for idx in repeated: the_list[idx] = None
    # Remove Nones from list
    removals = the_list.count(None)
    for removal in range(removals): the_list.remove(None)
    return the_list

def print_2Dmatrix(thematrix):
    nrows, ncols = thematrix.shape
    for row in range(nrows):
        col_string = ""
        for col in range(ncols):
            col_string = col_string + " %+6.3f "%float(thematrix[row,col])
        print col_string
#----------------------------------#



#-----------------------------------------#
#   Functions related to interpolations,  #
#   splines, extrema, and integration     #
#-----------------------------------------#
def localspline(xx,yy,idx,nps=(2,2),spl=3,local_data=False):
    '''
    idx may have any value if nps == "all"
    '''

    # Select points
    if nps == "all":
       local_xx = xx
       local_yy = yy
    else:
       idx_first = max(0,idx - nps[0])
       idx_last  = idx + nps[1] +1
       local_xx  = xx[idx_first : idx_last]
       local_yy  = yy[idx_first : idx_last]

    # Check spline order and get spline
    if spl >= len(local_xx): spl = len(local_xx) - 1
    the_spline = UnivariateSpline(local_xx,local_yy,k=spl,s=0,ext=3)

    if local_data is False: return the_spline
    if local_data is True:  return the_spline, local_xx, local_yy

def interpolate(xvalues,yvalues,x,nps=(3,2),spl=3,d=0):
    '''
    * x: value where interpolation has to be done
    * nps: number of point to create spline (left and right)
    * spl: order of the spline (spl <= 5)
    * d:  0 to get y(x), 1 to get y(x) and dy(x)
    '''
    if spl > 5: spl = 5

    target = None
    min_x  = +float("inf")
    max_x  = -float("inf")

    if x <= xvalues[0]:
       target = 0
       x      = xvalues[0]
    elif x >= xvalues[-1]:
       target = len(xvalues) -1
       x      = xvalues[-1]
    else:
       for idx in range(len(xvalues)):
           if xvalues[idx] < min_x: min_x = xvalues[idx]
           if xvalues[idx] > max_x: max_x = xvalues[idx]

           if xvalues[idx] < x or abs(xvalues[idx]-x) < cons.ZERO: target = idx
           else: break
       if target is None:
          errorstring = "Problems with interpolation! Asking for value at x=%f,"%x
          errorstring = errorstring + " but data points are in [%f,%f]"%(min_x,max_x)
          sys.exit(errorstring)

    # Create spline
    if nps == "all":
       function   = localspline(xvalues,yvalues,target,nps="all",spl=spl)
    else:
       function   = localspline(xvalues,yvalues,target,nps=(nps[0]-1,nps[1]),spl=spl)
    # Get values
    y  = function(x)
    if d == 0: return y
    # Get derivative
    derivative = function.derivative()
    dy = derivative(x)
    return y, dy

def interpolation_recalc(xvalues,yvalues,idx,nps=2,spl=3):
    '''
    Recalculates the value of a given point ignoring it
    from a list of points and interpolating it with a spline
    '''
    x_idx = xvalues[idx]
    xx = xvalues[idx-nps:idx] + xvalues[idx+1:idx+nps+1]
    yy = yvalues[idx-nps:idx] + yvalues[idx+1:idx+nps+1]
    the_spline = localspline(xx,yy,idx,nps="all",spl=3)
    new_y = the_spline(x_idx)
    return new_y

def obtain_extremum(xvalues , yvalues , xtr="min", full=False):
    '''
    * xtr: 'min" for minima, "max" for maxima
    '''
    assert len(xvalues) == len(yvalues), "x data and y data have different length"

    xx = np.array(xvalues)
    yy = np.array(yvalues)

    # If asked for maxima, multiply by -1
    if xtr=="max": yy = -1.0 * yy

    # Find "minima" in the discrete list
    targets = []
    if yy[0]  < yy[+1]: targets.append(0)
    for idx in range(1,len(xx)-1):
        y_a = yy[idx-1]
        y_b = yy[idx]
        y_c = yy[idx+1]
        if y_b < y_a and y_b < y_c: targets.append(idx)
    if yy[-1] < yy[-2]: targets.append( len(yy)-1 )

    #-------------------------------------#
    # Get points around and create spline #
    #-------------------------------------#
    
    crit_points = []
    for idx in targets:
        x_guess   = xx[idx]
        spline_22, xx_22, yy_22 = localspline(xx,yy,idx,nps=(2,2),spl=3,local_data=True)
        spline_33, xx_33, yy_33 = localspline(xx,yy,idx,nps=(3,3),spl=3,local_data=True)
        spline_44, xx_44, yy_44 = localspline(xx,yy,idx,nps=(4,4),spl=3,local_data=True)
        xm_22, ym_22, a,a,a = fmin(spline_22, x_guess, disp=False,full_output=True)
        xm_33, ym_33, a,a,a = fmin(spline_33, x_guess, disp=False,full_output=True)
        xm_44, ym_44, a,a,a = fmin(spline_44, x_guess, disp=False,full_output=True)

        if   ym_22 < ym_33 and ym_22 < ym_44:
           data_tuple = (ym_22,xm_22,xx_22,yy_22)
        elif ym_33 < ym_22 and ym_33 < ym_44:
           data_tuple = (ym_33,xm_33,xx_33,yy_33)
        else:
           data_tuple = (ym_44,xm_44,xx_44,yy_44)

        crit_points.append( data_tuple )

    crit_points.sort()
    x_absmin = crit_points[0][1]
    if xtr=="min": y_absmin = + crit_points[0][0]
    if xtr=="max": y_absmin = - crit_points[0][0]

    # Return data
    if full is False : return x_absmin, y_absmin


    # Correct data if asked for a maximum
    full_data = []
    for y_min, x_min, loc_xx, loc_yy in crit_points:
        x_alpha, x_omega = loc_xx[0], loc_xx[-1]
        if xtr=="min": thespline = localspline(loc_xx,     loc_yy,None,nps="all",spl=3)
        if xtr=="max": thespline = localspline(loc_xx,-1.0*loc_yy,None,nps="all",spl=3); y_min = - y_min
        full_data.append(  (x_min,y_min,thespline,x_alpha,x_omega)  )
    return full_data

def integrate(xx,yy,a,b):

    def function(x,xx,yy):
        return interpolate(xx,yy,x,nps=(3,2),spl=3,d=0)
        y_array = []
        for xi in x:
            y_array.append( interpolate(xx,yy,xi,nps=(3,2),spl=3,d=0) )
        return y_array
    int_value = quad(function,a,b, args=(xx,yy),full_output=1)[0]
    return int_value

def trapezoidal(f, a, b, n=100,args=None):

    if args is not None: args = tuple(args)
    else:                args = ()
    h = float(b - a) / n
    s = 0.0
    s += f( *( (a,)+args  )  ) /2.0
    for i in range(1, n):
        s += f( *((a + i*h,)+args) )
    s += f( *( (b,)+args  )  ) /2.0
    return s * h
#-----------------------------------------#



#-----------------------------#
#   Other type of functions   #
#-----------------------------#
def time2human(t, units):
    while True:
        if t < 1 and units == "secs":
           t = t * 1000
           units = "msecs"
        elif t > 60 and units == "secs":
           t = t / 60.0
           units = "mins"
        elif t > 60 and units == "mins":
           t = t / 60.0
           units = "hours"
        elif t > 24 and units == "hours":
           t = t / 24.0
           units = "days"
        else: break
    return (t, units)


def readfile(filename,hashtag=True,strip=True,ommitblank=True):
    '''
    read a file, excluding comments! (#)
    '''
    # Read file and exclude comments
    the_file = open(filename,'r')
    lines    = []
    for line in the_file:
        if line.endswith("\n"):
           line = line[:-1]
        if hashtag:
           line = line.split("#")[0]
        if strip:
           line = line.strip()
        if ommitblank and line == "":
           continue
        lines.append(line)
    the_file.close()

    return lines

def select_lines(lines,start,end,ignorecase=False):
    '''
    to use after readfile function
    '''
    selected_lines = []

    if ignorecase:
       start = start.lower()
       end   = end.lower()

    record = False
    for line in lines:
        if ignorecase: LINE = line.lower()
        else         : LINE = line
        # Check line and save
        if LINE.startswith(end): break
        if record: selected_lines.append(line)
        if LINE.startswith(start): record = True
    return selected_lines


def read_xyz(xyzfile):
    '''
    Read standard .xyz file
    Coordinates in angstroms
    '''
    xyz = open(xyzfile,'r')
    natoms   = int(xyz.readline())
    comment  = xyz.readline()
    xvector  = []
    symbols  = []
    masslist = []
    isotopic = False
    for line in range(natoms):
        data = xyz.readline().split()
        if len(data) == 4:
           symbol, x, y, z = data
           atonum = cons.dict_symbol2z[symbol]
           mass   = cons.dict_atomasses[atonum]
        if len(data) == 5:
           symbol, x, y, z, mass = data
           mass = float(mass) / cons.amu
           # Is an isotopic modification?
           atonum = cons.dict_symbol2z[symbol]
           mass2  = cons.dict_atomasses[atonum]
           diff = abs(mass-mass2)*cons.amu
           if diff > 0.01: isotopic = True
        masslist.append(mass)
        symbols.append(symbol)
        xvector += [float(x),float(y),float(z)]
    xyz.close()
    xvector = np.array(xvector)
    if "[BOHR]" in comment: xvector = xvector * cons.angstrom
    if not isotopic: masslist = None
    return xvector, symbols, masslist

def write_xyz(symbols,xvector,xyzfile,infoline="info line",mode="w"):
    '''
    Write standard .xyz file
    Coordinates in angstroms
    * mode = "a" or "w"
    '''
    natoms = len(symbols)

    # Prepare coordinate vector
    xyz_list = []
    if len(xvector) != natoms:
       for idx in range(0,len(xvector),3):
           x,y,z = xvector[idx:idx+3]
           xyz_list.append( np.array( (x,y,z) ) )
    else:
       xyz_list = xvector

    # Write file
    xyz = open(xyzfile,mode)
    xyz.write(" %i\n"%natoms)
    xyz.write("%s\n"%infoline)
    for idx in range(natoms):
        symbol = symbols[idx]
        x,y,z  = xyz_list[idx]
        xyz.write(" %2s   %+10.5f   %+10.5f   %+10.5f\n"%(symbol,x,y,z))
    xyz.close()

#-----------------#
# Wilson Matrices #
#-----------------#
'''
Some references:
 [1] D.F. McIntosh and K.H. Michaelian, The Wilson GF Matrix Method of Vibrational Analysis. Part I, II, and III
     Can.J.Spectros., 24, 1-10 (1979) ; Can.J.Spectros., 24, 35-40 (1979) ; Can.J.Spectros., 24, 65-74 (1979)
 [2] Ian H. Williams, Torsional Internal Coordinates in Normal Coordinate Calculations
     J. Mol. Spectros., 66, 288-301 (1977)
 [3] C.F. Jackels, Z. Gu, D.G. Truhlar, Reaction-path potential and vibrational frequencies in terms of curvilinear internal coordinates
     J. Chem. Phys. 102, 3188-3201 (1995)
 [4] Y-Y Chuang and D.G. Truhlar, Reaction-Path Dynamics in Redundant Internal Coordinates
     J. Phys. Chem. A, 102, 242-247 (1998)
 [5] V. Bakken and T. Helgaker, The efficient optimization of molecular geometries using redundant internal coordinates
     J. Chem. Phys., 117, 9160-9174 (2002)
'''
def get_BondVectors(x_cc, natoms):
    '''
    This functions calculates the distance
    between each pair of atoms (dij) and also
    the unit bond vector eij = (rj - ri)/dij
    '''
    bond_vectors = {}
    for i in range(natoms):
        xi, yi, zi = x_cc[3*i:3*i+3]
        for j in range(natoms):
            if i == j: continue
            xj, yj, zj = x_cc[3*j:3*j+3]
            eij = np.array( [float(xj-xi), float(yj-yi), float(zj-zi)] )
            dij = norm(eij)
            eij = np.array(eij) / norm(eij)
            bond_vectors[(i,j)] = (eij, dij)
    return bond_vectors

def angle_2vecs(vec1,vec2):
    '''
    This function calculates the angle between
    two vectors by using the dot product
    '''
    dot_prod = np.dot(vec1,vec2)
    mod_prod = norm(vec1) * norm(vec2)
    cos_angle = dot_prod / mod_prod
    # In case of numerical problems
    if   abs(cos_angle - 1.0) < cons.ZERO:
       angle = 0.0
    elif abs(cos_angle + 1.0) < cons.ZERO:
       angle = np.pi
    else:
       angle = np.arccos( cos_angle )
    return angle

def icBondStretch(bond_vectors,ij,natoms):
    '''
    Returns the row of the B matrix associated to
    the bond length between atom i and atom j and
    also the corresponding C matrix.
    Check ref [1] and [5].
    '''

    # Using nomenclature of reference [5]
    n = min(ij)
    m = max(ij)
    u, r = bond_vectors[(n,m)]

    ################################################
    # Calculating 1st derivatives: row of B matrix #
    #----------------------------------------------#
    B_row = np.zeros(3*natoms)
    #B_row = [0,0,0]*natoms
    for a in [m,n]:
        # Get zeta values
        if a == m: zeta_amn = +1.0
        if a == n: zeta_amn = -1.0
        for i in [0,1,2]:
            # Get B element
            dr_dai = zeta_amn * u[i]
            B_row[3*a+i] = dr_dai
    ################################################

    ######################################################
    # Calculating 2nd derivatives: 2D matrix of C tensor #
    #----------------------------------------------------#
    C_matrix = np.zeros( (3*natoms,3*natoms) )
    for a in [m,n]:
      for i in [0,1,2]:
          for b in [m,n]:
            for j in [0,1,2]:
               if C_matrix[3*a+i,3*b+j] != 0.0: continue
               # Get delta values
               if a == b: delta_ab = 1.0
               else:      delta_ab = 0.0
               if i == j: delta_ij = 1.0
               else:      delta_ij = 0.0
               # Get C element
               dr_daidbj = ((-1.0) ** delta_ab) * (u[i]*u[j] - delta_ij) / r
               # Save data in both positions
               C_matrix[3*a+i,3*b+j] = dr_daidbj
               C_matrix[3*b+j,3*a+i] = dr_daidbj
    ######################################################

    return [B_row], [C_matrix]

def icBondAngle(bond_vectors,ijk,natoms):
    '''
    Returns the row of the B matrix associated to the
    i-j-k angle bend and the corresponding C matrix.
    Check ref [1] and [5].
    '''

    # Using nomenclature of reference [5]
    m, o, n = ijk
    u, lambda_u = bond_vectors[(o,m)]
    v, lambda_v = bond_vectors[(o,n)]

    # Get internal coordinate: bond angle
    q = angle_2vecs(u,v)
    sinq = np.sin(q)
    cosq = np.cos(q)

    # Generation of w
    w = np.cross(u,v);  w = w / norm(w)

    uxw = np.cross(u,w)
    wxv = np.cross(w,v)

    ################################################
    # Calculating 1st derivatives: row of B matrix #
    #----------------------------------------------#
    B_row = np.zeros(3*natoms)
    for a in [m,o,n]:
        # Get zeta values
        if a == m: zeta_amo = +1.0; zeta_ano =  0.0
        if a == o: zeta_amo = -1.0; zeta_ano = -1.0
        if a == n: zeta_amo =  0.0; zeta_ano = +1.0
        for i in [0,1,2]:
            # Get B element
            dq_dai =  zeta_amo * uxw[i] / lambda_u + zeta_ano * wxv[i] / lambda_v
            B_row[3*a+i] = dq_dai
    ################################################

    ######################################################
    # Calculating 2nd derivatives: 2D matrix of C tensor #
    #----------------------------------------------------#
    C_matrix = np.zeros( (3*natoms,3*natoms) )
    if abs(sinq) < cons.ZERO: return B_row, C_matrix
    for a in [m,o,n]:
      for i in [0,1,2]:
        for b in [m,o,n]:
          for j in [0,1,2]:
            if C_matrix[3*a+i,3*b+j] != 0.0: continue
            # Define all delta and zeta values
            if a == m: zeta_amo = +1.0; zeta_ano =  0.0
            if a == o: zeta_amo = -1.0; zeta_ano = -1.0
            if a == n: zeta_amo =  0.0; zeta_ano = +1.0
            if b == m: zeta_bmo = +1.0; zeta_bno =  0.0
            if b == o: zeta_bmo = -1.0; zeta_bno = -1.0
            if b == n: zeta_bmo =  0.0; zeta_bno = +1.0
            if i == j: delta_ij = 1.0
            else:      delta_ij = 0.0
            # Get second derivative
            t1 = zeta_amo*zeta_bmo*(u[i]*v[j]+u[j]*v[i]-3*u[i]*u[j]*cosq+delta_ij*cosq)/(lambda_u**2 * sinq)
            t2 = zeta_ano*zeta_bno*(v[i]*u[j]+v[j]*u[i]-3*v[i]*v[j]*cosq+delta_ij*cosq)/(lambda_v**2 * sinq)
            t3 = zeta_amo*zeta_bno*(u[i]*u[j]+v[j]*v[i]-u[i]*v[j]*cosq-delta_ij)/(lambda_u*lambda_v*sinq)
            t4 = zeta_ano*zeta_bmo*(v[i]*v[j]+u[j]*u[i]-v[i]*u[j]*cosq-delta_ij)/(lambda_u*lambda_v*sinq)
            t5 = cosq / sinq * B_row[3*a+i] * B_row[3*b+j]
            dr_daidbj = t1 + t2 + t3 + t4 - t5
            C_matrix[3*a+i,3*b+j] = dr_daidbj
            C_matrix[3*b+j,3*a+i] = dr_daidbj
    return [B_row], [C_matrix]


def icLinealAngle(m,o,n,MON,natoms):
    '''
    '''
    def getB(m,o,n,k):
        om = m - o
        on = n - o
        dom = norm(om)
        don = norm(on)
        qk = ((n[2]-o[2])*(m[k]-o[k]) - (n[k]-o[k])*(m[2]-o[2])) / dom / don
        Bk = []
        Dk = []
        for a in ["m","o","n"]:
            for i in [0,1,2]:
                if a == "m": dam, dao, dan = 1.0, 0.0, 0.0
                if a == "o": dam, dao, dan = 0.0, 1.0, 0.0
                if a == "n": dam, dao, dan = 0.0, 0.0, 1.0
                if i ==  2 : di2 = 1.0
                else       : di2 = 0.0
                if i == k  : dik = 1.0
                else       : dik = 0.0
                # Numerator
                N_ai = dik*(dao-dan)*m[2] + di2*(dan-dao)*m[k] +\
                       dik*(dan-dam)*o[2] + di2*(dam-dan)*o[k] +\
                       dik*(dam-dao)*n[2] + di2*(dao-dam)*n[k]
                # Denominator
                D_ai = (dam-dao)*(m[i]-o[i]) * don / dom + \
                       (dan-dao)*(n[i]-o[i]) * dom / don
                Dk.append(D_ai)
                # Whole derivative
                B_ai = (N_ai - qk*D_ai) / (dom*don)
                Bk.append(B_ai)
        return np.array(Bk), np.array(Dk)

    def get_C(m,o,n,k,Bk=None,Dk=None):
        om = m - o
        on = n - o
        dom = norm(om)
        don = norm(on)
        if (Bk is None) or (Dk is None): Bk, Dk = get_B(m,o,n,k,True)
        Ck = np.zeros( (9,9) )
        for a in ["m","o","n"]:
            for i in [0,1,2]:
                if a == "m": dam, dao, dan = 1.0, 0.0, 0.0; A=0
                if a == "o": dam, dao, dan = 0.0, 1.0, 0.0; A=1
                if a == "n": dam, dao, dan = 0.0, 0.0, 1.0; A=2
                if i ==  2 : di2 = 1.0
                else       : di2 = 0.0
                if i == k  : dik = 1.0
                else       : dik = 0.0
                for b in ["m","o","n"]:
                    for j in [0,1,2]:
                        if b == "m": dbm, dbo, dbn = 1.0, 0.0, 0.0; B=0
                        if b == "o": dbm, dbo, dbn = 0.0, 1.0, 0.0; B=1
                        if b == "n": dbm, dbo, dbn = 0.0, 0.0, 1.0; B=2
                        if j ==  2 : dj2 = 1.0
                        else       : dj2 = 0.0
                        if j == k  : djk = 1.0
                        else       : djk = 0.0
                        # Term 1
                        term1 = (djk*di2-dik*dj2)*( (dan-dao)*dbm + (dam-dan)*dbo + (dao-dam)*dbn )
                        # Term 2
                        term2 = Bk[3*B+j] * Dk[3*A+i]
                        # Term 3
                        term3 = Bk[3*A+i] * Dk[3*B+j]
                        #
                        Ck[3*A+i,3*B+j] = (term1 - term2 - term3) / dom / don
        return np.matrix(Ck)

    M,O,N = MON
    #------------------------------------------#
    # Redefine new Cartesian coordinate system #
    #------------------------------------------#
    # Define new z axis
    new_z = n - m
    new_z = new_z / norm(new_z)

    # Define new y axis
    ref  = np.array([+1.0,+1.0,+1.0])
    prod = np.dot(new_z,ref) / norm(new_z) / norm(ref)
    if abs(prod) == 0.0:
       ref = np.array([+1.0,-1.0,0.0])
    new_y = np.cross(ref,new_z)
    new_y = new_y / norm(new_y)

    # Define new x axis
    new_x = np.cross(new_y,new_z)
    new_x = new_x / norm(new_x)

    # Define rotation matrix such as new_r = R * old_r,
    # considering new_r and old_r as column vectors
    R = np.matrix([new_x,new_y,new_z])
    Tbar = np.zeros( (9,9)  )
    Tbar[0:3,0:3] = R
    Tbar[3:6,3:6] = R
    Tbar[6:9,6:9] = R

    # Define position vectors with this new matrix
    m = R * np.matrix(m).transpose(); m = np.array( m.transpose().tolist()[0] )
    o = R * np.matrix(o).transpose(); o = np.array( o.transpose().tolist()[0] )
    n = R * np.matrix(n).transpose(); n = np.array( n.transpose().tolist()[0] )

    #-----------------------------------------------#
    # Get row of B matrix and 2D matrix of C tensor #
    #-----------------------------------------------#
    B_rows , C_matrices = [] , []
    for k in [0,1]:
        Bk, Dk = getB(m,o,n,k)
        Ck = get_C(m,o,n,k,Bk,Dk)
        # In old Cartesian coord. systems
        Bk = np.matrix(Bk) * Tbar.transpose()
        Ck = Tbar * Ck * Tbar.transpose()
        # Append data
        B_rows.append(np.array(Bk.tolist()[0]))
        C_matrices.append(Ck)
    
    #------------------------------------------#
    # Consider all the atoms to create B and C #
    #------------------------------------------#
    final_B = []
    for Bk in B_rows:
        fBk = np.zeros(3*natoms)
        for a in [M,O,N]:
            for i in [0,1,2]:
                if a == M: fBk[3*M+i] = Bk[0+i]
                if a == O: fBk[3*O+i] = Bk[3+i]
                if a == N: fBk[3*N+i] = Bk[6+i]
        final_B.append(fBk)

    final_C = []
    for Ck in C_matrices:
        fCk = np.zeros( (3*natoms,3*natoms) )
        for a in [M,O,N]:
            for i in [0,1,2]:
                # Select row in C
                if a == M: row = 0+i
                if a == O: row = 3+i
                if a == N: row = 6+i
                for b in [M,O,N]:
                    for j in [0,1,2]:
                        # Select col in C
                        if b == M: col = 0+j
                        if b == O: col = 3+j
                        if b == N: col = 6+j
                        fCk[3*a+i,3*b+j] = Ck[row,col]
        final_C.append(fCk)

    return final_B, final_C


def icDihedralAngle(bond_vectors,ijkl,natoms):
    '''
    Returns the row of the B matrix associated to the
    i-j-k-l torsion and the corresponding C matrix.
    Check ref [1] and [5].
    '''

    # Using nomenclature of reference [5]
    m, o, p, n = ijkl
    u, lambda_u = bond_vectors[(o,m)]
    v, lambda_v = bond_vectors[(p,n)]
    w, lambda_w = bond_vectors[(o,p)]

    uxw = np.cross(u,w)
    vxw = np.cross(v,w)

    cosPhi_u = np.dot(u,w) / norm(u) / norm(w)
    sinPhi_u = np.sqrt(1.0 - cosPhi_u**2)
    cosPhi_v = -np.dot(v,w) / norm(v) / norm(w)
    sinPhi_v = np.sqrt(1.0 - cosPhi_v**2)

    # Get internal coordinate: dihedral angle
    cosq = np.dot(uxw,vxw) / sinPhi_u / sinPhi_v
    if   abs(cosq - 1.0) < cons.ZERO:
       cosq = +1.0; q = 0.0
    elif abs(cosq + 1.0) < cons.ZERO:
       cosq = -1.0; q = np.pi
    else: q = np.arccos(cosq)

    ################################################
    # Calculating 1st derivatives: row of B matrix #
    #----------------------------------------------#
    B_row = np.zeros(3*natoms)
    for a in [m,o,p,n]:
        # Get zeta values
        if a == m: zeta_amo = +1.0; zeta_apn =  0.0; zeta_aop =  0.0
        if a == o: zeta_amo = -1.0; zeta_apn =  0.0; zeta_aop = +1.0
        if a == p: zeta_amo =  0.0; zeta_apn = +1.0; zeta_aop = -1.0
        if a == n: zeta_amo =  0.0; zeta_apn = -1.0; zeta_aop =  0.0
        for i in [0,1,2]:
            # Get B element
            dq_dai =  zeta_amo * uxw[i] / lambda_u / sinPhi_u / sinPhi_u + \
                      zeta_apn * vxw[i] / lambda_v / sinPhi_v / sinPhi_v + \
                      zeta_aop * uxw[i] * cosPhi_u / lambda_w / sinPhi_u / sinPhi_u + \
                      zeta_aop * vxw[i] * cosPhi_v / lambda_w / sinPhi_v / sinPhi_v
            B_row[3*a+i] = dq_dai
    ################################################

    ######################################################
    # Calculating 2nd derivatives: 2D matrix of C tensor #
    #----------------------------------------------------#
    C_matrix = np.zeros( (3*natoms,3*natoms) )
    for a in [m,o,p,n]:
      for i in [0,1,2]:
        for b in [m,o,p,n]:
          for j in [0,1,2]:
            if C_matrix[3*a+i,3*b+j] != 0.0: continue
            # Define all delta and zeta values
            if a == m: zeta_amo = +1.0; zeta_anp =  0.0; zeta_apo =  0.0; zeta_aop =  0.0; zeta_ano =  0.0
            if a == o: zeta_amo = -1.0; zeta_anp =  0.0; zeta_apo = -1.0; zeta_aop = +1.0; zeta_ano = -1.0
            if a == p: zeta_amo =  0.0; zeta_anp = -1.0; zeta_apo = +1.0; zeta_aop = -1.0; zeta_ano =  0.0
            if a == n: zeta_amo =  0.0; zeta_anp = +1.0; zeta_apo =  0.0; zeta_aop =  0.0; zeta_ano = +1.0

            if b == m: zeta_bom = -1.0; zeta_bmo = +1.0; zeta_bnp =  0.0; zeta_bop =  0.0; zeta_bpo =  0.0
            if b == o: zeta_bom = +1.0; zeta_bmo = -1.0; zeta_bnp =  0.0; zeta_bop = +1.0; zeta_bpo = -1.0
            if b == p: zeta_bom =  0.0; zeta_bmo =  0.0; zeta_bnp = -1.0; zeta_bop = -1.0; zeta_bpo = +1.0
            if b == n: zeta_bom =  0.0; zeta_bmo =  0.0; zeta_bnp = +1.0; zeta_bop =  0.0; zeta_bpo =  0.0
            zeta_bpn = - zeta_bnp

            if a == b: delta_ab = 1.0
            else:      delta_ab = 0.0
            # Get second derivative
            t01 = uxw[i]*(w[j]*cosPhi_u-u[j])/(lambda_u**2)/(sinPhi_u**4)
            t02 = uxw[j]*(w[i]*cosPhi_u-u[i])/(lambda_u**2)/(sinPhi_u**4)

            t03 = vxw[i]*(w[j]*cosPhi_v+v[j])/(lambda_v**2)/(sinPhi_v**4)
            t04 = vxw[j]*(w[i]*cosPhi_v+v[i])/(lambda_v**2)/(sinPhi_v**4)

            t05 = uxw[i]*(w[j]-2*u[j]*cosPhi_u+w[j]*cosPhi_u**2)/(2*lambda_u*lambda_w*sinPhi_u**4)
            t06 = uxw[j]*(w[i]-2*u[i]*cosPhi_u+w[i]*cosPhi_u**2)/(2*lambda_u*lambda_w*sinPhi_u**4)

            t07 = vxw[i]*(w[j]+2*v[j]*cosPhi_v+w[j]*cosPhi_v**2)/(2*lambda_v*lambda_w*sinPhi_v**4)
            t08 = vxw[j]*(w[i]+2*v[i]*cosPhi_v+w[i]*cosPhi_v**2)/(2*lambda_v*lambda_w*sinPhi_v**4)

            t09 = uxw[i]*(u[j]+u[j]*cosPhi_u**2-3*w[j]*cosPhi_u+w[j]*cosPhi_u**3) / (2*lambda_w**2*sinPhi_u**4)
            t10 = uxw[j]*(u[i]+u[i]*cosPhi_u**2-3*w[i]*cosPhi_u+w[i]*cosPhi_u**3) / (2*lambda_w**2*sinPhi_u**4)

            t11 = vxw[i]*(v[j]+v[j]*cosPhi_v**2+3*w[j]*cosPhi_v-w[j]*cosPhi_v**3) / (2*lambda_w**2*sinPhi_v**4)
            t12 = vxw[j]*(v[i]+v[i]*cosPhi_v**2+3*w[i]*cosPhi_v-w[i]*cosPhi_v**3) / (2*lambda_w**2*sinPhi_v**4)

            if i != j and a != b:
               k = [0,1,2]; k.remove(i); k.remove(j); k = k[0]
               t13 = (w[k]*cosPhi_u-u[k]) / lambda_u / lambda_w / sinPhi_u**2
               t14 = (w[k]*cosPhi_v+v[k]) / lambda_v / lambda_w / sinPhi_v**2
            else:
               t13 = 0.0
               t14 = 0.0

            dr_daidbj = zeta_amo*zeta_bmo*(t01 + t02) + \
                        zeta_anp*zeta_bnp*(t03 + t04) + \
                       (zeta_amo*zeta_bop+zeta_apo*zeta_bom)*(t05 + t06) + \
                       (zeta_anp*zeta_bpo+zeta_aop*zeta_bpn)*(t07 + t08) +\
                        zeta_aop*zeta_bpo*(t09 + t10) + \
                        zeta_aop*zeta_bop*(t11 + t12) + \
                       (zeta_amo*zeta_bop+zeta_apo*zeta_bom)*(1-delta_ab)*(j-i) / (-2.)**(abs(j-i))*t13 +\
                       (zeta_anp*zeta_bpo+zeta_aop*zeta_bpn)*(1-delta_ab)*(j-i) / (-2.)**(abs(j-i))*t14
            # Save data in both positions
            C_matrix[3*a+i,3*b+j] = dr_daidbj
            C_matrix[3*b+j,3*a+i] = dr_daidbj

    return [B_row], [C_matrix]
#########################################


def get_B_and_C(xvector, natoms, icoords):

    # Get all unit bond vectors
    bond_vectors = get_BondVectors(xvector,natoms)
    B, C = [], []

    # Get bond stretch rows of B
    stretch_list = [ic for kind,ic in icoords if kind == "1"]
    for ij_list in stretch_list:
        rowB, matrC = icBondStretch(bond_vectors,ij_list,natoms)
        B += rowB; C += matrC

    # Get angular bending rows of B
    abend_list   = [ic for kind,ic in icoords if kind == "2"]
    for ijk_list in abend_list:
        rowB, matrC = icBondAngle(bond_vectors,ijk_list,natoms)
        B += rowB; C += matrC

    # Get lineal bending rows of B
    lbend_list   = [ic for kind,ic in icoords if kind == "3"]
    for ijk_list in lbend_list:
        i,j,k = ijk_list
        r_m = np.array(xvector[3*i:3*i+3],copy=True); r_m = np.array(r_m)
        r_o = np.array(xvector[3*j:3*j+3],copy=True); r_o = np.array(r_o)
        r_n = np.array(xvector[3*k:3*k+3],copy=True); r_n = np.array(r_n)
        rowB, matrC = icLinealAngle(r_m,r_o,r_n,ijk_list,natoms)
        B += rowB; C += matrC
    
    # Get torsion rows of B
    torsion_list = [ic for kind,ic in icoords if kind == "4"]
    for ijkl_list in torsion_list:
        rowB,matrC = icDihedralAngle(bond_vectors,ijkl_list,natoms)
        B += rowB; C += matrC
 
    B = np.matrix(B)

    return B, C



#-------------#
# Point Group #
#-------------#
def get_pgs(atom_num,atom_mass,geom_xyz,toldist=0.05,tolsym=3e-2,epsilon=3e-2):
  '''
  This module finds the symmetry point
  group symmetry of a molecule 
  
  *---------------------------------------*
  | Main Author:  Antonio Fernandez-Ramos |
  | Last Update:  May 10th 2017 (by DFC)  |
  *---------------------------------------*
  '''
  
  #  -------Parameters----------
  # atomic numbers
  # atomic masses
  # Geometry as [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3] ....
  # toldist : tolerance of difference of distances
  # tolsym : tolerance in the symmetry operation
  # epsilon : tolerance in the cross products
  #           usually used to check if two C2 axis are parallel

  # /////  Matrix representation of symmetry operators \\\\\
  
  # in case geom_xyz is a 3N list 
  if len(geom_xyz) == 3*len(atom_num):
     xyz = []
     for idx in range(0,len(geom_xyz),3):
         x,y,z = geom_xyz[idx+0:idx+3]
         xyz.append( [x,y,z] )
     geom_xyz = xyz

  natom = len(atom_num)

  if natom == 1: return 'C1',1

  def inversion():
      i_center = ([[-1.,0.,0.],
                   [0.,-1.,0.],
                   [0.,0.,-1.]])
      return i_center
  
  def reflex(uvw,main_plane):
  # reflection about a plane perpendicular to uvw
    plane = np.zeros( (3,3) )
    if main_plane == 'xy': uvw = [0.,0.,1.]
    if main_plane == 'xz': uvw = [0.,1.,0.]
    if main_plane == 'yz': uvw = [1.,0.,0.]
    x = uvw[0]
    y = uvw[1]
    z = uvw[2]
    plane[0,0] = 1.-2.*x**2 
    plane[1,1] = 1.-2.*y**2 
    plane[2,2] = 1.-2.*z**2 
    plane[0,1] = -2.*x*y
    plane[0,2] = -2.*x*z
    plane[1,2] = -2.*y*z
    plane[1,0] = plane[0,1] 
    plane[2,0] = plane[0,2] 
    plane[2,1] = plane[1,2] 
    return plane
  
  def Cngen(n,u):
  # rotacion about the u vector
     Cngen = np.zeros( (3,3) )
     t=2.*np.pi/float(n)
     ux = u[0]
     uy = u[1]
     uz = u[2]
     cosu = np.cos(t)
     sinu = np.sin(t)
     Cngen[0,0] = cosu + ux**2*(1.-cosu)
     Cngen[1,1] = cosu + uy**2*(1.-cosu)
     Cngen[2,2] = cosu + uz**2*(1.-cosu)
     Cngen[0,1] = ux*uy*(1.-cosu)-uz*sinu
     Cngen[0,2] = ux*uz*(1.-cosu)+uy*sinu
     Cngen[1,2] = uy*uz*(1.-cosu)-ux*sinu
     Cngen[1,0] = ux*uy*(1.-cosu)+uz*sinu
     Cngen[2,0] = ux*uz*(1.-cosu)-uy*sinu
     Cngen[2,1] = uy*uz*(1.-cosu)+ux*sinu
     return Cngen

### ///// END \\\\
  
  def cdmass(mass,xyz,natom):
   # center of masses
     cdm = [0.,0.,0.]
     tot_mass = sum(mass)
     for i in range(0,natom):
       cdm[0] += mass[i]*xyz[i,0]/tot_mass
       cdm[1] += mass[i]*xyz[i,1]/tot_mass
       cdm[2] += mass[i]*xyz[i,2]/tot_mass
     for i in range(0,natom):
       xyz[i,0] += - cdm[0]
       xyz[i,1] += - cdm[1]
       xyz[i,2] += - cdm[2]
     return xyz
  
  def I_diff(fx,fy,tol):
      abs_rel = np.fabs(fx-fy)/fy*100.
      if abs_rel <= tol: 
        return True
      else:
        return False
  
  def calc_tensor(mass,xyz,natom):
    # tensor of inertia 
      tensor_iner = np.zeros( (3,3) )
      for i in range(0,natom):
        tensor_iner[0,0] += mass[i]*(xyz[i,1]**2+xyz[i,2]**2)
        tensor_iner[1,1] += mass[i]*(xyz[i,0]**2+xyz[i,2]**2)
        tensor_iner[2,2] += mass[i]*(xyz[i,0]**2+xyz[i,1]**2)
        tensor_iner[0,1] += -mass[i]*xyz[i,0]*xyz[i,1]
        tensor_iner[0,2] += -mass[i]*xyz[i,0]*xyz[i,2]
        tensor_iner[1,2] += -mass[i]*xyz[i,1]*xyz[i,2]
      tensor_iner[1,0] = tensor_iner[0,1]
      tensor_iner[2,0] = tensor_iner[0,2]
      tensor_iner[2,1] = tensor_iner[1,2]
      #eigvals, eigvecs = np.linalg.eigvalsh(tensor_iner)
      Im,Ivec =  np.linalg.eigh(tensor_iner)
      Ivec = Ivec.transpose()
      return Im,Ivec
  
  
  def dist_mat(xyz,natom):
  # Matrix of distances
     dmx = []
     for i in range(0,natom):
       vtmp = xyz[i]
       for j in range (0,natom):
         vtmp2 = xyz[j]
         vx=vtmp-vtmp2
  # Norm + storage
         dmx.append(np.linalg.norm(vx))
     dmx=np.reshape(dmx,(natom,natom))
     return dmx
  
  def dist_c(v1,v2):
     L_novec = False
     vxn = []
     vx = (v1+v2)/2
     vxn=np.linalg.norm(vx)
     if vxn == 0.: 
       L_novec = True
       return vx,L_novec
     vx /= vxn
     return vx,L_novec 
     
  def dist_v(v1,v2):
     L_nosig = False
     vxn = []
     vx = (v1-v2)
     vxn=np.linalg.norm(vx)
     if vxn == 0.: 
       L_nosig = True
       return vx,L_nosig
     vx /= vxn
     return vx,L_nosig
    
  
  def get_sea(dict_sea,mat_geo,idx):
  # Mass and geometry of each SEA
    coor_sea = []
    mass_sea = []
    idum=len(dict_sea[idx])
    for i in range(idum):
          k = dict_sea[idx][i]
          coor_sea.append(mat_geo[k])
          mass_sea.append(atom_mass[k])
    coor_sea=np.reshape(coor_sea,(idum,3))
    return idum,mass_sea,coor_sea
  
  def compare_geom(A,B):
    nrows, ncols = A.shape
    tot_error = 0.0
    indices_comparison = []

    for row_i in range(nrows):
        Arow = A[row_i,:]
        diff    = float("inf")
        idx     = None
        for row_j in range(nrows):
            Brow = B[row_j,:]
            error = np.linalg.norm(Arow - Brow)
            if error < diff:
               diff = error
               idx  = row_j
        tot_error += diff
        indices_comparison.append( (row_i,idx)  )
    tot_error = 1.0 * tot_error / nrows
    return tot_error 

  def get_sym(dict_sea,mat_geo,mat_sym,num_set,strsym,tolsym):
    # returns True or False, i.e. if the molecule is invariant after a given
    # symmetry operation (strsym) 
    total_error = 0.
    mat_tmp = mat_geo.T 
    mat_tmp = mat_sym*mat_tmp
    mat_tmp = mat_tmp.T
    L_sym = False
    for i in range(num_set):
      idx = i 
      isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx)
      isea,mass_tmp,coor_tmp =  get_sea(dict_sea,mat_tmp,idx)
      coor_sea = np.array(coor_sea)
      error_sym_oper = compare_geom(coor_sea,coor_tmp)
      total_error += error_sym_oper
    total_error /= num_set
#    print 'Error ',total_error
    if total_error < tolsym:
      L_sym = True

    return L_sym  
  
  def get_c2(dict_sea,mat_geo,num_set,uvw,L_round,epsilon):
    #
    # uvw is for all type of molecules except spheric
    # L_round = True for spherical molecules 
    # search for C2 axis: 1) passing through the middle of a bond
    #                     2) passing through atoms
    coor_sea = []
    mass_sea = []
    new_c2 = 0
    new_c2_red = 0  
    list_tot = []
    list_red = []
    dxc2 = []

    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx)
      if isea > 1:
        for i in range(isea):
          for j in range(i+1,isea):
            dx,L_novec = dist_c(coor_sea[i],coor_sea[j])
            if L_novec:
              continue
            elif L_round:
              dxm.append(dx)
              new_c2_red += 1
            else:
              vcros = np.cross(dx,uvw) 
              vcrosn=np.linalg.norm(vcros)
              if vcrosn > 1.-epsilon: 
                new_c2_red += 1
                dxm.append(dx)
          vnorm=np.linalg.norm(coor_sea[i])
          if vnorm > epsilon:
            dxm.append(coor_sea[i]/vnorm)
            new_c2_red += 1
  #
  # ---remove the redundancies from the C2 axes
  #
    for i in range(new_c2_red): list_tot.append(i)  
    set_tot = set(list_tot)
    if new_c2_red > 1:
      for i in range(new_c2_red):
        vdxi = dxm[i]
        for j in range(i+1,new_c2_red): 
          vdxj = dxm[j]
          vcros = np.cross(vdxi,vdxj) 
          vcrosn=np.linalg.norm(vcros)
          if j in list_red:
            continue
          else:
            if vcrosn < epsilon: 
              list_red.append(j)
    set_red = set(list_red)
    set_not_red = set_tot - set_red
    dxc2.append([dxm[i] for i in set_not_red])
    new_c2 = len(set_not_red)
    dxc2=np.reshape(dxc2,(new_c2,3))
    return new_c2,dxc2

  def get_sigma_d(dict_sea,mat_geo,num_set,mat_c2_axis,numc2,tolsym):
  #
  #  Diagonal planes: considers that are perpendicular to the
  #    vectors defined by the C2 axes and the coordinate of
  #    an atom of the SEA 
  #    Important for Td molecules
  #
    coor_sea = []
    mass_sea = []
    mat_sigma_d = [] 
    num_sigma_d = 0
    vdxi = []
    vdxj = []
    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx)
      if isea > 1:
        for i in range(isea):
          vdxi = coor_sea[i]
          for j in range(numc2):
            vdxj = mat_c2_axis[j]
            vnormal = np.cross(vdxi,vdxj)
            vnorm=np.linalg.norm(vnormal)
            if vnorm > tolsym:
              vnormal /= vnorm 
              mat_sigma_d = reflex(vnormal,' ')
              L_sigma_d = get_sym(dict_sea,mat_geo,mat_sigma_d,num_set,'sigma_d',tolsym) 
              if L_sigma_d:
                num_sigma_d += 1
            else:
              continue
    return num_sigma_d 
    
  
  #####initial variables##################
  
  L_cn = False
  L_c2_p = False
  L_c2_sp = False
  L_sigma_v = False
  L_sigma_h = False
  L_sigma_x = False
  L_i = False
  L_s2n = False
  L_cs = False
  cnrot = []
  cn_from_I = [1]
  sigma_h = []
  sigma_v = []
  sigma_x = []
  sn_mat = []
  i_center = inversion()
  udum = []
  uvw = [0.,0.,0.]
  axis = ['x','y','z']
  plane = ['yz','xz','xy']
  dxm = [] 
  
  
  mat_geo = np.matrix(geom_xyz)
  # Puts the molecule in the center of masses
  mat_geo = cdmass(atom_mass,mat_geo,natom)
  # Evaluates the tensor of inertia
  Ipmol,Ipvec = calc_tensor(atom_mass,mat_geo,natom)
  # Rotates the molecule so Ia, Ib and Ic coincide with i,j,k
  mat_geo = Ipvec * mat_geo.T
  mat_geo = mat_geo.T
  
  
  dist = dist_mat(mat_geo,natom)
  distv = []
  # sum the distances of each row
  for i in range(natom):
     rkk = sum(dist[i,:])/float(natom)
     distv.append(rkk)
  
  # now it checks for possible SEA
  # index of SEAs stored in the dictionary
  dict_sea = { }
  idx = -1
  sea_at_idx = []
  for i in range(natom):
     next = False
     if i not in sea_at_idx: 
       idx += 1
     for j in range(i,natom):
       
       if atom_mass[i] == atom_mass[j] and atom_num[i] == atom_num[j] and np.fabs(distv[i] - distv[j])< toldist and j not in sea_at_idx:
         sea_at_idx.append(j)
         dict_sea[idx]  = dict_sea.get(idx,[]) + [j]
  
  # number of SEAs
  # maximum number of atoms in SEA
  num_set = idx+1
  atoms_sea = []
  for i in range(num_set):
    atoms_sea.append(len(dict_sea[i])) 
  atom_sea_max = max(atoms_sea)
  
  
  # Is it linear?
  if np.fabs(Ipmol[0]) <= epsilon:
  # Linear molecule
     linear = True
     L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
     if L_i == True:
       return 'Dinfv',2
     else:
       return 'Cinfv',1
  
  # Non linear molecule
      # spheric symmetry (Ia = Ib = Ic)
      # achatado : oblate (Ia = Ib < Ic)
      # alargado : prolate (Ia < Ib = Ic)
      # asymmetric (Ia != Ib != Ic)
  L_sphera = False
  L_asym = False
  L_oblate = I_diff(Ipmol[0],Ipmol[1],epsilon*10)
  L_prolate = I_diff(Ipmol[1],Ipmol[2],epsilon*10)
  if L_oblate and L_prolate:
    L_sphera = True
    L_asym = False
  if L_sphera:
    L_prolate = False
    L_oblate = False
  if not L_oblate and not L_prolate and not L_sphera: 
    L_asym = True

  #
  # --- SPHERICAL MOLECULES
  #
  if L_sphera: 
    # ---C2 axes
    uvw_asym = [0.,0.,0]
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon)
    c2sp = []
    dxc2_save = []
    number_c2 =0
    for i in range(new_c2):
      c2sp = Cngen(2,dxc2[i])
      L_c2_sp=get_sym(dict_sea,mat_geo,c2sp,num_set,'C2',tolsym) 
      if L_c2_sp:
        number_c2 += 1
        dxc2_save.append(dxc2[i]) 
    # ---diagonal planes
    if number_c2 <= 3:
      num_sigma_d = get_sigma_d(dict_sea,mat_geo,num_set,dxc2_save,number_c2,tolsym)
    # ---inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
    if number_c2 == 3:
      if L_i:
        return 'Th',12
      elif num_sigma_d >= 6: 
        return 'Td',12 
      else:
        return 'T',12 
    if number_c2 == 9:
      if L_i:
        return 'Oh',24
      else:
        return 'O',24 
    if number_c2 > 9:
      if L_i:
        return 'Ih',60
      else:
        return 'I',60 
  
  #
  # --- SYMMETRIC MOLECULES (PROLATE or OBLATE)
  #
  #      Cn axis search
  #
  if L_oblate or L_prolate:
      L_cn = False
      L_cn_tmp = False
      L_cn2_tmp = False
      if L_oblate: 
         main_axis = axis[2]
         main_plane = plane[2]
         uvw = [0.,0.,1.]
      if L_prolate: 
         main_axis = axis[0]
         main_plane = plane[0]
         uvw = [1.,0.,0.]
      for nrot in range(6,1,-1):
        cn_axis = nrot
        cad_axis = 'c'+str(cn_axis)+main_axis
        cnrot =  Cngen(cn_axis,uvw)
        L_cn=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym) 
        if L_cn: 
          max_cn_axis = cn_axis
          break
      if not L_cn:
        L_asym = True
  
  #
  # --- ASYMMETRIC MOLECULES
  #
  #      C2 axis search
  #
  if L_asym:
    max_cn_axis = 1
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      if iaxis == 0: uvw = [1.,0.,0.]
      if iaxis == 1: uvw = [0.,1.,0.]
      if iaxis == 2: uvw = [0.,0.,1.]
      cn_axis = 2
      cad_axis = 'c'+str(cn_axis)+init_axis
      cnrot = []
      cnrot =  Cngen(cn_axis,uvw)
      L_cn_tmp=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym) 
      if L_cn_tmp:
        L_cn = L_cn_tmp
        max_cn_axis = 2
        main_axis = init_axis
        main_plane = init_plane
        cad_plane = 'sigma_'+main_plane
        dxm.append(uvw) 
        uvw_asym = uvw
  
  # --- C2 axis perpendicular to the main axis with redundancies
  # ... Check for vertical planes too
  # --- If not C2 then looks for a Cs plane

  #
  # ---finds C2 axis 
  #
  if L_cn:
    if L_prolate or L_oblate: uvw_asym = uvw 
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon)
    new_c2_p = 0
    if new_c2 > 0:
      for i in range(new_c2):
        cnrot = Cngen(2,dxc2[i])
        L_kk=get_sym(dict_sea,mat_geo,cnrot,num_set,'C2p',tolsym) 
        if L_kk: new_c2_p += 1
    if new_c2_p >= max_cn_axis : L_c2_p = True
  #
  # ---finds vertical planes sigma_v
  #
    coor_sea = []
    mass_sea = []
    num_sigma_v = 0
    for l in range(num_set):
     idx = l
     isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx)
     if isea > 1:
      for i in range(isea):
        for j in range(i+1,isea):
          dv,L_nosig = dist_v(coor_sea[i],coor_sea[j])
          if L_nosig:
            continue
          else: 
            sigma_v = reflex(dv,' ')
            L_sigma_tmp=get_sym(dict_sea,mat_geo,sigma_v,num_set,' ',tolsym) 
            if(L_sigma_tmp): num_sigma_v += 1
    
  
    if num_sigma_v > 0: L_sigma_v = True
  #
  #--- look for horizontal reflection
  #
    sigma_h = np.matrix(reflex(udum,main_plane))
    cad_plane = 'sigma_h'
    L_sigma_h=get_sym(dict_sea,mat_geo,sigma_h,num_set,cad_plane,tolsym) 
  #--- look for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
  # ---look for S2n symmetry operation
    sn_axis = 2*max_cn_axis
    cnrot =  Cngen(sn_axis,uvw)
    snrot = sigma_h * cnrot
    cad_saxis = 's'+str(sn_axis)
    L_s2n=get_sym(dict_sea,mat_geo,snrot,num_set,cad_saxis,tolsym) 
  else:
   #-- looking for a reflection plane
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      cad_cs = 'Cs in axis ',axis[iaxis]
      csplane = reflex(udum,init_plane)
      L_cs=get_sym(dict_sea,mat_geo,csplane,num_set,cad_cs,tolsym)
      if L_cs:
        break
    uvw = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    num_sigma_cs = 0
    num_sigma_cs = get_sigma_d(dict_sea,mat_geo,num_set,uvw,3,tolsym)
    if num_sigma_cs > 0: L_cs = True
   #-- looking for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)

# Any Cn axis (n >=2)?
#    YES
#       Any C2 axis perpendicular to the Cn? 
#         YES  
#            Horizontal plane?
#               YES -> Dnh
#               NO
#                  Improper rotation?
#                     YES -> Dnd
#                     NO  -> Dn
#         NO 
#            Horizontal plane?
#               YES -> Cnh
#               NO  
#                  Vertical planes?
#                     YES -> Cnv
#                     NO
#                        Improper rotation?
#                           YES -> S2n
#                           NO  -> Cn
#    NO  
#      Any plane?
#         YES -> Cs
#            Inversion?
#               YES -> Ci
#               NO  -> C1

  if L_cn:
    if L_c2_p:
      if L_sigma_h:
        return 'D'+str(max_cn_axis)+'h',2*max_cn_axis
      else:
        if L_s2n:
          return 'D'+str(max_cn_axis)+'d',2*max_cn_axis
        else:
          return 'D'+str(max_cn_axis),2*max_cn_axis
    else:
      if L_sigma_h:
        return 'C'+str(max_cn_axis)+'h',max_cn_axis
      else: 
        if L_sigma_v:
          return 'C'+str(max_cn_axis)+'v',max_cn_axis
        else:
          if L_s2n:
            return 'S'+str(sn_axis),max_cn_axis
          else:
            return 'C'+str(max_cn_axis),max_cn_axis
  else:
    if L_cs: 
      return 'Cs',1
    else:
      if L_i: 
        return 'Ci',1
      else:
        return 'C1',1

#for line in readfile("cag.py",hashtag=False,strip=False,ommitblank=False): print line
