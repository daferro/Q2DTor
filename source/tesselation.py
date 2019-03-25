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

 This is a module for the tesselation

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
''' 

import constants as cons
import numpy  as np
from   scipy.optimize import brentq

def cross_point(p1,v1,p2,v2):
    '''
    line 1: p1 with vector v1
    line 2: p2 with vector v2
    '''
    x1, y1 = p1
    x2, y2 = p2
    # Calculate x and y coordinates
    if   v1[0] == 0.0:
         xx = x1
         m2 = v2[1]/v2[0]
         yy = m2*(xx-x2)+y2
    elif v2[0] == 0.0:
         xx = x2
         m1 = v1[1]/v1[0]
         yy = m1*(xx-x1)+y1
    else:
         m1 = v1[1]/v1[0]
         m2 = v2[1]/v2[0]
         xx = ((m1*x1-y1)-(m2*x2-y2))/(m1-m2)
         yy = m1*(xx-x1)+y1
    return np.array( [xx,yy] )

def get_dist2edge(p,pA,pB):
    '''
    returns distance between p and A-B
    '''
    p  = np.array(p)
    pA = np.array(pA)
    pB = np.array(pB)
    # Line 1: A --> B
    p1, v1 = pA, pB-pA
    # Line 2: p --> A-B
    p2, v2 = p , np.array( [v1[1],-v1[0]] )
    # Cross point
    cp  = cross_point(p1,v1,p2,v2)
    # Distance to edge
    dist = np.linalg.norm(cp-p)

    return dist, cp

def get_closest_edge(p0,indices,delaunay):
    '''
    p0: the target point
    indices: the indices of the polygon which contains p0
    '''

    # Assert we have a loop
    if indices[0] != indices[-1]: indices = indices + [indices[0]]

    # Now check each edge and select the closest one
    dist = float("inf")
    for i in range(0,len(indices)-1):
        idxA = indices[i]
        idxB = indices[i+1]
        # Position of each vertex
        pA   = np.array(delaunay.points[idxA])
        pB   = np.array(delaunay.points[idxB])
        # Distance to edge
        distAB, cpAB = get_dist2edge(p0,pA,pB)
        # Update
        if distAB < dist:
           dist  = distAB
           edge  = (idxA,idxB)
           cross = cpAB

    # Return data
    return dist, edge, cross

def get_pointVmid(pA,pB,fPES):

    def edgepot(t,fPES,Vmid,p0,v0):
       return fPES(p0[0]+t*v0[0],p0[1]+t*v0[1]) - Vmid

    VA   = fPES(pA[0],pA[1])
    VB   = fPES(pB[0],pB[1])
    Vmid = (VA+VB)/2.0

    if abs(VA-VB) < 1E-12: return (pA+pB)/2.0, 0.5

    if VB > VA: p0, v0 = pA, pB-pA
    else      : p0, v0 = pB, pA-pB

    t = brentq(edgepot, 0.0, 1.0, args=(fPES,Vmid,p0,v0))
    p = p0 + t*v0
    return p, t

def get_Vmidpoints(delaunay,fPES):
    midpoints = {}
    # For each triangle
    for idxA, idxB, idxC in delaunay.simplices:
        pA = delaunay.points[idxA]
        pB = delaunay.points[idxB]
        pC = delaunay.points[idxC]
        # Point A inside region of interest?
        Ainside1 = pA[0] >= 0.0-cons.ZERO_RAD and pA[0] <= 2.0*np.pi+cons.ZERO_RAD
        Ainside2 = pA[1] >= 0.0-cons.ZERO_RAD and pA[1] <= 2.0*np.pi+cons.ZERO_RAD
        Ainside  = Ainside1 and Ainside2
        # Point B inside region of interest?
        Binside1 = pB[0] >= 0.0-cons.ZERO_RAD and pB[0] <= 2.0*np.pi+cons.ZERO_RAD
        Binside2 = pB[1] >= 0.0-cons.ZERO_RAD and pB[1] <= 2.0*np.pi+cons.ZERO_RAD
        Binside  = Binside1 and Binside2
        # Point C inside region of interest?
        Cinside1 = pC[0] >= 0.0-cons.ZERO_RAD and pC[0] <= 2.0*np.pi+cons.ZERO_RAD
        Cinside2 = pC[1] >= 0.0-cons.ZERO_RAD and pC[1] <= 2.0*np.pi+cons.ZERO_RAD
        Cinside  = Cinside1 and Cinside2
       #if 204 in [idxA,idxB,idxC] and 132 in [idxA,idxB,idxC]:
       #   print pA*cons.R2D, Ainside
       #   print pB*cons.R2D, Binside
       #   print pC*cons.R2D, Cinside
       #   print
        # Any point of triangle in interest region?
        if (not Ainside) and (not Binside) and (not Cinside): continue
        # Get mid points
        if (idxA,idxB) not in midpoints.keys():
           midAB,tAB = get_pointVmid(pA,pB,fPES)
           midpoints[(idxA,idxB)] = midAB,tAB
           midpoints[(idxB,idxA)] = midAB,tAB
        if (idxA,idxC) not in midpoints.keys():
           midAC,tAC = get_pointVmid(pA,pC,fPES)
           midpoints[(idxA,idxC)] = midAC,tAC
           midpoints[(idxC,idxA)] = midAC,tAC
        if (idxB,idxC) not in midpoints.keys():
           midBC,tBC = get_pointVmid(pB,pC,fPES)
           midpoints[(idxB,idxC)] = midBC,tBC
           midpoints[(idxC,idxB)] = midBC,tBC
    return midpoints

def check_edge(pA,pB,fPES):
    '''
    Checks if points between pA and pB
    present a energy below or above VA and VB
    PS:
      VA is the energy at pA
      VB is the energy at pB
    '''
    # Coordinates of each point and energy
    pA, VA  = np.array(pA), fPES.value(pA[0],pA[1])
    pB, VB  = np.array(pB), fPES.value(pB[0],pB[1])
    Vmin    = min(VA,VB)
    Vmax    = max(VA,VB)
    # Vector from A to B
    vAB = pB-pA
    # Check points between A and B
    for t in range(1,11,1):
        p = pA + (t/10.0)*vAB
        V = fPES.value(p[0],p[1])
        if V > Vmax or V < Vmin: return False
    return True


#----------------------#
# Method 1 : distances #
#----------------------#
def method1_idx(p0,triangle,delaunay):
    dist = float('inf')
    for idxA in triangle:
        # Get point
        pA = np.array(delaunay.points[idxA])
        # Get distance
        dA = np.linalg.norm( pA - p0 )
        # Compare
        if dA < dist: dist, idx = dA, idxA
    return idx

    idxA, idxB, idxC = triangle
    # Get points
    pA = np.array(delaunay.points[idxA])
    pB = np.array(delaunay.points[idxB])
    pC = np.array(delaunay.points[idxC])
    # Get distances
    dB  = np.linalg.norm( pB - p0 )
    dC  = np.linalg.norm( pC - p0 )
    # Select vertex
    dist,idx = min( [(dA,idxA),(dB,idxB),(dC,idxC)] )
    return idx

#----------------------#
# Method 2 : potential #
#----------------------#
def method2_idx(p0,triangle,delaunay,fPES):
    V0 = fPES.value(p0[0],p0[1])
    dV = float('inf')
    for idxA in triangle:
        # Get point
        pA = np.array(delaunay.points[idxA])
        # Get energy
        dVA = abs(fPES.value(pA[0],pA[1]) - V0)
        # Compare
        if dVA < dV: dV, idx = dVA, idxA
    return idx

#-----------------------------------#
# Method 3 : potential and distance #
#-----------------------------------#
def method3_idx(p0,triangle,delaunay,fPES):
    idxA, idxB, idxC = triangle
    # Get points
    pA = np.array(delaunay.points[idxA])
    pB = np.array(delaunay.points[idxB])
    pC = np.array(delaunay.points[idxC])
    # Get energies
    V0  = fPES.value(p0[0],p0[1])
    VA  = fPES.value(pA[0],pA[1])
    VB  = fPES.value(pB[0],pB[1])
    VC  = fPES.value(pC[0],pC[1])
    # Vertex with intermediate energy
    VI,idxI = sorted([(VA,idxA),(VB,idxB),(VC,idxC)])[1]
    # Select vertex
    dV, idx = min( [(abs(VA-V0),idxA),(abs(VB-V0),idxB),(abs(VC-V0),idxC)] )
    if idx == idxI:
       dA  = np.linalg.norm( pA - p0 )
       dB  = np.linalg.norm( pB - p0 )
       dC  = np.linalg.norm( pC - p0 )
       dist,idx = min( [(dA,idxA),(dB,idxB),(dC,idxC)] )
    return idx

#--------------------#
# Method 4 : heights #
#--------------------#
def method4_idx(p0,indices,delaunay,fPES,Vmidpoints):

    # Get closest edge
    dist, edge, cross = get_closest_edge(p0,indices,delaunay)

    # V(idx1) < V(idx2)
    idx1, idx2 = edge
    p1 = np.array(delaunay.points[idx1])
    p2 = np.array(delaunay.points[idx2])
    V1 = fPES(p1[0],p1[1])
    V2 = fPES(p2[0],p2[1])
    if V1 > V2:
       idx1, p1, V1, idx2, p2, V2 = idx2, p2, V2, idx1, p1, V1

    # Get middle point (in V) for the edge
    if (idx1,idx2) in Vmidpoints.keys():
       pmid, tmid = Vmidpoints[(idx1,idx2)]
    else:
       pmid, tmid = get_pointVmid(p1,p2,fPES)

    # Value of t for the cross point
    v12 = p2 - p1
    if abs(v12[0]) > abs(v12[1]): t_cross = (cross[0]-p1[0])/v12[0]
    else                        : t_cross = (cross[1]-p1[1])/v12[1]

    if t_cross <= tmid: return idx1
    if t_cross >  tmid: return idx2

#-----------------------------------#
# Method 5 : remove edge + method 4 #
#-----------------------------------#
def method5_polygon(p0,simplex,delaunay,fPES):

    # Indices of the triangle
    idxA, idxB, idxC = list(delaunay.simplices[simplex])

    # Initialize 'indices' list
    indices = [idxA]

    # Get points
    pA = np.array(delaunay.points[idxA])
    pB = np.array(delaunay.points[idxB])
    pC = np.array(delaunay.points[idxC])

    # check edge A-B
    eAB  = check_edge(pA,pB,fPES)
    edge = set([idxA,idxB])
    if eAB is False:
       for neighbor in delaunay.neighbors[simplex]:
           node = set(delaunay.simplices[neighbor]).difference(edge)
           if len(node) == 1: break
       idxAB = node.pop()
       indices += [idxAB,idxB]
    else:
       indices += [idxB]

    # check edge B-C
    eBC  = check_edge(pB,pC,fPES)
    edge = set([idxB,idxC])
    if eBC is False:
       for neighbor in delaunay.neighbors[simplex]:
           node = set(delaunay.simplices[neighbor]).difference(edge)
           if len(node) == 1: break
       idxBC = node.pop()
       indices += [idxBC,idxC]
    else:
       indices += [idxC]

    # check edge A-C
    eAC  = check_edge(pA,pC,fPES)
    edge = set([idxA,idxC])
    if eAC is False:
       for neighbor in delaunay.neighbors[simplex]:
           node = set(delaunay.simplices[neighbor]).difference(edge)
           if len(node) == 1: break
       idxAC = node.pop()
       indices += [idxAC,idxA]
    else:
       indices += [idxA]

    # Return new polygon
    return list(indices)

def method5_idx(p0,simplex,delaunay,fPES):
    # Get polygon
    polygon = method5_polygon(p0,simplex,delaunay,fPES)
    # Get idx with method 2
    idx = method2_idx(p0,polygon,delaunay,fPES)
    # Get value of function
    return idx

#-----------------------------------#
# Method 6 : method 4 + tanh interp #
#-----------------------------------#
def method6_value(p0,triangle,delaunay,fPES,Vmidpoints,data4interp):
    # Get edge
    dist, edge, cross = get_closest_edge(p0,triangle,delaunay)
    # Define 1 and 2
    idx1, idx2 = edge
    p1 = np.array(delaunay.points[idx1])
    p2 = np.array(delaunay.points[idx2])
    # Get point with V=(VA+VB)/2
    if (idx1,idx2) in Vmidpoints.keys():
       pmid, dummy = Vmidpoints[(idx1,idx2)]
    else:
       pmid, dummy = get_pointVmid(p1,p2,fPES)
    # Get value of t for the middle point (tm) and for the cross point (t)
    v12 = p2-p1
    if abs(v12[0]) > abs(v12[1]):
       tm = ( pmid[0]-p1[0])/v12[0]
       t  = (cross[0]-p1[0])/v12[0]
    else:
       tm = ( pmid[1]-p1[1])/v12[1]
       t  = (cross[1]-p1[1])/v12[1]
    if t> tm: idx = idx2
    else    : idx = idx1

    #--------------------#
    # tanh interpolation #
    #--------------------#
    Q1 = data4interp[idx1]
    Q2 = data4interp[idx2]
    # Assert we go in increasing direction
    if Q1 > Q2:
       idx1, Q1, idx2, Q2 = idx2, Q2, idx1, Q1
       tm = 1.0-tm
       t  = 1.0-t
    # Get gamma value
    relerror = 0.5 # in percentage
    eps1 = abs(relerror*Q1/100.0)
    eps2 = abs(relerror*Q2/100.0)
    if eps1 == 0.0: eps1 = eps2
    if eps2 == 0.0: eps2 = eps1
    DELTA = (Q2-Q1)/2.0
    if DELTA <= eps1: gamma1 = 0.0
    else            : gamma1 = np.arctanh(eps1/DELTA-1.0)/(-tm)
    if DELTA <= eps2: gamma2 = 0.0
    else            : gamma2 = np.arctanh(1.0-eps2/DELTA)/(1.0-tm)
    gamma  = max(gamma1,gamma2)
    #gamma = 100.0 (with high value, we have method 4)
    # Value of the partition function
    Q = ((Q1+Q2)+(Q2-Q1)*np.tanh(gamma*(t-tm)))/2.0
    return Q, idx








#-------------------#
# Values for Q2DTor #
#-------------------#
def q2dtor_interpolation(p0,delaunay,spdata,method=6,args=[None,None]):
    '''
    '''

    fPES,Vmidpoints = args

    # Delaunay triangle
    p0 = np.array(p0)
    simplex  = delaunay.find_simplex( (p0[0],p0[1]) )
    triangle = list(delaunay.simplices[simplex])
    idxA, idxB, idxC = triangle

    # Really close to a vertex
    pA = np.array(delaunay.points[idxA]); distA = np.linalg.norm(pA-p0)
    pB = np.array(delaunay.points[idxB]); distB = np.linalg.norm(pA-p0)
    pC = np.array(delaunay.points[idxC]); distC = np.linalg.norm(pA-p0)
    mdist, midx = min([(distA,idxA),(distB,idxB),(distC,idxC)])
    if mdist < 2.0 * cons.D2R: idx = midx
    else                     : idx = None

    # Methods 1-5
    if method in [1,2,3,4,5]:
       if   method == 1 and idx is None: idx = method1_idx(p0,triangle,delaunay)
       elif method == 2 and idx is None: idx = method2_idx(p0,triangle,delaunay,fPES)
       elif method == 3 and idx is None: idx = method3_idx(p0,triangle,delaunay,fPES)
       elif method == 4 and idx is None: idx = method4_idx(p0,triangle,delaunay,fPES,Vmidpoints)
       elif method == 5 and idx is None: idx = method5_idx(p0,simplex ,delaunay,fPES)
       return list(spdata[:,idx]), idx

    if method == 6:
       data2return = []
       indices = []
       nrows, ncols = spdata.shape
       for row in range(nrows):
           data = spdata[row,:]
           value, idx = method6_value(p0,triangle,delaunay,fPES,Vmidpoints,data)
           data2return.append(value)
           indices.append(idx)
       return data2return, indices[0]

#---------------------#
# Selection of method #
#---------------------#

def which_statpoint(p0,delaunay,method=1,args=None):

    fPES,Vmidpoints = args

    # Delaunay triangle
    p0 = np.array(p0)
    simplex  = delaunay.find_simplex( (p0[0],p0[1]) )
    triangle = list(delaunay.simplices[simplex])
    idxA, idxB, idxC = triangle

    # Really close to a vertex
    pA = np.array(delaunay.points[idxA]); distA = np.linalg.norm(pA-p0)
    pB = np.array(delaunay.points[idxB]); distB = np.linalg.norm(pA-p0)
    pC = np.array(delaunay.points[idxC]); distC = np.linalg.norm(pA-p0)
    mdist, midx = min([(distA,idxA),(distB,idxB),(distC,idxC)])
    if mdist < 2.0 * cons.D2R: idx = midx
    else                     : idx = None

    # Now, select idx
    if   method == 1 and idx is None: idx = method1_idx(p0,triangle,delaunay)
    elif method == 2 and idx is None: idx = method2_idx(p0,triangle,delaunay,fPES)
    elif method == 3 and idx is None: idx = method3_idx(p0,triangle,delaunay,fPES)
    elif method == 4 and idx is None: idx = method4_idx(p0,triangle,delaunay,fPES,Vmidpoints)
    elif method == 5 and idx is None: idx = method5_idx(p0,simplex ,delaunay,fPES)

    return idx



#------------------------------------#
# Replication of points for Delaunay #
#------------------------------------#

def replicate_points(molname,dict_points):
    '''
    dict_points[point_name] = (phi1, phi2, etc...)
    point name is obtained using molname, phi1, and phi2
    '''
    def get_ppname(phi1,phi2,molname="",units="rad"):
        '''
        This function generates a name for each point in the 2D-torsional PES
        * units: indicates the units of phi1 and phi2 ('rad' or 'deg')
        PS: points that are less than a degree apart cannot be differentiated by the generated name
        '''
        if units == "rad": int_phi1, int_phi2 = int(round(phi1*cons.R2D)), int(round(phi2*cons.R2D))
        if units == "deg": int_phi1, int_phi2 = int(round(phi1    )), int(round(phi2    ))
        return molname+"_%03i_%03i"%(int_phi1,int_phi2)

    list_points = dict_points.keys()

    for point_name in list_points:
        phi1 = dict_points[point_name][0]
        phi2 = dict_points[point_name][1]

        phi1_a = phi1 - 2.0*np.pi
        phi2_a = phi2 - 2.0*np.pi
        name_a = get_ppname(phi1_a,phi2_a,molname,"rad")
        dict_points[name_a] = [phi1_a,phi2_a]+dict_points[point_name][2:]

        phi1_b = phi1 - 2.0*np.pi
        phi2_b = phi2
        name_b = get_ppname(phi1_b,phi2_b,molname,"rad")
        dict_points[name_b] = [phi1_b,phi2_b]+dict_points[point_name][2:]

        phi1_c = phi1 - 2.0*np.pi
        phi2_c = phi2 + 2.0*np.pi
        name_c = get_ppname(phi1_c,phi2_c,molname,"rad")
        dict_points[name_c] = [phi1_c,phi2_c]+dict_points[point_name][2:]

        phi1_d = phi1
        phi2_d = phi2 - 2.0*np.pi
        name_d = get_ppname(phi1_d,phi2_d,molname,"rad")
        dict_points[name_d] = [phi1_d,phi2_d]+dict_points[point_name][2:]

        phi1_e = phi1
        phi2_e = phi2 + 2.0*np.pi
        name_e = get_ppname(phi1_e,phi2_e,molname,"rad")
        dict_points[name_e] = [phi1_e,phi2_e]+dict_points[point_name][2:]

        phi1_f = phi1 + 2.0*np.pi
        phi2_f = phi2 - 2.0*np.pi
        name_f = get_ppname(phi1_f,phi2_f,molname,"rad")
        dict_points[name_f] = [phi1_f,phi2_f]+dict_points[point_name][2:]

        phi1_g = phi1 + 2.0*np.pi
        phi2_g = phi2
        name_g = get_ppname(phi1_g,phi2_g,molname,"rad")
        dict_points[name_g] = [phi1_g,phi2_g]+dict_points[point_name][2:]

        phi1_h = phi1 + 2.0*np.pi
        phi2_h = phi2 + 2.0*np.pi
        name_h = get_ppname(phi1_h,phi2_h,molname,"rad")
        dict_points[name_h] = [phi1_h,phi2_h]+dict_points[point_name][2:]

    return dict_points
#-------------------------------------------------------#


