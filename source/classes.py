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

 This module constains different Classes

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
'''

#------------------------------------------------#
# >>           Importation  section           << #
#------------------------------------------------#
import math, cmath
import sys
import time
#------------------------------------------------#
import numpy          as     np
import random
from   scipy.optimize import curve_fit
from   scipy.optimize import newton
#random.seed(1111)
#------------------------------------------------#
import constants      as     cons
import helpfns        as     hf
from   gtsfile        import read_gtsfile
# >>>>>>>>>>>>>>>>>>>>> ## <<<<<<<<<<<<<<<<<<<<< #




#>>>>>>>>>>>>>>>>>>*
# CLASS: InfoStr   *
#>>>>>>>>>>>>>>>>>>*
class InfoStr:
    def __init__(self):
        self.text = ""
    def add(self, string, n=0, b=0):
        self.text = self.text + b*" " + string + n*"\n"
    def __str__(self):
        return self.text
    def blanklines(self, n=1):
        self.text = self.text + n*"\n"
    def endinblank(self):
        if not self.text.endswith("\n"): self.text = self.text + "\n"
    def clear(self):
        self.text = ""
    def removeblank(self):
        if self.text.endswith("\n\n"): self.text = self.text[:-1]


#>>>>>>>>>>>>>>>>>>*
# CLASS: Logger    *
#>>>>>>>>>>>>>>>>>>*
class Logger(object):
    '''
    Class used to save in a file
    what is printed in terminal
    Requirements:
      sys library
    Use:
      sys.stdout = Logger(f_out,tprint)
    '''
    def __init__(self,output,tprint,mode="a"):
        self.terminal = sys.stdout
        self.log = open(output, mode)
        self.tprint = tprint

    def write(self, message):
        if self.tprint: self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    


#>>>>>>>>>>>>>>>>>>*
# CLASS: Fourier2D *
#>>>>>>>>>>>>>>>>>>*
class Fourier2D(object):
      def __init__(self, terms, coefs=None, imag=True):
          # All terms or only not imaginary?
          if imag:
             self._terms = list(terms)
          else:
             self._terms = []
             for term, idx1, idx2 in terms:
                 if term == "sin"   : continue
                 if term == "sincos": continue
                 if term == "cossin": continue
                 self._terms.append( (term, idx1, idx2) )
          self.NP      = len(self._terms)
          self._coefs  = coefs


      def set_coefs(self,coefs):
          self._coefs = coefs

      def __call__(self, Phi1, Phi2):
         return self.value(Phi1,Phi2)

      def value(self,Phi1,Phi2):
          '''
          Returns the value of the functions
          at the given point
          '''

          value = 0.0
          for term,coef in zip(self._terms,self._coefs):
              term_type, idx1, idx2 = term
              if term_type == "const": value += coef
              elif term_type == "cos" and coef != 0.0:
                   if idx1 == "-":     value += coef * np.cos(idx2*Phi2)
                   if idx2 == "-":     value += coef * np.cos(idx1*Phi1)
              elif term_type == "sin" and coef != 0.0:
                   if idx1 == "-":     value += coef * np.sin(idx2*Phi2)
                   if idx2 == "-":     value += coef * np.sin(idx1*Phi1)
              elif term_type == "coscos" and coef != 0.0:
                   value += coef * np.cos(idx1*Phi1)*np.cos(idx2*Phi2)
              elif term_type == "sincos" and coef != 0.0:
                   value += coef * np.sin(idx1*Phi1)*np.cos(idx2*Phi2)
              elif term_type == "sinsin" and coef != 0.0:
                   value += coef * np.sin(idx1*Phi1)*np.sin(idx2*Phi2)
              elif term_type == "cossin" and coef != 0.0:
                   value += coef * np.cos(idx1*Phi1)*np.sin(idx2*Phi2)
          return value

      def derphi1(self,Phi1,Phi2):
          '''
          '''

          value = 0.0
          for term,coef in zip(self._terms,self._coefs):
              term_type, idx1, idx2 = term
              if   term_type == "cos" and coef != 0.0:
                   if idx2 == "-":     value += -idx1 * coef * np.sin(idx1*Phi1)
              elif term_type == "sin" and coef != 0.0:
                   if idx2 == "-":     value +=  idx1 * coef * np.cos(idx1*Phi1)
              elif term_type == "coscos" and coef != 0.0:
                   value += -idx1 * coef * np.sin(idx1*Phi1)*np.cos(idx2*Phi2)
              elif term_type == "sincos" and coef != 0.0:
                   value +=  idx1 * coef * np.cos(idx1*Phi1)*np.cos(idx2*Phi2)
              elif term_type == "sinsin" and coef != 0.0:
                   value +=  idx1 * coef * np.cos(idx1*Phi1)*np.sin(idx2*Phi2)
              elif term_type == "cossin" and coef != 0.0:
                   value += -idx1 * coef * np.sin(idx1*Phi1)*np.sin(idx2*Phi2)
          return value

      def derphi2(self,Phi1,Phi2):
          '''
          '''
          value = 0.0
          for term,coef in zip(self._terms,self._coefs):
              term_type, idx1, idx2 = term
              if   term_type == "const": continue
              elif term_type == "cos" and coef != 0.0:
                   if idx1 == "-":     value += - idx2 * coef * np.sin(idx2*Phi2)
              elif term_type == "sin" and coef != 0.0:
                   if idx1 == "-":     value +=   idx2 * coef * np.cos(idx2*Phi2)
              elif term_type == "coscos" and coef != 0.0:
                   value += -idx2 * coef * np.cos(idx1*Phi1)*np.sin(idx2*Phi2)
              elif term_type == "sincos" and coef != 0.0:
                   value += -idx2 * coef * np.sin(idx1*Phi1)*np.sin(idx2*Phi2)
              elif term_type == "sinsin" and coef != 0.0:
                   value +=  idx2 * coef * np.sin(idx1*Phi1)*np.cos(idx2*Phi2)
              elif term_type == "cossin" and coef != 0.0:
                   value +=  idx2 * coef * np.cos(idx1*Phi1)*np.cos(idx2*Phi2)
          return value

      def value2(self,XY_array,*coefs):
          '''
          Just created to be used in 'fit' method
          '''
          # Set new coefs
          self._coefs = coefs
          # Get output
          output = []
          for (Phi1,Phi2) in XY_array:
             fcn = self.value(Phi1,Phi2)
             output.append(fcn)
          return np.array(output)

      def fit(self,xdata,ydata,weight=0.0,guess=None):
          '''
          A function created in order to fit
          '''
          t0 = time.time()

          # Min, average and max values in ydata
          y_min  = min(ydata)
          y_mean = sum(ydata) / len(ydata)
          y_max  = max(ydata)

          # Define weights
          if weight in [0.0,0,None,False]: weights = None
          else                           : weights = [ (y-y_min+0.01)**float(weight) for y in ydata ]

          # Define guess values
          if guess in [None,False,[]]: guess = [y_mean] + [0.] * (self.NP-1)

          # Define bounds
          low_bound = [y_min]+[-np.inf]*(self.NP-1)
          upp_bound = [y_max]+[+np.inf]*(self.NP-1)
          bounds    = (low_bound,upp_bound)

          # Perform fitting
          try:
             popt, pcov = curve_fit(self.value2, xdata, ydata, p0=guess, sigma=weights, bounds=bounds)
          except:
             popt, pcov = curve_fit(self.value2, xdata, ydata, p0=guess, sigma=weights)

          # Set coefficients
          self.set_coefs(popt)

          # Check
          ydata_fit = self.value2(xdata,*popt)
          averAbsErr = 0.0
          averAbsErr_small = 0.0
          n_small = 0.0
          SStot = 0.0
          SSres = 0.0
          for y,yf in zip(ydata,ydata_fit):
              SSres += (y-yf)**2
              SStot += (y-y_mean)**2
              averAbsErr += abs(y-yf)
              if y < y_mean:
                 n_small += 1
                 averAbsErr_small += abs(y-yf)
          rsquare = 1.0- SSres/SStot
          averAbsErr /= len(ydata)
          averAbsErr_small /= n_small

          t1 = time.time()
          fitting_time = t1 - t0
          return self._terms, popt, (rsquare,averAbsErr,averAbsErr_small), fitting_time


#>>>>>>>>>>>>>>>>>>*
# CLASS: splineVaG *
#>>>>>>>>>>>>>>>>>>*
class SplineVaG:

      def __init__(self,xx,yy):
          self.ZERO = 1e-6
          self.xx = xx
          self.yy = yy
          self.s_alpha = xx[0]
          self.V_alpha = yy[0]
          self.s_omega = xx[-1]
          self.V_omega = yy[-1]
          for idx in range(len(xx)):
              if xx[idx] == 0.0: self.V_saddle = yy[idx]

          # a) Obtain maxima inside xx
          full_maxdata = hf.obtain_extremum(xx,yy,xtr="max",full=True)

          # Data of "global" maximum
          self.smax = full_maxdata[0][0]
          self.VAG  = full_maxdata[0][1]

          # b) Remove coincident ranges in maximum splines!
          self.splines = []
          ranges       = []
          for xxx,xxx,the_spline,si,sj in full_maxdata:
              for sa_prime,sb_prime in ranges:
                  if abs(si - sb_prime) < self.ZERO: continue
                  if abs(sj - sa_prime) < self.ZERO: continue
                  if si + self.ZERO > sa_prime and sb_prime + self.ZERO > si: si = sb_prime
                  if sj + self.ZERO > sa_prime and sb_prime + self.ZERO > sj: sj = sa_prime
                  if abs(sj-si) < self.ZERO: break
              if abs(sj-si) < self.ZERO: continue
              if sj < si: continue
              self.splines += [ (si,sj,the_spline)  ]
              ranges.append( (si,sj)  )
          self.splines.sort()
          ranges.sort()

          # c) Add new splines (if no coincidence in ranges)
          for idx in range(0,len(xx)-1):
              si = xx[idx]
              sj = xx[idx+1]
              for sa_prime,sb_prime in ranges:
                  if abs(si - sb_prime) < self.ZERO: continue
                  if abs(sj - sa_prime) < self.ZERO: continue
                  if si + self.ZERO > sa_prime and sb_prime + self.ZERO > si: si = sb_prime
                  if sj + self.ZERO > sa_prime and sb_prime + self.ZERO > sj: sj = sa_prime
                  if abs(sj-si) < self.ZERO: break
              if abs(sj-si) < self.ZERO: continue
              if sj < si: continue
              the_spline = hf.localspline(xx,yy,idx,nps=(2,2),spl=3)
              self.splines.append( (si, sj, the_spline) )
              ranges.append( (si,sj)  )
          self.splines.sort()
          ranges.sort()

      def __call__(self, svalue):
         return self.interpolate(svalue)

      def get_bwvalues(self):
          return self.s_alpha, self.V_alpha

      def get_fwvalues(self):
          return self.s_omega, self.V_omega

      def get_saddle(self): return self.V_saddle

      def get_max(self):
          return self.smax, self.VAG

      def interpolate(self,si,E=0.0):
          if si <= self.s_alpha: return self.V_alpha - E
          if si >= self.s_omega: return self.V_omega - E

          for idx in range(len(self.xx)):
              if abs(self.xx[idx]-si) < self.ZERO: return self.yy[idx] - E

          for s1, s2, spline in self.splines:
              if s1 < si and si < s2: return spline(si) - E

      def reduceMEP(self):
          self.s_alpha = self.xx[ 1]
          self.V_alpha = self.yy[ 1]
          self.s_omega = self.xx[-2]
          self.V_omega = self.yy[-2]

      def undoreduceMEP(self):
          self.s_alpha = self.xx[ 0]
          self.V_alpha = self.yy[ 0]
          self.s_omega = self.xx[-1]
          self.V_omega = self.yy[-1]

      def returnpoints(self,E,ds):
          s_current = self.s_alpha
          rtn_points = []

          while s_current <= self.s_omega:
                V1 = self.interpolate(s_current)
                V2 = self.interpolate( min(s_current+ds,self.s_omega) )

                if   V1 > V2: rpoint_type = "-"
                elif V1 < V2: rpoint_type = "+"
                else: s_current = s_current+ds; continue

                if (V1 <= E and E < V2) or (V1 >= E and E > V2):
                   s_guess = s_current + 0.5*ds
                   rpoint = newton(self.interpolate,s_guess,args=(E,))
                   rtn_points.append( (rpoint,rpoint_type) )

                s_current = s_current+ds

          if E == V2: rtn_points.append( (self.s_omega,rpoint_type) )

          return rtn_points


#>>>>>>>>>>>>>>>>>>*
# CLASS: Freq      *
#>>>>>>>>>>>>>>>>>>*
class Freq:
    def __init__(self , scaling=1.0, mu=1.0/cons.amu, rmode=1):
       '''
       if rmode == 1: for imag frequencies, pfn = 1 and its derivates are 0
       if rmode == 2: for imag frequencies, pfn = 1E10 and its derivates are 0
       '''
       self.scaling = scaling 
       self.mu      = mu
       # Basic ones
       self.evalue  = None
       self.angfreq = None
       self.wavenum = None
       self.evector = None
       # Derivated ones
       self.vibT    = None
       self.zpe     = None
       self.tpoint  = None
       self.real    = None
       # rmode
       self._rmode   = rmode

    def __str__(self):
       return self.str(f="%6.1f")

    def __float__(self):
       if self.real: return  self.wavenum.real
       else:         return -self.wavenum.imag

    #####################
    # Setting variables #
    #####################
    def set_scaling(self,scaling):
        self.scaling = scaling

    def set_mu(self,mu):
        self.mu = mu

    def set_evalue(self , evalue):
        self.evalue = evalue
        # a) Calculate angular frequency
        self.angfreq  = cmath.sqrt(self.evalue/self.mu)
        self.angfreq *= self.scaling
        # If real, keep only real part
        if self.angfreq.imag == 0.0:
           self.angfreq = self.angfreq.real
           self.real = True
        else:
           self.real = False
        # b) Calculate wavenum
        self.wavenum  = self.angfreq / (2.0*math.pi*cons.c0)

    def set_wavenum(self , wavenum):
        if wavenum.imag == 0.0:
           if wavenum.real < 0.0:
              self.wavenum = complex(0,abs(wavenum.real))
              self.real = False
           else:
              self.wavenum = wavenum.real
              self.real = True
        else:
           if wavenum.real == 0.0:
              self.wavenum = wavenum
              self.real = False
           else:
              sys.exit("This frequency has both real and imaginary terms...")
        # a) Calculate angular frequency
        self.angfreq = self.wavenum * 2.0 * math.pi * cons.c0
        # b) Calculate eigenvalue
        self.evalue = (self.angfreq / self.scaling)**2 * self.mu
        if self.evalue.imag == 0.0: self.evalue = self.evalue.real

    def set_evector(self , evector):
        self.evector = evector

    def str(self,f="%4i"):
        real = self.wavenum.real / cons.cm
        imag = self.wavenum.imag / cons.cm
        if   imag == 0.0: string = f%real+" "
        elif real == 0.0: string = f%imag+"i"
        else:             string = "%4i%+4ii"%(real,imag)
        return string

    ###############################
    # Calculate rest of variables #
    ###############################
    def calc_derivated_magnitudes(self):
        # (a) The frequency is imaginary
        if not self.real:
           self.zpe    = 0.0
           self.vibT   = None
           self.tpoint = 1e10
        # (b) The frequency is real
        else:
           self.zpe    = cons.hbar * self.angfreq / 2.0
           self.vibT   = cons.hbar * self.angfreq / cons.kB
           self.tpoint = math.sqrt( cons.hbar / self.angfreq / self.mu )

    ###############
    # Return data #
    ###############
    def copy(self):
        copy_freq = Freq(scaling = self.scaling, mu = self.mu, rmode= self._rmode)
        copy_freq.set_evalue(self.evalue)
        copy_freq.set_evector(self.evector)
        if self.zpe is not None: copy_freq.calc_derivated_magnitudes()
        return copy_freq
  
    def isItImag(self):    return not self.real

    def get(self,variable):
        if variable == "wavenum": return self.wavenum
        if variable == "angfreq": return self.angfreq
        if variable == "mu"     : return self.mu
        if variable == "evalue" : return self.evalue
        if variable == "evector": return self.evector
        if variable == "zpe"    : return self.zpe
        if variable == "vibT"   : return self.vibT
        if variable == "tpoint" : return self.tpoint

    def get_qvib(self,T):
        '''
        Calculates vibrational partition function,
        taking as zero the energy of level n=0
        (i.e. the zero point energy)
        '''

        if (not self.real) and (self._rmode == 1): return 1.0 , self.zpe
        if (not self.real) and (self._rmode == 2): return 1E10, self.zpe

        exp  = np.exp(-self.vibT/T)
        qvib = 1.0/(1.0-exp)

        return qvib, self.zpe

    def get_fdln(self,T):
        '''
        '''
        if not self.real: return 0.0
        hw   = cons.hbar * self.angfreq
        bhw  = (1.0/cons.kB/T) * hw
        exp  = np.exp(-bhw)
        fdln = - hw * exp / (1.0-exp)
        return fdln

    def get_sdln(self,T):
        '''
        '''
        if not self.real: return 0.0
        hw = cons.hbar * self.angfreq
        bhw = (1.0/cons.kB/T) * hw
        # Exponential part
        exp_part = 1.0 / ((np.exp(bhw)-1.0)*(1.0-np.exp(-bhw)))
        # second derivative of natural log
        sdln = hw*hw * exp_part
        return sdln


#>>>>>>>>>>>>>>>>>>*
# CLASS: Queue     *
#>>>>>>>>>>>>>>>>>>*
class Queue:
    """
    A simple implementation of a FIFO queue.
    """
    def __init__(self):
        self._items = []
    def __len__(self):
        return len(self._items)
    def __iter__(self):
        for item in self._items:
            yield item
    def __str__(self):
        return str(self._items)
    def enqueue(self, item):
        self._items.append(item)
    def dequeue(self):
        return self._items.pop(0)
    def clear(self):
        self._items = []


#>>>>>>>>>>>>>>>>>>*
# CLASS: Stack     *
#>>>>>>>>>>>>>>>>>>*
class Stack:
    """
    A simple implementation of a LIFO stack
    """
    def __init__(self):
        self._items = []
    def __len__(self):
        return len(self._items)
    def __iter__(self):
        for item in self._items:
            yield item
    def __str__(self):
        return str(self._items)
    def push(self, item):
        self._items = [item] + self._items
    def pop(self):
        return self._items.pop(0)
    def clear(self):
        self._items = []


#>>>>>>>>>>>>>>>>>>*
# CLASS: UGRAPH    *
#>>>>>>>>>>>>>>>>>>*
class UGRAPH:
      """
      A simple implementation of a undirected graph
      """

      def __init__(self):
         self._ugdict = {}

      def __str__(self):
          num_nodes = self.get_nnodes()
          num_edges = self.get_nedges()
          return "(n,e)=(%i,%i)"%(num_nodes,num_edges)

      def get_nnodes(self):
          '''
          Returns number of nodes in the ugraph
          '''
          return len(self._ugdict.keys())

      def get_edges(self):
          '''
          Returns the edges in the ugraph
          '''
          edges = set([])
          for node1 in self._ugdict.keys():
              for node2 in self._ugdict[node1]:
                  edge = tuple(sorted((node1,node2)))
                  edges.add(edge)
          return edges

      def get_nedges(self):
          '''
          Returns number of edges in the ugraph
          '''
          edges = set([])
          for node1 in self._ugdict.keys():
              for node2 in self._ugdict[node1]:
                  edge = tuple(sorted((node1,node2)))
                  edges.add(edge)
          return len(edges)

      def neighbors(self,node):
          return self._ugdict[node].copy()

      #----------------------------#
      # Add/Remove nodes and edges #
      #----------------------------#
      def add_node(self,node):
          if node not in self._ugdict.keys():
             self._ugdict[node] = set([])

      def add_edge(self,node1,node2):
          self.add_node(node1)
          self.add_node(node2)
          self._ugdict[node1].add(node2)
          self._ugdict[node2].add(node1)

      def remove_node(self,node1):
          # Remove node
          self._ugdict.pop(node1)
          # Remove edges with that node
          for node2 in self._ugdict.keys():
              self._ugdict[node2].discard(node1)
          
      def remove_edge(self,node1,node2):
          self._ugdict[node1].discard(node2)
          self._ugdict[node2].discard(node1)
      #----------------------------#

      def bfsearch(self,start_idx):
          '''
          Breadth First Search for undirected graph
          Input:
            * graph_dict: a dict of the graph representing the
                          adjacency list
                          - key  : integer
                          - value: list of integers
            * start_idx : the index where to start the BFS
          '''
          # Initialize queue
          queue   = Queue()
          visited = [start_idx]
          queue.enqueue(start_idx)
          # Start BFS
          while len(queue) != 0:
               # Take node out of queue
               target_idx = queue.dequeue()
               # Get neighbors
               neighbors  = self._ugdict[target_idx]
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   queue.enqueue(neighbor)
          return visited

      def dfsearch(self,start_idx):
          '''
          Depth First Search
          Breadth First Search for undirected graph
          Input:
            * graph_dict: a dict of the graph representing the
                          adjacency list
                          - key  : integer
                          - value: list of integers
            * start_idx : the index where to start the BFS
          '''
          # Initialize queue
          stack   = Stack()
          visited = [start_idx]
          stack.push(start_idx)
          # Start BFS
          while len(stack) != 0:
               # Take node out of queue
               target_idx = stack.pop()
               # Get neighbors
               neighbors  = self._ugdict[target_idx]
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   stack.push(neighbor)
          return visited

      def bfsearch1d(self,idx1,idx2):
          '''
          Using a BFS algorithm, goes through
          the graph.
          However, it does it in the idx1-->idx2
          directions.
          '''
          # Initialize queue
          queue      = Queue()
          neighbors1 = self._ugdict[idx1]
          old2       = None

          # idx1 and idx2 are not bonded, there is a node in the middle (idx1--idxJ--idx2)
          if idx2 not in neighbors1:
             neighbors2 = self._ugdict[idx2]
             idxJ = list(neighbors1.intersection(neighbors2))
             if idxJ == []:
                return None
             else:
                old2 = idx2
                idx2 = idxJ[0]

          visited = [idx2]
          queue.enqueue(idx2)
          # Start BFS
          while len(queue) != 0:
               # Take node out of queue
               target_idx = queue.dequeue()
               # Get neighbors
               neighbors  = list(self._ugdict[target_idx])
               if target_idx == idx2:
                  neighbors.remove(idx1)
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   queue.enqueue(neighbor)
          visited.remove(idx2)
          if old2 is not None: visited.remove(old2)
          return visited

      def get_fragments(self):
          fragments = []
          nodes = self._ugdict.keys()

          visited_nodes = set([])
          for node in nodes:
              if node in visited_nodes: continue
              fragment = self.bfsearch(node)
              visited_nodes = visited_nodes.union(fragment)
              fragments.append(fragment)
          return fragments

      def longest_path(self,start_idx,visited=[]):
          '''
          Naive algorithm to explore the graph, starting at start_idx,
          and return the longest path
          Limitations:
            * start_idx has to be a terminal node
              if not, part of the path may be ommitted
          '''
          # Get neighbors, excluding previously visited ones
          neighbors = [node for node in self._ugdict[start_idx] if node not in visited]

          if len(neighbors) == 0: return [start_idx]
          #if visited == []: assert len(neighbors) == 1

          # Get longest from non-visited neighbors
          length = - float("inf")
          for neighbor in neighbors:
              visited_i = visited + [start_idx,neighbor]
              path_i    = self.longest_path(neighbor,visited=visited_i)
              if len(path_i) > length:
                 length = len(path_i)
                 the_path = path_i
          return [start_idx] + the_path

      def get_layers(self,center):
          '''
           returns a list of layers for the node center
              * 1st layer: neighbors of node center
              * 2nd layer: neighbors of neighbors of center (excluding repetitions of previous layers)
          '''
          layers  = [set([center])]
          current = [center]
          visited = set([center])
          nnodes  = len(self._ugdict.keys())

          while len(visited) != nnodes:
                layer = []
                for node in current:
                    neighbors = self._ugdict[node]
                    layer     = layer + list(neighbors)
                layer = set(layer).difference(visited)
                visited = visited.union(layer)
                layers.append(layer)
                current = list(layer)
         #for layer in layers: print layer
         #raw_input("")
          return layers
                
      #----------------------------#
      # Get matrix representations #
      #----------------------------#
      def gen_laplacian(self):
          num_nodes = self.get_nnodes()
          laplacian = np.zeros((num_nodes,num_nodes))
          for node in self._ugdict.keys():
              neighbors = self._ugdict[node]
              for neighbor in neighbors:
                  laplacian[node,node] = laplacian[node,node] + 1
                  laplacian[node,neighbor] = -1

          # Eigenvalues
          vals, vecs = np.linalg.eigh(laplacian)

          # Degenerancies?
          degs = [0]*len(vals)
          for i in range(len(vals)):
              val_i = vals[i]
              for j in range(len(vals)):
                 val_j = vals[j]
                 if abs(val_i-val_j) < 1e-3: degs[i] = degs[i]+1

          # Data for each node
          dict_vecs = {}
          for node in self._ugdict.keys():
              vector = [ abs(float(vecs[node,idx])) for idx in range(len(vals)) if degs[idx] == 1]
              dict_vecs[node] = vector
          return dict_vecs


#>>>>>>>>>>>>>>>>>>*
# CLASS: Struct    *
#>>>>>>>>>>>>>>>>>>*
class Struct(object):
    '''
    The Struct class - includes graph theory implementations

    Mandatory variables:
       * name
         a string to give a name to the structure

       * x_cc
         the non-scaled cartesian coordinates;
         both 3Nx1 and Nx3 formats are allowed

       * atonums_or_symbols
         a list of atomic numbers or a list of atomic symbols

    Non-mandatory variables:
       * masslist
         a list of atomic masses or a list of (atom index,mass) tuples
         if None, masses are taken from a dictionary in cons module

       * stype
          0 for a stationary point with 0 imaginary frequencies
          1 for a stationary point with 1 imaginary frequencies
          2 for a stationary point with 2 imaginary frequencies
          ...
         -1 for a non-stationary point
    '''

    def __init__(self, name, x_cc, atonums_or_symbols, masslist=None, stype=0):

       self._name    = name
       self._type    = stype
       self._natoms  = len(atonums_or_symbols)
       self._xcc     = hf.xvecformat(x_cc,self._natoms,'3Nx1')

       #---------------------------#
       # Deal with atonums/symbols #
       #---------------------------#
       # atonums_or_symbols is atonums
       try:
          self._atonums = [int(atonum) for atonum in atonums_or_symbols]
          self._symbols = [cons.dict_z2symbol[atonum].strip() for atonum in self._atonums]
       # atonums_or_symbols is symbols
       except:
          self._symbols = [symbol.strip() for symbol in atonums_or_symbols]
          self._atonums = [cons.dict_symbol2z[symbol] for symbol in self._symbols]
       # Get molecular formula
       self._molformula = hf.get_molformula(self._symbols)

       #-----------------------------------------------#
       # Deal with masslist and isotopic modifications #
       #-----------------------------------------------#
       self._masslist = [cons.dict_atomasses[atonum] for atonum in self._atonums]
       if masslist is not None:
          # masslist is a list of (idx,mass) tuples
          try:
             for idx,mass in masslist:
                 # mass may be a string, like "2H"
                 if mass in cons.dict_isomasses:
                    mass = cons.dict_isomasses[mass]
                 self._masslist[idx] = float(mass)
          # masslist is a list of masses
          except:
             self._masslist = masslist
       self._totmass  = sum(self._masslist)
       self._masslist = np.array(self._masslist)

       #---------------------------#
       # Generate undirected graph #
       #---------------------------#
       self._ugraph     = UGRAPH()
       for idx in range(self._natoms): self._ugraph.add_node(idx)

       #-----------------------------#
       # Basic analysis of structure #
       #-----------------------------#
       self._linear = hf.isitlinear(self._xcc)
       if   self._natoms == 1:
            self._ntra = 3
            self._nrot = 0
            self._nvib = 0
            self._kind = "Atom"
       elif self._natoms == 2:
            self._ntra = 3
            self._nrot = 2
            self._nvib = 1
            self._kind = "Diatomic Molecule"
       elif self._linear:
            self._ntra = 3
            self._nrot = 2
            self._nvib = 3*self._natoms - self._ntra - self._nrot
            self._kind = "Linear Molecule"
       else:
            self._ntra = 3
            self._nrot = 3
            self._nvib = 3*self._natoms - self._ntra - self._nrot
            self._kind = "Non-linear Molecule"

       #--------------------------------#
       # Initializing rest of variables #
       #--------------------------------#
       self._freqscal  = 1.0
       self._rmode     = 1
       self._mu        = 1.0 / cons.amu
       self._ch        = None
       self._mtp       = None
       self._elstates  = None
       self._gcc       = None
       self._Fcc       = None
       self._Etot      = None
       self._pgroup    = None
       self._rotsigma  = None

       self._imoments  = None
       self._xms       = None
       self._gms       = None
       self._Fms       = None
       self._ccfreqs   = None
       self._icfreqs   = None
       self._zpe       = None
       self._Vadi      = None
       self._rotT      = None
       self._vibT      = None

       self._v0        = None # in mass-scaled
       self._v1        = None # in mass-scaled
       self._meppoint  = None

    def __str__(self):
        return self._molformula

    def __repr__(self):
        return "struct"

    #-----------------------------------#
    # Set/Get variable                  #
    #-----------------------------------#
    def set(self,variables,values):

        if type(variables) == type(""):
           variables = [variables]
           values    = [values]

        for variable,value in zip(variables,values):
            if   variable == "freqscal" : self._freqscal  = value
            elif variable == "masslist" : self._masslist  = value; self._totmass = sum(self._masslist)
            elif variable == "mu"       : self._mu        = value
            elif variable == "gcc"      : self._gcc       = hf.xvecformat(value,self._natoms,'3Nx1')
            elif variable == "Fcc"      :
                 if value is not None   : self._Fcc       = hf.hessianformat(value,self._natoms)
            elif variable == "ch"       : self._ch        = float(value)
            elif variable == "mtp"      : self._mtp       = float(value)
            elif variable == "Etot"     : self._Etot      = float(value)
            elif variable == "pgroup"   : self._pgroup    = str(value)
            elif variable == "rotsigma" : self._rotsigma  = int(value)
            elif variable == "elstates" : self._elstates  = value
            elif variable == "meppoint" : self._meppoint  = value
            elif variable == "v0"       : self._v0        = value
            elif variable == "v1"       : self._v1        = value
            else                        : sys.exit("Trying to set an unknown variable (%s)..."%variable)

    def get(self,variable):
        # Variables defined in __init__
        if   variable == "name"      : return self._name
        elif variable == "type"      : return self._type
        elif variable == "xcc"       :
             if self._xcc is None    : return None
             else                    : return np.array(self._xcc,copy=True)
        elif variable == "natoms"    : return self._natoms
        elif variable == "symbols"   :
             if self._symbols is None: return None
             else                    : return list(self._symbols)
        elif variable == "atonums"   :
             if self._atonums is None: return None
             else                    : return list(self._atonums)
        elif variable == "molformula": return str(self._molformula)
        elif variable == "masslist"  : return list(self._masslist)
        elif variable == "totmass"   : return float(self._totmass)
       #elif variable == "ugraph"    : return self._ugraph
        elif variable == "lineal"    : return self._linear
        elif variable == "ntra"      : return self._ntra
        elif variable == "nrot"      : return self._nrot
        elif variable == "nvib"      : return self._nvib
        elif variable == "kind"      : return self._kind
        # Variables that can be defined through set function
        elif variable == "freqscal"  : return float(self._freqscal)
        elif variable == "mu"        : return float(self._mu)
        elif variable == "gcc"       :
             if self._gcc is None    : return None
             else                    : return np.array(self._gcc,copy=True)
        elif variable == "Fcc"       :
             if self._Fcc is None    : return None
             else                    : return np.matrix(self._Fcc,copy=True)
        elif variable == "ch"        : return int(self._ch)
        elif variable == "mtp"       : return int(self._mtp)
        elif variable == "Etot"      : return float(self._Etot)
        elif variable == "pgroup"    : return str(self._pgroup)
        elif variable == "rotsigma"  : return int(self._rotsigma)
        elif variable == "elstates"  : return list(self._elstates)
        elif variable == "meppoint"  : return self._meppoint
        elif variable == "v0"        :
             if self._v0 is None     : return None
             else                    : return np.array(self._v0,copy=True)
        elif variable == "v1"        :
             if self._v1 is None     : return None
             else                    : return np.array(self._v1,copy=True)
        # Derivated variables that cannot be defined using set function
        elif variable == "imoments"  : return list(self._imoments)
        elif variable == "xms"       :
             if self._xms is None    : return None
             else                    : return np.array(self._xms,copy=True)
        elif variable == "gms"       :
             if self._gms is None    : return None
             else                    : return np.array(self._gms,copy=True)
        elif variable == "Fms"       :
             if self._Fcc is None    : return None
             else                    : return np.matrix(self._Fms,copy=True)
        elif variable == "ccfreqs"   :
             if self._ccfreqs is None: return None
             else                    : return list(self._ccfreqs)
        elif variable == "icfreqs"   : return list(self._icfreqs)
        elif variable == "zpe"       : return float(self._zpe)
        elif variable == "Vadi"      : return float(self._Vadi)
        elif variable == "rotT"      : return list(self._rotT)
        elif variable == "vibT"      : return list(self._vibT)

        else                     : sys.exit("Requesting an unknown variable (%s)..."%variable)


    def get_ifreq(self):
        if self._ccfreqs is None: return None
        for freq in self._ccfreqs:
            if freq.isItImag(): return freq.copy()

    #-----------------------------------#
    # Basic setup of structure          #
    #-----------------------------------#
    def basic_setups(self,setup="all"):
        '''
        '''
        if type(setup) == type(1): setup = [setup]
        if      setup  == "all"  : setup = range(5)


        for stype in setup:

            # Origin in center of mass
            if stype == 0:
               self._xcc = hf.shift2cm(self._xcc,self._masslist)

            # Point Group and rotsigma
            if stype == 1 and self._rotsigma is None:
               self._pgroup, self._rotsigma = hf.get_pgs(self._atonums,self._masslist,self._xcc)

            # elstates
            if stype == 2 and self._elstates is None:
                  self._elstates = [ (0.0,self._mtp) ]

            # Mass-scale
            if stype == 3:
               self._xms = hf.cc2ms_x(self._xcc,self._masslist,self._mu)
               if self._gcc is not None:
                  self._gms = hf.cc2ms_g(self._gcc,self._masslist,self._mu)
               if self._Fcc is not None:
                  self._Fms = hf.cc2ms_F(self._Fcc,self._masslist,self._mu)

            # Rotational temperatures
            if stype == 4 and self._natoms > 1:
               # (a) Calculation of the inertia tensor (a.u.)
               inertia = np.zeros((3,3))
               for i in range(self._natoms):
                   # Diagonal elements
                   inertia[0][0] += self._masslist[i] * (self._xcc[i*3+1]**2 + self._xcc[i*3+2]**2)
                   inertia[1][1] += self._masslist[i] * (self._xcc[i*3+0]**2 + self._xcc[i*3+2]**2)
                   inertia[2][2] += self._masslist[i] * (self._xcc[i*3+0]**2 + self._xcc[i*3+1]**2)
                   # Offdiagonal elements
                   inertia[0][1] -= self._masslist[i] * self._xcc[i*3+0] * self._xcc[i*3+1]
                   inertia[0][2] -= self._masslist[i] * self._xcc[i*3+0] * self._xcc[i*3+2]
                   inertia[1][2] -= self._masslist[i] * self._xcc[i*3+1] * self._xcc[i*3+2]
               inertia[1][0] = inertia[0][1]; inertia[2][0] = inertia[0][2]; inertia[2][1] = inertia[1][2]
               # (b) Get its eigenvalues
               evalsI, evecsI = np.linalg.eigh(inertia)
               self._imoments = []
               if self._nrot == 2:
                   evalsI = np.sort(evalsI)
                   cocient = evalsI[1] / evalsI[2]
                   if abs(cocient - 1.0) < cons.ZERO: self._imoments = [evalsI[1]]
                   else: sys.exit("ERROR! Molecule is not linear!")
               else:
                   self._imoments = evalsI
               # (c) Get rotational Temperatures
               self._rotT   = [cons.hbar**2 / (2*I_i*cons.kB) for I_i in self._imoments]

    #-----------------------------------#
    #  Cartesian-Coordinate FREQUENCIES #
    #-----------------------------------#
    def calc_ccfreqs(self):
        '''
        info about projection matrix: Appendix D of JChemPhys 1988, 88, 922-935
        '''
        if self._natoms < 2         : return None
        if self._ccfreqs is not None: return None

        #-----------------------#
        # Get projection matrix #
        #-----------------------#
        # (a) Translation vectors (b1,b2,b3)
        T_vecs = []
        for i in range(3):
            T = np.zeros(3*self._natoms)
            for j in range(self._natoms):
                T[3*j+i] = math.sqrt(self._masslist[j])
            T = T / np.linalg.norm(T)
            T_vecs.append(T)
        # (b) Rotation vectors (b4,b5,b6)
        R1 = np.zeros(3*self._natoms)
        R2 = np.zeros(3*self._natoms)
        R3 = np.zeros(3*self._natoms)
        for i in range(self._natoms):
            R1[3*i + 1] =   math.sqrt(self._masslist[i]) * self._xcc[3*i + 2]
            R1[3*i + 2] = - math.sqrt(self._masslist[i]) * self._xcc[3*i + 1]
            R2[3*i + 0] = - math.sqrt(self._masslist[i]) * self._xcc[3*i + 2]
            R2[3*i + 2] =   math.sqrt(self._masslist[i]) * self._xcc[3*i + 0]
            R3[3*i + 0] =   math.sqrt(self._masslist[i]) * self._xcc[3*i + 1]
            R3[3*i + 1] = - math.sqrt(self._masslist[i]) * self._xcc[3*i + 0]
        R_vecs = []
        for R in (R1,R2,R3):
            normR = np.linalg.norm(R)
            if normR > 1e-7: R = R / normR; R_vecs.append(R)
        # (c) Apply Gram Schmidt method (v1 to v6 vectors)
        X = np.matrix(T_vecs+R_vecs).transpose() # each column is a vector
        X_gs, R = np.linalg.qr(X)
        # (d) Get Translation-Rotation matrix
        R_matrix  = X_gs * X_gs.H
        # (e) Get projection matrix (eq D1)
        if self._type == -1: proj_matrix = R_matrix + np.matrix(self._v0).transpose() * np.matrix(self._v0)
        else               : proj_matrix = R_matrix

        #-----------------------#
        # Get projected hessian #
        #-----------------------#
        # Get identity matrix
        I = np.identity(3*self._natoms)
        # Calculate projected hessian matrix
        proj_F = (I - proj_matrix) * self._Fms * (I - proj_matrix)

        # Diagonalize
        evalsF, evecsF = np.linalg.eigh(proj_F)

        # Remove translation and rotation degrees of freedom
        removals = self._ntra+self._nrot
        if self._type == -1: removals += 1
        for removal in range(removals):
            ref_value = float("inf")
            for idx in range(len(evalsF)):
                absval = abs(evalsF[idx])
                if absval < ref_value:
                   ref_value = absval
                   target    = idx
            evalsF = np.delete(evalsF,target)
            evecsF = np.delete(evecsF,target,1)

        # Generate list of cc vibrational frequencies
        self._ccfreqs = []
        for idx in range(len(evalsF)):
            FreqInst = Freq(scaling=self._freqscal , mu=self._mu, rmode= self._rmode)
            FreqInst.set_evalue(  evalsF[idx]   )
            # Get evector as a numpy array
            evector = np.array(evecsF[:,idx].transpose().tolist()[0])
            FreqInst.set_evector( evector )
            FreqInst.calc_derivated_magnitudes()
            self._ccfreqs.append(  FreqInst  )

    #-----------------------------------#
    #  Internal-Coordinate FREQUENCIES  #
    #-----------------------------------#
    def get_Dmatrix(self,nricoords):
        '''
        Get the D matrix for two torsions
        Last two torsions are the selected torsions
        '''
        B_wilson, C_wilson = hf.get_B_and_C(self._xcc, self._natoms, nricoords)

        mass_array = []
        for m in self._masslist: mass_array += [m,m,m]
        u = [ 1.0/mass for mass in mass_array]
        u = np.diag(u)

        # Calculate h matrix (h* in [3], cursive G in [4])
        h = B_wilson * u * B_wilson.transpose()

        # Calculate D matrix
        Dmatrix = np.linalg.inv(h)

        # Units of each element of Dmatrix is [distance]*[mass]^2 (in a.u.)
        Dmatrix = Dmatrix[-2:,-2:]

        return Dmatrix

    def calc_icfreqs(self,icoords):
        '''
        '''

        self._icfreqs = []
    
        #---------------------------#
        # Get B matrix and C tensor #
        #---------------------------#
        B_wilson, C_wilson = hf.get_B_and_C(self._xcc, self._natoms, icoords)
        num_ric, num_3N = B_wilson.shape
        num_nric = self._nvib
    
        #--------------------------#
        # Calculate h*, h^-1 and A #
        #--------------------------#
        mass_array = []
        for m in self._masslist: mass_array += [m,m,m]
        u = [ 1.0/mass for mass in mass_array]
        u = np.diag(u)
        # Calculate h matrix (h* in [3], cursive G in [4])
        h = B_wilson * u * B_wilson.transpose()
        # Calculate inverse of h*
        if num_ric <= num_nric:
           h_inv = np.linalg.inv(h)
        else:
           h_evals, h_evecs = np.linalg.eigh(h)
           K , Kprime = [], []
           Gamma = []
           Gamma_zeros = []
           for idx in range(len(h_evals)):
               h_eval = h_evals[idx]
               h_evec = [float(i) for i in h_evecs[:,idx]]
               # Non-zero eigenvalue
               if abs(h_eval) > 1e-15:
                  K.append( h_evec )
                  Gamma.append( h_eval )
               # Zero eigenvalue - save as zero
               else:
                  Kprime.append( h_evec )
                  Gamma_zeros.append( 0.0 )
           h_evecs_sorted = np.matrix(K + Kprime)
           h_evecs_sorted = h_evecs_sorted.transpose()
           h_evals_inv    = [1.0/evalue for evalue in Gamma] + Gamma_zeros
           h_evals_inv    = np.diag(h_evals_inv)
           h_inv = h_evecs_sorted * h_evals_inv * h_evecs_sorted.transpose()
        # Get A matrix
        A = u * B_wilson.transpose() * h_inv
    
        #--------------------------------------------------------#
        # Gradiend and hessian in redundant internal coordinates #
        #--------------------------------------------------------#
        F_ric = A.transpose() * self._Fcc * A
        # A MEP point is being analyzed
        if self._type == -1:
           g_cc  = np.array(self._gcc,copy=True)
           g_cc  = np.matrix(g_cc).transpose() # save as column matrix
           g_ric = A.transpose() * g_cc
           for idx in range(g_ric.shape[0]):
               F_ric -= float(g_ric[idx]) * A.transpose() * C_wilson[idx] * A
           proj_rc = True
        # A critical point of the PES is being analyzed
        else:
           g_cc  = np.zeros((3*self._natoms,3*self._natoms))
           g_ric = A.transpose() * g_cc
           proj_rc = False
    
        #------------------------------------------------------------#
        # Gradiend and hessian in non-redundant internal coordinates #
        # and with the reaction coordinate projected out             #
        #------------------------------------------------------------#
        # Project from ric to non-redundant internal coordinates (nric)
        P = h * h_inv
        f_nric = P * F_ric * P
        # Project out the reaction coordinate
        if proj_rc:
           g_nric = P * g_ric
           p      = g_nric * g_nric.transpose() / (g_nric.transpose() * h * g_nric)
           proj_F = (np.identity(num_ric) - p*h) * f_nric * (np.identity(num_ric)-h*p)
        else:
           proj_F = f_nric
        # Get eigenvalues and eigenvectors
        # PS: as h*proj_F isn't symmetric, pF_evals and pF_evecs have imaginary components
        pF_evals, pF_evecs = np.linalg.eig(h*proj_F)

        #---------------------------------#
        # Get indices of zero eigenvalues #
        #---------------------------------#
        indices_zero = []
        for idx in range(len(pF_evals)):
            F_evalue = pF_evals[idx]
            F_evalue = F_evalue*self._mu
            freq = Freq(scaling=self._freqscal,mu=self._mu,rmode= self._rmode)
            freq.set_evalue(F_evalue)
            if abs(freq.get("wavenum")) < cons.ZERO_wavenum: indices_zero.append(idx)
    
        #--------------------------------------------------------------#
        # Get eigenvectors in mass-scaled cartesian coordinates        #
        # PS: single value decomposition (svd) used to get pF_evecs^-1 #
        #--------------------------------------------------------------#
        u_svd, s_svd, v_svd = np.linalg.svd(pF_evecs,full_matrices=True,compute_uv=True)
        # Pseudoinverse of s_svd
        s_svd_inv = []
        for idx in range(len(s_svd)):
            if abs(s_svd[idx]) < cons.ZERO2: s_svd_inv.append( s_svd[idx] )
            else:                       s_svd_inv.append( s_svd[idx]**-1.0 )
        pF_evecs_inv = np.dot(v_svd.transpose(), np.dot(np.diag(s_svd_inv),u_svd.transpose()))
        # Get normalized vectors (but only those of non-zero eigenvalue)
        pF_norm_evecs = []
        nrows, ncols = pF_evecs.shape
        zero_evec    = [0.0] *  nrows
        C = pF_evecs_inv * h * pF_evecs_inv.transpose()
        for idx in range(len(pF_evals)):
          # Only normalize those of non-zero evalue
          if abs(pF_evals[idx].real) > cons.ZERO2:
             # remove imag part of C element
             if abs( C[idx,idx].imag ) < cons.ZERO: Cii = float(C[idx,idx].real)
             else:                                  Cii = C[idx,idx]
             # Get normalized L vector
             norm_L = pF_evecs[:,idx] * np.sqrt(Cii)
             pF_norm_evecs.append( norm_L.transpose().tolist()[0] )
          else:
             pF_norm_evecs.append( zero_evec )
        pF_norm_evecs = np.matrix(pF_norm_evecs).transpose()
    
        # Get chi matrix
        chi_matrix = A * pF_norm_evecs
    
        # Get normal-mode eigenvectors in mass-scaled cartesian
        self._icfreqs = []
        for j in range(num_ric):
            if j in indices_zero: continue
            F_evalue  = pF_evals[j]
            F_evector = chi_matrix[:,j]
            # mass-scalde eigenvalue
            F_evalue = F_evalue*self._mu
            # mass-scaled eigenvector
            for i in range(3*self._natoms):
                m_i = mass_array[i]
                F_evector[i] = np.sqrt(m_i) * F_evector[i]
            if np.linalg.norm(F_evector) != 0.0: F_evector = F_evector / np.linalg.norm(F_evector)
            # Evector (remove zero imaginary component)
            F_evector = np.array([float(Fi.real) for Fi in F_evector.transpose().tolist()[0]])
            # Generate frequency
            freq = Freq(scaling=self._freqscal,mu=self._mu,rmode= self._rmode)
            freq.set_evalue(F_evalue)
            freq.set_evector(F_evector)
            freq.calc_derivated_magnitudes()
            self._icfreqs.append( (F_evalue,freq) )
        # Sort freqs
        self._icfreqs.sort()
        self._icfreqs = [j for (i,j) in self._icfreqs]

    def check_icoords(self,icoords):
        '''
        Compare frequencies obtained with
        cartesian coordinates and with icoords
        '''
        # cc-freqs
        ccfreqs = self.get("ccfreqs")
        if ccfreqs is None or ccfreqs == []:
           print "ERROR: cc-freqs are not calculated!"
           sys.exit()
        # ic-freqs
        self.calc_icfreqs(icoords)
        icfreqs = self.get("icfreqs")
        # Compare number of freqs
        if len(ccfreqs) != len(icfreqs): return False
        # Compare frequency value
        for idx in range(len(ccfreqs)):
            ccfreq = ccfreqs[idx].get("wavenum")/cons.cm
            icfreq = icfreqs[idx].get("wavenum")/cons.cm
            diff = abs(ccfreq-icfreq)
            if diff > 0.5: return False
        return True

    #-----------------------------------#
    #  Adiabatic potential              #
    #-----------------------------------#
    def calc_Vadi(self,k="cc",nimag=None):

        if nimag is None: nimag = self._type
        if nimag < 0.0:   nimag = 0
        self._zpe  = 0.0
        if k == "cc" and self._natoms>1: 
           for ccfreq in self._ccfreqs[nimag:]:
               self._zpe += ccfreq.get("zpe")
        if k == "ic" and self._natoms>1: 
           for icfreq in self._icfreqs[nimag:]:
               self._zpe += icfreq.get("zpe")
        self._Vadi = self._Etot + self._zpe

    #-----------------------------------#
    #  Partition Functions              #
    #-----------------------------------#
    def get_partition_functions(self,T,k="cc"):
        '''
        Input:
          * T: the temperature, in Kelvin degrees
          * k: to indicate if Cartesian-coordinate frequencies (cc)
               or internal coordinate frequencies are to be used
        Returns:
          * ptra (phi_tra) : per volume unit (bohr^3)
          * Qrot
          * Qvib
          * Qele
          * Vadi - the reference energy in hartree
        '''

        # Number of vibrational freqs to consider
        if self._type >= 0: nimag = self._type
        else:               nimag = 0

        # Sort elstates
        self._elstates.sort()

        # In case of asking for zero kelvin
        if T == 0.0:
           ptra = 0.0
           qrot = 0.0
           qvib = 1.0
           qele = self._elstates[0][1]
           if k == "cc": zpe = sum([freq.get("zpe") for freq in self._ccfreqs[nimag:]])
           if k == "ic": zpe = sum([freq.get("zpe") for freq in self._icfreqs[nimag:]])

        # At other temperature
        else:
           beta = 1.0 / (cons.kB * T)

           # (a) Translational partition function
           ptra = ( 2 * np.pi * self._totmass * cons.kB * T )**(3./2.) / (cons.h**3)

           # (b) Rotational partition function
           qrot = 1.0
           if   self._nrot == 2:
                qrot = 1.0 * (T / self._rotT[0])
           if self._nrot == 3:
                product = self._rotT[0] * self._rotT[1] * self._rotT[2]
                qrot = math.sqrt(np.pi * T**3 / product)
           qrot = qrot/self._rotsigma

           # (c) Vibrational partition function
           qvib = 1.0
           zpe  = 0.0
           if self._nvib != 0:
              if k=="cc": 
                 for freq in self._ccfreqs[nimag:]:
                     qvib_i, zpe_i = freq.get_qvib(T)
                     qvib *= qvib_i
                     zpe  += zpe_i
              if k=="ic": 
                 for freq in self._icfreqs[nimag:]:
                     qvib_i, zpe_i = freq.get_qvib(T)
                     qvib *= qvib_i
                     zpe  += zpe_i

           # (d) Electronic partition function
           qele = sum( [mtp*np.exp(-beta*relE) for relE,mtp in self._elstates]  )

        # Reference for molecule and external reference
        Vadi = self._Etot + zpe

        return ptra, qrot, qvib, qele, Vadi

    def get_pfns(self,T,k="cc"):
        # Single temperature or a list of them
        singleT = (type(T) == type(0)) or (type(T) == type(0.0))
        # Return data for single T
        if singleT:
           return self.get_partition_functions(T,k)
        # Return data for list of T
        else:
           vec_ptra = []
           vec_qrot = []
           vec_qvib = []
           vec_qele = []
           for idx in range(len(T)):
               ptra, qrot, qvib, qele, Vadi = self.get_partition_functions(T[idx],k)
               vec_ptra.append(ptra)
               vec_qrot.append(qrot)
               vec_qvib.append(qvib)
               vec_qele.append(qele)
           vec_ptra = np.array(vec_ptra)
           vec_qrot = np.array(vec_qrot)
           vec_qvib = np.array(vec_qvib)
           vec_qele = np.array(vec_qele)
           return vec_ptra,vec_qrot,vec_qvib,vec_qele, Vadi

    #-----------------------------------#
    # Functions related to the atoms    #
    #-----------------------------------#
    def atom_xcc(self,idx):
        return np.array(self._xcc[3*idx:3*idx+3],copy=True)


    #-----------------------------------#
    # Distances/angles in the molecule  #
    #-----------------------------------#
    def icvalue(self,indices):
        # Asking for distance
        if len(indices) == 2:
           idx1, idx2 = indices
           x1 = self.atom_xcc(idx1)
           x2 = self.atom_xcc(idx2)
           return float(hf.calc_distance(x1,x2))
        # Asking for angle
        if len(indices) == 3:
           idx1, idx2, idx3 = indices
           x1 = self.atom_xcc(idx1)
           x2 = self.atom_xcc(idx2)
           x3 = self.atom_xcc(idx3)
           return float(hf.calc_angle(x1,x2,x3))
        # Asking for dihedral
        if len(indices) == 4:
           idx1, idx2, idx3, idx4 = indices
           x1 = self.atom_xcc(idx1)
           x2 = self.atom_xcc(idx2)
           x3 = self.atom_xcc(idx3)
           x4 = self.atom_xcc(idx4)
           return float(hf.calc_dihedral(x1,x2,x3,x4))

    def dist_between_frags(self, fragmentA, fragmentB):
        min_distance = + float("inf")
        for atom_i in fragmentA:
            for atom_j in fragmentB:
                d_ij = self.icvalue( (atom_i,atom_j) )
                if d_ij < min_distance:
                   min_distance = d_ij
                   idxA, idxB   = atom_i, atom_j
        return min_distance, idxA, idxB


    #-----------------------------------#
    # Functions related to Graph Theory #
    #-----------------------------------#
    def graph_getbonds(self):
        bonds = self._ugraph.get_edges()
        return bonds

    def graph_addbond(self, idx1, idx2):
        self._ugraph.add_edge(idx1,idx2)

    def graph_removebond(self, idx1, idx2):
        self._ugraph.remove_edge(idx1,idx2)

    def graph_nbonds(self):
        return self._ugraph.get_nedges()

    def graph_autoconnect(self,lengthfactor=1.3,bonds=[]):
        for idx1 in range(self._natoms):
            atonum1 = self._atonums[idx1]
            for idx2 in range(idx1+1,self._natoms):
                atonum2 = self._atonums[idx2]
                if (idx1,idx2) in bonds or (idx2,idx1) in bonds:
                   connected = True
                else:
                   dist = self.icvalue( (idx1,idx2) )
                   dref = cons.dict_covradii[atonum1] + cons.dict_covradii[atonum2]
                   dref = lengthfactor*dref
                   if dist <= dref: connected = True
                   else:            connected = False
                # Connected atoms
                if connected: self.graph_addbond(idx1,idx2)

    def graph_fragconnect(self,maxfrags=1):
        '''
        connect defined fragments until you have,
        as much, a total of maxfrags
        '''

        fragments = self._ugraph.get_fragments()
        nf        = len(fragments)
        if nf <= maxfrags: return []

        # Calculate min distances between fragments
        distances = []
        for fA in range(nf):
            fragmentA = fragments[fA]
            for fB in range(fA+1,nf):
                fragmentB = fragments[fB]
                # Get distance
                dist_AB, idxA, idxB = self.dist_between_frags(fragmentA, fragmentB)
                distances.append( (dist_AB,fA,fB,idxA,idxB)  )

        # Sort by distances
        distances.sort()

        # Create bonds
        added_edges = []
        connections = []
        for dist_AB,fA,fB,idx1,idx2 in distances:
            connect = True
            append  = True
            for connection in connections:
                if   fA in connection and fB in connection:
                    connect = False
                    append  = False
                elif fA in connection and fB not in connection:
                    connection += [fB]
                    append  = False
                elif fA not in connection and fB in connection:
                    connection += [fA]
                    append  = False
            if append: connections.append( [fA,fB] )
            if connect:
               self.graph_addbond(idx1,idx2)
               added_edges.append( (idx1,idx2) )
               nf -= 1
            if nf <= maxfrags: return added_edges

    def graph_neighbors(self,idx):
        return list(self._ugraph.neighbors(idx))

    def graph_valence(self,idx):
        return len(self._ugraph.neighbors(idx))

    def graph_atostretchings(self,idx):
        neighbors = self.graph_neighbors(idx)
        neighbors.sort()
        N = len(neighbors)

        stretchings = []
        for neighbor in neighbors:
            if idx > neighbor: stretchings.append( (neighbor,idx) )
            else             : stretchings.append( (idx,neighbor) )
        return stretchings

    def graph_atobendings(self,idx):
        neighbors = self.graph_neighbors(idx)
        neighbors.sort()
        N = len(neighbors)

        bendings = []
        if N > 1:
           for idx1 in range(N):
               IDX1 = neighbors[idx1]
               for idx2 in range(idx1+1,N):
                   IDX2 = neighbors[idx2]
                   bending = (IDX1,idx,IDX2)
                   bendings.append(bending)
        return bendings

    def graph_atowaggings(self,idx):
        neighbors = self.graph_neighbors(idx)
        neighbors.sort()
        N = len(neighbors)

        waggings = []
        if N > 2:
           for idx1 in range(N):
               IDX1 = neighbors[idx1]
               for idx2 in range(idx1+1,N):
                   IDX2 = neighbors[idx2]
                   for idx3 in range(idx2+1,N):
                       IDX3 = neighbors[idx3]
                       wagging = (IDX1,IDX2,IDX3,idx)
                       waggings.append(wagging)
        return waggings

    def graph_nolinealangle(self,M_idx,N_idx,visited=[]):
        '''
        Given the selected bond in the molecule (M-N)
        this functions returns the first angle that,
        starting at M and going in the direction of N,
        presents a non-lineal value
    
                       ...-M-N-->...
    
        * Returns a tuple of three indices starting at M_idx
        * Returns None if no angle was found
        '''
    
        # Get nodes and corresponding neighbors
        M_neighbors = self.graph_neighbors(M_idx)
        N_neighbors = self.graph_neighbors(N_idx)
    
        # Remove M from N neighbors and vice versa
        if M_idx in N_neighbors: N_neighbors.remove(M_idx)
        if N_idx in M_neighbors: M_neighbors.remove(N_idx)
    
        # Remove common nodes, if any
        common_nodes = list(set(M_neighbors).intersection(set(N_neighbors)))
        for node in common_nodes:
            M_neighbors.remove(node)
            N_neighbors.remove(node)
    
        # Sort neighbors of N according to connectivity and to atomic number
        sorted_Nneighbors = []
        for neighbor in N_neighbors:
            atonum  = self._atonums[neighbor]
            valence = self.graph_valence(neighbor)
            sorted_Nneighbors += [(valence,atonum,neighbor)]
        sorted_Nneighbors.sort(reverse=True)
        N_neighbors = [neighbor for nbonds,atonum,neighbor in sorted_Nneighbors]
    
        # Define all visited neighbors
        visited = visited + [M_idx,N_idx] + M_neighbors + common_nodes
    
        # No neighbors to visit??
        if N_neighbors == []: return None
    
        next_N = None
        for L_idx in N_neighbors:
            if L_idx in visited: continue
            else: visited.append(L_idx)
            next_N = L_idx
            theta = self.icvalue( (M_idx,N_idx,L_idx) )
            if theta < cons.LINEAR: return (M_idx,N_idx,L_idx)
    
        if next_N is None: return None
        return self.graph_nolinealangle(M_idx,next_N,visited=visited)


    def graph_torsion(self,bond):
        '''
        Returns the atoms involved in the
        torsion coordinate around the bond
        '''

        M_idx, N_idx = bond

        # Valid L neighbors of M
        NML = self.graph_nolinealangle(N_idx,M_idx)

        # Valid O neighbors of N
        MNO = self.graph_nolinealangle(M_idx,N_idx)

        if NML is None: return None
        if MNO is None: return None

        L, M, N, O = NML[2],NML[1],MNO[1],MNO[2]
        if   L > O: return (L,M,N,O)
        elif L < O: return (O,N,M,L)
        else      : return None
 
    def graph_irotation(self,bonds,thetas):
        '''
        To generate a structure after internal rotations around bonds
        * bonds: a list of bonds to rotate around
        * thetas: the corresponding angle of rotation around each bond
        Generates the rotation around each bond 
        '''

        if type(thetas) != type(list()) and type(thetas) != type(tuple()):
           bonds  = [bonds]
           thetas = [thetas]

        # The geometry
        xvector = np.array(self._xcc,copy=True)
        failed = False

        for bond , theta in zip(bonds,thetas):

            # Get axis vector and origin
            idxA, idxB = bond

            # Get two fragments
            A_frag = set(self._ugraph.bfsearch1d(bond[1],bond[0]))
            B_frag = set(self._ugraph.bfsearch1d(bond[0],bond[1]))

            # Compare fragments. They may be equal in case of cyclic systems
            if (B_frag is None) or (A_frag is None) or (B_frag == A_frag):
               failed = True
               break

            # Choose smaller fragment
            if len(A_frag) > len(B_frag):
               # if B_frag, rotation around A-->B
               x0   = xvector[3*idxA:3*idxA+3]
               axis = xvector[3*idxB:3*idxB+3] - x0
               target_fragment = B_frag.copy()
            else:
               # if A_frag, rotation around B-->A
               x0   = xvector[3*idxB:3*idxB+3]
               axis = xvector[3*idxA:3*idxA+3] - x0
               target_fragment = A_frag.copy()
            axis = axis / np.linalg.norm(axis)

            # Remove indices of the bond
            target_fragment.discard(idxA)
            target_fragment.discard(idxB)

            # Get rotation matrix
            R, rot_theta = hf.gen_rotmatrix(axis,theta)

            # Rotate atoms in fragment
            rotated_xyz = []
            for idx in range(self._natoms):
                xyz    = xvector[3*idx:3*idx+3]
                symbol = self._symbols[idx]
                if idx in target_fragment:
                    xyz = xyz - x0
                    xyz = R * np.matrix(xyz).transpose()
                    xyz = np.array((xyz.transpose()).tolist()[0])
                    xyz = xyz + x0
                rotated_xyz += xyz.tolist()
            xvector = np.array(rotated_xyz)

        if failed: return None

        return xvector

#    def graph_nricoords(self, utorsions = []):
#        '''
#        utorsions: a list of torsions defined by the user
#        '''
#        def helper(dihedral,deftorsions,waggings_k):
#            '''
#            dihedral: i-j-k-l, l being the target
#            '''
#            l = dihedral[3]
#            sent = "torsion"
#            dihedral = tuple(dihedral)
#            for deftorsion in deftorsions:
#                same_bond    = sorted(list(deftorsion[1:3])) == sorted(list(dihedral[1:3]))
#                same_torsion = sorted(list(deftorsion))      == sorted(list(dihedral))
#                # The dihedral is one of the defined by the user
#                if   same_torsion:
#                   dihedral = tuple(utorsion)
#                   break
#                # The dihedral cannot be used and has to be included as a wagging
#                elif same_bond:
#                   twag = None
#                   for wag in range(len(waggings_k)):
#                       if l in waggings_k[wag]:
#                          twag = wag
#                          break
#                   if twag is None:
#                      dihedral = None
#                      break
#                   else: dihedral = tuple(waggings_k.pop(twag))
#                   sent = "wagging"
#                else:
#                   continue
#            return dihedral, waggings_k, sent
#
#        utorsions = [tuple(utorsion) for utorsion in utorsions]
#
#        stretchings, bendings, waggings = [] , [] , []
#        torsions  = list(utorsions)
#
#        #---------------------#
#        # (a) Find main chain #
#        #---------------------#
#        max_length = - float("inf")
#        main_chain = None
#        # Start at terminal atoms
#        for i in range(self._natoms):
#            if self.graph_valence(i) == 1:
#               chain = self._ugraph.longest_path(i)
#               if len(chain) > max_length:
#                  max_length = len(chain)
#                  main_chain = chain
#        # In case of having no terminal atoms
#        if main_chain is None:
#           main_chain = self._ugraph.longest_path(0)
#
#        #----------------------------------#
#        # (b) A dictionary of waggings     #
#        #     and a dictionary of valences #
#        #----------------------------------#
#        dict_waggings = {}
#        dict_valences = {}
#        for atom in range(self._natoms):
#            dict_valences[atom] = self.graph_valence(atom)
#            dict_waggings[atom] = self.graph_atowaggings(atom)
#
#        #----------------------#
#        # (c) Visit main chain #
#        #----------------------#
#        ramified_nodes = []
#        for idx in range(1,len(main_chain)):
#
#            atom_i = main_chain[idx]
#
#            # Add nodes to ramified_nodes for next step
#            if dict_valences[atom_i] > 2: ramified_nodes += [atom_i]
#
#            # Define streching with previous atom
#            atom_h = main_chain[idx-1]
#            bond   = (atom_h,atom_i)
#            stretchings.append(bond)
#            if idx == 1: continue
#
#            # Define bending with previous atoms
#            atom_g = main_chain[idx-2]
#            angle  = (atom_g,atom_h,atom_i)
#            bendings.append(angle)
#            if idx == 2: continue
#
#            # Define dihedral with previous atoms
#            atom_f   = main_chain[idx-3]
#            dihedral = (atom_f,atom_g,atom_h,atom_i)
#            dihedral, dict_waggings[atom_h], sent = helper(dihedral,utorsions,dict_waggings[atom_h])
#            if dihedral is None: continue
#            if sent == "torsion": torsions.append(dihedral)
#            if sent == "wagging": waggings.append(dihedral)
#
#        #-------------------------#
#        # (d) Visit ramifications #
#        #-------------------------#
#        visited       = set(main_chain)
#        while ramified_nodes != []:
#              atom_i      = ramified_nodes.pop(0)
#              wagg_list   = dict_waggings[atom_i]
#              neighbors_i = set(self.graph_neighbors(atom_i))
#              for node in visited: neighbors_i.discard(node)
#
#              for atom_j in neighbors_i:
#
#                  # Add atom to ramified_nodes?
#                  neighbors_j = set(self.graph_neighbors(atom_j))
#                  for node in visited: neighbors_j.discard(node)
#                  if len(neighbors_j) > 0:
#                     ramified_nodes.append(atom_j)
#
#                  # The stretching
#                  stretching = (atom_i,atom_j)
#
#                  # The bending
#                  for icoord in stretchings:
#                      a,b = icoord
#                      if a == atom_i: bending  = (atom_j,a,b)
#                      if b == atom_i: bending  = (atom_j,b,a)
#
#                  # The torsion
#                  # a) was atom_j included previously in a (user-defined) dihedral?
#                  needtorsion = True
#                  for utorsion in utorsions:
#                      if utorsion[0] == atom_j or utorsion[3] == atom_j:
#                         needtorsion = False
#                         break
#                  # b) Find a torsion for this atom
#                  if needtorsion:
#                      for icoord in bendings:
#                          a,b,c = icoord
#                          if a == atom_i: dihedral = (c,b,a,atom_j)
#                          if c == atom_i: dihedral = (a,b,c,atom_j)
#                      dihedral, dict_waggings[atom_i], sent = helper(dihedral,torsions,dict_waggings[atom_i])
#                  else:
#                      dihedral = None
#
#                  # Add icoords
#                  stretchings.append( stretching )
#                  bendings.append( bending )
#                  if dihedral is None : continue
#                  if sent == "torsion": torsions.append( dihedral )
#                  if sent == "wagging": waggings.append( dihedral )
#                  visited.add(atom_j)
#
#        #-------------------------------------------#
#        # (e) Defined torsions should go at the end #
#        #-------------------------------------------#
#        final_torsions = []
#        for torsion in list(set(torsions)):
#            if torsion in utorsions: continue
#            final_torsions.append(torsion)
#        final_torsions = final_torsions + utorsions
#
#        return stretchings, bendings, list(set(waggings)), final_torsions

    def gen_ricoords(self,torsions=[],check=False,cleaned=True):
        '''
        '''
        if torsions != []:
           for idx in range(len(torsions)):
               torsion = torsions[idx]
               if torsion[0] > torsion[3]:
                  torsions[idx] = tuple(torsion[::-1])
               else:
                  torsions[idx] = tuple(torsion)

        linear  = 178.0 * cons.D2R
        angular = 168.0 * cons.D2R

        #-----------------#
        # Get stretchings #
        #-----------------#
        ics_stretch = self.graph_getbonds()
        ics_stretch = [ set(bond) for bond in ics_stretch]

        #--------------#
        # Get bendings #
        #--------------#
        ics_abend = []
        ics_lbend = []
        for idxA in range(len(ics_stretch)):
            for idxB in range(idxA+1,len(ics_stretch)):
                bondA = ics_stretch[idxA]
                bondB = ics_stretch[idxB]
                central = bondA.intersection( bondB )
                if len(central) == 1:
                   at1 = bondA.difference(central)
                   at3 = bondB.difference(central)
                   if list(at1)[0] < list(at3)[0]: bending = list(at1) + list(central) + list(at3)
                   else                          : bending = list(at3) + list(central) + list(at1)
                   angle = self.icvalue(bending)
                   if   angle > linear:
                        ics_lbend.append( tuple(bending) )
                   elif angle < angular:
                        ics_abend.append( tuple(bending) )
                   else:
                        ics_lbend.append( tuple(bending) )
                        ics_abend.append( tuple(bending) )
    
        #--------------------------------------#
        # Get torsions associated to lbendings #
        #--------------------------------------#
        ics_ltors = []
        for lbend in ics_lbend:
            atI = lbend[0]
            atJ = lbend[1]
            atK = lbend[2]
            # atoms bonded to I
            bondedI = self.graph_neighbors(atI)
            bondedI = [at for at in bondedI if at != atJ]
            if len(bondedI) == 0: continue
            # atoms bonded to K
            bondedK = self.graph_neighbors(atK)
            bondedK = [at for at in bondedK if at != atJ]
            if len(bondedK) == 0: continue
            # Angles H-I-J
            partHI = None
            while len(bondedI) > 0:
                  atH = bondedI.pop()
                  angleHIJ = self.icvalue([atH,atI,atJ])
                  if angleHIJ < linear:
                     partHI = [atH,atI]
                     break
            if partHI is None: continue
            # Angles J-K-L
            partKL = None
            while len(bondedK) > 0:
                  atL = bondedK.pop()
                  angleJKL = self.icvalue([atJ,atK,atL])
                  if angleJKL < linear:
                     partKL = [atK,atL]
                     break
            if partKL is None: continue
            # Save data
            ltorsion = partHI+partKL
            ics_ltors.append(ltorsion)

        #-----------------------#
        # Get waggings/torsions #
        #-----------------------#
        ics_wagg = []
        ics_tors = []
        for idxA in range(len(ics_abend)):
            for idxB in range(idxA+1,len(ics_abend)):
                angleA = ics_abend[idxA]
                angleB = ics_abend[idxB]
                common = set(angleA).intersection(set(angleB))
                if len(common) == 2:
                   # Wagging:
                   if angleA[1] == angleB[1]:
                      around = sorted(list(set([angleA[0],angleA[2],angleB[0],angleB[2]])))
                      wagging = around + [angleA[1]]
                      if wagging not in ics_wagg:
                         ics_wagg.append(wagging)
                   # Torsion:
                   else:
                      cycle1 = set([angleA[0],angleA[2]]) == set([angleB[1],angleB[2]])
                      cycle2 = set([angleA[0],angleA[2]]) == set([angleB[0],angleB[1]])
                      cycle3 = set([angleB[0],angleB[2]]) == set([angleA[1],angleA[2]])
                      cycle4 = set([angleB[0],angleB[2]]) == set([angleA[0],angleA[1]])
                      cycle5 = set([angleA[0],angleA[2]]) == set([angleB[0],angleB[2]])
                      # cycle?
                      if cycle1 or cycle2 or cycle3 or cycle4 or cycle5:
                         if cycle1: torsion = [angleB[0]]+list(angleA)
                         if cycle2: torsion = [angleB[2]]+list(angleA)
                         if cycle3: torsion = [angleA[0]]+list(angleB)
                         if cycle4: torsion = [angleA[2]]+list(angleB)
                         if cycle5: torsion = [angleB[1]]+list(angleA)
                         if torsion[0] > torsion[3]: torsion = torsion[::-1]
                      # Normal torsion
                      else:
                         non_common = sorted(list((set(angleA).union(set(angleB))).difference(common)))
                         left  = non_common[0]
                         if   left  == angleA[0]: left  = [left , angleA[1]]
                         elif left  == angleA[2]: left  = [left , angleA[1]]
                         elif left  == angleB[0]: left  = [left , angleB[1]]
                         elif left  == angleB[2]: left  = [left , angleB[1]]
                         right = non_common[1]
                         if   right == angleA[0]: right = [angleA[1], right]
                         elif right == angleA[2]: right = [angleA[1], right]
                         elif right == angleB[0]: right = [angleB[1], right]
                         elif right == angleB[2]: right = [angleB[1], right]
                         torsion = left+right
                         if torsion not in ics_tors: ics_tors.append(torsion)
    
        ics_stretch = [tuple(sorted(list(ic))) for ic  in ics_stretch]
        ics_abend   = [tuple(ic) for ic  in ics_abend  ]
        ics_lbend   = [tuple(ic) for ic  in ics_lbend  ]
        ics_ltors   = [tuple(ic) for ic  in ics_ltors  ]
        ics_wagg    = [tuple(ic) for ic  in ics_wagg   ]
        ics_tors    = [tuple(ic) for ic  in ics_tors   ]

        if cleaned:
           #--------------------------------#
           # Keep only one wagging per atom #
           #--------------------------------#
           toremove = []
           newlist  = []
           for atom in range(self._natoms):
               waggs = [ic for ic in ics_wagg if ic[3] == atom]
               if waggs != []:
                  ic = random.choice(waggs)
                  newlist = newlist + [ic]
                  # Bending to remove
                  bending = (ic[0],ic[3],ic[1])
                  toremove.append(bending)
           ics_wagg = newlist
   
           #--------------------------------------#
           # Remove bendings included in waggings #
           #--------------------------------------#
           ics_abend = [ic for ic in ics_abend if ic not in toremove]

           #---------------------------#
           # Remove redundant torsions #
           #---------------------------#
           newlist  = torsions
           included = [sorted(torsion[1:3]) for torsion in newlist]
           for torsion in ics_tors:
               bond = sorted(torsion[1:3])
               if torsion in newlist:
                  included.append(bond)
               elif bond not in included:
                  included.append(bond)
                  newlist = [torsion] + newlist
           ics_tors = newlist

        #------------------#
        # Put all together #
        #------------------#
        nics = len(ics_stretch)+len(ics_abend)+2*len(ics_lbend)+len(ics_wagg)+len(ics_tors)
        ricoords = []
        for stretching in ics_stretch: ricoords.append( ("1",stretching) )
        for abending   in ics_abend  : ricoords.append( ("2",abending  ) )
        for lbending   in ics_lbend  : ricoords.append( ("3",lbending  ) )
        for ltorsion   in ics_ltors  : ricoords.append( ("4",ltorsion  ) )
        for wagging    in ics_wagg   : ricoords.append( ("4",wagging   ) )
        for torsion    in ics_tors   : ricoords.append( ("4",torsion   ) )

        if check:
          ccfreqs = self.get("ccfreqs")
          if ccfreqs is None or ccfreqs == []:
             print "ERROR: cc-freqs are not calculated!"
             sys.exit()
          self.calc_icfreqs(ricoords)
          icfreqs = self.get("icfreqs")
          if len(ccfreqs) != len(icfreqs):
             print "ERROR: redundant internal coordinates are not adequate!"
             sys.exit()
          else:
             for idx in range(len(ccfreqs)):
                 ccfreq = ccfreqs[idx].get("wavenum")/cons.cm
                 icfreq = icfreqs[idx].get("wavenum")/cons.cm
                 diff = abs(ccfreq-icfreq)
                 if diff > 0.5:
                    print "ERROR: cc-freqs and ic-freqs differ more than 0.5 cm^-1"
                    sys.exit()

        return ricoords, nics

    def purify_ricoords(self,ricoords,torsions=[],show=False,rbonds=False):
        '''
        '''

        def get_nics(icoords):
            nics = 0
            for kind,ic in icoords:
                if kind == "3": nics += 2
                else          : nics += 1
            return nics

        nn = get_nics(ricoords)
        if nn < self._nvib:
           print "ERROR: System requires %i internal coordinates, but only %i are given"%(self._nvib,nn)
           sys.exit()

        ccfreqs = self.get("ccfreqs")
        if ccfreqs is None or ccfreqs == []:
           print "ERROR: cc-freqs are not calculated!"
           sys.exit()
        ncc       = len(self._ccfreqs)
        if rbonds:
           nricoords = []
           toremove  = [(kind,tuple(ic)) for kind,ic in ricoords]
        else:
          nricoords = [(kind,tuple(ic)) for kind,ic in ricoords if kind=="1"]
          toremove  = [(kind,tuple(ic)) for kind,ic in ricoords if kind!="1"]

        # Anyone to keep?
        keep = []
        for torsion in torsions:
            ttuple = ("4",tuple(torsion))
            if ttuple in toremove:
               idx = toremove.index(ttuple)
               keep.append( toremove.pop(idx) )

        # Purify icoords
        while get_nics(nricoords) < ncc-len(keep):
              idx      = random.choice(range(len(toremove)))
              kind, ic = toremove.pop(idx)
        
              try:
                 nricoords2 = nricoords + toremove + keep
                 # Calc ic freqs
                 self.calc_icfreqs(nricoords2)
                 icfreqs = self.get("icfreqs")
                 remove = True
                 if len(icfreqs) != len(ccfreqs):
                    remove = False
                 else:
                    for idx in range(len(ccfreqs)):
                        ccfreq = ccfreqs[idx].get("wavenum")/cons.cm
                        icfreq = icfreqs[idx].get("wavenum")/cons.cm
                        diff = abs(ccfreq-icfreq)
                        if diff > 0.5:
                           remove = False
                           break
              except:
                 remove = False

              if not remove: nricoords = nricoords + [(kind,ic)]
              elif show:
                   if kind != "3": print "             removing IC %s..."%("-".join( [str(atom+1) for atom in ic]))
                   else          : print "             removing IC %s..."%("=".join( [str(atom+1) for atom in ic]))

        nricoords = nricoords + keep
        return nricoords, get_nics(nricoords)

    #-----------------------------------#
    # Functions related to the MEP      #
    #-----------------------------------#
    def imag_dir(self,icoords=None):
        '''
        Get the pair of atoms whose distance
        varies the most when following the
        imaginary frequency
        '''
        assert self._type == 1, "Problems in imag_dir"
        for freq in self._ccfreqs:
            if freq.isItImag():
               ievec = freq.get("evector")
               break
        
        if icoords is None:
           coordinates = [(i,j) for i in range(self._natoms) for j in range(i+1,self._natoms)]
        else:
           coordinates = [ic for ictype,ic in icoords if ictype != "3"]
        # Initial geom and new geom
        x_initial = np.array(self._xms,copy=True)
        x_final   = hf.ms2cc_x(x_initial + ievec,self._masslist,self._mu)
        newstrut  = Struct("", x_final, self._atonums, masslist=self._masslist)
        maxvar    = -float("inf")
        for coord in coordinates:
            value_initial =     self.icvalue(coord)
            value_final   = newstrut.icvalue(coord)
            diff          = (value_final-value_initial)
            if abs(diff) > maxvar:
               maxvar = abs(diff)
               if diff < 0.0: idir = (coord , "-")
               if diff > 0.0: idir = (coord , "+")
        return idir

    def v0dir_ifreq(self):
        for freq in self._ccfreqs:
            if freq.isItImag():
               self._v0 = np.array(freq.get("evector"),copy=True)
               ifreq = freq.copy()
               break
        return ifreq

    def v0dir_grad(self):
        assert self._gms is not None, "ERROR: gms is None"
        self._v0 = - self._gms / np.linalg.norm(self._gms)

    def v0dir_check(self,dir_v0):
        coord, effect = dir_v0
        # Prepare geometries
        x_initial = np.array(self._xms,copy=True)
        x_final   = hf.ms2cc_x(x_initial + self._v0,self._masslist,self._mu)
        newstrut  = Struct("", x_final, self._atonums, masslist=self._masslist)
        # Calculate differentce
        value_initial =     self.icvalue(coord)
        value_final   = newstrut.icvalue(coord)
        diff          = (value_final-value_initial)
        # Is it correct?
        if diff > 0.0 and effect == "-": return False
        if diff < 0.0 and effect == "+": return False
        return True

    def v0dir_invert(self):
        if (self._v0 is not None): self._v0 = - self._v0
        if (self._v1 is not None): self._v1 = - self._v1

    def v1dir_hess(self):
        v0Fv0 = float( np.matrix(self._v0) * self._Fms * np.matrix(self._v0).transpose() )
        component_A = np.array(self._Fms * np.matrix(self._v0).transpose()).transpose()[0]
        component_B = v0Fv0*self._v0
        self._v1 = (component_A - component_B) / np.linalg.norm(self._gms)

    def nextTaylor(self,ds,bw=False,qt=False):
        '''
        * Uses the quadratic Taylor expansion of the path:
                   x(s0+ds) = x(s0) + v0*ds + 0.5 *v1 * ds^2

        * If bw is True, it uses -v^(0) and -v^(1)

        * If qt is False, the quadratic term is ommitted:
                   x(s0+ds) = x(s0) + v0*ds
        '''
        ds = abs(ds)

        v0 = self.get("v0")
        v1 = self.get("v1")

        if bw:        v0 = -v0
        if bw and qt: v1 = -v1

        x_next = self._xms + v0 * ds
        if qt: x_next += 0.5 * v1 * (ds**2)

        # Get it also in non-scaled coordinates
        x_next_cc = hf.ms2cc_x(x_next,self._masslist,self._mu)

        return x_next_cc, x_next


    def write_molden(self,moldenfile):

        molden = open(moldenfile,'w')
        molden.write("[Molden Format]\n")
        molden.write("[FR-COORD] # Coordinates in bohr\n")
        for idx in range(self._natoms):
            symbol = self._symbols[idx]
            x, y, z = self._xcc[3*idx:3*idx+3]
            molden.write(" %2s  %+11.6f  %+11.6f  %+11.6f \n"%(symbol,x,y,z))

        molden.write("[FREQ] # Frequencies in cm^-1\n")
        evecs = []
        for vibfreq in self._ccfreqs:
            f, v = vibfreq.get("wavenum"), vibfreq.get("evector")
            if f.imag != 0.0: f = - f.imag
            f = f/cons.cm
            molden.write(" %9.4f\n"%f)
            # As they are displacements, ms2cc_x has to be used
            if v is not None:
               evec = hf.ms2cc_x(v,self._masslist,self._mu)
               evec = evec / np.linalg.norm(evec)
               evecs.append( evec )

        molden.write("[FR-NORM-COORD] # Displacements in bohr\n")
        nv = 1
        for evector in evecs:
            molden.write("vibration  %i\n"%nv )
            for j in range(0,len(evector),3):
                vx, vy, vz = evector[j:j+3]
                molden.write("   %+9.3f  %+9.3f  %+9.3f\n"%(vx,vy,vz))
            nv += 1

        molden.close()


#--------------------------------------#
# Some functions associated to classes #
#--------------------------------------#
def basic2Struct(name,xcc,atonums,ch,mtp,Etot,gcc,Fcc,stype=-1,masslist=None):
    variables = ["ch","mtp","Etot"]
    values    = [ ch , mtp , Etot ]
    if gcc is not None:  variables.append("gcc"); values.append(gcc)
    if Fcc is not None:  variables.append("Fcc"); values.append(Fcc)
    structure = Struct(name,xcc,atonums,masslist=masslist,stype=stype)
    structure.set(variables,values)
    return structure
#--------------------------------------#
def xyz2Struct(xyzfile,name="",stype=0):
    xvector, symbols, masslist = hf.read_xyz(xyzfile)
    xvector = xvector / cons.angstrom
    structure = Struct(name,xvector,symbols,masslist=masslist,stype=stype)
    return structure, masslist
#--------------------------------------#
def gts2Struct(gtsfile,name="",masslist=None,stype=0):
    xyz_list , atonum_list , ch, mtp, Etot, pgroup, rotsigma, gcc , Fcc, freqs = read_gtsfile(gtsfile)
    variables = ["ch","mtp","Etot"]
    values    = [ ch , mtp , Etot ]
    if masslist is None: variables += ["pgroup","rotsigma"]; values += [pgroup , rotsigma ]
    if gcc is not None:  variables.append("gcc"); values.append(gcc)
    if Fcc is not None:  variables.append("Fcc"); values.append(Fcc)
    structure = Struct(name,xyz_list,atonum_list,masslist=masslist,stype=stype)
    structure.set(variables,values)
    # No hessian but freqs
    if Fcc is None and freqs is not None:
       ccfreqs = []
       for freq in freqs:
           instance = Freq()
           instance.set_wavenum(freq * cons.cm)
           instance.calc_derivated_magnitudes()
           ccfreqs.append(instance)
       structure._ccfreqs = ccfreqs
    # Return structure
    return structure
#--------------------------------------#


def test():
    structure = gts2Struct("XTR_gts/S19_000_000.gts")
    structure.basic_setups()
    structure.graph_autoconnect()

    torsion1 = (13-1,12-1, 1-1, 2-1)
    torsion2 = (16-1,13-1,12-1, 1-1)
    utorsions = [torsion1,torsion2]
    #utorsions = []
    r_stretchings , r_bendings , r_waggings , r_torsions  = structure.graph_ricoords()
    nr_stretchings, nr_bendings, nr_waggings, nr_torsions = structure.graph_nricoords(utorsions)

    for icoord in nr_stretchings: print icoord
    print "--"
    for icoord in nr_bendings   : print icoord
    print "--"
    for icoord in nr_waggings   : print icoord
    print "--"
    for icoord in nr_torsions   : print icoord

    structure.calc_ccfreqs()
    ccfreqs = structure.get("ccfreqs")

    structure.calc_icfreqs(r_stretchings, r_bendings, [], r_waggings+r_torsions)
    ricfreqs = structure.get("icfreqs")

    structure.calc_icfreqs(nr_stretchings, nr_bendings, [], nr_waggings+nr_torsions)
    nricfreqs = structure.get("icfreqs")

    print "----"
    for idx in range(len(ccfreqs)):
        print ccfreqs[idx], ricfreqs[idx], nricfreqs[idx]
    print "----"

    print structure.graph_atostretchings(11)
    print structure.graph_atobendings(11)
    print structure.graph_atowaggings(11)
    print structure.graph_nolinealangle(0,11)
    print structure.graph_nolinealangle(11,0)
    print structure.graph_torsion((0,11))

#test()

