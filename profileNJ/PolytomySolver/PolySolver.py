"""
Author: Manuel Lafond
Date: 07/2015
Functions for the the dynamic programming scheme of polytomy resolution 
"""
import os
import sys

from pprint import pprint
import copy

from TreeLib import *

class PolytomySolver:
  
  """
  Usage example : 
  lcamap = TreeUtils.lcaMapping(genetree, specietree, multspeciename=False)
            
  ps = PolytomySolver(genetree, specietree, lcamap)
  ps.computeCostsTable()
  
  print "COST=", ps.getTableValue(specietree, 1)
  
  r = ps.getResolution()
  print r
  """
  
  
  def __init__(self, polytomy, speciestree, lcaMapping):
    self.polytomy = polytomy
    self.speciestree = speciestree
    self.lcaMapping = lcaMapping
    self.multiplicities = {}
    self.cup_values = {}
  
  def computeMultiplicities(self):
    
    """ Compute multiplicities for each species s, i.e. the number of children x of self.polytomy such that lcaMapping[x] = s
    """
    
    self.multiplicities = {}
    for g in self.polytomy.children:
      s = self.lcaMapping[g]
      
      if not s in self.multiplicities:
        self.multiplicities[s] = 0
      
      self.multiplicities[s] += 1
        

  def computeCostsTable(self):
      """ Compute costs table, in time O(S)
      """
      
      #The stuff in comments is the classical table-filling dynamic programming algorithm - left here for testing purposes
      #dpsize = 20
      #dp_values = {}
      
      self.computeMultiplicities()
      
      self.cup_values = {}   #for a species s, cup_values[s] = ( bottom plateau value,  breakpt left, breakpt right )
      
      for s in self.speciestree.traverse("postorder"):
        
        mult = 0
        if s in self.multiplicities:
          mult = self.multiplicities[s]
        
        if(s.is_leaf()):
          
            if mult == 0:
               self.cup_values[s] = ( 1, 1, 1 )
            else:
               self.cup_values[s] = ( 0, mult, mult )
            
            """
            dp_values[s] = dpsize * [0]
            for k in range(0, dpsize):
              dp_values[s][k] = abs(k - mult)
            print self.cup_values[s]
            print dp_values[s]
            """
        
        else:
            
            cl = self.cup_values[s.children[0]]
            cr = self.cup_values[s.children[1]]
            
            lmin = cl[0]
            l1 = cl[1]
            l2 = cl[2]
            
            rmin = cr[0]
            r1 = cr[1]
            r2 = cr[2]
            
            #here we go, magic
            if l1 < r1 and l2 < r1:
              smin = lmin + rmin + r1 - l2
              b1 = l2
              b2 = r1
            elif l1 < r1 and r1 <= l2 and l2 <= r2:
              smin = lmin + rmin
              b1 = r1
              b2 = l2
            elif l1 < r1 and l2 > r2:
              smin = lmin + rmin
              b1 = r1
              b2 = r2
            elif r1 <= l1 and l1 <= r2 and r1 <= l2 and l2 <= r2:
              smin = lmin + rmin
              b1 = l1
              b2 = l2
            elif r1 <= l1 and l1 <= r2 and l2 > r2:
              smin = lmin + rmin
              b1 = l1
              b2 = r2
            elif l1 > r2 and l2 > r2:
              smin = lmin + rmin + l1 - r2
              b1 = r2
              b2 = l1
            else:
              print "CASE NOT COVERED"
              
            b1 += mult  #todo : verify this
            b2 += mult
              
            self.cup_values[s] = ( smin, b1, b2 )
            
            """
            dp_values[s] = dpsize * [0]
            bb1 = 9999
            bb2 = 9999
            bmin = 9999
            s1 = s.children[0]
            s2 = s.children[1]
            for k in range(0, dpsize):
              dp_values[s][k] = dp_values[s1][k] + dp_values[s2][k]
              if dp_values[s][k] < bmin:
                bmin = dp_values[s][k]

            
            for k in range(0, dpsize):
              if bb1 == 9999:
                if dp_values[s][k] == bmin:
                  bb1 = k
              elif bb1 != 9999 and bb2 == 9999:
                if dp_values[s][k] != bmin:
                  bb2 = k - 1
            
            for k in range(0, dpsize):
              if k < bb1:
                dp_values[s][k] = bmin + (bb1 - k)
              elif k > bb2:
                dp_values[s][k] = bmin + k - bb2
            
            #print self.cup_values[s]
            #print dp_values[s]
            """
            
      #print self.cup_values[self.speciestree]     
      #print dp_values[self.speciestree]

      
  def getTableValue(self, s, k):
    """ Returns the equivalent of table[s, k], assuming that computeCostsTable has been called previously
    """
    
    c = self.cup_values[s]   #format is a list with c = (min, b1, b2)
    
    if k < c[1]:
      return c[0] + c[1] - k
    elif k > c[2]:
      return c[0] + k - c[2]
    else:
      return c[0]


  def getResolution(self):
    """ Returns the first resolution found in Newick format, assuming that computeCostsTable has been called previously
    """
    return self.getResolutions(self.speciestree, 1)[0]

  def getResolutions(self, s, k):
    """ Returns one set of k resolutions rooted at s
    """
        
    v = self.getTableValue(s, k)
    vright = self.getTableValue(s, k + 1)
    vleft = self.getTableValue(s, k - 1)
    
    
    #------------------------------------------
    # the leaf case
    #------------------------------------------
    if s.is_leaf():
      
      if v == 0:
        return [s.name] * k
        
      elif v == vright + 1:
        
        resz = self.getResolutions(s, k + 1)
        r1 = '(' + resz[0] + ',' + resz[1] + ')'
        del resz[0]
        resz[0] = r1
        return resz
        
      else:       #v == vleft + 1
        resz = self.getResolutions(s, k - 1)
        resz.append(s.name + '_LOSS')
        return resz
    #------------------------------------------
    # the non-leaf case
    #------------------------------------------
    else:  
      
      mult = 0
      if s in self.multiplicities:
        mult = self.multiplicities[s]
      
      s1 = s.children[0]
      s2 = s.children[1]
      vup1 = self.getTableValue(s1, k - mult)
      vup2 = self.getTableValue(s2, k - mult)
      
      if k - mult > 0 and v == vup1 + vup2:     #speciation
      
        resz_s1 = self.getResolutions(s1, k - mult)
        resz_s2 = self.getResolutions(s2, k - mult)
        
        resz = []
        #we just make the obvious joins 
        for i in range(0, k - mult):
          if i > k - mult - 1:
            resz.append('LOSS')
          else:
            res = '(' + resz_s1[i] + ',' + resz_s2[i] + ')'
            resz.append(res)
        return resz
        
      #TODO : this is actually the same as above (the leaf case)  
      elif v == vright + 1:
      
        resz = self.getResolutions(s, k + 1)
        r1 = '(' + resz[0] + ',' + resz[1] + ')'
        del resz[0]
        resz[0] = r1
        return resz
      else:  #hopefully if we here, v == vleft + 1 (otherwise there's something I don't understand)
      
        #print "k=", k, " v=",v, " vup1=", vup1, " vup2=", vup2, " vright=", vright, " vleft=",vleft
      
        resz = self.getResolutions(s, k - 1)
        resz.append('LOSS')
        return resz