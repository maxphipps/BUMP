""" 
NEB package: class storage
"""

import os
from math import floor

VERBOSE = False

def mean(arr):
  return sum(arr)/len(arr)

def printv(string):
  # verbose print
  if (VERBOSE): print(string)

class XYZ:
  '''XYZ data class'''

  def __init__(self, atdat=[], xyzdat=[], coord=0, localfile=False, readVelo=True):
    self.at = atdat
    self.xyz = xyzdat
    # set meta data
    # number of atoms
    self.nat = len(self.xyz[0])
    # replica energy pathway coordinate
    self.coord_val = coord
    self.coord_str = '{:3.2f}'.format( coord )

    # If not a 'local' xyz, then create directory metadata
    if (localfile):
      self.datdir =      ''
    else:
      self.datdir =      'im_'+self.coord_str

      # make directory for running onetep files
      self.makedatdir()

    # construct filename for this replica
    self.filenamexyz = 'im_'+self.coord_str+'_init.xyz'
    self.filenamedat = self.datdir+\
                       '/im_'+self.coord_str+'.dat'
    self.filenamevelo = self.datdir+\
                       '/im_'+self.coord_str+'.velo.xyz'
    self.filenameout = self.datdir+\
                       '/onetep.out'

    # xyz velocity data, used for velocity verlet algorithm
    if (readVelo):
      # read in velocity data from file
      self.readVelo(self.filenamevelo, self.coord_val)
    else:
      # default: zeroes
      self.xyzvelo = [[[0.0 for i in range(3)] for j in range(self.nat)]]

  def makedatdir(self):
    try:
      os.stat(self.datdir)
    except:
      os.mkdir(self.datdir) 

  @classmethod
  def float2str( self, coord ):
    # Parses coord float to string
    coordstr = '{:3.2f}'.format( coord )
    return 'im_'+coordstr

  @classmethod
  #def fromdata(self, atdat, xyzdat, coord):
  def fromsinglexyzdata(self, atdat, xyzdat, readVelo):
    # from single xyz frame data:
    coord = -1
    tmpXYZ = XYZ(atdat, xyzdat, coord, localfile=True, readVelo=readVelo)
    #tmpXYZ.datdir =     '' 
    return tmpXYZ

  def applypositionchange(self,deltapos):
    for ii in range(self.nat):
      #print self.xyz[ii][0]
      self.xyz[-1][ii][0] += deltapos[ii][0]
      self.xyz[-1][ii][1] += deltapos[ii][1]
      self.xyz[-1][ii][2] += deltapos[ii][2]
      #print self.xyz[-1][ii][0]

  @classmethod
  def readXYZdata(self,filename):
    ''' returns filename's xyz data
    (returns the final xyz data if more than one xyz structure is present 
    in the file.) '''
  
    # read lines to list
    with open(filename) as f:
      lines = list(f)
  
    # find number of XYZ structure entries in the file
    nat = int(lines[0])
    nentries = int(floor( float(len(lines))/float(nat+2) ))
    printv(''.join([filename, '-> Reading in structures up to final entry number ',str(nentries)]))
  
    # init xyz data
    xyzdat = [[None]*nat for i in xrange(nentries)]
    atdat = ['']*nat
  
    for entry in range(0,nentries):
  
      # calculate line number this xyz structure entry
      # begins at
      linestart = int((entry)*(nat+2))
       
      # read xyz data
      ii = 0
      for line in lines[linestart+2:linestart+nat+2]:
        [at,x,y,z] = line.split()
        xyzdat[entry][ii] = [float(x),float(y),float(z)]
        atdat[ii] = at
        ii += 1
  
    return [atdat,xyzdat,nat]
    #return XYZ(atdat,xyzdat,coord)

  @classmethod
  def fromXYZ(self,filename,coord,readVelo=False):
    [atdat,xyzdat,nat] = XYZ.readXYZdata(filename)
    return XYZ(atdat,xyzdat,coord,readVelo=readVelo)

  def readVelo(self,filename,coord):
    # Reads in xyz velocity data
    [atdat,xyzvelodat,nat] = XYZ.readXYZdata(self.filenamevelo)
    self.xyzvelo = xyzvelodat

  def writeXYZ(self, ixyz=-1, appendfile=False, filename=None):
    print 'writing replica at pathway coordinate ',str(self.coord_str), \
        ' to file ', self.filenamexyz
    # Get filename from argument, else use default filename 
    filename = filename or self.filenamexyz
    if (appendfile):
      f = open(filename, 'a')
    else:
      f = open(filename, 'w')
    # header: number of atoms and comment line:
    f.write(str(self.nat)+'\n\n')
    # xyz data:
    for ii in range(self.nat):
      f.write( '{}  {:4f}  {:4f}  {:4f}\n'.format( self.at[ii], self.xyz[ixyz][ii][0], self.xyz[ixyz][ii][1], self.xyz[ixyz][ii][2] ) )
    f.close()

  def writeVelo(self, ixyz=-1, appendfile=False):
    if (appendfile):
      f = open(self.filenamevelo, 'a')
    else:
      f = open(self.filenamevelo, 'w')
    # header: number of atoms and comment line:
    f.write(str(self.nat)+'\nCART_VELO\n')
    # xyz data:
    for ii in range(self.nat):
      f.write( '{}  {:4f}  {:4f}  {:4f}\n'.format( self.at[ii], self.xyzvelo[ixyz][ii][0], self.xyzvelo[ixyz][ii][1], self.xyzvelo[ixyz][ii][2] ) )
    f.close()

  def centerXYZ(self,cellParams):
    # Translates the (last) xyz to the cell center
    # Calc. current center of xyz
    cent_xyz = [0.0]*3
    cent_xyz[0] = mean( [self.xyz[-1][ii][0] for ii in range(self.nat)] )
    cent_xyz[1] = mean( [self.xyz[-1][ii][1] for ii in range(self.nat)] )
    cent_xyz[2] = mean( [self.xyz[-1][ii][2] for ii in range(self.nat)] )
  
    # Calc. cell center
    cent_cell = [0.0]*3
    cent_cell[0] = cellParams[0]/2.
    cent_cell[1] = cellParams[1]/2.
    cent_cell[2] = cellParams[2]/2.
  
    # Apply translation vector
    for ii in range(self.nat):
      self.xyz[-1][ii][0] += (cent_cell[0] - cent_xyz[0])
      self.xyz[-1][ii][1] += (cent_cell[1] - cent_xyz[1])
      self.xyz[-1][ii][2] += (cent_cell[2] - cent_xyz[2])
  
    return self

  def getboundbox(self):
    # Returns the box that bounds the (last) xyz system
    xmin = +1E8; ymin = +1E8; zmin = +1E8
    xmax = -1E8; ymax = -1E8; zmax = -1E8
    for ii in range(self.nat):
      xmin = min(xmin,self.xyz[-1][ii][0])
      ymin = min(ymin,self.xyz[-1][ii][1])
      zmin = min(zmin,self.xyz[-1][ii][2])
      xmax = max(xmax,self.xyz[-1][ii][0])
      ymax = max(ymax,self.xyz[-1][ii][1])
      zmax = max(zmax,self.xyz[-1][ii][2])
    return [[xmin, xmax], [ymin, ymax], [zmin, zmax]]

  def getboundbox_lengths(self):
    # Returns the box that bounds the (last) xyz system
    xmin = +1E8; ymin = +1E8; zmin = +1E8
    xmax = -1E8; ymax = -1E8; zmax = -1E8
    for ii in range(self.nat):
      xmin = min(xmin,self.xyz[-1][ii][0])
      ymin = min(ymin,self.xyz[-1][ii][1])
      zmin = min(zmin,self.xyz[-1][ii][2])
      xmax = max(xmax,self.xyz[-1][ii][0])
      ymax = max(ymax,self.xyz[-1][ii][1])
      zmax = max(zmax,self.xyz[-1][ii][2])
    return [abs(xmin - xmax), abs(ymin - ymax), abs(zmin - zmax)]
    
