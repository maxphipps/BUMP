""" 
NEB package: class storage
"""

import os

class XYZ:
  '''XYZ data class'''

  def __init__(self, atdat=[], xyzdat=[], coord=0):
    self.at = atdat
    self.xyz = xyzdat
    # set meta data
    # number of atoms
    self.nat = len(self.xyz[0])
    # replica energy pathway coordinate
    self.coord_val = coord
    self.coord_str = '{:3.2f}'.format( coord )
    # construct filename for this replica
    self.datdir =      'im_'+self.coord_str
    self.filenamexyz = 'im_'+self.coord_str+'_init.xyz'
    self.filenamedat = self.datdir+\
                       '/im_'+self.coord_str+'.dat'
    self.filenameout = self.datdir+\
                       '/onetep.out'

    # make directory for running onetep files
    self.makedatdir()

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
  def fromdata(self, atdat, xyzdat, coord):
    # from xyz data:
    return XYZ(atdat, xyzdat, coord)

  def applyforces(self,forces):
    for ii in range(self.nat):
      #print self.xyz[ii][0]
      self.xyz[-1][ii][0] += forces[ii][0]
      self.xyz[-1][ii][1] += forces[ii][1]
      self.xyz[-1][ii][2] += forces[ii][2]
      #print self.xyz[ii][0]

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
    f.write(str(glob_nat)+'\n\n')
    # xyz data:
    for ii in range(self.nat):
      f.write( '{}  {:4f}  {:4f}  {:4f}\n'.format( self.at[ii], self.xyz[ixyz][ii][0], self.xyz[ixyz][ii][1], self.xyz[ixyz][ii][2] ) )
    f.close()

