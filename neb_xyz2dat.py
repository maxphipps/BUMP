#!/usr/bin/env python

"""xyz2dat.py: Converts xyz geometry to ONETEP dat file. """

import sys
from math import floor

datlines1="threads_max : 1 \n\
maxit_palser_mano : 0 \n\
write_density_plot : F\n\
kernel_cutoff 1000 bohr\n\
\n\
write_forces : T\n\
\n\
cutoff_energy       : 600 eV\n\
ngwf_threshold_orig : 0.000002\n\
k_zero              : 3.5\n\
write_xyz true\n\
\n\
! elec_cg_max 0\n\
occ_mix 1.0\n\
\n\
minit_lnv 10 \n\
maxit_lnv 15 \n\
\n\
write_denskern F\n\
write_tightbox_ngwfs F\n\
\n\
output_detail BRIEF\n\
\n\
xc_functional PBE\n\
TASK SINGLEPOINT\n\
! dispersion 1\n\
\n"

datlines1B="coulomb_cutoff_type SPHERE\n\
charge 0\n\
spin_polarized : F\n\
\n\
maxit_ngwf_cg 25\n\
\n\
%block lattice_cart\n\
ang\n"

#pbc_correction_cutoff 3.70426 bohr ! dimensionless\n\

datlines2="%endblock lattice_cart\n\
\n\
\n\
%block positions_abs\n\
ang\n"

datlines3="%endblock positions_abs\n\
\n\
\n\
%block species\n\
H   H   1 1 3.70426\n\
O   O   8 4 3.70426\n\
%endblock species\n\
\n\
\n\
%block species_pot\n\
H   '../hpc_files/h-optgga1.recpot'\n\
O   '../hpc_files/o-optgga1.recpot'\n\
%endblock species_pot\n\
"
#C   C   6 4 3.70426\n\
#C   '../hpc_files/c-optgga1.recpot'\n\

BOHR2ANG = 0.529
NGWF_RADII = 7. * BOHR2ANG


from neb_classes import XYZ
#class XYZ:
#  '''XYZ data class: struct'''
#
#  def __init__(self):
#    self.at = None
#    self.xyz = None
#    self.nat = None


def mean(arr):
  return sum(arr)/len(arr)

def secondLargest(arr):
  # Returns the second largest element of an array
  return sorted(arr, reverse=True)[1]

# TODO:
# TODO: This definition is duplicated from neb.py
# TODO:
def getXYZdata(filename):
  global nat
  ''' returns filename's xyz data as atdat and xyzdat arrays.
  (returns the final xyz data if more than one xyz structure is present 
  in the file.) '''

  # read lines to list
  with open(filename) as f:
    lines = list(f)

  # find number of XYZ structure entries in the file
  nat = int(lines[0])
  nentries = floor( float(len(lines))/float(nat+2) )
  print 'Reading in structure entry number ',str(int(nentries))

  # calculate line number the final xyz structure entry
  # begins at
  linestart = int((nentries-1)*(nat+2))

  # init xyz data
  xyzdat = [[[0.0,0.0,0.0]]*nat]
  atdat = ['']*nat

  # read xyz data
  ii = 0
  for line in lines[linestart+2:linestart+nat+2]:
    [at,x,y,z] = line.split()
    xyzdat[-1][ii] = [float(x),float(y),float(z)]
    atdat[ii] = at
    ii += 1

  return [atdat,xyzdat,nat]


def centerXYZ(dat,cellParams):
  # Translates the xyz to the cell center
  # Calc. current center of xyz
  cent_xyz = [0.0]*3
  cent_xyz[0] = mean( [dat.xyz[-1][ii][0] for ii in range(dat.nat)] )
  cent_xyz[1] = mean( [dat.xyz[-1][ii][1] for ii in range(dat.nat)] )
  cent_xyz[2] = mean( [dat.xyz[-1][ii][2] for ii in range(dat.nat)] )

  # Calc. cell center
  cent_cell = [0.0]*3
  cent_cell[0] = cellParams[0]/2.
  cent_cell[1] = cellParams[1]/2.
  cent_cell[2] = cellParams[2]/2.

  # Apply translation vector
  for ii in range(dat.nat):
    dat.xyz[-1][ii][0] += (cent_cell[0] - cent_xyz[0])
    dat.xyz[-1][ii][1] += (cent_cell[1] - cent_xyz[1])
    dat.xyz[-1][ii][2] += (cent_cell[2] - cent_xyz[2])

  return dat.xyz

def getCellParams(dat,molSize):
  # Calculates the cell parameters
  # Calculate maximum distance between atoms (brute force)
  Rmax = 0.0
  for ii in range(dat.nat):
    for ij in range(ii+1,dat.nat):
      Rtest = ( \
        ( dat.xyz[-1][ii][0] - dat.xyz[-1][ij][0])**2. + \
        ( dat.xyz[-1][ii][1] - dat.xyz[-1][ij][1])**2. + \
        ( dat.xyz[-1][ii][2] - dat.xyz[-1][ij][2])**2. )**0.5
      if (Rtest > Rmax): Rmax = Rtest

  margin = Rmax + NGWF_RADII*2.
  xmax = 2.*margin + molSize[0]
  ymax = 2.*margin + molSize[1]
  zmax = 2.*margin + molSize[2]

  return [xmax,ymax,zmax]

def calcMolSize(dat):
  xSize = (max( [dat.xyz[-1][ii][0] for ii in range(dat.nat)] ) \
     - min( [dat.xyz[-1][ii][0] for ii in range(dat.nat)] ))
  ySize = (max( [dat.xyz[-1][ii][1] for ii in range(dat.nat)] ) \
     - min( [dat.xyz[-1][ii][1] for ii in range(dat.nat)] ))
  zSize = (max( [dat.xyz[-1][ii][2] for ii in range(dat.nat)] ) \
     - min( [dat.xyz[-1][ii][2] for ii in range(dat.nat)] ))
  return [xSize,ySize,zSize]

def fromXYZ(fnameIn,fnameOut):
  # Read xyz file
  # Parse to XYZ object
  # Use fromdata() method with XYZ object to create dat file
  xyzobj = getXYZdata(fnameIn)

  fromdata(xyzobj,fnameOut)

def fromdata(xyzobj,fnameOut):
  # Use XYZ data to create dat file

  #print xyzobj.xyz
  molSize = calcMolSize(xyzobj)

  #centerXYZ(xyzobj)

  # Read cellline from cell_ref.dat file
  #cellline = read_cellref()
  # Infer cellline from xyz data
  cellParams = getCellParams(xyzobj,molSize)

  xmax = cellParams[0]
  ymax = cellParams[1]
  zmax = cellParams[2]

  cellline = [str(xmax)+'  0.00  0.00'+'\n', \
    '0.00  '+str(ymax)+'  0.00'+'\n', \
    '0.00  0.00  '+str(zmax)+'\n']
  
  # Center the xyz in the cell
  xyz = centerXYZ(xyzobj,cellParams)

  # Calculate coulomb_cutoff_radius
  # = the maximum diameter of the molecule,
  # plus twice the radius of the NGWFs, 
  # plus a little bit extra for good luck (~2Bohr)
  # TODO: This is approximated as the maximum cuboid container for the molecule
  max1 = max(cellParams)
  max2 = secondLargest(cellParams)
  coulomb_cutoff_radius = (max1*max1 + max2*max2)**0.5 + \
    + NGWF_RADII*2. \
    + BOHR2ANG*2.

  main(xyzobj,cellline,coulomb_cutoff_radius,fnameOut)


#def read_cellref():
#  # Read ONETEP cell reference data
#  fInCell = open('cell_ref.dat','r')
#  cellline = ['']*3
#  null = fInCell.readline()
#  cellline[0] = fInCell.readline()
#  cellline[1] = fInCell.readline()
#  cellline[2] = fInCell.readline()
#  fInCell.close()
#
#  return cellline
  
  

def main(xyzobj,cellline,coulomb_cutoff_radius,fnameOut):

  #################################################################

  # Writeout the dat file
  fOut = open(fnameOut,'w')
  
  # Dump template lines
  fOut.write(datlines1)

  # Write cell lines
  fOut.write("coulomb_cutoff_radius "+str(coulomb_cutoff_radius)+" ang\n")

  # Dump template lines
  fOut.write(datlines1B)
  
  # Write cell lines
  fOut.write(cellline[0])
  fOut.write(cellline[1])
  fOut.write(cellline[2])
  
  # Dump template lines
  fOut.write(datlines2)
  
  # Write geometry lines
  for ii in range(xyzobj.nat):
    fOut.write( '{}   {:4f}   {:4f}   {:4f}\n'.format( xyzobj.at[ii], xyzobj.xyz[-1][ii][0], xyzobj.xyz[-1][ii][1], xyzobj.xyz[-1][ii][2] ) )
  
  # Dump template lines
  fOut.write(datlines3)
  
  fOut.close()
  
  #################################################################

if __name__ == "__main__":
  fnameIn = sys.argv[1]
  fnameOut = sys.argv[2]
  fromXYZ(fnameIn,fnameOut)


