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
\n\
pbc_correction_cutoff 3.70426 bohr ! dimensionless\n\
charge 0\n\
spin_polarized : F\n\
\n\
maxit_ngwf_cg 25\n\
\n\
%block lattice_cart\n\
ang\n"


datlines2="%endblock lattice_cart\n\
\n\
\n\
%block positions_abs\n\
ang\n"

datlines3="%endblock positions_abs\n\
\n\
\n\
%block species\n\
C   C   6 4 3.70426\n\
H   H   1 1 3.70426\n\
O   O   8 4 3.70426\n\
%endblock species\n\
\n\
\n\
%block species_pot\n\
C   '../hpc_files/c-optgga1.recpot'\n\
H   '../hpc_files/h-optgga1.recpot'\n\
O   '../hpc_files/o-optgga1.recpot'\n\
%endblock species_pot\n\
"


# TODO:
# TODO: This definition is duplicated from neb.py
# TODO:
def getXYZdata(filename):
  """ returns filename's xyz data as atdat and xyzdat arrays.
  (returns the final xyz data if more than one xyz structure is present 
  in the file.) """

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

  return [atdat,xyzdat]

#def getXYZdata(filename):
#  """ Returns filename's XYZ data as atdat and xyzdat arrays """
#  f = open(filename, 'r')
#  # Read number of lines
#  nat = int(f.readline())
#  comment = f.readline()
#  # Init xyz data
#  xyzdat = [[0.0,0.0,0.0]]*nat
#  atdat = ['']*nat
#  # Read xyz data
#  for ii in range(nat):
#    [at,x,y,z] = f.readline().split()
#    xyzdat[ii] = [float(x),float(y),float(z)]
#    atdat[ii] = at
#  f.close()
#  return [atdat,xyzdat]

def fromXYZ(fnameIn,fnameOut):
  xyzdata = getXYZdata(fnameIn)
  at = xyzdata[0]
  xyz = xyzdata[1]

  main(at,xyz,fnameOut)


def fromdata(at,xyz,fnameOut):
  main(at,xyz,fnameOut)


def main(at,xyz,fnameOut):

  nat = len(at)

  # Read ONETEP cell reference data
  fInCell = open('cell_ref.dat','r')
  cellline = ['']*3
  null = fInCell.readline()
  cellline[0] = fInCell.readline()
  cellline[1] = fInCell.readline()
  cellline[2] = fInCell.readline()
  fInCell.close()
  
  
  #################################################################
  
  # Writeout the dat file
  fOut = open(fnameOut,'w')
  
  # Dump template lines
  fOut.write(datlines1)
  
  # Write cell lines
  fOut.write(cellline[0])
  fOut.write(cellline[1])
  fOut.write(cellline[2])
  
  # Dump template lines
  fOut.write(datlines2)
  
  # Write geometry lines
  for ii in range(nat):
    fOut.write( '{}   {:4f}   {:4f}   {:4f}\n'.format( at[ii], xyz[-1][ii][0], xyz[-1][ii][1], xyz[-1][ii][2] ) )
  
  # Dump template lines
  fOut.write(datlines3)
  
  fOut.close()
  
  #################################################################

if __name__ == "__main__":
  fnameIn = sys.argv[1]
  fnameOut = sys.argv[2]
  fromXYZ(fnameIn,fnameOut)


