#!/usr/bin/env python

""" neb.py: A nudged elastic band interface for ONETEP v4.5.10.4.

---OVERVIEW---
Implementation of NEB without kink-fixing as per:
http://theory.cm.utexas.edu/henkelman/pubs/jonsson98_385.pdf

Rpath with constraint:
For a recent review of this method compared with other chain-of-states methods,
see:
'comparison of three chain-of-states methods: nudged elastic band
and replica path with restraints or constraints', J. Chem. Theory Comput.
2012, 8, 5035-5051.

A brief overview of the code's usage follows, and assumes the user to be working
on the intended HPC platform that the calculations will be ran on:
1. End point geometries are optimised, aligned and shifted to the cell centre.
2. The following files are then created/updated by the user:
   - Optimised reactant xyz renamed and located to ./R.xyz .
   - Optimised product xyz renamed and located to ./P.xyz .
   - ./cell_ref.dat updated to set the cell parameters.
   - ./neb_xyz2dat.py updated to set the XC, energy cutoff and pseudopotential file 
     locations (default is within ./hpc_files/).
   - Onetep binary and appropriate recpot files placed in ./hpc_files/ .
   - Total number of images optionally set.  Currently this is set using the
     glob_nim variable (default = 6).  This value includes the reactant and product images, 
     and so should be adjusted accordingly.  Please note that modifying this should be performed 
     with caution, and may require updating of the ./sub_neb_forces.sh file.
3. neb.py is ran with the -init flag, providing the optimised reactant and product .xyz files,
   e.g. python2.7 neb.py -init R.xyz P.xyz .
   This will result in the creation of a directory structure containing multiple ONETEP input files,
   and a series of (linearly interpolated) images that represent the transit from reactants to
   products.
4. run_neb.sh is ran, which iterates (default number of steps = 3) the optimisation.
   Specifically, this results in the following jobs being submitted:
   (A) A batch of (parallel) singlepoint and force energy calculations for each of the images.
   (B) Once the above batch is fully complete, a (serial) call to the neb.py python script is made
     to update the geometries in the search direction of the minimum energy pathway. 
     Specifically, forces are obtained from the ONETEP singlepoint calculations, 
     and a steepest descent geometry step taken.
   (C) The batch (A) is resubmitted once the neb.py task has completed.
5. Step 4 is repeated until convergence of the NEB minimum energy pathway. (No optimisation 
   convergence criteria has so far been implemented - this ranks high on my TODO list!)

The status of the execution can be confirmed using the -plot and -xyz flags, which are further detailed 
below in the flags section. It is recommended that calls be be made after completion of a complete 
run_neb.sh call (step 4), as issues can arise when calculations are half complete.  

An example of a minimum energy pathway produced by calling neb.py with the -plot flag can be seen in the
./mep.png file included.


--- EXAMPLE USAGE ---
  pythonw neb.py -init react_shift.xyz prod_shift.xyz
  pythonw neb.py -iter
  pythonw neb.py -plot
  pythonw neb.py -xyz


--- FLAGS ---
  [ -init ] Initialise replicas by linear interpolation from reactants (.xyz) to 
            products (.xyz)
  [ -iter ] Iterates to optimise geometry by calculating nudged forces from 
            singlepoint calculation and applying these to the xyz and dat files
  [ -plot ] Plots the energy pathway
  [ -xyz  ]  Extracts the most recent path -> path.xyz


--- NOTES ---
This driver has been implemented as an external python script for a number of reasons, namely:
  - To avoid British Crown and Accelrys/Biovia licensing issues should this have otherwise
    been implemented directly within ONETEP.
  - To (easily) allow multiple images to be ran independently of one another (parallelisation).
  - To possibly support future interfacing of ONETEP with the POLYRATE package.


"""

__author__      = "Maximillian J. S. Phipps"
__copyright__   = "GNU General Public License v3"
__credits__ = ["Maximillian J. S. Phipps"]
__version__ = "0.1"
__maintainer__ = "Max Phipps"
__email__ = "m.phipps@soton.ac.uk"
__status__ = "Development"


import sys
import os
from shutil import copyfile 
from math import floor

import neb_xyz2dat
import neb_getForces


# Local pointer to submission scripts and recpot files:
hpc_files_dir='hpc_files'

# Lambda for steepest descent
lambdaval = 0.5


class xyz:
  """XYZ data class"""

  def __init__(self, atdat=[], xyzdat=[], coord=0):
    self.at = atdat
    self.xyz = xyzdat
    # set meta data
    # number of atoms
    self.nat = len(self.xyz[0])
    # replica energy pathway coordinate
    self.coord = '{:3.2f}'.format( coord )
    # construct filename for this replica
    self.datdir =      'im_'+self.coord
    self.filenamexyz = 'im_'+self.coord+'_init.xyz'
    self.filenamedat = self.datdir+\
                       '/im_'+self.coord+'.dat'
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
  def coord2str( self, coord ):
    # Parses coord float to string
    coordstr = '{:3.2f}'.format( coord )
    return 'im_'+coordstr

  @classmethod
  def fromdata(self, atdat, xyzdat, coord):
    # from xyz data:
    return xyz(atdat, xyzdat, coord)

  @classmethod
  def fromXYZ(self, filename, coord):
    # from file:
    [atdat,xyzdat] = getXYZdata(filename)
    return xyz(atdat, xyzdat, coord)

  def applyforces(self,forces):
    for ii in range(self.nat):
      #print self.xyz[ii][0]
      self.xyz[-1][ii][0] += forces[ii][0]
      self.xyz[-1][ii][1] += forces[ii][1]
      self.xyz[-1][ii][2] += forces[ii][2]
      #print self.xyz[ii][0]

  def writeXYZ(self, ixyz=-1, appendfile=False, filename=None):
    print 'writing replica at pathway coordinate ',str(self.coord), \
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


def getXYZdata(filename):
  """ returns filename's xyz data as atdat and xyzdat arrays.
  (returns the final xyz data if more than one xyz structure is present 
  in the file.) """

  # read lines to list
  with open(filename) as f:
    lines = list(f)

  # find number of XYZ structure entries in the file
  nat = int(lines[0])
  nentries = int(floor( float(len(lines))/float(nat+2) ))
  print 'Reading in structures up to final entry number ',str(nentries)

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

  return [atdat,xyzdat]

def rmsdist(im1,im2,ixyz1,ixyz2):
  """ returns rms distance between two replicas """
  sqdist = [0.0]*glob_nat
  # calculate atom pair square distances (same atom indices only)
  for ii in range(glob_nat):
    sqdist[ii] = ((im2.xyz[ixyz2][ii][0] - im1.xyz[ixyz1][ii][0])**2.0 + \
      (im2.xyz[ixyz2][ii][1] - im1.xyz[ixyz1][ii][1])**2.0 + \
      (im2.xyz[ixyz2][ii][2] - im1.xyz[ixyz1][ii][2])**2.0)
  # mean square distances
  ms = sum(sqdist)/glob_nat
  # root mean squared 
  rms = ms**0.5
  return rms

def tangvec(im1,im2,im3):
  ''' returns the unit vector tangent to the path for im2 '''
  r1 = im1.xyz[-1]
  r2 = im2.xyz[-1]
  r3 = im3.xyz[-1]
  rc = [[0.0,0.0,0.0]]*glob_nat
  for ii in range(glob_nat):
    # calculate ra = r1-r2 and rb = r3-r2
    ra = [ r1[ii][xx] - r2[ii][xx] for xx in [0,1,2]]
    rb = [ r3[ii][xx] - r2[ii][xx] for xx in [0,1,2]]
    #print "ra,rb:",str(rb),str(ra)
    # bisect the vectors (=rc)
    rc[ii] = [ (ra[xx] + rb[xx])/2. for xx in [0,1,2]]
    # normalise the vector
    norm = sum( [ rc[ii][xx]**2.0 for xx in [0,1,2] ] )**0.5
    # divide by zero error handling:
    if (norm == 0):
      # unit vector in 1,1,1
      tmp = 1.0/(3.**0.5)
      rc[ii] = [tmp,tmp,tmp]
    else:
      rc[ii] = [ rc[ii][xx] / norm for xx in [0,1,2] ]
  #print "print: ",str(rc),'\n'
  return rc

#########################################################################

# mode parameters
mode_initialise = 0
mode_iterate = 1
mode_plot = 2
mode_xyz = 3

# number of replicas
glob_nim = 6

MODE = -999
if ('-init' in sys.argv[1:]):
  MODE = mode_initialise
elif ('-iter' in sys.argv[1:]):
  MODE = mode_iterate
elif ('-plot' in sys.argv[1:]):
  MODE = mode_plot
elif ('-xyz' in sys.argv[1:]):
  MODE = mode_xyz
else:
  sys.exit('\nERROR: No mode flag given. Flag options:\n\
[ -init ] Initialise replicas by linear interpolation from reactants [.xyz] to products [.xyz]\n\
[ -iter ] Iterate by calculating nudged forces from singlepoint calculation and applying these to the xyz and dat files\n\
[ -plot ] Plots the energy pathway\n\
[ -xyz  ] Extracts the most recent path -> path.xyz\
\n')

if (MODE == mode_initialise):
  ''' Initialise the NEB optimisation by constructing guess 
  replicas by linear interpolation from reactant geometry to product geometry. '''

  # import optimised reactants and products
  react = xyz.fromXYZ(sys.argv[2],0)
  prod = xyz.fromXYZ(sys.argv[3],1)
  
  # Take shared atom names and number of atoms from reactants
  glob_nat = react.nat
  glob_at = react.at
  
  # initialise replicas
  ims = [[]]*(glob_nim+2)
  ims[0] = react
  ims[glob_nim] = prod
  for im in range(0,glob_nim):
    # linear interpolation between reactants and products
    ratio = float(im)/float(glob_nim-1)
    print 'initialising image at energy pathway coordinate ', ratio
  
    xyzdat = [[[0.0,0.0,0.0]]*glob_nat]
    for ii in range(glob_nat):
      xyzdat[0][ii] = [react.xyz[-1][ii][0] * (1.-ratio) + prod.xyz[-1][ii][0] * (ratio), \
          react.xyz[-1][ii][1] * (1.-ratio) + prod.xyz[-1][ii][1] * (ratio), \
          react.xyz[-1][ii][2] * (1.-ratio) + prod.xyz[-1][ii][2] * (ratio)]
  
    ims[im] = xyz.fromdata(glob_at,xyzdat,ratio)

  # writeout the initial xyzs to this directory
  for im in range(0,glob_nim):
    ims[im].writeXYZ()
  
  # convert our shifted images to dat files
  for im in range(glob_nim):
    print '-> energy pathway coordinate ', ims[im].coord
    print 'writing initial dat file to ', ims[im].filenamedat
    neb_xyz2dat.fromXYZ(ims[im].filenamexyz, ims[im].filenamedat)
    #print 'duplicating onetep submission script to ', ims[im].datdir
    #copyfile(hpc_files_dir+'/sub.sh', ims[im].datdir+'/sub.sh')

  #for im in range(0,glob_nim-1):
  #  print rmsdist(ims[im],ims[im+1])

else:

  # Read in the previous replica geometries (including reactants and products)
  # from the replica directories
  ims = [[]]*(glob_nim+2)
  for im in range(0,glob_nim):
    ratio = float(im)/float(glob_nim-1)
    #print ratio
    coord = xyz.coord2str(ratio)
    filename = coord+'/'+coord+'.xyz'
    ims[im] = xyz.fromXYZ(filename,ratio)

  # Take shared atom names and number of atoms from reactants
  glob_nat = ims[0].nat
  glob_at = ims[0].at

  if (MODE == mode_iterate):
    ''' Iterate a NEB optimisation step by extracting the forces from
    the previous singlepoint forces calculation for each replica and 
    steepest descent minimising the ion coordinates. '''
  
    # for each replica,
    # calculate and apply nudged forces
    for im in range(1,glob_nim-1):
    
      # Read forces from singlepoint forces calculation
      fion = neb_getForces.main(ims[im].filenameout, glob_nat)

      # Scale forces for steepest descent
      for ii in range(glob_nat):
         fion[ii][0] *= lambdaval
         fion[ii][1] *= lambdaval
         fion[ii][2] *= lambdaval

      #print "\nIMAGE"
      #print fion
      #print '' 
  
      # calculate path tangent unit vectors
      tv = tangvec(ims[im-1], ims[im], ims[im+1])
      
      # 1. calculate nudged forces
      ffinal = [[0.0,0.0,0.0]]*glob_nat
      for ii in range(glob_nat):
        # fs: sprung forces
        k3=0.05
        k1=0.05
        r1 = ims[im-1].xyz[-1][ii]
        r2 = ims[im].xyz[-1][ii]
        r3 = ims[im+1].xyz[-1][ii]
        fs = [ k3*(r3[xx]-r2[xx]) - k1*(r2[xx]-r1[xx]) for xx in [0,1,2] ] 
    
        # fs: project with unit vector tangent to path
        fsproj = [ fs[xx] * (tv[ii][xx]*tv[ii][xx]) for xx in [0,1,2] ]
    
        # v: forces with projection to path tangent
        fproj = [ fion[ii][xx] - fion[ii][xx]*(tv[ii][xx]*tv[ii][xx]) for xx in [0,1,2] ]
    
        # final neb forces:
        ffinal[ii] = [ fproj[xx] + fsproj[xx] for xx in [0,1,2] ]
  
      #print ffinal
    
      # 2. apply nudged forces
      ims[im].applyforces(ffinal)
  
      # convert our shifted images to dat files
      print '-> energy pathway coordinate ', ims[im].coord
      print 'backing up previous dat file to ', ims[im].filenamedat+'.bak'
      copyfile(ims[im].filenamedat, \
               ims[im].filenamedat+'.bak')
      print 'writing steepest descent updated dat file to ', ims[im].filenamedat
      neb_xyz2dat.fromdata(ims[im].at, ims[im].xyz, ims[im].filenamedat)
    
  
  elif (MODE == mode_plot):
    ''' Plot of the NEB minimisation '''
    import matplotlib.pyplot as plt
    from numpy import arange, concatenate, cumsum
  
    # Distances calculated using replica neighbour RMS distances, and normalised
    plt.xlabel('Minimum energy pathway coordinate')
    #plt.ylabel('Energy (a.u.)')
    plt.ylabel('Energy (kcal/mol)')


    # TODO: Calculate MEP coordinate by calling rmsdist() for the neighbouring replicas
    #x = arange(0.0,glob_nim,1.0)
    # TODO: x using RMS distances, but using the final set of xyz coordinates
    x = [[0.0]]
    num_paths = len(ims[1].xyz)

    # For each energy pathway...
    for ipath in range(num_paths):

      # Calculate the RMS distance between the neighouring replicas
      tmpx = [0.0]
      for im in range(0,glob_nim-1):
        # xyz entries:
        # for intermediate replicas, use the current xyz entry iteration
        ixyz1 = ipath
        ixyz2 = ipath
        # for reactants and products, there is only one entry, so override
        if (im == 0):          ixyz1 = 0
        if (im == glob_nim-2): ixyz2 = 0
        # calculate RMS distance 
        tmpx.append( rmsdist(ims[im],ims[im+1],ixyz1,ixyz2) )
      # distribute the distances by cumulative sum
      x.append( cumsum(tmpx) )
      # normalise to 1
      x[-1] = [xval/x[-1][-1] for xval in x[-1]]
    x = x[1:]
  

    # conversion ratio for y values (from Hartree to kcal/mol)
    conversion_ratio = 627.509

    # iterate intermediate replicas
    y = [[0.0]]*glob_nim
    shift = 0.0
    #for im in range(1,glob_nim-1):
    for im in range(glob_nim):
  
      # read lines to list
      filename = ims[im].filenameout
      with open(filename) as f:
        lines = list(f)
  
      # grep for '<-- CG'
      num_found = 0
      for line in lines:
        if '<-- CG' in line:
          E = float(line.split()[2]) * conversion_ratio + shift
          if (num_found > 0):
            y[im].append(E)
          else:
            y[im] = [E]
          num_found += 1

      # shift for y values so that reactants (im=0) is at energy=0
      if (im == 0):
        shift = -E
        y[0][0] += shift

    # # calculate number of paths (= maximum sub-array length)
    # num_paths = max( [ len(vec) for vec in y ] ) 

    # plot each path
    for ipath in range(num_paths):
      # NOTE: IF THIS FAILS, CHECK THAT ALL THE CALCULATIONS WERE FINISHED 
      # CLEANLY (i.e. WITH A '<-- CG' VALUE )
      ypathmid = [ y[xx][ipath] for xx in range(1,glob_nim-1) ]
      ypath = concatenate((y[0],ypathmid,y[-1]))
      xpath = x[ipath]
      color = str(1.0-float(ipath+1)/float(num_paths+1))
      plt.plot(xpath,ypath,'x-',color=color)
  
    # plt.title('Nudged Elastic Band Minimum Energy Path')
    plt.grid(True)
    plt.savefig("mep.png")
    plt.show()

  elif (MODE == mode_xyz):
    ''' Extract the most recent NEB pathway '''

    # writeout the latest xyzs to path.xyz
    # reactant's xyz
    ims[0].writeXYZ(appendfile=False, filename='path.xyz')
    # intermediate replicas' xyz
    for im in range(1,glob_nim-1):
      print im
      ims[im].writeXYZ(appendfile=True, ixyz=-1, filename='path.xyz')
    # product's xyz
    ims[glob_nim-1].writeXYZ(appendfile=True, filename='path.xyz')



