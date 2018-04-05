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
   Input geometry files can be easily produced using the command:
     python neb_xyz2dat.py in.xyz out.dat
   Be sure to modify the 'TASK' keyword setting to 'GEOMETRYOPTIMIZATION'.
2. The following files are then created/updated by the user:
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
     and a velocity verlet or steepest descent geometry step taken.
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

#from numpy import *
import numpy as np
from numpy import matrix

import neb_xyz2dat
import neb_getForces
from neb_classes import XYZ

# Debug flags
NO_FORWARD_X = False # True

# Energy in kJ/mol or kcal/mol
EinkJ = True

# Steepest descent of velocity verlet (without half step) 
MINMODE_SD = 0
MINMODE_VV = 1
MINMODE = MINMODE_SD

# Time step/step size for optimisation
deltat = 1.2

# Spring constants
SPRING_K = 1.0 # 0.10

# mode parameters
mode_initialise = 0
mode_iterate = 1
mode_plot = 2
mode_xyz = 3

# number of replicas
glob_nim = 5
  
class prettyfloat(float):
  def __repr__(self):
    return "%0.1f" % self
    #return "%0.2f" % self

def mat2arr(M):
  # numpy utility tool
  return np.squeeze(np.asarray(M))

def rmsdist(im1,im2,ixyz1,ixyz2):
  ''' returns rms distance between two replicas '''
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

def norm_vec3(vec):
  return sum( [ vec[xx]**2.0 for xx in [0,1,2] ] )**0.5

def unit_vec3(vec):
  norm = sum( [ vec[xx]**2.0 for xx in [0,1,2] ] )**0.5
  # divide by zero error handling:
  if (norm == 0):
    # unit vector in 1,1,1
    print "DIV0 ERROR!"
    #tmp = 1.0/(3.**0.5)
    tmp = 0.0
    vec = [tmp,tmp,tmp]
  else:
    vec = [ vec[xx] / norm for xx in [0,1,2] ]
  return vec

def tangvec(im1,im2,im3):
  ''' returns the unit vector tangent to the path for im2 '''
  r1 = im1.xyz[-1]
  r2 = im2.xyz[-1]
  r3 = im3.xyz[-1]
  rtang = [[0.0 for i in range(3)] for j in range(glob_nat)]
  for ii in range(glob_nat):
    # calculate vectors along the path
    #ra = r1-r2 and rb = r3-r2
    ra = [ r2[ii][xx] - r1[ii][xx] for xx in [0,1,2]]
    rb = [ r3[ii][xx] - r2[ii][xx] for xx in [0,1,2]]

    # normalise the vectors

    norma = sum( [ ra[xx]**2.0 for xx in [0,1,2] ] )**0.5
    # divide by zero error handling:
    if (norma == 0):
      # unit vector in 1,1,1
      print "DIV0 ERROR!"
      #tmp = 1.0/(3.**0.5)
      tmp = 0.0
      ra = [tmp,tmp,tmp]
    else:
      ra = [ ra[xx] / norma for xx in [0,1,2] ]

    normb = sum( [ rb[xx]**2.0 for xx in [0,1,2] ] )**0.5
    # divide by zero error handling:
    if (normb == 0):
      # unit vector in 1,1,1
      print "DIV0 ERROR!"
      #tmp = 1.0/(3.**0.5)
      tmp = 0.0
      rb = [tmp,tmp,tmp]
    else:
      rb = [ rb[xx] / normb for xx in [0,1,2] ]

    # norm of ra+rb
    normtang = sum( [ (ra[xx]+rb[xx])**2.0 for xx in [0,1,2] ] )**0.5
    # divide by zero error handling:
    if (normtang == 0):
      # unit vector in 1,1,1
      print "DIV0 ERROR!"
      #tmp = 1.0/(3.**0.5)
      tmp = 0.0
      rtang[ii] = [tmp,tmp,tmp]
    else:
      rtang[ii] = [ (ra[xx]+rb[xx]) / normtang for xx in [0,1,2] ]

    ##print norma, normb
    #print "RB:", rb
    #print im2.coord_val, im3.coord_val
    ##print [ (ra[xx]+rb[xx])**2.0 for xx in [0,1,2] ]
    ##print normtang

  return rtang

#########################################################################

if __name__ == "__main__":

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
    react = XYZ.fromXYZ(sys.argv[2],0,readVelo=False)
    prod = XYZ.fromXYZ(sys.argv[3],1,readVelo=False)
    
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
    
      xyzdat = [[[0.0 for i in range(3)] for j in range(glob_nat)]]
      for ii in range(glob_nat):
        xyzdat[0][ii] = [react.xyz[-1][ii][0] * (1.-ratio) + prod.xyz[-1][ii][0] * (ratio), \
            react.xyz[-1][ii][1] * (1.-ratio) + prod.xyz[-1][ii][1] * (ratio), \
            react.xyz[-1][ii][2] * (1.-ratio) + prod.xyz[-1][ii][2] * (ratio)]
    
      ims[im] = XYZ(glob_at,xyzdat,ratio,readVelo=False)
  
    # writeout the initial xyzs to this directory
    # and the initial (zero) velocities
    for im in range(0,glob_nim):
      ims[im].writeXYZ()
      ims[im].writeVelo()
    
    # convert our shifted images to dat files
    for im in range(glob_nim):
      print '-> energy pathway coordinate ', ims[im].coord_str
      print 'writing initial dat file to ', ims[im].filenamedat
      neb_xyz2dat.fromXYZ(ims[im], ims[im].filenamedat)
  
    #for im in range(0,glob_nim-1):
    #  print rmsdist(ims[im],ims[im+1])
  
  else:
  
    # Read in the previous replica geometries (including reactants and products)
    # from the replica directories
    ims = [[]]*(glob_nim+2)
    for im in range(0,glob_nim):
      ratio = float(im)/float(glob_nim-1)
      #print ratio
      coord_str = XYZ.float2str(ratio)
      filename = coord_str+'/'+coord_str+'.xyz'
      ims[im] = XYZ.fromXYZ(filename,ratio,readVelo=True)
  
    # Take shared atom names and number of atoms from reactants
    glob_nat = ims[0].nat
    glob_at = ims[0].at
  
    if (MODE == mode_iterate):
      ''' Iterate a NEB optimisation step by extracting the forces from
      the previous singlepoint forces calculation for each replica and 
      taking a minimisation step of the ion coordinates. '''

      deltapos = [[[0.0 for i in range(3)] for j in range(glob_nat)] for k in range(glob_nim)]
  
      # for each replica,
      # calculate and apply nudged forces
      for im in range(1,glob_nim-1):
      
        # Read forces from singlepoint forces calculation
        fion = neb_getForces.main(ims[im].filenameout, glob_nat)
  
        #print "\nIMAGE"
        #print fion
        #print '' 
    
        # calculate path tangent unit vectors
        tv = tangvec(ims[im-1], ims[im], ims[im+1])

        # 1. calculate nudged forces
        ffinal = [[0.0 for i in range(3)] for j in range(glob_nat)]
        for ii in range(glob_nat):

          # calculate path tangent unit vector projector for the atom ii
          tv_mat = matrix(tv[ii])
          #pprint.pprint( tv_mat )
          tvproj = mat2arr( tv_mat.T.dot( tv_mat ) )
          #pprint.pprint( tv_mat.T.dot( tv_mat ) )

          # fs: spring forces
          k3=SPRING_K
          k1=SPRING_K
          r1 = ims[im-1].xyz[-1][ii]
          r2 = ims[im].xyz[-1][ii]
          r3 = ims[im+1].xyz[-1][ii]

          r32 = norm_vec3([r3[xx]-r2[xx] for xx in [0,1,2]])
          r21 = norm_vec3([r2[xx]-r1[xx] for xx in [0,1,2]])
          # NEB springs: project with unit vector tangent to path
          fsproj = [ (k3*r32 - k1*r21) * (tv[ii][xx]) for xx in [0,1,2] ]
      
          # true potential: forces with projection to path normal
          fprojdottt = [ sum(fion[ii][kk] * tvproj[kk][xx] for kk in [0,1,2]) for xx in [0,1,2]] 
          fproj = [ fion[ii][xx] - fprojdottt[xx] for xx in [0,1,2] ]
          #fproj = [ fion[ii][xx] * (1.-(tv[ii][xx])) for xx in [0,1,2] ] # ORIGINAL
      
          # final neb forces:
          ffinal[ii] = [ ( fproj[xx] + fsproj[xx] ) for xx in [0,1,2] ]
	  #print "TEST IM",im, map(prettyfloat,ffinal[ii][:]), map(prettyfloat,tv[ii][:])
	  #print "TEST IM",im, map(prettyfloat,ffinal[ii][:]), map(prettyfloat,fproj[:]), map(prettyfloat,fsproj[:])
	  #print "TEST IM",im, map(prettyfloat,ffinal[ii][:]), map(prettyfloat,fion[ii][:]), map(prettyfloat,fprojdottt[:])
	  #print "TEST IM",im, tv[ii][:]
    
        #print "FFINAL",ffinal
        #print "FFINAL AV im=", im, ((sum(ffinal[0])+sum(ffinal[1])+sum(ffinal[2]))/(3*len(ffinal)))
      

        # 2. apply nudged forces
        if (MINMODE == MINMODE_SD):
          # Scale forces by timestep -> deltapos = 0.5*a0*dt**2
          for ii in range(glob_nat):
            deltapos[im][ii][0] = 0.5 * ffinal[ii][0] * deltat**2
            deltapos[im][ii][1] = 0.5 * ffinal[ii][1] * deltat**2
            deltapos[im][ii][2] = 0.5 * ffinal[ii][2] * deltat**2
            #print "DELTAPOS IM ", im, " AT ", ii, "=", deltapos[im][ii]

        elif (MINMODE == MINMODE_VV):
          # Velocity verlet:

          # x1 = x0 + 0.5 a0 dt**2 + v0 dt 
          for ii in range(glob_nat):
            deltapos[im][ii][0] = 0.5 * ffinal[ii][0] * deltat**2 - ims[im].xyzvelo[-1][ii][0]*deltat
            deltapos[im][ii][1] = 0.5 * ffinal[ii][1] * deltat**2 - ims[im].xyzvelo[-1][ii][1]*deltat
            deltapos[im][ii][2] = 0.5 * ffinal[ii][2] * deltat**2 - ims[im].xyzvelo[-1][ii][2]*deltat


          #u =  [[0.0 for i in range(3)] for j in range(glob_nat)]
          #uF =  [[0.0 for i in range(3)] for j in range(glob_nat)]
          #uFF = [[0.0 for i in range(3)] for j in range(glob_nat)]
          ## uF
          #for ii in range(glob_nat):
          #  for xx in range(2):
          #    u[ii*3+xx] = ims[im].xyzvelo[-1][ii][xx]
          #    # TODO make u into unit vector
          ## uF
          #for ii in range(glob_nat):
          #  for xx in range(2):
          #    uF[ii*3+xx] = ims[im].xyzvelo[-1][ii][xx]*ffinal[ii][xx]
          ## uFF
          #for ii in range(glob_nat):
          #  for xx in range(2):
          #    uFF[ii*3+xx] = ims[im].xyzvelo[-1][ii][xx]*ffinal[ii][xx]
            

          # update velocities
          # v1 = v0 + (a0+a1)/2 * dt 
          # NOTE: this is not the formal velocity verlet algo, as we do not average the forces.
          # Instead, we use the set of forces associated with the previous coordinates
          # (i.e. a0, associated with x0) to update the velocities.
          # TODO: calculation of a1
          for ii in range(glob_nat):
            # To prevent overshooting the minimum:
            # If (velocity component * force component < 0) then we
            # have overshot the minimum, so set velocity component to 0
            # TODO: modify to use sign, rather than doing the full multiplication (speedup)


            if (ims[im].xyzvelo[-1][ii][0] * ffinal[ii][0] < 0):
              # Overshoot -> zero the velocity component
              ims[im].xyzvelo[-1][ii][0] = 0
            else:
              # Add the force component to the velocity component
              ims[im].xyzvelo[-1][ii][0] += ffinal[ii][0] * deltat

            if (ims[im].xyzvelo[-1][ii][1] * ffinal[ii][1] < 0):
              ims[im].xyzvelo[-1][ii][1] = 0
            else:
              ims[im].xyzvelo[-1][ii][1] += ffinal[ii][1] * deltat

            if (ims[im].xyzvelo[-1][ii][2] * ffinal[ii][2] < 0):
              ims[im].xyzvelo[-1][ii][2] = 0
            else:
              ims[im].xyzvelo[-1][ii][2] += ffinal[ii][2] * deltat

          # Write the updated velocities
          ims[im].writeVelo(appendfile=True)

      #print "--C-H-E-C-K-I-N-G--" # DEBUG

      for im in range(1,glob_nim-1):

        #for ii in range(glob_nat): # DEBUG
        #  print "DELTAPOS IM ", im, " AT ", ii, "=", deltapos[im][ii] # DEBUG

        # Add the position delta
        ims[im].applypositionchange(deltapos[im])
    
        # convert our shifted images to dat files
        print '-> energy pathway coordinate ', ims[im].coord_str
        print 'backing up previous dat file to ', ims[im].filenamedat+'.bak'
        copyfile(ims[im].filenamedat, \
                 ims[im].filenamedat+'.bak')
        print 'writing gradient updated dat file to ', ims[im].filenamedat
        neb_xyz2dat.fromXYZ(ims[im], ims[im].filenamedat)
      
    
    elif (MODE == mode_plot):
      ''' Plot of the NEB minimisation '''
  
      # Set backend to Agg
      # NOTE: This will turn OFF to-screen rendering!
      if (NO_FORWARD_X):
        import matplotlib as mpl
        mpl.use('Agg')
  
      import matplotlib.pyplot as plt
      from numpy import arange, concatenate, cumsum
    
      if (EinkJ):
        # conversion ratio for y values (from Hartree to kJ/mol)
        conversion_ratio = 627.509438736/0.238845896627
      else:
        # conversion ratio for y values (from Hartree to kcal/mol)
        conversion_ratio = 627.509438736
  
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
    
  
      # iterate intermediate replicas
      y = [0.0]*glob_nim
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
  
      def plot_path(output=False):
        # Energy pathway plotting function

        # Distances calculated using replica neighbour RMS distances, and normalised
        plt.xlabel('Minimum energy pathway coordinate')
        #plt.ylabel('Energy (a.u.)')
        if (EinkJ):
          plt.ylabel('Energy (kJ/mol)')
        else:
          plt.ylabel('Energy (kcal/mol)')
  
        for ipath in plot_range:
          # NOTE: IF THIS FAILS, CHECK THAT ALL THE CALCULATIONS WERE FINISHED 
          # CLEANLY (i.e. WITH A '<-- CG' VALUE )
          ypathmid = [ y[xx][ipath] for xx in range(1,glob_nim-1) ]
          ypath = concatenate((y[0],ypathmid,y[-1]))
          xpath = x[ipath]
          color = str(1.0-float(ipath+1)/float(num_paths+1))
          plt.plot(xpath,ypath,'x-',color=color)

          # Print values to screen
          if (output):
            print "NEB pathway:"
            if (EinkJ):
              print "coord_arb energy_kj_mol"
            else:
              print "coord_arb energy_kcal_mol"
            for ii in range(len(xpath)):
              print xpath[ii], ypath[ii]
         
        # plt.title('Nudged Elastic Band Minimum Energy Path')
        plt.grid(True)
        plt.savefig(filename_mep)
        plt.show()
        plt.clf()
  
      #filename_mep = ""
      #if (PLOT_ALL):
      #  # plot each path
      #  plot_range = list(range(num_paths))
      #  filename_mep = "mep.png"
      #else:
      #  # plot final path
      #  plot_range = [num_paths-1]
      #  filename_mep = "mep_final.png"
      
      # plot all paths
      plot_range = list(range(num_paths))
      filename_mep = "mep.png"
      plot_path(output=False)
  
      # plot final path
      plot_range = [num_paths-1]
      filename_mep = "mep_final.png"
      plot_path(output=True)
  
  
    elif (MODE == mode_xyz):
      ''' Extract the most recent NEB pathway '''
  
      # writeout the latest xyzs to path.xyz
      # reactant's xyz
      ims[0].writeXYZ(appendfile=False, filename='path.xyz')
      # intermediate replicas' xyz
      for im in range(1,glob_nim-1):
        ims[im].writeXYZ(appendfile=True, ixyz=-1, filename='path.xyz')
      # product's xyz
      ims[glob_nim-1].writeXYZ(appendfile=True, filename='path.xyz')
  
  
