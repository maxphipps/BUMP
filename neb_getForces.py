#!/usr/bin/env python

"""getForces.py: Extracts forces from ONETEP output file """

import sys

def main(fnameIn,nat):

  # Read ONETEP output file
  with open(fnameIn) as f:
    lines = list(f)

  # Find line number forces begin at
  readline = 0
  for linenum in range(len(lines)):
    if '****************** Unconstrained *******************' in lines[linenum]:
      if '********************** Forces **********************' in lines[linenum+1]:
        readline = linenum + 7
        # We do not break here, as this file may contain a series of 
        # calculation outputs, of which we wish to use the final output

  if (readline == 0):
    sys.exit('\nERROR: Could not find forces in the output file. Either:\n\
- A force write-out was not requested in the calculation.\n\
- The file is not a ONETEP output file.\n\
- You are running a new version of ONETEP which this script doesn\'t properly parse.\n\
- This script is broken.\n')

  # Construct forces array
  forces = [[]]*nat
  ii = 0
  for line in lines[readline:readline+nat]:
    forces_str = line.split()[3:6]
    forces[ii] = [float(forces_str[xx]) for xx in [0,1,2]]
    ii += 1
      
  return forces
  
if __name__ == "__main__":
  fnameIn = sys.argv[1]
  nat = int(sys.argv[2])
  print main(fnameIn,nat)

