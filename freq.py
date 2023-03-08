import sys
import numpy as np


xyz_file = sys.argv[1]

# Read in xyz file
with open(xyz_file, 'r') as f:
    lines = f.readlines()

# Get number of atoms
num_atoms = int(lines[0])

# Get atoms and coordinates
atoms, positions = [], np.empty((num_atoms, 3))
for i in range(2, num_atoms + 2):
    atom, x, y, z = lines[i].split()
    atoms.append(atom)
    positions[i - 2, :] = float(x), float(y), float(z)

# Write orca job
with open(f'{xyz_file[:-4]}_freq.inp', 'w') as f:
    f.write('! FREQ ωB97X def2-TZVPPD D4\n')
    f.write('*xyz 0 1\n')
    for atom, position in zip(atoms, positions):
        f.write(f'{atom} {position[0]} {position[1]} {position[2]} \n')
    f.write('*')

