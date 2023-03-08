from argparse import ArgumentParser 
import numpy as np



parser = ArgumentParser()
parser.add_argument('--xyz_file', help='xyz file')
parser.add_argument('--gbw', default=None,help='gbw file')
parser.add_argument('--charge',default=0, help='charge')
parser.add_argument('--mult',default=1, help='multiplicity')
parser.add_argument('--job', help='orca job')
parser.add_argument('--procs',default=1, help='number of processors')

args = parser.parse_args()


xyz_file = args.xyz_file

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
with open(f'{args.job}_freq.inp', 'w') as f:
    if not args.gbw is None:
        f.write('! FREQ ωB97X def2-TZVPPD D4 MORead\n')
        f.write(f'%moinp "{args.gbw}"\n')
    else:
        f.write('! FREQ ωB97X def2-TZVPPD D4\n')
    f.write(f'%pal nprocs {args.procs} end\n')
    f.write(f'*xyz {args.charge} {args.mult}\n')
    for atom, position in zip(atoms, positions):
        f.write(f'{atom} {position[0]} {position[1]} {position[2]} \n')
    f.write('*')

