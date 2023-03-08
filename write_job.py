from fragment import Fragment 
import sys


frag = Fragment(xyz_file=sys.argv[1], charge=0, mult=2)
frag.set_job(job='Opt PBE0 def2-TZVP D4')
frag.add_fragment(lst=[i for i in range(0,18)])
frag.add_fragment(lst=[i for i in range(18, 25)])
frag.set_procs(32)
frag.connect_fragments(frag_one=1, frag_two=2,atom_one=1,atom_two=18)
frag.constrain_atoms(1, 18, dist=sys.argv[3])
frag.write_orca_file(filename=sys.argv[2])