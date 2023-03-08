import os 
import numpy as np
from warnings import warn



# orca script to automatically generate fragment orca.inp file

class Fragment:
  
    # Initialize fragment class via xyz file
    def __init__(self, xyz_file):
        with open(xyz_file) as f:
            lines=f.readlines()
        self.num_atoms=int(lines[0])
        atoms=np.empty(shape=(self.num_atoms,1),dtype='str')
        positions=np.empty(shape=(self.num_atoms,3),dtype='float')
        for i in range(self.num_atoms):
            atom, x, y, z = lines[i+2].split() #Skip first two lines
            atoms[i] = atom 
            positions[i,:] = float(x), float(y), float(z)
        self.atoms=atoms 
        self.positions=positions 
        self.charge=0
        self.mult=1
        self.procs=1
        self.fragments=None
        self.bond_constraints=None
        self.angle_constraints=None
        self.dihedral_constraints=None
        self.frag_constraints=None
        self.frag_connections=None

    def add_fragment(self, start=None, end=None, lst=None):
        if start is None and end is None and lst is None:
            raise Exception('Either starting and ending indexes are provided or a list of atom indexes')
        elif start is not None and end is not None and list is not None:
            warn('Provided starting and ending indexes and a list, just using list')
        elif start is not None and end is not None and lst is None: 
            lst = [i for i in range(start, end + 1)]
        if self.fragments is None:
            self.fragments=[lst] 
        else:
            self.fragments.append(lst)
    
    def connect_fragments(self, frag_one, frag_two, atom_one, atom_two):
        connection = [frag_one, frag_two, atom_one, atom_two]
        if self.frag_connections is None:
            self.frag_connections=[connection]
        else:
            self.frag_connections.append(connection)

    def set_charge(self, charge):
        self.charge=charge

    def set_mult(self,mult):
        self.mult=mult
    
    
    def set_procs(self, procs):
      self.procs=procs
        
    def set_job(self, job):
        # Set the first line of the orca input file 
        self.job=job

    def constrain_fragments(self,frag):
        if self.frag_constraints is None:
            self.frag_constraints=[frag]
        else:
            self.frag_constraints.append(frag)
    
    def constrain_atoms(self, atom_one, atom_two, dist=None):
        #If dist is none, then defaults to current distance
        if dist is None:
            dist = np.linalg.norm(self.positions[atom_one,:]-self.positions[atom_two,:])
        if self.bond_constraints is None:
            self.bond_constraints=[[atom_one, atom_two, dist]]
        else:
            self.bond_constraints.append([atom_one, atom_two, dist])

    def constrain_angles(self, atom_one, atom_two, atom_three, angle=None):
        # If angle is none, then defaults to current angle
        if angle is None:
            vec_one = self.positions[atom_one,:] - self.positions[atom_two,:]
            vec_two = self.positions[atom_three,:] - self.positions[atom_two,:]
            angle = np.arccos(np.dot(vec_one, vec_two)/(np.linalg.norm(vec_one)*np.linalg.norm(vec_two)))
        if self.angle_constraints is None:
            self.angle_constraints=[[atom_one, atom_two, atom_three, angle]]
        else:
            self.angle_constraints.append([atom_one, atom_two, atom_three, angle])
        
    def constrain_dihedrals(self, atom_one, atom_two, atom_three, atom_four, dihedral=None):
        # If dihedral is none, then defaults to current dihedral
        if dihedral is None:
            vec_one = self.positions[atom_one,:] - self.positions[atom_two,:]
            vec_two = self.positions[atom_three,:] - self.positions[atom_two,:]
            vec_three = self.positions[atom_four,:] - self.positions[atom_three,:]
            vec_four = self.positions[atom_four,:] - self.positions[atom_three,:]
            dihedral = np.arccos(np.dot(vec_one, vec_two)/(np.linalg.norm(vec_one)*np.linalg.norm(vec_two)))
        if self.dihedral_constraints is None:
            self.dihedral_constraints=[[atom_one, atom_two, atom_three, atom_four, dihedral]]
        else:
            self.dihedral_constraints.append([atom_one, atom_two, atom_three, atom_four, dihedral])

        

    
    def write_orca_file(self, filename=None):
        if filename is None:
            filename='orca.inp'
        
        with open(filename,'w') as orca_file:
            orca_file.write(f'! {self.job}\n')
            orca_file.write(f'%PAL NPROCS {self.procs} END')
            

            # Add geom constraints for the fragment 
            orca_file.write('%geom\n')
            if not self.bond_constraints is None and not self.angle_constraints is None and not self.dihedral_constraints is None:
                orca_file.write('\tConstraints\n')
            if not self.bond_constraints is None:
                for bond_constraint in self.bond_constraints:
                    orca_file.write('\t\t{')
                    orca_file.write(f'B {bond_constraint[0]} {bond_constraint[1]} {bond_constraint[2]} C')
                    orca_file.write(' }\n')
            
            if not self.angle_constraints is None:
                for angle_constraint in self.angle_constraints:
                    orca_file.write('\t\t{')
                    orca_file.write(f'A {angle_constraint[0]} {angle_constraint[1]} {angle_constraint[2]} {angle_constraint[3]} C')
                    orca_file.write(' }\n')
            
            if not self.dihedral_constraints is None:
                for dihedral_constraint in self.dihedral_constraints:
                    orca_file.write('\t\t{')
                    orca_file.write(f'D {dihedral_constraint[0]} {dihedral_constraint[1]} {dihedral_constraint[2]} {dihedral_constraint[3]} {dihedral_constraint[4]} C')
                    orca_file.write(' }\n')

            orca_file.write('\t\tend\n')

            if not self.frag_constraints is None:
                constraints=''
                orca_file.write('\tConstrainFragments { ')
                for i in self.frag_constraints:
                    orca_file.write(f'{i} ')
                orca_file.write('} end\n')

            orca_file.write('\tConnectFragments\n')
            for connection in self.frag_connections:
                orca_file.write(f'\t\t{{{connection[0]} {connection[1]} C {connection[2]} {connection[3]}}}\n')
            orca_file.write('\tend\n')
            orca_file.write('end\n')



            orca_file.write(f'* xyz {self.charge} {self.mult}\n')
            for n, atom in enumerate(self.atoms):
                for j, frag in enumerate(self.fragments):
                    if n in frag:
                        orca_file.write(f'{atom[0]}({j+1})\t{self.positions[n,0]}\t{self.positions[n,1]}\t{self.positions[n,2]}\n')
                        break 
            orca_file.write('*')
