import numpy as np 
import ase
import ase.io
from chemiscope import write_input
import plumed
import sys
import argparse
import warnings


warnings.filterwarnings("ignore")

# Create the parser
my_parser = argparse.ArgumentParser(description='Generate a chemiscope json.gz file from a trajectory')

# Add the arguments
my_parser.add_argument('Trajectory',
                       metavar='traj',
                       type=str,
                       help='the pdb trajectory file (make sure to check atom name so they all start with a letter.)')

my_parser.add_argument('CVfile',
                      metavar='cv_file',
                      type=str,
                      help='The file listing the CVs to calculate. They should be written in `name: [CV definition]` style.')
my_parser.add_argument('Output',
                      metavar='output_file',
                      type=str,
                      help='The full name of the chemiscope file (it should end in `.json.gz`)')
my_parser.add_argument('-p','--plumed', required=True
                      )

# Execute the parse_args() method
args = my_parser.parse_args()

traj_file = args.Trajectory

cv_file = args.CVfile

out_file = args.Output

with open(cv_file) as f:
    cvs = f.readlines()
CVs = [i for i in cvs if i !='\n']
  

nCVs = len(CVs)

for cv in cvs:
   if 'PATH' in cv:
       nCVs += 1


print(f'Making a Chemiscope file with {nCVs} CVs. They are :\n{CVs}')


cv_np = np.zeros((nCVs,1))


# Read in trajectory using ase
traj = ase.io.read(traj_file,':')

# Setup plumed object to do calculation
p = plumed.Plumed()
p.cmd("setMDEngine","python")
p.cmd("setTimestep", 1.)
p.cmd("setKbT", 1.)
natoms = len(traj[0].positions)
p.cmd("setNatoms",natoms)
p.cmd("setLogFile","test.log")
p.cmd("init")

# Read plumed input 
if args.plumed:
    p.cmd("readInputLine",f'INCLUDE FILE={args.plumed}')
else:
    p.cmd("readInputLine",f"MOLINFO STRUCTURE={traj_file}") 

names=[]
n=0
for CV in CVs:
    p.cmd(f"readInputLine",CV)
    if 'PATH' in CV:
        name = CV.split()[0].strip(':')
        
        names.append(name+'.s')
        shape = np.zeros(1, dtype=np.int_)
        p.cmd(f"getDataRank {name}.s",shape)
        p.cmd(f"setMemoryForData {name}.s",cv_np[n])
        n += 1
        names.append(name+'.z')
        p.cmd(f"getDataRank {name}.z", shape)
        p.cmd(f"setMemoryForData {name}.z", cv_np[n])
        n+=1
      
    else:
        names.append(name)
        p.cmd("readInputLine", CV)
        # Get the shape of the value calculated in plumed
        shape = np.zeros( 1, dtype=np.int_ )
        p.cmd(f"getDataRank {name}", shape )

        # Now setup some memory to hold the variable that is shared 
        # between plumed and the underlying code
        p.cmd(f'setMemoryForData {name}', cv_np[n])
        n += 1

# # Loop over trajectory and get data from plumed
nfram, tt, vs, box = 0, [], [[] for x in range(nCVs)], np.array([[100.,0,0],[0,100.,0],[0,0,100]])
charges, forces, virial = np.zeros(natoms,dtype=np.float64), np.zeros([natoms,3]), np.zeros((3,3),dtype=np.float64)

for ts in traj :
    p.cmd("setStep",nfram)
    p.cmd("setBox",box )
    p.cmd("setMasses", ts.get_masses() )
    p.cmd("setCharges", charges )
    pos = np.array(ts.get_positions(), dtype=np.float64 )
    p.cmd("setPositions", pos )
    p.cmd("setForces", forces )
    p.cmd("setVirial", virial )
    p.cmd("calc")
    tt.append(nfram)

    for cv in range(nCVs):
        if cv_np[cv,0] > 0:
            vs[cv].append(cv_np[cv,0])
        else:
            vs[cv].append(2*np.pi + cv_np[cv,0])
        
    nfram = nfram + 1
    
atom_dict= {'H':1, 'C':6, 'N':7, 'O':8, 'S':16, '1':1, '2':1,'3':1}


for frame in traj:
    frame.numbers = np.array(
        [
            atom_dict[am[0]] 
            for am in frame.arrays["atomtypes"]
        ]
    )

    
# This constructs the dictionary of properties for chemiscope
properties = {}  
properties["time"]={"target": "structure","values": tt,"description": "Simulation step number"}
for i,name in enumerate(names):

    properties[name]= {'target':'structure', "values":vs[i], "description":'fix this late'} #CVs[i]}


# # This generates our chemiscope output
write_input(out_file, frames=traj, properties=properties )
