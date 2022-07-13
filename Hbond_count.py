import numpy as np
from ase import io, Atoms, neighborlist
from ase.visualize import view
from scipy import sparse
import argparse

parser = argparse.ArgumentParser()
#list of arguments
parser.add_argument('structure', metavar='G', help='Structure to hollow. In any ase structure format.', type=str)
parser.add_argument('-minH', '--min_Hlength', help='Minimun Hbond length (in A).', type=float, default='1.5')
parser.add_argument('-maxH', '--max_Hlength', help='Maximum Hbond length (in A).', type=float, default='2.2')
parser.add_argument('-ac', '--angle_cutoff', help='Maximun angle between two molecules for a Hbond to be compted (in degrees).', type=float, default='40')
parser.add_argument('-nonbr', '--no_number', help='Do not show the total number of Hbond H--O and O--H.', default=False, action='store_true')
parser.add_argument('-nol', '--no_length', help='Do not print the distance for each Hbond H--O.', default=False, action='store_true')
parser.add_argument('-con', '--connection', help='Print the molecules id for each Hbond length.', default=False, action='store_true')

args = parser.parse_args()

min_Hbond = args.min_Hlength
max_Hbond = args.max_Hlength
angle_cutoff = 180 - args.angle_cutoff
no_number = args.no_number
no_length = args.no_length
connection = args.connection

def molecule_Id(atoms):
    ''' 
    Takes the atoms object and associates each atoms with its corresponding molecule
    return the number of molecules and an array the size of atoms where the values component_list[i] gives the molecule id of the atom i
    '''
    cutOff = neighborlist.natural_cutoffs(atoms, mult=0.95)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=False)
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    #or: matrix = neighborlist.get_connectivity_matrix(neighborList.nl)
    n_components, component_list = sparse.csgraph.connected_components(matrix)
    return n_components, component_list

def second_minimun(distances):
    '''
    Takes an array of distances and select the second minimum distance
    '''

    first_minimum_arg = np.argmin(distances)
    second_minimum_arg = np.argmin(np.array(distances)[distances != np.min(distances)])

    if second_minimum_arg >= first_minimum_arg:
        second_minimum_arg = int(second_minimum_arg + 1)

    return second_minimum_arg


grain = io.read(args.structure) #select the grain to stucture
grain_points = grain.get_positions()
n_mol, list_atoms = molecule_Id(grain) #Search for the molecules inside the structure and id them

All_distances = grain.get_all_distances() #matrix with the distances of every atoms between every atoms

nbr_Hbond_H = 0
nbr_Hbond_other = 0
for i in range(len(grain)):
    #if i !=0: continue
    if grain[i].symbol =='H': #computes Hbonds from H atoms
        connected_atoms = [k for k, length in enumerate(All_distances[i]) if length > min_Hbond and length < max_Hbond and grain[k].symbol !='H' and list_atoms[i] != list_atoms[k] and (grain.get_angle(second_minimun(All_distances[i]), i, k) >= angle_cutoff) ]
        nbr_Hbond_H += len(connected_atoms)
        if no_length == False and len(connected_atoms) > 0:
            for atom in connected_atoms:
                if connection == False:
                    print(All_distances[i,atom])
                else:
                    print(i, atom, All_distances[i,atom])
    else: #computes Hbonds from O atoms
        connected_atoms = [k for k, length in enumerate(All_distances[i]) if length > min_Hbond and length < max_Hbond and grain[k].symbol =='H' and list_atoms[i] != list_atoms[k] and (grain.get_angle(second_minimun(All_distances[k]), k, i) >= angle_cutoff) ]
        nbr_Hbond_other += len(connected_atoms)

if no_number == False:
    print(nbr_Hbond_H, nbr_Hbond_other)
    
    