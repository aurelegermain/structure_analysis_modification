import matplotlib.pyplot as plt
from ase import io, Atoms, neighborlist, geometry
from ase.visualize import view
from scipy import sparse
from scipy.spatial.transform import Rotation as R
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
#list of arguments
parser.add_argument('structure', metavar='G', help='Structure to hollow. In any ase structure format.', type=str)
parser.add_argument('-ml', '--max_length', help='Each molecules at an inferior distance from the center of the structure will be deleted.', type=float, default='10')
parser.add_argument('-nolist', '--no_list', help='Do not print the list of atoms kept.', default=False, action='store_true')
parser.add_argument('-nos', '--no_structure', help='Do not print the hollowed out structure.', default=False, action='store_true')
parser.add_argument('-vs', '--view_structure', help='View the hollowed out structure.', default=False, action='store_true')

args = parser.parse_args()
file_name = args.structure
max_len = args.max_length
no_list = args.no_list
no_structure = args.no_structure
view_structure = args.view_structure

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

split_file_name = os.path.splitext(file_name) #split the file name into name and extension

#select the grain to study
grain = io.read(file_name)
grain_points = grain.get_positions()
grain_distances = geometry.get_distances([0,0,0], grain_points)[1][0]
n_mol, list_atoms = molecule_Id(grain) #Search for the molecules inside the structure and id them

list_atoms_inside  = [i for i,length in enumerate(grain_distances) if length <= max_len] #list of atoms inside the max length radius
list_molecules_inside = list(dict.fromkeys([i for atom_inside in list_atoms_inside for i in list(np.where(list_atoms == list_atoms[atom_inside])[0])])) #some atoms of the same molecules could not be include in the previous list so this correct this error
list_molecules_outside = [i for i in range(len(grain)) if i not in list_molecules_inside] #makes a list of molecules outside of the radius

if len(list_molecules_outside) == 0:
    print('No molecules in the structure. Try smaller length.')
    exit()

if no_list == False: #print list of atoms not deleted
    print(list_molecules_outside)

if no_structure == False:
    new_file_name = split_file_name[0] + '_hollow' + split_file_name[1] #name of the hollowed out new structure, with the same extension as the original structure file
    io.write(new_file_name, grain[list_molecules_outside]) #write the hollowed out file structure

if view_structure == True: #view the structure
    view(grain[list_molecules_outside])