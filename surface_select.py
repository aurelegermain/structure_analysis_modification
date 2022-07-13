import numpy as np
from ase import io, Atoms, neighborlist
from ase.visualize import view
from scipy import sparse
import argparse

parser = argparse.ArgumentParser()
#list of arguments
parser.add_argument('structure', metavar='G', help='structure to sample select the outer layer in any ase structure format.', type=str)
parser.add_argument('-d', '--distance', help='Radius of molecules selected around the molecule studied (in A).', type=float, default="5")
parser.add_argument('-ma', '--max_angle', help='Angle of the cone on top of the molecule studied. If a molecule is in this cone, it is considered on top of the studied one (in degrees).', type=float, default='140')
parser.add_argument('-mt', '--max_top', help='RMax number of molecules allowed to be on top of the one studied to be considered as on the outer layer.', type=int, default='1')
parser.add_argument('-xtbin', '--xtb_input', help='Print an xtb input with the top molecules constrained. The xtb input will be names "xtb.inp"', type=bool, default=False)
parser.add_argument('-nolist', '--no_list', help='Do not print the list of atoms of the outer layer.', type=bool, default=False)
parser.add_argument('-vs', '--view_structure', help='View the structure with the selected surface molecules.', type=bool, default=False)

args = parser.parse_args()

distance_select_mol = args.distance #in A, radius of the sphere of molecules we look at
max_angle = args.max_angle #in degrees, max angle so that a molecule is considered on top of the studied one
max_list_len = args.max_top #max list lentgh of number of molecules on top of the studied one

#print or not the outputs
xtb_input = args.xtb_input
no_list = args.no_list
view_structure = args.view_structure

def molecule_Id(atoms):
    ''' 
    takes the atoms object and associates each atoms with its corresponding molecule
    return the number of molecules and an array the size of atoms where the values component_list[i] gives the molecule id of the atom i
    '''
    cutOff = neighborlist.natural_cutoffs(atoms, mult=0.95)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=False)
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    #or: matrix = neighborlist.get_connectivity_matrix(neighborList.nl)
    n_components, component_list = sparse.csgraph.connected_components(matrix)
    return n_components, component_list

def barycentre(positions):
    ''' 
    Takes the x y z coordinates (as in atoms.get_coordinates()) and return the barycentre of the given atoms
    '''
    return np.sum(positions, axis = 0)/len(positions)

def centre(point, positions):
    ''' 
    Takes a point in x y z coordinates and a list of positions (as in atoms.get_positions()) a return the list of positions centered around the point 
    '''
    return positions - point

def spherical_coordinates(positions):
    ''' 
    Takes a list of positions in x y z coordinates (as in atoms.get_positions() ) and return the spherical coordinates as r theta phi (theta 0 to pi, phi 0 to 2pi)
    '''
    r = np.sqrt(np.sum(np.square(positions), axis = 1))
    theta = np.arctan2(np.sqrt(np.square(positions[:,0]) + np.square(positions[:,1])), positions[:,2]) 
    phi = np.arctan2(positions[:,1], positions[:,0])

    return np.r_['1,2,0', r, theta, phi]

def cartesian_coordinates(positions):
    '''
    Takes a list of posistions in spherical coordinates (as in the sphercal coordinates function) and return the cartesion coordinates x y z (as in atoms.get_coordinates() )
    '''
    x = np.cos(positions[:,2])*np.sin(positions[:,1])
    y = np.sin(positions[:,2])*np.sin(positions[:,1])
    z = np.cos(positions[:,1])                                       
    return np.r_['1,2,0', x, y, z]

grain = io.read(args.structure) #select the grain to stucture
grain_points = grain.get_positions()
n_mol, list_atoms = molecule_Id(grain) #Search for the molecules inside the structure and id them

#Compute the barycentre of each molecule to make a new structure out of those barycentre 
positions_barycentre = np.zeros([n_mol, 3])
for i in range(n_mol):
    list_i = [index for index in range(len(list_atoms)) if list_atoms[index] == i]
    positions_barycentre[i,:] = barycentre(grain_points[list_i,:])
barycentre_structure = Atoms(str(n_mol) + 'N', positions=positions_barycentre)
barycentre_structure += Atoms('C', positions=[[0,0,0]])

#Construct the list of molecules (not atoms) present at the structure surface
list_mol_surface = []
for i in range(n_mol):
    distances = barycentre_structure.get_distances(i, range(n_mol)) #get the distances between the selected molecule i and all the others
    molecule_select = [i for i in range(len(distances)) if distances[i] <= distance_select_mol] #select only the molecules inside a certain radius of the studied molecule i
    angles = [barycentre_structure.get_angle(n_mol, i, j) if j != i else 0 for j in molecule_select] #compute the angle between the centre of the structure and the molecule selected i for each molecules in the molecule_select list
    list_mol_up = [molecule_select[j] for j in range(len(molecule_select)) if angles[j] > max_angle and molecule_select[j]] #make a list of every molecules with an angle of more than max_angle (suposedly on top of the molecule selected)
    
    if len(list_mol_up) <= max_list_len: #if len of list is inferior or equal to max_list_len then it means the molecule selected is at the surface
        list_mol_surface += [i] #adds i to the list of surface molecules

list_atoms_surface = list(dict.fromkeys([i for atom_inside in range(len(list_mol_surface)) for i in list(np.where(list_atoms == list_mol_surface[atom_inside])[0])])) #convert list_mol_up_final in a list of atoms instead of molecules

if no_list == False:
    print(list_atoms_surface)

if view_structure == True: #show the structure if the user wants
    list_atoms_not_surface = [i for i in range(len(grain)) if i not in list_atoms_surface] #list of atoms not at the surface
    chemical_symbols = grain.get_chemical_symbols()
    chemical_symbols_not_up_final = [chemical_symbols[i] for i in range(len(chemical_symbols)) if i in list_atoms_not_surface]
    final_structure_atoms = Atoms(chemical_symbols_not_up_final, positions=grain_points[list_atoms_not_surface]) + Atoms(str(len(list_atoms_surface)) + 'C', positions=grain_points[list_atoms_surface])
    view(final_structure_atoms)

if xtb_input == True: #Print an xtb.inp constraining the outer layer atoms if xtb_input is True
    file_xtb_unfixed_input = open('./xtb.inp','w')
    print('$constrain', file=file_xtb_unfixed_input)
    print('    atoms: ', end='', file=file_xtb_unfixed_input)
    list_atoms_surface = np.atleast_1d(list_atoms_surface)
    for k in range(len(list_atoms_surface)):
        if k!=0:
            if k==len(list_atoms_surface)-1:
                if last_fix == k - 1:
                    print('-' + str(list_atoms_surface[k]+1), end='', file=file_xtb_unfixed_input)
                    j = j + 1
                else:
                    print(',' + str(list_atoms_surface[k]+1), end='', file=file_xtb_unfixed_input)
            else:
                if list_atoms_surface[last_fix] == list_atoms_surface[k] - 1:
                    last_fix = k
                else:
                    print('-' + str(list_atoms_surface[last_fix]+1) + ',' + str(list_atoms_surface[k]+1), end='', file=file_xtb_unfixed_input)
                    last_fix = k
                    j = j + 1
        elif k==0:
            j = 0
            last_fix = k
            print(list_atoms_surface[k]+1, end='', file=file_xtb_unfixed_input)
    print('\n$end', file=file_xtb_unfixed_input)
    file_xtb_unfixed_input.close()