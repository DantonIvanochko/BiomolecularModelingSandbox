# modules

# pyrosetta
from pyrosetta import * 
pyrosetta.init('-include_sugars -write_pdb_link_records -ignore_zero_occupancy false') # within parenthesis is required for working with glycans ('-include_sugars -write_pdb_link_records')
# pyrosetta - clean
from pyrosetta.toolbox import cleanATOM
# pyrosetta - add glycan
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover, GlycanTreeModeler
# pyrosetta - relax
from pyrosetta.rosetta.protocols.relax import FastRelax
# pyrosetta - build disulfides
from pyrosetta.rosetta.protocols.denovo_design import DisulfidizeMover
from rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.protocols.backrub import BackrubMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
# misc
import os
import math
import random
# for mutate_residue
from pyrosetta.teaching import *
import sys, os
import random
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory

def clean_structure(input_pdb):
    print("cleaning structure")
    cleanATOM(input_pdb, input_pdb.replace("..", "").replace(".pdb", "") + ".clean.pdb")



def add_glycan(input_pdb, residue_index, glycan_chain):

    print("adding glycan to pose")
    test_pose = pose_from_pdb(input_pdb)
    
    select = ResidueIndexSelector()
    select.set_index(residue_index)
    
    glycosylate = SimpleGlycosylateMover()
    glycosylate.set_glycosylation(glycan_chain)
    glycosylate.set_residue_selector(select)
    glycosylate.apply(test_pose)
    
    print("fixing glycan on pose")
    glycan_modeler = GlycanTreeModeler()
    scorefxn = get_fa_scorefxn()  
    glycan_modeler.set_scorefxn(scorefxn)
    glycan_modeler.set_hybrid_protocol(True)
    glycan_modeler.set_use_shear(True)
    glycan_modeler.set_use_gaussian_sampling(True)
    glycan_modeler.apply(test_pose)
    
    
    print(test_pose.residue(residue_index))
    print(test_pose.residue_type(residue_index).carbohydrate_info())
    test_pose.dump_pdb(input_pdb.replace(".pdb", "") + "-" + str(residue_index) +"_" + glycan_chain + ".pdb")



def relax_structure(input_pose_to_relax, constrain_to_input=False, cycles=5, cartesian=True):

    print(f'preparing {input_pose_to_relax} structure with FastRelax()')
    
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_branches(setting=True) # sugars
    mmf.all_nu(setting=True) # sugars

    
    relax = FastRelax(standard_repeats=cycles)
    
    if cartesian == True:
        mmf.set_cartesian(setting=True)
        scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        relax.set_scorefxn(scorefxn)
        relax.set_movemap_factory(mmf)
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
    else:
        mmf.set_cartesian(setting=False)
        relax.set_movemap_factory(mmf)
        scorefxn = get_fa_scorefxn()
        relax.set_scorefxn(scorefxn)
    
    relax.constrain_relax_to_start_coords(constrain_to_input)

    print(relax)
    relax.apply(input_pose_to_relax)
    
    return input_pose_to_relax


    




def build_disulfide(input_pose_for_disulfide, ds_res1, ds_res2, cycles=1000, distance_cutoff = 5):

    input_pose_for_disulfide = pose_from_pdb(input_pose_for_disulfide)
    
   
    def resi_within(pose, angstrom, the_sun):
        
        def resi_atom_coords(res_in_pose): # input = the pose.residue(number) object
            xyz_coords_list = [float(i) for i in str(res_in_pose.atoms())[31:-1].split(',')] # All atom coords in a residue
            atoms_coords_list = [xyz_coords_list[j:j+3] for j in range(0, len(xyz_coords_list), 3)] # List of lists. xyz of individual residue
            return atoms_coords_list
        
        def distance_between_two_xyz_coords(a,b):
            distance = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
            return distance
        
        resi_within = []
        center = resi_atom_coords(pose.residue(the_sun))
        for resi in range(1, pose.total_residue()+1): # loop over protein
            planets = resi_atom_coords(pose.residue(resi))
            for coords in center:
                for other_coords in planets:
                    distance = distance_between_two_xyz_coords(coords,other_coords)
                    if distance <= angstrom:
                        resi_within.append(resi)
                    
        return sorted(list(set(resi_within)))[0:-1]
    
    
    def disulfide_SG_distance(pose_for_ds_dist,ds_residue_1_for_ds_dist,ds_residue_2_for_ds_dist):
        ds_dist_initial_vector_for_ds_dist = pose_for_ds_dist.residue(ds_residue_1_for_ds_dist).atom('SG').xyz() - pose_for_ds_dist.residue(ds_residue_2_for_ds_dist).atom('SG').xyz()
        ds_dist_initial_for_ds_dist = ds_dist_initial_vector_for_ds_dist.norm()
        return ds_dist_initial_for_ds_dist
        
    def delta_disulfide_SG_distance(input_disulfide_SG_distance, ideal_disulfide_SG_distance):
        disulfide_SG_distance_value = abs(ideal_disulfide_SG_distance - input_disulfide_SG_distance)
        return disulfide_SG_distance_value


    
    # build 2 cys residues
    scorefxn = get_fa_scorefxn()
    #scorefxn.set_weight(ScoreType.dslf_fa13, 1000.0) # to encourage disulfide formation, default values is 1.250
    print(scorefxn)
    DisulfidizeMover().make_disulfide(input_pose_for_disulfide, ds_res1, ds_res2, True, scorefxn)


    # minimize distance between SG atoms using backrub and minimization movers
    
    mmap = MoveMap()
    backrub = BackrubMover()
    min_mover = MinMover()
    min_mover.score_function(scorefxn)


    reference_pose = input_pose_for_disulfide.clone()
    working_pose = reference_pose.clone()

    print("the STARTING reference difference is " + str(delta_disulfide_SG_distance(disulfide_SG_distance(reference_pose, ds_res1, ds_res2), 2.04)))


    for n in range(0,cycles):

        print("START OF LOOP ROUND " + str(n))
        
        delta_reference_distance = delta_disulfide_SG_distance(disulfide_SG_distance(reference_pose, ds_res1, ds_res2), 2.04)
        
        if  delta_reference_distance > 0.01:
            
            print("before backrub ds distance is " + str(disulfide_SG_distance(reference_pose, ds_res1, ds_res2)))
            
            for i in resi_within(working_pose, distance_cutoff, ds_res1):
                mmap.set_bb(i, True)
                mmap.set_chi(i,True)
                mmap.set_jump(i,True)
            backrub.set_movemap(mmap)
            backrub.apply(working_pose)
            mmap.set_bb(False)
            mmap.set_chi(False)
            mmap.set_jump(False)
 
            for j in resi_within(working_pose, distance_cutoff, ds_res2):
                mmap.set_bb(j, True)
                mmap.set_chi(j,True)
                mmap.set_jump(j,True)
            backrub.set_movemap(mmap)
            backrub.apply(working_pose)
            mmap.set_bb(False)
            mmap.set_chi(False)
            mmap.set_jump(False)
            
            print("after backrub ds distance is " + str(disulfide_SG_distance(working_pose, ds_res1, ds_res2)))
 
            delta_working_distance = delta_disulfide_SG_distance(disulfide_SG_distance(working_pose, ds_res1, ds_res2), 2.04)
            print("working difference is " + str(delta_working_distance))
  
 
            if delta_working_distance <= 0.01: 
                for i in resi_within(working_pose, distance_cutoff, ds_res1):
                    mmap.set_bb(i, True)
                    mmap.set_chi(i,True)
                    mmap.set_jump(i,True)
                min_mover.movemap(mmap)
                min_mover.apply(working_pose)
                mmap.set_bb(False)
                mmap.set_chi(False)
                mmap.set_jump(False)
                for j in resi_within(working_pose, distance_cutoff, ds_res2):
                    mmap.set_bb(j, True)
                    mmap.set_chi(j,True)
                    mmap.set_jump(j,True)
                min_mover.movemap(mmap)
                min_mover.apply(working_pose)
                mmap.set_bb(False)
                mmap.set_chi(False)
                mmap.set_jump(False)
                
                print("after minimization ds distance is " + str(disulfide_SG_distance(working_pose, ds_res1, ds_res2)))
                reference_pose = working_pose.clone()

                
            else:
                if delta_working_distance < delta_reference_distance:
                    reference_pose = working_pose.clone()
                else:
                    working_pose = reference_pose.clone()
 

 
        else:
            print("SG distance is ~2.04 A")
            break
    
    print("")
    print("FINISHED")
    print("INITIAL Di-S distance was " + str(disulfide_SG_distance(input_pose, ds_res1, ds_res2)))
    print(" FINAL  Di-S distance  is " + str(disulfide_SG_distance(working_pose, ds_res1, ds_res2)))
    print("")
    
    if 2.03 <= disulfide_SG_distance(working_pose, ds_res1, ds_res2) <= 2.05:
        working_pose.dump_pdb(input_pdb.replace("pdb", "") + str(ds_res1) + "_" + str(ds_res2) + "_DISULFIDE.pdb")
    elif 2.05 < disulfide_SG_distance(working_pose, ds_res1, ds_res2) < 2.50 or disulfide_SG_distance(working_pose, ds_res1, ds_res2) < 2.03:
        working_pose.dump_pdb(input_pdb.replace("pdb", "") + str(ds_res1) + "_" + str(ds_res2) + "_MAYBE_disulfide.pdb")
    else:
        print("NO NEW DISULFIDES COULD BE MODELLED...")
        working_pose.dump_pdb(input_pdb.replace("pdb", "")+ str(ds_res1) + "_" + str(ds_res2) + "_NO_disulfide.pdb")
    
    print("")
    print("FIN")
    print("")





def point_mutation(pose_to_mutate, pdb_chain_to_mutate, pdb_position_to_mutate, mutant_aa, repack_radius=10, scorefxn=pyrosetta.create_score_function("ref2015_cart.wts")):
    # Mutate the residue at the given position
    pose_position_to_mutate = pose_to_mutate.pdb_info().pdb2pose(pdb_chain_to_mutate, pdb_position_to_mutate)
    pyrosetta.toolbox.mutants.mutate_residue(pose_to_mutate, pose_position_to_mutate, mutant_aa)

    # Create residue selector for the mutated residue and its neighbors within repack_radius
    center_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(str(pose_position_to_mutate))
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(center_residue_selector, repack_radius, True)
    

    # Create a packer task and restrict it to repacking only the mutated residue and its neighbors
    task_factory = pyrosetta.rosetta.core.pack.task.TaskFactory()
    task_factory.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    task_factory.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    task_factory.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    task_factory.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    task_factory.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), nbr_selector, True))
    
    # Repack the mutated residue and its neighbors
    pack_rotamers_mover = PackRotamersMover()
    pack_rotamers_mover.nloop(10)
    pack_rotamers_mover.score_function(scorefxn)
    pack_rotamers_mover.task_factory(task_factory)
    pack_rotamers_mover.apply(pose_to_mutate)
    
    # Create a movemap to enable minimization of only the mutated residue and its neighbors
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(False)
    move_map.set_jump(False)
    for res in pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose_to_mutate)):
        move_map.set_bb(res, True)
        move_map.set_chi(res, True)
        move_map.set_jump(res, True)
    
    # Minimize only the mutated residue and its neighbors
    min_mover = MinMover()
    min_mover.cartesian(True)
    min_mover.score_function(scorefxn)
    min_mover.movemap(move_map)
    min_mover.tolerance(0.00001)
    min_mover.min_type('lbfgs_armijo_nonmonotone')

    min_mover.apply(pose_to_mutate)
    
    return pose_to_mutate




def simple_fold(input_pose):
    
    # This is a work in progress!

    
    scorefxn = get_fa_scorefxn() 
    
    kT = 1.0
    
    mc = MonteCarlo(input_pose, scorefxn, kT)
    
    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(True)
    mmap.set_jump(True)
    
    min_mover = MinMover()
    min_mover.score_function(scorefxn)
    
    small_mover = SmallMover(mmap, kT, 5)
    
    ncycles=20
    
    # run simulatioin
    
    for i in range(1, ncycles):
        small_mover.apply(input_pose)
        min_mover.movemap(mmap)
        min_mover.apply(input_pose)
        mc.boltzmann(input_pose)
        print(scorefxn.score(input_pose))
        
    mc.recover_low(input_pose)
    
    return input_pose




def fast_relax_subset(pose, residue_indices, constrain_to_input=False, cycles=5, cartesian=True):
    """
    Perform fast relaxation on a subset of residues in a PyRosetta pose.

    Parameters:
        pose (pyrosetta.Pose): The PyRosetta pose object.
        residue_indices (list): List of residue indices to relax.

    Returns:
        pyrosetta.Pose: The relaxed pose.
    """

    # Create a MoveMap to control which residues are allowed to move
    movemap = rosetta.core.kinematics.MoveMap()
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    
    for i in range(1, pose.total_residue() + 1):
        if i in residue_indices:
            # Allow backbone and side-chain movement for the specified residue
            movemap.set_bb(i, True)
            movemap.set_chi(i, True)
            movemap.set_jump(i, True)
            movemap.set_nu(i, True) # this only applies to glycans
            movemap.set_branches(i, True) # this only applies to glycans
        else:
            # Prevent movement for other residues
            movemap.set_bb(i, False)
            movemap.set_chi(i, False)
            movemap.set_jump(i, False)
            movemap.set_nu(i, False) # this only applies to glycans
            movemap.set_branches(i, False) # this only applies to glycans


    # Create a FastRelax object
    relax = rosetta.protocols.relax.FastRelax(standard_repeats=cycles)

    # Set FastRelax options
    if cartesian == True:
        mmf.set_cartesian(setting=True)
        scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        relax.set_scorefxn(scorefxn)
        relax.set_movemap_factory(mmf)
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
    else:
        mmf.set_cartesian(setting=False)
        relax.set_movemap_factory(mmf)
        scorefxn = get_fa_scorefxn()
        relax.set_scorefxn(scorefxn)
    
    relax.constrain_relax_to_start_coords(constrain_to_input)
    
    # Apply MoveMap to FastRelax
    relax.set_movemap(movemap)

    # Perform relaxation
    relax.apply(pose)

    return pose








def make_pdb2pose_list(pose, pdb_chains, pdb_residue_indices): 
    """
    Convert PDB residue indices to corresponding PyRosetta residue indices.
    
    Parameters:
        pose (pyrosetta.Pose): The PyRosetta pose object.
        pdb_chains (list): A list of PDB chain identifiers.
        pdb_residue_indices (list): A list of residue indices within the PDB chains.
    
    Returns:
        list: A list of PyRosetta residue indices corresponding to the input PDB indices.
    """
    pyrosetta_residue_indices = []
    for pdb_chain in pdb_chains: 
        for pdb_residue_index in pdb_residue_indices:
            pyrosetta_index = pose.pdb_info().pdb2pose(pdb_chain, pdb_residue_index)
            #print(f'{pdb_chain},{pdb_residue_index},{pyrosetta_index}')
            pyrosetta_residue_indices.append(pyrosetta_index)
    return pyrosetta_residue_indices



def count_clashes(pose, selection1, selection2=None, overlap=0.4):
    """
    Count the number of clashes in a given pose
    [1] either within a selection 
    [2] or between two selections of residues,
    based on a 0.4 Å overlap of van der Waals radii.
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose): The pose object containing the protein structure.
        selection1 (list): List of pyrosetta residue indices.
        selection2 (list; OPTIONAL): List of pyrosetta residue indices for the second selection.
        overlap (float): The overlap distance for bad clashes (default is 0.4 Å).
        
    Returns:
        int: The number of (bad) clashes.
    """
    
    vdw_radii = {
    "H": 1.2,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    # Add more elements as needed
    }
    
    def count_atom_clashes(residue1, residue2):
        clashes = 0
        for atom1 in range(1, residue1.natoms() + 1):
            for atom2 in range(1, residue2.natoms() + 1):
                atom1_xyz = residue1.xyz(atom1)
                atom2_xyz = residue2.xyz(atom2)
                element1 = residue1.atom_type(atom1).element()
                element2 = residue2.atom_type(atom2).element()
                vdw_radius1 = vdw_radii.get(element1, 1.7) # Default to 1.7 Å (C) if element not found
                vdw_radius2 = vdw_radii.get(element2, 1.7) # Default to 1.7 Å (C) if element not found
                distance = (atom1_xyz - atom2_xyz).norm() # same as sqrt((ax-bx)^2 + (ay-by)^2 + (az-bz)^2)
                clash_threshold = vdw_radius1 + vdw_radius2 - overlap
                if distance < clash_threshold:
                    clashes += 1
        return clashes

    bad_clashes = 0
    if selection2 is None:
        for i, res1 in enumerate(selection1):
            for res2 in selection1[i + 1:]:
                bad_clashes += count_atom_clashes(pose.residue(res1), pose.residue(res2))
    else:
        for res1 in selection1:
            for res2 in selection2:
                if res1 != res2:
                    bad_clashes += count_atom_clashes(pose.residue(res1), pose.residue(res2))

    return bad_clashes













### simple threading of a sequences into a pose and pack mutant rotamers
def thread_and_pack(pose, sequence):
    """
    Thread a sequence onto a pose and pack rotamers of mutations.
    """
    three_letter = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN',
                    'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                    'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET',
                    'G': 'GLY', 'P': 'PRO', 'C': 'CYS', 'Y': 'TYR', 'A': 'ALA'}

    if len(sequence) != pose.total_residue():
        raise ValueError("Sequence length must match the pose residue count")

    mutator = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
    mutated_residues = []

    for i in range(1, pose.total_residue() + 1):  # PyRosetta is 1-indexed
        current_res = pose.residue(i).name1()
        target_res = sequence[i - 1]
        target_res_3 = three_letter[target_res]

        if current_res != target_res:
            print(f"Mutating residue {i}: {current_res} -> {target_res}")
            mutator.set_res_name(target_res_3)
            mutator.set_target(i)
            mutator.apply(pose)
            mutated_residues.append(i)

    # Apply packer to all mutated residues at once
    if mutated_residues:
        
        # set up task factory for only mutated residues
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

        #prevent all residues except selected from design and repacking
        nonmutated_residues = [x for x in range(1, pose.total_residue() + 1) if x not in mutated_residues]
        selected_NO_packing = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for nonmut in nonmutated_residues:
            selected_NO_packing.append_index(nonmut)
        disallow = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(disallow, selected_NO_packing, False )
        tf.push_back(prevent_subset_repacking)
        
        # allow selected residues only repacking (=switch off design)
        selected_packing = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for mut in mutated_residues:
            selected_packing.append_index(mut)
        print(selected_packing)
        allow = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
        restrict_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(allow , selected_packing, False )
        tf.push_back(restrict_subset_repacking)

        #print(f"task factory : {tf.create_task_and_apply_taskoperations(pose)}")

        #MinPacker = pyrosetta.rosetta.protocols.minimization_packing.MinPackMover()
        Packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
        scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        Packer.task_factory(tf)        
        Packer.score_function(scorefxn)
        initial_score = scorefxn(pose)        
        Packer.apply(pose)
        final_score = scorefxn(pose)

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(f" starting score : {initial_score} REU")
        print(f"    final score : {final_score} REU")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")















