# modules

# pyrosetta
from pyrosetta import * 
pyrosetta.init('-include_sugars -write_pdb_link_records') # within parenthesis is required for working with glycans ('-include_sugars -write_pdb_link_records')
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
    cleanATOM(input_pdb, input_pdb.replace("../models/", "").replace(".pdb", "") + ".clean.pdb")



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



def relax_structure(input_pdb_to_relax, contrain_to_input=False, standard_repeats=5, cartesian=True):

    print("prepare starting structure with FastRelax()")
    poseRelaxed = pose_from_pdb(input_pdb_to_relax)
    relax = FastRelax(standard_repeats=5)
    
    if cartesian == True:
        scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
    else:
        scorefxn = get_fa_scorefxn()                
    
    relax.set_scorefxn(scorefxn)
    relax.constrain_relax_to_start_coords(contrain_to_input)

    print(relax)
    relax.apply(poseRelaxed)
    
    return poseRelaxed

    




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





def point_mutation(pose_to_mutate, position_to_mutate, mutant_aa, pack_radius, pack_scorefxn, relax_cycles=1):

    pyrosetta.toolbox.mutants.mutate_residue(pose_to_mutate, position_to_mutate, mutant_aa, pack_radius=5, pack_scorefxn=get_fa_scorefxn())
    
    focus_pos = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    focus_pos.set_index(position_to_mutate)
    nbr = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr.set_distance(pack_radius)
    nbr.set_focus_selector(focus_pos)
    nbr.set_include_focus_in_subset(True)
    pose_to_mutate.update_residue_neighbors()
    nbr.apply(pose_to_mutate)
    
    rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    
    repack_tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(rlt, nbr, True))
    
    #jump_num = mutant_pose.num_jump()
    #jump_selector = pyrosetta.rosetta.core.select.jump_selector.JumpIndexSelector()
    #jump_selector.jump(jump_num)
    
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.add_chi_action(pyrosetta.rosetta.core.select.movemap.mm_enable, nbr)
    mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.mm_enable, nbr)
    #mmf.add_jump_action(pyrosetta.rosetta.core.select.movemap.mm_enable, nbr)
    #mmf.add_jump_action(pyrosetta.rosetta.core.select.movemap.mm_enable, jump_selector)
    
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(pack_scorefxn, standard_repeats = relax_cycles)
    relax.constrain_relax_to_start_coords(True) 
    relax.set_task_factory(repack_tf)
    relax.set_movemap_factory(mmf)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    #relax.max_iter(300)
    relax.apply(pose_to_mutate)
    
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
