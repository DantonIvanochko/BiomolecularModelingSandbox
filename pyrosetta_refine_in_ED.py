# pyrosetta
from pyrosetta import *
pyrosetta.init('-include_sugars -write_pdb_link_records -ignore_zero_occupancy false') # within parenthesis is required for working with glycans ('-include_sugars -write_pdb_link_records')
# pyrosetta - clean
from pyrosetta.toolbox import cleanATOM
# pyrosetta - relax/minimize
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from rosetta.core.scoring import ScoreType
# pyrosetta - ss
from pyrosetta.rosetta.core.simple_metrics.metrics import SecondaryStructureMetric

# misc
import os
import math
import random
from glob import glob






def clean_structure(input_pdb):
    print("cleaning structure")
    #cleanATOM(input_pdb, input_pdb.replace("..", "").replace(".pdb", "") + ".clean.pdb")
    cleanATOM(input_pdb, input_pdb.replace(".pdb", "") + ".clean.pdb")










def fast_relax_ED(pose, map_file, weight):

    """
    map_file must be .ccp4 or .mrc
    weight should be large (>30) for higher res
    """
    # score function
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

    # set electron density info
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_file)
    initial_fit = ED.matchPose(pose)

    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)

    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxn.set_weight(elec_dens_fast, weight)
    # Create a MoveMap to control which residues are allowed to move
    movemap = rosetta.core.kinematics.MoveMap()

    for i in range(1, pose.total_residue() + 1):
        # allow backbone + side-chain movement
        movemap.set_bb(i, True)
        movemap.set_chi(i, True)
        movemap.set_jump(i, True)

    # Create a FastRelax object
    relax = rosetta.protocols.relax.FastRelax(standard_repeats=1)
    relax.set_scorefxn(scorefxn)
    relax.set_movemap(movemap)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)

    # Perform relaxation
    relax.apply(pose)
    final_fit = ED.matchPose(pose)

    # RESULTS
    print("\n\nRESULTS -- FASTRELAX IN ELECTRON DENSITY")
    print(f"map_file : {map_file} weight : {weight}")
    print(f"initial_fit : {initial_fit}")
    print(f"final_fit   : {final_fit}")
    print("RESULTS END\n\n")

    return pose








def min_mover_ED(pose, map_file, weight):

    """
    map_file must be .ccp4 or .mrc
    weight should be large (>30) for higher res
    residue_indices_BB_SC is for backbone and sidechain movements
    residue_indices_BB_SC is for only sidechain movements

    """

    # set electron density info
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_file)
    initial_fit = ED.matchPose(pose)

    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)

    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    scorefxn.set_weight(elec_dens_fast, weight)


    # Create a MoveMap to control which residues are allowed to move
    movemap = rosetta.core.kinematics.MoveMap()
    for i in range(1, pose.total_residue() + 1):
        # allow backbone + side-chain movement
        movemap.set_bb(i, True)
        movemap.set_chi(i, True)
        movemap.set_jump(i, True)

    # create a MinMover 
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.cartesian(True)
    minmover.score_function(scorefxn)
    minmover.set_movemap(movemap)

    # minimize
    minmover.apply(pose)
    final_fit = ED.matchPose(pose)

    # RESULTS
    print("\n\nRESULTS -- MINMOVER IN ELECTRON DENSITY")
    print(f"map_file : {map_file} weight : {weight}")
    print(f"initial_fit : {initial_fit}")
    print(f"final_fit   : {final_fit}")
    print("RESULTS END\n\n")

    return pose







def min_mover_ED_2subsets(pose, map_file, weight, residue_indices_BB_SC, residue_indices_SC):

    """
    map_file must be .ccp4 or .mrc
    weight should be large (>30) for higher res
    residue_indices_BB_SC is for backbone and sidechain movements
    residue_indices_BB_SC is for only sidechain movements

    """

    # set electron density info
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_file)
    initial_fit = ED.matchPose(pose)

    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)

    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    scorefxn.set_weight(elec_dens_fast, weight)


    # Create a MoveMap to control which residues are allowed to move
    movemap = rosetta.core.kinematics.MoveMap()
    for i in range(1, pose.total_residue() + 1):
        if i in residue_indices_BB_SC:
            # allow backbone + side-chain movement
            movemap.set_bb(i, True)
            movemap.set_chi(i, True)
            movemap.set_jump(i, True)
        elif i in residue_indices_SC:
            # allow side-chain movement
            movemap.set_bb(i, False)
            movemap.set_chi(i, True)
            movemap.set_jump(i, True)
        else:
            # no movement for other residues
            movemap.set_bb(i, False)
            movemap.set_chi(i, False)
            movemap.set_jump(i, False)

    # Create a Min mover object
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.cartesian(True)
    minmover.score_function(scorefxn)
    minmover.set_movemap(movemap)

    # Perform relaxation
    minmover.apply(pose)
    final_fit = ED.matchPose(pose)

    # RESULTS
    print("\n\nRESULTS -- MINMOVER IN ELECTRON DENSITY")
    print(f"map_file : {map_file} weight : {weight}")
    print(f"initial_fit : {initial_fit}")
    print(f"final_fit   : {final_fit}")
    print("RESULTS END\n\n")

    return pose







def fast_relax_ED_2subsets(pose, map_file, weight, residue_indices_BB_SC, residue_indices_SC):

    """
    map_file must be .ccp4 or .mrc
    weight should be large (>30) for higher res
    residue_indices_BB_SC is for backbone and sidechain movements
    residue_indices_BB_SC is for only sidechain movements
    """

    # score function
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

    # set electron density info
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_file)
    initial_fit = ED.matchPose(pose)

    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)

    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxn.set_weight(elec_dens_fast, weight)

    # Create a MoveMap to control which residues are allowed to move
    movemap = rosetta.core.kinematics.MoveMap()

    for i in range(1, pose.total_residue() + 1):
        if i in residue_indices_BB_SC:
            # allow backbone + side-chain movement
            movemap.set_bb(i, True)
            movemap.set_chi(i, True)
            movemap.set_jump(i, True)
        elif i in residue_indices_SC:
            # allow side-chain movement
            movemap.set_bb(i, False)
            movemap.set_chi(i, True)
            movemap.set_jump(i, True)
        else:
            # no movement for other residues
            movemap.set_bb(i, False)
            movemap.set_chi(i, False)
            movemap.set_jump(i, False)

    # Create a FastRelax object
    relax = rosetta.protocols.relax.FastRelax(standard_repeats=1)
    relax.set_scorefxn(scorefxn)
    relax.set_movemap(movemap)
    relax.cartesian(True)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)

    # Perform relaxation
    relax.apply(pose)
    final_fit = ED.matchPose(pose)

    # RESULTS
    print("\n\nRESULTS -- FASTRELAX IN ELECTRON DENSITY")
    print(f"map_file : {map_file} weight : {weight}")
    print(f"initial_fit : {initial_fit}")
    print(f"final_fit   : {final_fit}")
    print("RESULTS END\n\n")

    return pose











# get a list from specific PDB formated chain and indices
def make_pdb2pose_list(pose, pdb_chains, pdb_residue_indices):
    """
    Convert PDB residue indices to list of rosetta residue indices.
    inputs : pose, [pdb_chains], [pdb_residue_indices]
    outputs : list of pyrosetta residue indices corresponding to the input pdb indices.
    """
    pyrosetta_residue_indices = []
    for pdb_chain in pdb_chains:
        for pdb_residue_index in pdb_residue_indices:
            pyrosetta_index = pose.pdb_info().pdb2pose(pdb_chain, pdb_residue_index)
            pyrosetta_residue_indices.append(pyrosetta_index)
    return pyrosetta_residue_indices






# get list of lists indicating ss in pose
def get_secondary_structure_regions(input_pose):
    # use dssp to assign ss to pose
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    DSSP.apply(input_pose)
    # set ss string
    ss_string = input_pose.secstruct()
    # lisst to hold residues for each secondary structure element
    ss_helices, ss_sheets, ss_loops, ss_unassigned = [], [], [], []

    # loop through DSSP assignments and store residue ranges
    for i, ss in enumerate(ss_string):
        if ss == "H":
            ss_helices.append(i + 1)  # because rosetta is base-1 indexed
        elif ss == "E":
            ss_sheets.append(i + 1)
        elif ss == "L":
            ss_loops.append(i + 1)
        else:
            ss_unassigned.append(i + 1)

    return ss_helices, ss_sheets, ss_loops, ss_unassigned




























# EXMAPLE BELOW

def main():

    print("MAIN")

    input_pdb_list = ["D1D6-RUPA97-LMIV01-2A2-20250106_real_space_refined_001.pdb"]

    for input_pdb_unclean in input_pdb_list:

        ccp4_map = "../J684_004_volume_map_sharp.mrc"

        clean_structure(input_pdb_unclean)

        input_pdb_cleaned = input_pdb_unclean.replace(".pdb", "") + ".clean.pdb"

        input_pose = pose_from_pdb(input_pdb_cleaned)

        relax_pose = input_pose.clone()

        wt = 50

        print("\n\n    relaxing pose -- 1\n")
        fast_relax_ED(relax_pose, ccp4_map, wt)
        print("\n\n    saving fast relaxed pose\n")
        relax_pose.dump_pdb(f"D1D6-RUPA97-LMIV01-2A2-3pt5-20250106.fr.{wt}.1.pdb")
        print("\n\n    relaxing pose -- 2\n")
        fast_relax_ED(relax_pose, ccp4_map, wt)
        print("\n\n    saving fast relaxed pose\n")
        relax_pose.dump_pdb(f"D1D6-RUPA97-LMIV01-2A2-3pt5-20250106.fr.{wt}.2.pdb")
        print("\n\n    relaxing pose -- 3\n")
        fast_relax_ED(relax_pose, ccp4_map, wt)




if __name__ == "__main__":
    main()
