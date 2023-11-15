# modules


# pymol
import pymol
pymol.pymol_argv = ['pymol','-qc']
pymol.finish_launching()
from pymol import cmd, stored, util
# misc
import os
import math


    


def Cbeta_distance(input_pdb, res_range1, res_range2, fluctuation=7):

    # Cb-Cb distance cutoff 
    cutoff = fluctuation # better value TBD

    # prepare a csv file for results
    Cbeta_csv = open(str(input_pdb)+'.'+str(cutoff) + 'A_cutoff_Cbeta_distances.csv','w')
    Cbeta_csv.write('model,chain_1,residue_1,index_1,chain_2,residue_2,index_2,CB-CB_distance \n') # header 
    
    # get interface residues (this allows a distance cutoff so all interdomain Cbeta distance do not have to be measured)
    
    print("")
    print("Working on " + str(input_pdb))
    cmd.load(input_pdb)
    input_pdb_object = cmd.get_names()
    
    print("")
    print("not including residues which are already forming disulfide bonds")
    stored.list_of_native_disulfides=[]
    cmd.select("native_disulfides", " resn CYS and name SG bound_to resn CYS and name SG")
    cmd.iterate("native_disulfides","stored.list_of_native_disulfides.append((model,chain,segi,resi,resn))")
    stored.list_of_native_disulfides = sorted(set(stored.list_of_native_disulfides))
    print(stored.list_of_native_disulfides)
    
    print("")
    print("1st interface")
    stored.list_of_1st_interface_residues=[]
    cmd.select("1st_interface", "br. {} and  {} within {} of {} and  {}".format(input_pdb_object[0],res_range1,cutoff,input_pdb_object[0],res_range2))
    cmd.iterate("1st_interface","stored.list_of_1st_interface_residues.append((model,chain,segi,resi,resn))")
    stored.list_of_1st_interface_residues = sorted(set(stored.list_of_1st_interface_residues))
    stored.list_of_1st_interface_residues = [x for x in stored.list_of_1st_interface_residues if x not in stored.list_of_native_disulfides] # dont count disulfided cys
    stored.list_of_1st_interface_residues = [ x for x in stored.list_of_1st_interface_residues if "GLY" not in x[-1] ] # dont count gly
    stored.list_of_1st_interface_residues = [ x for x in stored.list_of_1st_interface_residues if "PRO" not in x[-1] ] # dont count pro
    print(stored.list_of_1st_interface_residues)
    
    print("")
    print("2nd interface")
    stored.list_of_2nd_interface_residues=[]
    cmd.select("2nd_interface", "br. {} and  {} within {} of {} and  {}".format(input_pdb_object[0],res_range2,cutoff,input_pdb_object[0],res_range1))
    cmd.iterate("2nd_interface","stored.list_of_2nd_interface_residues.append((model,chain,segi,resi,resn))")
    stored.list_of_2nd_interface_residues = sorted(set(stored.list_of_2nd_interface_residues))
    stored.list_of_2nd_interface_residues = [x for x in stored.list_of_2nd_interface_residues if x not in stored.list_of_native_disulfides] # dont count disulfided cys
    stored.list_of_2nd_interface_residues = [ x for x in stored.list_of_2nd_interface_residues if "GLY" not in x[-1] ] # dont count gly
    stored.list_of_2nd_interface_residues = [ x for x in stored.list_of_2nd_interface_residues if "PRO" not in x[-1] ] # dont count pro
    print(stored.list_of_2nd_interface_residues)
    
    # get distance of each Cbeta from interface 1 with each Cbeta from interface 2
    
    for i in range(0, len(stored.list_of_1st_interface_residues)):
        for j in range(0, len(stored.list_of_2nd_interface_residues)):
            SGSGdist = cmd.get_distance(
            "chain {} and segi {} and resi {} and name CB".format(stored.list_of_1st_interface_residues[i][1],stored.list_of_1st_interface_residues[i][2],stored.list_of_1st_interface_residues[i][3]),
            "chain {} and segi {} and resi {} and name CB".format(stored.list_of_2nd_interface_residues[j][1],stored.list_of_2nd_interface_residues[j][2],stored.list_of_2nd_interface_residues[j][3]))
            if SGSGdist <= cutoff:
                Cbeta_csv.write("{},{},{},{},{},{},{},{} \n".format(
                stored.list_of_1st_interface_residues[i][0],
                stored.list_of_1st_interface_residues[i][1],
                stored.list_of_1st_interface_residues[i][4],
                stored.list_of_1st_interface_residues[i][3],
                stored.list_of_2nd_interface_residues[j][1],
                stored.list_of_2nd_interface_residues[j][4],
                stored.list_of_2nd_interface_residues[j][3],
                SGSGdist))
            else:
                continue
    
    # delete object
    cmd.delete('all')
    print("")
    print("FIN")
    print("")





def average_b(selection):
	stored.tempfactor = 0
	stored.atomnumber = 0
	cmd.iterate(selection, "stored.tempfactor = stored.tempfactor + b")
	cmd.iterate(selection, "stored.atomnumber = stored.atomnumber + 1")
	print("Your selection: %s" % selection)
	print("sum of B factors: %s" % stored.tempfactor)
	print("number of atoms: %s" % stored.atomnumber)
	averagetempfactor = stored.tempfactor / stored.atomnumber
	print("average B of '%s': %s" % (selection, averagetempfactor))
cmd.extend("average_b", average_b)



def pLDDT(selection):
	stored.atom_pLDDT = 0
	stored.atom_number = 0
	cmd.iterate(selection, "stored.atom_pLDDT = stored.atom_pLDDT + b")
	cmd.iterate(selection, "stored.atom_number = stored.atom_number + 1")
	average_pLDDT = stored.atom_pLDDT / stored.atom_number
	print(f"pLDDT of {selection} : {average_pLDDT}")
	return average_pLDDT
cmd.extend("pLDDT", pLDDT)




def RMSD(query_structure, reference_structure):
        rmsd = cmd.align(f"{query_structure} and n. CA+C+CN", f"{reference_structure} and n. CA+C+CN", cycles=0)[0]
        print(f"backbone RMSD of {query_structure} and {reference_structure} : {rmsd}")
        return rmsd
cmd.extend("RMSD", RMSD)

def count_interchain_clashes(chain1, chain2, clash_cutoff=0.4):

        """
        Compute the number of bad clashes between two selections.
        Bad clashes (MolProbity) occur when the van der Waals radii overlap by ≥ 0.4 Å.
        """

        object_sele1 = pymol.cmd.get_names(selection=chain1)[0]
        object_sele2 = pymol.cmd.get_names(selection=chain2)[0]

        atom_ids_sele1 = pymol.cmd.identify(chain1)
        atom_ids_sele2 = pymol.cmd.identify(chain2)

        atom_ids_sele1 = "+".join(str(x) for x in atom_ids_sele1)
        atom_ids_sele2 = "+".join(str(x) for x in atom_ids_sele2)

        pymol.stored.list_sele1 = []
        pymol.stored.list_sele2 = []

        pymol.cmd.iterate( f"({object_sele1} and ID {atom_ids_sele1})", "stored.list_sele1.append([model, segi, chain, resn, resi, ID, name, elem, vdw])" )
        pymol.cmd.iterate( f"({object_sele2} and ID {atom_ids_sele2})", "stored.list_sele2.append([model, segi, chain, resn, resi, ID, name, elem, vdw])" )

        clash_count = 0

        for i in stored.list_sele1:
            for j in stored.list_sele2:
                clashing_distance = i[8] + j[8] - 0.4 
                measured_distance = pymol.cmd.get_distance(f"(/{i[0]}/{i[1]}/{i[2]}/{i[4]}/{i[6]})",f"(/{j[0]}/{j[1]}/{j[2]}/{j[4]}/{j[6]})")
                if measured_distance < clashing_distance:
                    clash_count += 1

        return clash_count

pymol.cmd.extend("count_interchain_clashes", count_interchain_clashes)


def get_aa_seq(pymol_object): # only one chain at a time
        fastastr = cmd.get_fastastr(pymol_object).splitlines()[1:]
        aa_seq = ''.join(fastastr)
        print(aa_seq)
	return aa_seq
cmd.extend("get_aa_seq", get_aa_seq)






