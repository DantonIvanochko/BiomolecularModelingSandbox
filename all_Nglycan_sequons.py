
"""
Generate all Nglycan sequons in sequence at defined positions
Nglycan sequon == Nx[S/T] where x is any residue except Proline
Inputs:
  name of sequence
  aa sequence
  list of positions
Output:
  fasta file of N glycan sequon containing sequences
"""


from itertools import product
import pandas as pd



def all_Nglycan_sequons(id, seq, position):
    sequon_N  = ["N"]
    sequon_X  = ["A","C","D","E","F","G","H","I","K","L","M","N","Q","R","S","T","V","W","Y"]
    sequon_ST = ["S","T"]
    sequons = [''.join(list(i)) for i in product(sequon_N,sequon_X,sequon_ST)]
    index_python_format = position - 1 
    sequon_dict = {}
    for n in sequons:
        sequon_dict.update({f"{id}_{position}_{n}":seq[:index_python_format] + n + seq[index_python_format + 3:]})
    sequon_df = pd.DataFrame(sequon_dict.items())
    return sequon_df



def main():

    # USE INPUTS HERE (basename, testSeq, list_of_positions)
    basename = "TEST"
    testSeq = "AAACCCDDDEEEFFFGGGHHHIIIKKKLLLMMMNNNPPPQQQ"
    list_of_positions = [4, 28]
    
    
    testSequons_list = []
    with open(f"{basename}_all_Nglycan_sequons.fasta", 'w') as the_file:
        pass
    for pos in list_of_positions:
        testSequons = all_Nglycan_sequons(basename, testSeq, pos)
        for index, row in testSequons.iterrows():
            print(f'>{row[0]}\n{row[1]}')
            with open(f"{basename}_all_Nglycan_sequons.fasta", 'a') as the_file:
                the_file.write(f'>{row[0]}\n{row[1]}\n')
    print("\n\n~END~\n\n")



if __name__ == "__main__":
    main()
