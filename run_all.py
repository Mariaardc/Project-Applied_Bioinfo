#!/usr/bin/env python3
# -*- coding: utf-8 -*

#MIT License

#Copyright (c) 2017 Maria Araceli Diaz Cruz


# First Part alignment

from glob import glob
for filename in glob('PDIA3_species1.fa'):
    with open(filename) as f:
	       output = str(filename)
	       output += 'a'
	       in_file = str(filename)
from Bio.Align.Applications import MafftCommandline
mafft_cline = MafftCommandline(input = in_file)
print(mafft_cline)
stdout, stderr = mafft_cline()
with open(output, 'w') as handle:
	    handle.write(stdout.upper())

# Second part: Calculate consensus sequence and position score matrix:

from Bio import AlignIO
alignment= AlignIO.read('PDIA3_species1.faa','fasta')
#print(alignment)

#Summary of the sequence

from Bio.Align import AlignInfo
summary_align = AlignInfo.SummaryInfo(alignment)
#print (summary_align)

#Consensus sequence

consensus = summary_align.dumb_consensus(ambiguous='-')
print(consensus)

#position specific score matrix (representation of the motifs in our sequence, counting nucleotides in each column)

# Matrix based in the consensus sequence but excluding N residues:

my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['N'])
with open('output.txt', 'w') as output:
    print (my_pssm, file=output) #prints to output file
    print (my_pssm)

#Specific position of nucleotide:
#print (my_pssm[1]["A"])


#Third part: Search for the conserved nucleotide sequences: FIND PATTERNS!
#Storing the variables into two lists: line in lines and nucleotides in nucleotides.

lines=[]
a =0
nucleotides=[]


f = open('sequences1.fa', 'wt')               #output file where the sequences selected are going to be stored.
with open('output.txt', 'r') as InFile:
    InFile.readline()                        #skip first line
    for line in InFile:
        a=a+1                                #line counter
        line = line.strip('\n')
        lines=line                           #store each line in lines (list)
        if '-' not in line and '5.0' in line: #condition 1

            nucleotides.append(lines[0])      # append nucleotides which are in the first column.

        else:
            if len(nucleotides)>15:           #Second condition. If the intron is shorter change the condition to less restrictive (>10).
                b=a-len(nucleotides)+1
                print('')
                print('',file=f)
                print('>',b,'-',a)
                print('>',b,'-',a,file=f)
                print(''.join(nucleotides))
                print(''.join(nucleotides),file=f)
                del nucleotides[:]                    #Delete the list when sequences are stored.
            else:

                del nucleotides[:]                    # No sequences that reach the condition. Delete the list.

#Fourth part: Verification step : Important to only use one sequence at a time which format 'string' or FASTA file.
#Hits will be extracted in a XML file in the same folder as the script.

from Bio import SeqIO
from Bio.Blast import NCBIWWW
my_query = open("sequence.fasta").read()
result_handle = NCBIWWW.qblast("blastn", "nr", my_query, word_size=7,
                                gapcosts='5 2',nucl_reward=1,nucl_penalty='-3',expect=1000)  #Parameters optimized for query short sequences
blast_result = open("my_blast.xml", "w")
blast_result.write(result_handle.read())
blast_result.close()
result_handle.close()

#We can also parse through the result_handle by creating the blast_records object: 
#from Bio.Blast import NCBIXML
#blast_records = NCBIXML.parse(result_handle)
