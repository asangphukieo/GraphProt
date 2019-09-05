import re
import sys
from Bio import SeqIO

pattern = 'GGA[ACGT]{4,70}GGA[ACGT]{2,12}'
#pattern = 'GGA'
fasta_file = sys.argv[1]
output_file = sys.argv[2]

w2file=open(output_file,'w')
for seq_rec in SeqIO.parse(fasta_file,"fasta"):
	test_string=str(seq_rec.seq)
	#print(test_string)
	result = re.search(pattern, test_string)
	if result:
		#print(str(result))
		#print(result.group())
		#print("Pattern found!")
		w2file.write(seq_rec.description+'\t'+'Y'+'\n')
		#selected_seq.write('>'+seq_rec.description+'\n'+str(seq_rec.seq)+'\n')
	else:
		#print("Pattern not found!")	
		w2file.write(seq_rec.description+'\t'+'N'+'\n')

w2file.close()

##usage ##
#python GGA_regex.py input_fasta output_name
