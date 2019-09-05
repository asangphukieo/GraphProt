from Bio import SeqIO
import sys
import random as rd

#use conda of build_tree
#allow to get franking region

## requirement 
#[1] Genome in FASTA
#[2] Start,Stop position,Strand of CLIP peak region (in 3 columns separated by tab) 
#[3] include 'Franking' or 'NotFranking' region 
#[4] determine total length of sequence if use 'Franking' 
#[5] 'Random' or 'NoRandom' start and stop position (keep seq length identical to original coordinate) 
#[6] output file 
#[7] excluding longer sequences (CLIP peak region)

#usage
#python extract_sequence_by_coordinate2.py FQ312003.fasta 01_coordinate.coo 'Franking' 150 'NoRandom' 01_pos.fasta 100

inFile = open(sys.argv[1],'r') #input fasta file
coordinate=sys.argv[2]
franking=sys.argv[3]
total_len=int(sys.argv[4])  #including franking region
random=sys.argv[5]
fw=open(sys.argv[6],'w') #output fasta file
length_filter=int(sys.argv[7]) #filter out very long binding site ,default 75 folling paper

for record in SeqIO.parse(inFile,'fasta'):
	count=1
	for line in open(coordinate): #input coordinate file form blastout

		line= line.rstrip()
		col=line.split()
		start=int(col[0])-1
		end=int(col[1])
		strandS=str(col[2])

		if (end - start) <length_filter:

			if random == 'Random':
				org_len=end-start
				start=rd.randrange(1, len(record))
				end=start+(org_len)

			if franking =='Franking':
				
				franking_region=int((total_len-(end-start)) /2)
				#print(franking_region)
				F_start=start-franking_region
				F_end=end+franking_region	
					
				fw.write(">seq_" + str(count) +'_'+str(F_start)+':'+str(F_end)+';'+str(strandS) +"\n")
				#print(record.seq)
				if strandS=='+':
					fww = (str(record.seq[F_start:start]).lower()+str(record.seq[start:end])+str(record.seq[end:F_end]).lower() + '\n')
				else:
					fww = (str(record.seq[F_start:start].reverse_complement()).lower() + str(record.seq[start:end].reverse_complement())+str(record.seq[end:F_end].reverse_complement()).lower() + '\n')

			else:
				fw.write(">seq_" + str(count) +'_'+str(start)+':'+str(end)+';'+str(strandS) +"\n")
				#print(record.seq)
				if strandS=='+':
					fww = (str(record.seq[start:end]) + '\n')
				else:
					fww = (str(record.seq[start:end].reverse_complement()) + '\n')
			
			fw.write(fww)
			count+=1
fw.close()
inFile.close()

