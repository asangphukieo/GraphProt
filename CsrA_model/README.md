CsrA GraphProt model 

## Training data ##
Trained from CLIP-seq data from 2 different sources
1. Holmqvist, Erik, et al. "Global RNA recognition patterns of post‐transcriptional regulators Hfq and CsrA revealed by UV crosslinking in vivo." The EMBO journal 35.9 (2016): 991-1011.
2. Potts, Anastasia H., et al. "Global role of the bacterial post-transcriptional regulator CsrA revealed by integrated transcriptomics." Nature communications 8.1 (2017): 1596.

(Salmonella_Typhimurium and Escherichia_coli)

## Data preprocess ##
1. download FASTA files of full genome
	#01 https://www.ebi.ac.uk/ena/data/view/FQ312003
	#02 https://www.ncbi.nlm.nih.gov/genome/?term=txid511145[Organism:exp]

2. extract coordinate from CLIP-seq files and remove redundance
	awk -F',' '{print $3,$4,$5}' 01_EMBJ-35-991-s008.csv |sort -u > 01_coordinate.coo
	awk -F',' '{print $10,$11,$8}' 02_41467_2017_1613_MOESM3_ESM.csv |sort -u > 02_coordinate.coo

3. extract FASTA sequence using coordinate files (below script require biopython)
	## requirement [1] Genome in FASTA, [2] Start,Stop position,Strand of CLIP peak region (in 3 columns separated by tab) [3] include 'Franking' or 'NotFranking' region , [4] determine total length of sequence if use 'Franking' , [5] 'Random' or 'NoRandom' start and stop position (keep seq length identical to original coordinate) , [6] output file , [7] excluding longer sequences (CLIP peak region)

	python extract_sequence_by_coordinate2.py FQ312003.fasta 01_coordinate.coo 'Franking' 150 'NoRandom' 01_pos.fasta 100
	python extract_sequence_by_coordinate2.py GCF_000005845.2_ASM584v2_genomic.fna 02_coordinate.coo 'Franking' 150 'NoRandom' 02_pos.fasta 100

4. extract FASTA sequence by random coordinate but keep sequence length similar to the original coordinate file
	python extract_sequence_by_coordinate2.py FQ312003.fasta 01_coordinate.coo 'Franking' 150 'Random' 01_neg.fasta 100
	python extract_sequence_by_coordinate2.py GCF_000005845.2_ASM584v2_genomic.fna 02_coordinate.coo 'Franking' 150 'Random' 02_neg.fasta 100

5. sequence combination
cat 01_pos.fasta 02_pos.fasta > csrA_pos.fasta
cat 01_neg.fasta 02_neg.fasta > csrA_neg.fasta

## Model performance ##
1. use 10% of data to tune parameter
2. use 10-fold cross validation to observe model performance on 90% of data
	ROC	0.7976
	APR	0.7888

## Use model to predict CsrA in Bacillus Subtilis ##
1. run trained model [csrA.model] in folder GraphProtCsrA

GraphProt.pl \
--action predict \
--model  GraphProtCsrA/csrA.model \
--fasta    GraphProtCsrA/NC_000964_upfromstartpos_200_down_100.fa \
--prefix   GraphProtCsrA/csrA_B_subtilis \
--params   GraphProtCsrA/csrA.params


