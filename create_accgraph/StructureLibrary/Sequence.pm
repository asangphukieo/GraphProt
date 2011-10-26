package StructureLibrary::Sequence;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
			calculate_sequence_contents
			divide_into_x_equal_sequences
			extract_motif_with_flanks
			find_motif
			read_fasta_file
			read_geneIDs_from_fasta_without_version_number
			word_frequencies
			write_hash_to_fasta
			);
@EXPORT_OK = qw(
			words
);

##################################################################################
##################################################################################
## Package Sequences.pm 	AUTHOR(S) = Sita Lange
## This package includes a collection of methods that are useful when handling 
## RNA and DNA sequences.
##################################################################################
##################################################################################

##################################################################################
# This method parses a standard fasta file and returns the sequences in a hash
# with the name as the key. The name of the sequence is taken from the header,
# which is the first word after the '>' symbol until the first space. 
# The first word may contain any symbol (except spaces of course). 
# Furthermore, the method deals with multiple lines, and returns a single sequence.
# Input: 
#		file		The name of the fasta file
# Output: 
#		An array with 
#	(1)	A hash reference where the key is the item id and the value is the
#		sequence as a string.
#	(2)	An array reference including the ids in the order they are given in the
#		input file, $file. This information is necessary if you need the exact 
#		order, which is not given in the hash.
##################################################################################
sub read_fasta_file{
	my($file) = @_;
	my $FUNCTION = "read_fasta_file in Sequences.pm";
	
	my $id 			= "";
	my $seqstring 	= "";
	my %fasta 		= ();
	my %header		= ();
	my @order 		= ();
	my $line 		= "";
	open(IN_HANDLE, "<$file") || die "ERROR in $FUNCTION:\nCouldn't open the following file in package Tool,".
									 " sub read_fasta_file: $file/n";
	
	while($line = <IN_HANDLE>){
		chomp($line);
		
		# header (can contain one space after > symbol)
		if($line =~ /^\>\s?(\S+)\s*([\S*\s*]*)/){
			if($id){
				$fasta{$id} = $seqstring;
				$seqstring = "";
			}
			$id = $1;
			$header{$id} = $2;
			push(@order, $id);
		} else {
				$seqstring .= $line if ($id);
		}
	}
	
	if($id){
				$fasta{$id} = $seqstring;
				$seqstring = "";
	}
	my @return = (\%fasta, \@order, \%header);
	return \@return;
}

##################################################################################
# This method parses a standard fasta file and returns the sequences in a hash
# with the name as the key. The name of the sequence is taken from the header,
# which is the first word after the '>' symbol until the first dot or space. 
# The first word may contain any symbol (except dots and spaces of course). 
# Furthermore, the method deals with multiple lines, and returns a single sequence.
# Input: 
#		An array with 
#	(1)	A hash reference where the key is the item id and the value is the
#		sequence as a string.
#	(2)	An array reference including the ids in the order they are given in the
#		input file, $file. This information is necessary if you need the exact 
#		order, which is not given in the hash.
##################################################################################
sub read_fasta_without_version_number{
	my($file) = @_;
	
	print STDERR "reading $file...\n";
	
	my $id 			= "";
	my $seqstring 	= "";
	my %fasta 		= ();
	my @order		= ();
	my $line 		= "";
	open(IN_HANDLE, "<$file") || die "couldn't open the following file in package Tool,".
									 " sub read_fasta_file: $file/n";
	
	while($line = <IN_HANDLE>){
		
		# header (can contain one space after > symbol)
		if($line =~ /^\>\s?(\S+)\.?\d?\s*/){
			if($id){
				$fasta{$id} = $seqstring;
				$seqstring = "";
			}
			$id = $1;
			push (@order, $id);
		} else {
				$seqstring .= $line if ($id);
		}
	}
	
	if($id){
				$fasta{$id} = $seqstring;
				$seqstring = "";
	}
	my @return = (\%fasta, \@order);
	return \@return;
}


##################################################################################
# This method writes the data in a hash that has sequences as a value to a fasta file.
# NOTE: This method can be used in combination with the read_fasta_file* methods to
# convert a multi-lined fasta file to a single lined fasta file. First read the 
# fasta file to a hash and then call this method to write the hash to either the
# same file or a different one.
#
# Input: 
#		filename	The name of the fasta file with path
#		data_href	The reference to the hash containing the fasta data
# Output: none
##################################################################################
sub write_hash_to_fasta{
	my($filename, $data_href) = @_;
	
	open(OUT_HANDLE, ">$filename") || die "couldn't open the file for writing in package Tools,".
								"sub write_file_from_array: $filename\n";
		
	foreach my $key (sort(keys (%{$data_href}))){
		
		print OUT_HANDLE ">$key\n";
		print OUT_HANDLE $data_href->{$key}."\n\n";
	}
	close(OUT_HANDLE);
	return 1;
}

##################################################################################
# Calculates the sequence contents of a given sequence. The contents calculated are:
# A, G, C, U, A+G, A+C, A+U, G+C, G+U, and C+U. Here U=T.
# NOTE: Only works for RNA sequences, and DNA is translated to RNA. No other
# alphabets are supported.
# INPUT:
# sequencestring The sequence on which to perform calculations
# OUTPUT:
# The sequence contents in a hash, where the name of the content is given as the key
# and the percentage of each given as the value.
##################################################################################
sub calculate_sequence_contents{
	my ($sequencestring) = @_;
	
	my ($a, $c, $g, $u, $length) = 0;
	$length = length($sequencestring);
	my @seq_array = split('', $sequencestring);
	
	for ($i = 0 ; $i < $length ; $i ++){
		my $base = uc($seq_array[$i]);
		
		if($base eq 'A'){
			++$a;
		} elsif($base eq 'C'){
			++$c;
		} elsif($base eq 'G'){
			++$g;
		} elsif($base eq 'T'||$base eq 'U'){
			++$u;
		} else {
			print STDERR "Warning this is neither DNA or RNA: $base\n";
		}
	}
	
	my %contents = ();

	$contents{A} = $a/$length;
	$contents{C} = $c/$length;
	$contents{G} = $g/$length;
	$contents{U} = $u/$length;
	$contents{GC} = ($g + $c)/$length;
	$contents{AU} = ($a + $u)/$length;
	$contents{AG} = ($a + $g)/$length;
	$contents{AC} = ($a + $c)/$length;
	$contents{CU} = ($c + $u)/$length;
	$contents{GU} = ($g + $u)/$length;
	
	return \%contents;
}

##################################################################################
# Calculates all possible n-mer words for the given alphabet.
# INPUT: 
#	$n			The size of the words
#	$alphabet	The alphabet of the words
# RETURN: 
#	The array of all words
##################################################################################
sub words {
  my $n 		= shift;
  my $alphabet	= shift;
  
  my @n = split ('', $alphabet);
  if ($n == 1) {
    return \@n;
  } else {
    my $short = words($n - 1, $alphabet);
    my @long;
    for my $seq (@$short) {
      for my $n (@n) {
        push @long, $seq . $n;
      }
    }
    return \@long;
  }
}


##################################################################################
# Calculates word frequencies within given sequence seq and a fixed
# alphabet for words of size k.
#
# INPUT:
#		seq 	The sequence to calculate the frequencies on
#		k		The size of the words to calculate
#		aphabet	The alphabet the words and sequence are made up of
# OUTPUT:
#		The hash of each word as the key and the frequency as the value 
#		(includes 0 frequences)
##################################################################################
sub word_frequencies {
  	my $seq 		= shift;
  	my $k 		= shift;
  	my $alphabet	= shift;

  	my $words = words($k,$alphabet);

  	my %freq;
  	for my $word (@$words) {
    	$freq{$word} = 0;
  	}

  	for my $i (0 .. length($seq)-$k) {
    	my $word = substr($seq,$i,$k);
    	$freq{$word}++;
  	}
  	
  	return \%freq;

#	my @seq_freqs;
#	map {push(@seq_freqs,$freq{$_})} sort keys %freq;
#
#  	return \@seq_freqs;
}

##################################################################################
# This method divides the input sequence into x equal-sized subsequences.
# If the sequence does not divide evenly into x subsequences (as is the norm),
# the remainder is divided into two, so that the outer flanks of the sequences
# are discarded equally and the subsequences that are returned are of equal 
# length. The returned hash only contains the equal-sized subsequences in the
# order of the keys and does not contain the discarded flanks. This is good
# because the miRNA should not target the very ends anyway.
#
# INPUT: 
#		seqstring	The input sequence string
#		x			The number of equal-sized subsequences to cut the input seq into
# OUTPUT:
#		The hash with the equal-sized subsequences as the value and the order of
# 		the sequence as the key.
##################################################################################
sub divide_into_x_equal_sequences{
	my ($seqstring, $x) = @_;
	
	my $len = length($seqstring);
	print STDERR "length: $len\n";
	my $portionsize = int($len/$x);
	my $remain = $len - ($portionsize * $x);
	my $flanks = int($remain/2);
	my $lflank = $flanks;
	my $rflank = $remain - $lflank;
	
	my %portion_startpos_0to_nminus1 = (); 
	$portion_startpos_0to_nminus1{PORTIONSIZE} = $portionsize;
	
	for(my $i = 0 ; $i <$x ; $i++){
		my $offset = $lflank + $i*$portionsize;
		$portion_startpos_0to_nminus1{$offset} = substr($seqstring, $offset, $portionsize);
	}
	return \%portion_startpos_0to_nminus1;
}

##################################################################################
# Extract a sequence from a larger sequence with the given amount of context
# to the left and right of a motif sequence. Checks whether the flanks exist. 
# Returns the extracted sequence with the starting position of the motif in the 
# extract. This method needs the exact position(s) of the motif and does no checks.
# If you have the motif sequence and not the positions, use the method find_motif().
# NOTE: Indices are counted from 0 to n-1 for this method.
# 
# INPUT:
# 		seq 		The original sequence
#		motif_aref	An arrayref of the motif's starting position(s), 0-(n-1).
#					(Use method find_motif(seq, motif_seq))
#		motif_len	The length of the motif
# 		lflank		The length of the left flank to be extracted
# 		rflank		The length of the right flank to be extracted
#		seq_len		(Optional) The length of seq. Can be passed on if already 
#					calculated.
# OUTPUT:
#		An reference of an array of tuples. Each tuple contains
# 		the information about one extracted sequence in the following order
#		(position of motive in extract, extracted sequence)
##################################################################################
sub extract_motif_with_flanks{
	my ($seq,$motif_href, $motif_len, $lflank, $rflank, $seq_len) = @_;
	my $FUNCTION = "extract_motif_with_flanks in Sequences.pm";
	
	# length of sequence
	$seq_len = length($seq) unless ($seq_len);
	
	my @extracts = ();
	my @extract = ();
	my $ext_pos = $lflank + 1;
	
	foreach my $mpos (@{$motif_href}){
		my $lborder = $mpos - $lflank;
		my $extract_len = $lflank + $motif_len + $rflank;
		my $rborder = $lborder + $extract_len - 1;
		
		# check to see whether the extraction is possible
		if($lborder < 0){
			$lborder = 0;
			$ext_pos = $mpos;
		}
		if($rborder >= $seq_len){
			$rborder = $seq_len -1;
		}
		# create extract tuple
		push(@extract, $ext_pos);
		push(@extract, (substr($seq, $lborder, ($rborder - $lborder + 1))));
		# add to all extracts
		push(@extracts, \@extract);
		
	}
	return \@extracts;
}

##################################################################################
# Finds all occurences of a motif within a longer sequence and returns all 
# starting positons counting from 0 to n-1.
#
# INPUT:
#		seq			The original sequence
#		motif_seq	The motif sequence
#		seq_len		(Optional) The length of the sequence in seq, so it 
#					doesn't need to be computed.
# OUTPUT:
#		An array ref of all starting positions of the motif in the sequence
#		The array is empty if the motif is not found in the sequence.
#		Counting: 1 to n-1.
#
# STATUS: checked!
##################################################################################
sub find_motif{
	my ($seq, $motif_seq, $seq_len) = @_;
	my $FUNCTION = "find_motif in Sequences.pm";
	
	# get sequence length if necessary
	$seq_len = length($seq) unless($seq_len);
	# change both sequences to upper case
	$seq = uc($seq);
	$motif_seq = uc($motif_seq);
	
	my $curr_index = 0;
	my $index = 0;
	my @starts = ();
	
	# continue to search for motif until the end of the sequence is reached
	while($curr_index < $seq_len){
		
		# find motif from position in $curr_index
		$index = index($seq, $motif_seq, $curr_index);
		
		# the motif is not found
		if($index == -1){
			$curr_index = $seq_len;
		# the motif is found and continue searching to the right of the occurence
		} else {
			
			push(@starts, $index);
			$curr_index = $index +1; # set current index to the right of occurrence
		}
	}
	return \@starts;
}