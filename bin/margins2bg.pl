#!/usr/local/perl/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;

=head1 NAME

margins2bg.pl

=head1 SYNOPSIS

margins2bg.pl -bed COORDINATES.bed12 < margins

converts nucleotide-wise margins to bedGraph format
reads margins (or summarized margins) from stdtin or file, prints bedGraph
coordinates are taken from a bed file that should correspond to the
entries in the margins file

margins (summarized):
sequence id, sequence position, margin, (min, max, mean, median, sum)

Options:

    -bed        coordinates of margin entries
    -aggregate  which aggregate measure to use. choose one of
                (min, max, mean, median, sum)
                (default: mean)
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $bed;
my $aggregate;
my $help;
my $man;
my $debug;
my $result = GetOptions( "help" => \$help,
  "man"   => \$man,
  "debug" => \$debug,
  "bed=s" => \$bed,
  "aggregate=s" => \$aggregate );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);

( defined $bed ) or pod2usage("error: parameter bed mandatory");
( defined $aggregate ) or $aggregate = "mean";

# parse choice of aggregate measure
my %aggregates = ("min" => 0, "max" => 1, "mean" => 2, "median" => 3, "sum" => 4);
my $aggregate_choice = $aggregates{$aggregate};
if (not defined $aggregate_choice) {
	die("error: don't know about aggregate '$aggregate'. choose one of: min, max, mean, median, sum");
}

###############################################################################
# functions
###############################################################################

# convert margins to bedGraph format
sub cmp_bedgraph {
  $debug and say STDERR 'enter sub cmp_bedgraph';

  my ( $margins_aref, $bedline ) = @_;

  # parse bed entry
  my @bed = split( "\t", $bedline );
  if ( @bed != 12 and @bed != 6 ) {
    say STDERR "warning: skipping non bed6/12 entry: '$bedline'";
    return;
  }
  my ( $chrom, $chromStart, $chromEnd, undef, undef, $strand, undef, undef,
    undef, $blockCount, $blockSizes, $blockStarts ) = @bed;

  my @blockSizes;
  my @blockStarts;
  # parse bed12 entry
  if (@bed == 12) {
    @blockSizes  = split ",", $blockSizes;
    @blockStarts = split ",", $blockStarts;
    # sanity check
    if ( $blockCount != @blockSizes or $blockCount != @blockStarts ) {
      say STDERR "warning: skipping entry because of invalid number of blockSizes or blockStarts: '$bedline'";
      return;
    }
  } else {
    @blockSizes = ($chromEnd-$chromStart);
    @blockStarts = (0);
    $blockCount = 1;
  }

  # compute sequence of genome coordinates
  my @genome_coords;
  foreach my $blockStart (@blockStarts) {
    my $blockSize = shift @blockSizes;

    # compute coords relative to 0
    my @rel_coords = 0 .. $blockSize - 1;

    # compute coores relative to blockStart
    my @abs_coords = map { $_ + $blockStart + $chromStart } @rel_coords;

    # save sequence
    push @genome_coords, @abs_coords;
}

  # extract array of margins, use as chosen by user
  my @margins;
  foreach my $marginline (@$margins_aref) {

    # parse margins entry
    my ( $sequence_id, $sequence_position, $margin,
      @aggregates ) = split "\t", $marginline;
    
    # save value of chosen aggregate measure
    $margins[ $sequence_position - 1 ] = $aggregates[$aggregate_choice];
  }

  # reverse order of margins if on negative strand
  if ( $strand eq '-' ) {
    @margins = reverse @margins;
  }

  # print bedGraph lines
  if ( length(@margins) != length(@genome_coords) ) {
    say STDERR "warning: skipping entry because somehow we ended up with a " .
      "different number of coordinates and margins: " .
      scalar @genome_coords . " != " . scalar @margins;
    return;
  }
  $debug and say STDERR "iterating over " . scalar @genome_coords . " coordinates:";
  $debug and say STDERR join( "\t", @genome_coords );
  $debug and say STDERR "first coordinate: " . $genome_coords[0];
  $debug and say STDERR "last coordinate: " . $genome_coords[-1];
  foreach my $coord (@genome_coords) {
    $debug and say "printing coordinate $coord";
    my $margin = shift @margins;
    my $bedGraph_output = join( "\t", $chrom, $coord, $coord + 1, $margin );
    say $bedGraph_output;
  }
}

###############################################################################
# main
###############################################################################

# read all bed entries
my @bed;
open BED, '<', $bed or die "error opening file '$bed'\n$!";
while ( my $bedline = <BED> ) {
  chomp $bedline;
  push @bed, $bedline;
}
my $nbed = length(@bed);    # number of bed entries
close BED;

# parse margins; data for each sequence is forwarded to the conversion function
my $current_seqid = 1;      # fist sequence id is expected to be 1
my @linestack;
while ( my $line = <> ) {
  chomp $line;
  my ($seqid) = split( "\t", $line );
  if ( $seqid == $current_seqid ) {

    # if sequence id not changed, just record line
    push @linestack, $line;
  } else {

    # push margins and coordinates to conversion function
    my $bedline = shift @bed;
    cmp_bedgraph( \@linestack, $bedline );

    # start new linestack using current line
    @linestack = ($line);

    # update current_seqid
    $current_seqid = $seqid;
  }
}

# push margins and coordinates to conversion function
my $bedline = shift @bed;
cmp_bedgraph( \@linestack, $bedline );

# final sanity check
if ( $nbed != $current_seqid ) {
  die "error: read $nbed bed entries, but got " . $current_seqid . " sequences";
}
