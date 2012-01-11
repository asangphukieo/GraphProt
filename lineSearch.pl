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
use File::Copy;
use Cwd;

=head1 NAME

lineSearch.pl

=head1 SYNOPSIS

lineSearch.pl -gspan [gspan file]

Options:

    -gspan      optimize parameters for these graphs
    -affy       affinities for graphs
    -mf         makefile to use for crossvalidation
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# create temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
my $tmp_template = 'lineSearch-XXXXXX';
my $tmp_prefix = '/var/tmp/';
my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);
$SIG{'INT'} = 'end_handler';
$SIG{'TERM'} = 'end_handler';
$SIG{'ABRT'} = 'end_handler';
$SIG{'USR1'} = 'end_handler';
$SIG{'USR2'} = 'end_handler';
sub end_handler {
	print STDERR "signal ", $_[0], " caught, cleaning up temporary files\n";
	# change into home directory. deletion of the temporary directory will
	# fail if it is the current working directory
	chdir();
	File::Temp::cleanup();
	die();
}

###############################################################################
# parse command line options
###############################################################################
my $gspan;
my $affy;
my $mf;
my $help;
my $man;
my $debug;
my $result = GetOptions (	"gspan=s"	=> \$gspan,
							"affy=s"    => \$affy,
							"mf=s"		=> \$mf,
							"help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);
(defined $gspan) or pod2usage("error: -gspan parameter mandatory");
(-f $gspan) or die "error: could not find file '$gspan'";
(-f $mf) or die "error: could not find file '$mf'";
(-f $affy) or die "error: could not find file '$affy'";

###############################################################################
# main
###############################################################################

# global variables
my $CURRDIR = cwd();
($debug) and say STDERR "cwd: $CURRDIR";

# binaries
my $libsvm = '~/src/libsvm-3.0/svm-train';
my $libsvm_options = ' -v 5 -s 3 -t 0';

# we optimize these parameters
my @parameters = qw/ e c R D /;

# valid values for parameters
my %parameters;
$parameters{'e'}{default} = 0.1;
$parameters{'e'}{values} = [0.001, 0.01, 0.1, 1, 10, 100];
$parameters{'c'}{default} = 1;
$parameters{'c'}{values} = [0.001, 0.01, 0.1, 1, 10, 100];
$parameters{'R'}{default} = 1;
$parameters{'R'}{values} = [1, 2, 3, 4];
$parameters{'D'}{default} = 4;
$parameters{'D'}{values} = [1, 2, 3, 4, 5, 6];

for my $par (@parameters) {
	$parameters{$par}{current}=$parameters{$par}{default};
}

# print important variables
if ($debug) {
	say STDERR 'parameters to optimize: ', join(', ', @parameters);
	say STDERR 'keys in hash: ', join(', ', keys %{$parameters{'epsilon'}});
	while (my ($param, $param_h) = each %parameters) {
		while (my ($key, $values) = each %{$param_h}) {
			say STDERR join('/', $param, $key), ': ', $values;
		}
	}
}

# main loop: do until finished
my $optimization_finished = 0;
do {
	# optimize each parameter
	for my $par (@parameters) {
		for my $try_this ($parameters{$par}{values}) {
			say STDERR "optimizing $gspan";
			$tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);
			
			# check if parameter combination is valid
			# TODO
			
			# copy relevant files into tmp
			copy($gspan, $tmpdir);
			copy($affy, $tmpdir);
			copy($mf, $tmpdir);
			
			# test parameter combination / get value from previous run
			chdir($tmpdir);
			($debug) and say STDERR "changed directory: ", cwd();
			# create parameter file
			my $param_file = $gspan;
			$param_file =~ s/.gspan//;
			$param_file .= '.pars';
			$debug and say STDERR "param_file: $param_file";
			open PARS, '>', $param_file;
			for my $par (@parameters) {
				($debug) and say STDERR $par, ' ', $parameters{$par}{current};
				say PARS $par, ' ', $parameters{$par}{current};
			}
			say PARS 'b 14';
			close PARS;
			# TODO
			# call Makefile for cv
			# TODO
			# parse result
			# TODO
			# run test			
			system('ls -l');
			
			# exit temp directory
			chdir($CURRDIR);
			($debug) and say STDERR "changed directory: ", cwd();
			File::Temp::cleanup();
			
			# save result
			# TODO
		}
		# set current to the best parameter combination
	}
	$optimization_finished = 1;
} while (not $optimization_finished);

chdir();
