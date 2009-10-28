#!@PERL@
# -*- perl -*-
# @configure_input@
# Last modified: 2009-10-16.10


# Copyright (C) 2009 Stefan Bienert
#
# This file is part of CoRB.
#
# CoRB is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoRB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoRB.  If not, see <http://www.gnu.org/licenses/>.


# SETTINGS         - BEGIN
BEGIN
{
  # try to get bash compatible shell
  $ENV{'SHELL'} = '@SHELL@' if exists $ENV{'DJGPP'};
  $|++;
}
# SETTINGS         - END


# PERL packages    - BEGIN
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
# PERL packages    - END


# PRIVATE packages - BEGIN
use lib '@CORB_PERL5LIB@';
use CorbIO qw(:all);
use PBar;
# PRIVATE packages - END


# CONSTANTS        - BEGIN
# CONSTANTS        - END


# GLOBALS          - BEGIN
our $CmdCorb        = 'corb';
our $CmdRnaFold     = 'RNAfold';
our $CmdRnaDistance = 'RNAdistance';
our $Cpus           = 2;
# GLOBALS          - END


# FUNCTIONS        - BEGIN
# parse the arguments of the script and store them in a hash
#   parseargs(argument_hashref)
sub parse_args(\%)
{
    my ($argument_hashref) = @_;

    # wrap the signal handler for warnings to copy them to our message space
    local $SIG{__WARN__} = sub { msg_error_and_die("@_") };

    # set defaults
    $argument_hashref->{brot_args} = '';
    $argument_hashref->{outfile}   = '-';

    # parse @ARGV
    GetOptions(
        # brot
        'brot-args=s'           => \$argument_hashref->{brot_args},
        # results
        'outfile=s'             => \$argument_hashref->{outfile},
        # script control
        'corb-command=s'        => \$CmdCorb,
        'RNAfold-command=s'     => \$CmdRnaFold,
        'RNAdistance-command=s' => \$CmdRnaDistance,
        'cpus=i'                => \$Cpus,
        # info
        'verbose!'           => sub { enable_verbose; PBar::enable(); },
        'help'               => sub { pod2usage(-exitval => 0, -verbose => 1) },
        'man'                => sub { pod2usage(-exitval => 0, -verbose => 2) }
              );

    if ((is_verbose()) && ($argument_hashref->{outfile} eq '-'))
    {
        msg_error_and_die('Option verbose requires option outfile to be ',
                          "set.\n");
    }

    # get command line argument (an argument given without option)
    # fetch file name
    if ($#ARGV != 0)
    {
        msg_error_and_die("Structure file missing.\nTry \"-help\" or \"-man\" ",
                          "for more information.\n");
    }
    if ( ! -r $ARGV[0])
    {
        msg_error_and_die("Structure file does not exist or is not readable: ",
                          "$ARGV[0]\n");
    }
    $argument_hashref->{structfile} = $ARGV[0];
}

# predict sequences for a bunch of structures given. The arrays
# @structure_lengths and @structures have to be of same size.
#   predict_sequence($arg_string, @structure_lengths, @structures)
sub predict_sequence(\$ \@ \@ @)
{
    my ($arg_ref, $info_ary_ref, $lens_ary_ref, @structures) = @_;
    my (@runs, $fh);
    my $i;
    my @res;
    my @err_run;

    # start all prediction
    foreach (@structures)
    {
        local *FILE;
        
        unless(open(FILE, "$CmdCorb \"brot $$arg_ref $_\" |"))
        {
            msg_error_and_die ("Could not start $CmdCorb brot: $!\n");
        }
        #print "Started $_\n";
        push(@runs, *FILE);
    }

    # collect results
    for ($i = 0; $i <= $#runs; $i++)
    {
        @err_run = ();

        $fh = $runs[$i];
        foreach (<$fh>)
        {
            push(@err_run, $_);
            if ($_ =~ /ERROR/)
            {
                next;
            }           
            if ($_ =~ /(^[ACGUacgu]{${$lens_ary_ref}[$i]})/)
            {
                close($fh);
                $res[$i] = $1;
                ${${$info_ary_ref}[$i]}{'seq'} = $res[$i];
                next;
            }
        }
        close($fh);

        if (! defined($res[$i]))
        {
            msg_error_and_die("Running brot failed, output of\n"
                           ."\`$CmdCorb \"brot $$arg_ref $structures[$i]\"\`: ",
                              @err_run);
        }

        #print "finished: $res[$i]\n";
    }

    return @res;
}

# Calculate the 'energy' of a bunch of RNA 2D structures. Of course, the
# sequence array must provide a sequence for each structure and vice versa.
#   calc_en_structure(@sequences, @structures)
sub calc_en_structure(\@ \@ @)
{
    my ($info_ary_ref, $seq_ary_ref, @structures) = @_;
    my $i;
    my (@runs, $fh);
    my @res;
    my @err_run;

    # start calculation for each seq/ struct pair
    foreach ($i = 0; $i <= $#structures; $i++)
    {
        local *FH;

        unless(open(FH,
               "$CmdCorb \"er2de ${$seq_ary_ref}[$i] $structures[$i]\" 2>&1 |"))
        {
            msg_error_and_die ("Could not start $CmdCorb er2de: $!\n");
        }
        push(@runs, *FH);
    }

    # collect energies
    for ($i = 0; $i <= $#runs; $i++ )
    {
        @err_run = ();
        $fh = $runs[$i];

        foreach (<$fh>)
        {
            if ($_ =~ /G\s*\=\s*(\-?\d+\.?\d*)/)
            {
                $res[$i] = $1;
                ${${$info_ary_ref}[$i]}{'p_en'} = $res[$i];
                next;
            }
            push(@err_run, $_);
        }
        close($fh);

        if (! defined($res[$i]))
        {
            msg_error_and_die("Running er2de failed, output of\n"
                ."\`$CmdCorb \"er2de ${$seq_ary_ref}[$i] $structures[$i]\"\`: ",
                              @err_run);
        }
    }

    return @res;
}

# Predict the minimum free energy structure for a bunch of sequences. Returns
# an array of hashes, holding the structure and its energy.
#   predict_mfe_structure (@sequences)
sub predict_mfe_structure (\@ \@ \@)
{
    my ($info_ary_ref, $seq_ary_ref, $lengths_ary_ref) = @_;
    my $i;
    my (@runs, $fh);
    my @err_run;
    my %result;
    my @res;

    # start folding for all seqs
    foreach (@{$seq_ary_ref})
    {
        local *FH;

        unless(open(FH, "echo \"$_\" | $CmdRnaFold -d2 -p0 2>&1 |"))
        {
            msg_error_and_die ("Could not start $CmdRnaFold: $!\n");
        }

        push(@runs, *FH);
    }

    # collect structures and energies
    for ($i = 0; $i <= $#runs; $i++)
    {
        @err_run = ();
        %result = ();

        $fh = $runs[$i];

        foreach (<$fh>)
        {
            if ($_ =~ 
                /([\(\.\)]{${$lengths_ary_ref}[$i]})\s+\(\s*(\-?\d+\.\d+)\)/)
            {
                #$result_count++;
                $result{structure} = $1;
                $result{mfe_en} = $2;

                ${${$info_ary_ref}[$i]}{'mfe_structure'} = $result{structure};
                ${${$info_ary_ref}[$i]}{'mfe_en'} = $result{mfe_en};
            }
            elsif ($_ =~ /free\s+energy\s+of\s+ensemble\s+\=\s+(\-\d+\.\d+)/)
            {
                #$result_count++;
                
                $result{ens_en} = $1;
                ${${$info_ary_ref}[$i]}{'ens_en'} = $result{ens_en};
            }

            push(@err_run, $_);
        }
        
        close($fh);
        
        if (! defined ($result{structure}) || ! defined ($result{ens_en}))
        {
            msg_error_and_die("Running RNAfold failed, output of\n"
                  ."\`echo \"${$seq_ary_ref}[$i]\" | $CmdRnaFold -d2 -p0\"\`: ",
                              @err_run);
        }

        $res[$i] = { %result };
    }

    return @res;
}

# Compare two sets of structures. Only structures with the same index are
# compared.
#   compare_structures (@set1, @set2)
sub compare_structures (\@ \@ @)
{
    my ($info_ary_ref, $pred_ary_ref, @ref) = @_;
    my $i;
    my (@runs, $fh);
    my @err_run;
    my @res;

    # start RNAdistance for all index pairs
    for ($i = 0; $i <= $#ref; $i++)
    {
        local *FH;

        unless(open(FH,
"echo \"$ref[$i]\n${$pred_ary_ref}[$i]{structure}\" | $CmdRnaDistance -DP 2>&1 |"))
        {
            msg_error_and_die ("Could not start $CmdRnaDistance: $!\n");
        }

        push(@runs, *FH);
    }

    # collect distances
    for ($i = 0; $i <= $#runs; $i++)
    {
        @err_run = ();

        $fh = $runs[$i];

        foreach (<$fh>)
        {
            if ($_ =~ /P\:\s*(\d+)/)
            {
                $res[$i] = $1;
                ${${$info_ary_ref}[$i]}{'dist'} = $res[$i];
                next;
            }

            push(@err_run, $_);
        }
        
        close($fh);
        
        if (!defined($res[$i]))
        {
            msg_error_and_die("Running RNAdistance failed, output of\n"
."\`echo \"$ref[$i]\n${$pred_ary_ref}[$i]{structure}\" | $CmdRnaDistance -DP\`: ",
                              @err_run);
        }
    }

    return @res;
}
# FUNCTIONS        - END


# MAIN             - BEGIN
# open input file
my %arg_hash;
my $fh;
my $outfh;
my ($i, $j);
my @en_diff;
my $sum_en_diff = 0;
my $length;
my @diff;
my $sum_diff = 0;
my $correct = 0;
my @structures;
my $jobs;
my @struct_lens;
my @predicted_seqs;
my @pseq_ens;
my @mfe_predictions;
my %progress;
my @info;

# fetch options/ arguments
parse_args(%arg_hash);

$jobs = $Cpus;

# load structures and sort them according to length
$fh = open_or_die($arg_hash{structfile}, "<");

@structures = <$fh>;
close($fh);
for ($i = 0; $i <= $#structures; $i++)
{
    if ($structures[$i] !~ /^[\(\)\.]+$/)
    {
        splice(@structures, $i, 1);
        $i--;
    }
    else
    {
        chomp($structures[$i]);
    }
}

@structures = map  { $_->[0] }               # unmap
              sort { $b->[1] <=> $a->[1] }   # sort according to length
              map  { [$_, length($_)] }      # map each struct. with its length
              @structures;                   # take structures

# run through all structures
$outfh = open_or_die($arg_hash{outfile}, '>');
PBar::start(%progress, $#structures); # do we need $#structures + 1?
for ($i = 0; $i <= $#structures; $i += $jobs)
{
    if (($i + $jobs) > $#structures) { $jobs = ($#structures - $i) + 1 }

    @struct_lens = ();
    @struct_lens = map { length($_) } @structures[$i..($i + $jobs - 1)];

    # prepare output information
    for ($j = 0; $j < $jobs; $j++)
    {
        ${$info[$j]}{'structure'} = $structures[$i + $j];
    }

    # predict sequence
    @predicted_seqs = ();
    @predicted_seqs = predict_sequence($arg_hash{brot_args},
                                       @info,
                                       @struct_lens,
                                       @structures[$i..($i + $jobs - 1)]);

    # calc. energy of pred.seq. in des.struct.
    @pseq_ens = ();
    @pseq_ens = calc_en_structure(@info,
                                  @predicted_seqs,
                                  @structures[$i..($i + $jobs - 1)]);

    # predict structure from pred.seq
    @mfe_predictions = ();
    @mfe_predictions = predict_mfe_structure(@info,
                                             @predicted_seqs,
                                             @struct_lens);

    @en_diff = ();
    for ($j = 0; $j < $jobs; $j++)
    {
        $en_diff[$j] = $pseq_ens[$j] - $mfe_predictions[$j]{mfe_en};
    }

    # calc. distance between pred.seq.struct. and mfe.struct.
    @diff = ();
    @diff = compare_structures (@info, @mfe_predictions,
                                @structures[$i..($i + $jobs - 1)]);

    for ($j = 0; $j < $jobs; $j++)
    {
        $sum_en_diff += $en_diff[$j];
        $sum_diff += $diff[$j];

        if ($diff[$j] == 0)
        {
            $correct++;
        }

        #print (($i+$j+1).": $sum_diff $correct $sum_en_diff \n");

        print $outfh 'Structure ', ($i + $j), ": ${$info[$j]}{'structure'}\n";
        print $outfh 'Pred. seq.:  ', "${$info[$j]}{'seq'}\n";
        print $outfh 'G(Structure, Predicted seq): ', "${$info[$j]}{'p_en'}\n";
        print $outfh 'MFE struct.: ', "${$info[$j]}{'mfe_structure'}\n";
        print $outfh 'G(MFE Structure, Pred. seq): ',"${$info[$j]}{'mfe_en'}\n";
        print $outfh 'G(Ensemble, Pred. seq):      ',"${$info[$j]}{'ens_en'}\n";
        print $outfh 'Diff(G_MFE, G_sps): ', $en_diff[$j], "\n";
        print $outfh 'Dist(Structure, Pred. struct.): ',
                     "${$info[$j]}{'dist'}\n";
        print $outfh '================================================', "\n";
        print $outfh 'Globals update: Energy diff. = ',
                     "$sum_en_diff BPD = $sum_diff Correct = $correct\n";
        print $outfh '================================================', "\n\n";
    }

    PBar::update(%progress, $i);
}
PBar::finish(%progress);
close($outfh);

print "\nSums: Energy = $sum_en_diff BPD = $sum_diff Correct = $correct\n";
# MAIN             - END


__END__


=head1 NAME

eval_structure_file - Run B<brot> on a list of structures, predict and compare
structures and output the sums of correct predictions in total and energy difference.

=head1 SYNOPSIS

B<eval_structure_file> [options] FILE

=head1 DESCRIPTION

B<eval_structure_file> is intended to evaluate a set of B<brot> parameters
agains a whole list of structures. The structures to be tested have to be
provided via file in Vienna notation. The result is presented as total number
of good predictions (sequences folding up into the right shape) and the sum of
the base pair distances/ energy differences in case of wrong sequence
prediction.

For structure prediction and comparison, the Vienna RNA package is used
(RNAfold, RNAdistance).

Please note: All options also work if invoked with an unique prefix. The
default values of the options are the same as for the B<brot> binary.

=head1 OPTIONS

=over 8

=item B<-brot-args <STRING>>

Instead of reproducing all options of B<brot> here, with all the maintenance
effort, we just provide you the possibility to pass a string full of parameters
to B<brot>.

=item B<-outfile <FILENAME>>

Write the results of all the tests to <FILENAME> rather than to standard out.
<FILENAME> must not exist before the script runs.

This option is required if B<-verbose> is invoked.

=item B<-verbose>

Be verbose.

=item B<-corb-command>

This defines the command to invoke B<CoRB>. By default the script just tries to
call `corb`. If you have not installed B<CoRB> in your $PATH, just give the
path and name of the binary here. The name is required to provide the
possibility to invoke B<CoRB> binaries with differing names.

=item B<-RNAfold-command>

This defines the command to invoke B<RNAfold>. By default the script just tries
to call `RNAfold`. If you have not installed B<RNAfold> in your $PATH, just
give the path and name of the binary here. The name is required to provide the
possibility to invoke B<RNAfold> binaries with differing names.

=item B<-RNAdistance-command>

This defines the command to invoke B<RNAdistance>. By default the script just
tries to call `RNAdistance`. If you have not installed B<RNAdistance> in your
$PATH, just give the path and name of the binary here. The name is required to
provide the possibility to invoke B<RNAdistance> binaries with differing names.

=item B<-cpus <INT>>

Define how many instances of an external tool should be run in parallel.

=item B<-man>

Print a man page for the program. This is a more detailed description of the
script.

=item B<-help>

Print a help message and exit. This only contains the synopsis and the options
description.

=cut


# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
