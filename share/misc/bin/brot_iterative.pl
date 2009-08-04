#!@PERL@
# -*- perl -*-
# @configure_input@
# Last modified: 2009-07-20.15


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
# PRIVATE packages - END


# CONSTANTS        - BEGIN
# CONSTANTS        - END


# GLOBALS          - BEGIN
our $CmdCorb        = 'corb';
our $CmdRnaFold     = 'RNAfold';
our $CmdRnaDistance = 'RNAdistance';
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
    # we do not set defaults for brot here, since it comes with its own

    # parse @ARGV
    GetOptions(
        # brot
        'scoring=s'                   => \$argument_hashref->{scoring},
        'fixed-nuc=s'                 => \@{$argument_hashref->{fix}},
        'steps=i'                     => \$argument_hashref->{steps},
        'temp=f'                      => \$argument_hashref->{T},
        'seed=i'                      => \$argument_hashref->{seed},
        'negative-design-scaling=f'   => \$argument_hashref->{neg},
        'heterogenity-term-scaling=f' => \$argument_hashref->{het},
        'entropy-output'              => \$argument_hashref->{efile},
        'window-size=i'               => \$argument_hashref->{win},
        'sm-entropy=f'                => \$argument_hashref->{e_drop},
        'lambda=f'                    => \$argument_hashref->{lambda},
        'beta-long=f'                 => \$argument_hashref->{b_long},
        'beta-short=f'                => \$argument_hashref->{b_short},
        'speedup-threshold=f'         => \$argument_hashref->{s_thresh},
        'min-cool=f'                  => \$argument_hashref->{min_cool},
        'scale-cool=f'                => \$argument_hashref->{scale_cool},
        # script control
        'corb-command=s'              => \$CmdCorb,
        'RNAfold-command=s'           => \$CmdRnaFold,
        'RNAdistance-command=s'       => \$CmdRnaDistance,
        # info
        'verbose!'           => \&enable_verbose,
        'help'               => sub { pod2usage(-exitval => 0, -verbose => 1) },
        'man'                => sub { pod2usage(-exitval => 0, -verbose => 2) }
              );

    # get command line argument (an argument given without option)
    # fetch structure
    if ( ! defined($ARGV[0]))
    {
        msg_error_and_die("Structure to be designed missing.\nTry \"-help\" "
                         ."or \"-man\" for more information.\n");
    }
    $argument_hashref->{structure} = $ARGV[0];
}


# from options, assemble a string which can be directly passed to brot
#   assemble_brot_param_string(argument_hashref)
sub assemble_brot_param_string (\%)
{
    my ($arg_hashref) = @_;
    my $param_string = "";
    
    if (defined($arg_hashref->{scoring}))
    {
        $param_string .= " -c ".$arg_hashref->{scoring};
    }

    if (defined($arg_hashref->{steps}))
    {
        $param_string .= " -s ".$arg_hashref->{steps};
    }

    if (defined($arg_hashref->{T}))
    {
        $param_string .= " -t ".$arg_hashref->{T};
    }

    if (defined($arg_hashref->{seed}))
    {
        $param_string .= " -r ".$arg_hashref->{seed};
    }

    if (defined($arg_hashref->{neg}))
    {
        $param_string .= " -d ".$arg_hashref->{neg};
    }

    if (defined($arg_hashref->{het}))
    {
        $param_string .= " -h ".$arg_hashref->{het};
    }

    if (defined($arg_hashref->{efile}))
    {
        $param_string .= " -p ".$arg_hashref->{efile};
    }

    if (defined($arg_hashref->{win}))
    {
        $param_string .= " -w ".$arg_hashref->{win};
    }

    if (defined($arg_hashref->{e_drop}))
    {
        $param_string .= " -e ".$arg_hashref->{e_drop};
    }

    if (defined($arg_hashref->{lambda}))
    {
        $param_string .= " -l ".$arg_hashref->{lambda};
    }

    if (defined($arg_hashref->{b_long}))
    {
        $param_string .= " -o ".$arg_hashref->{b_long};
    }

    if (defined($arg_hashref->{b_short}))
    {
        $param_string .= " -i ".$arg_hashref->{b_short};
    }

    if (defined($arg_hashref->{s_thresh}))
    {
        $param_string .= " -u ".$arg_hashref->{s_thresh};
    }

    if (defined($arg_hashref->{min_cool}))
    {
        $param_string .= " -j ".$arg_hashref->{min_cool};
    }

    if (defined($arg_hashref->{scale_cool}))
    {
        $param_string .= " -q ".$arg_hashref->{scale_cool};
    }


    return $param_string;
}


# predict a sequence using brot
#   predict_sequence (parameters, fixed_sites, structure, structure_length)
sub predict_sequence (\$ \$ \$ $)
{
    my ($params_ref, $fixed_sites_ref, $struct_ref, $struct_len) = @_;
    my $seq;
    my @err_run;

    # call brot
    unless(open(FH,
       "$CmdCorb \"brot $$params_ref $$fixed_sites_ref $$struct_ref\" 2>&1 |"))
    {
        msg_error_and_die ("Could not start $CmdCorb brot: $!\n");
    }

    foreach (<FH>)
    {
        if ($_ =~ /(^[ACGUacgu]{$struct_len})/)
        {
            close(FH);
            return $1;
        }

        push(@err_run, $_);
    }
    close(FH);

    msg_error_and_die("Running brot failed, output of\n"
         ."\`$CmdCorb \"brot $$params_ref $$fixed_sites_ref $$struct_ref\"\`: ",
                      @err_run);
}


# predict a structure using RNAfold
#   predict_structure (sequence, sequence_len)
sub predict_structure (\$ $)
{
    my ($seq_ref, $seq_len) = @_;
    my $structure;
    my $mfe_en;
    my $ens_en;
    my @err_run;

    #call RNAfold
    unless(open(FH, "echo \"$$seq_ref\" | $CmdRnaFold -d2 -p0 2>&1 |"))
    {
        msg_error_and_die ("Could not start $CmdRnaFold: $!\n");
    }

    foreach (<FH>)
    {
        #print $_;
        if ($_ =~ /([\(\.\)]{$seq_len})\s+\(\s*(\-?\d+\.\d+)\)/)
        {
            $structure = $1;
            $mfe_en = $2;
        }
        elsif ($_ =~ /free\s+energy\s+of\s+ensemble\s+\=\s+(\-\d+\.\d+)/)
        {
            $ens_en = $1;
        }
        push(@err_run, $_);
    }
    close(FH);

    if (!defined($structure))
    {
        msg_error_and_die ("No structure found found for sequence "
                           ."\"$$seq_ref\" with command\n"
                           ."\`echo \"$$seq_ref\" | $CmdRnaFold -d2 -p0 \`: ",
                           @err_run);
    }
    elsif (!defined($ens_en))
    {
        msg_error_and_die ("No ensemble energy found found for sequence "
                          ."\"$$seq_ref\" with command\n"
                           ."\`echo \"$$seq_ref\" | $CmdRnaFold -d2 -p0 \`: ",
                           @err_run)
    }

    return ($structure, $mfe_en, $ens_en);
}


sub compare_structures (\$ \$)
{
    my ($ref_ref, $pred_ref) = @_;
    my (@ref_ary, @pred_ary);
    my @matches;
    my $i;
    # first idea: just compare structures, return array of matching positions

    # create arrays from structures
    @ref_ary = split(//, $$ref_ref);
    @pred_ary = split(//, $$pred_ref);

    # loop over structures and compare, store matching pos for fixing
    for ($i = 0; $i <= $#ref_ary; $i++)
    {
        #print $ref_ary[$i],":",$pred_ary[$i],"\n";
        if ($ref_ary[$i] eq $pred_ary[$i])
        {
            push(@matches, $i);
        }
    }

    #print join(":", @matches),"\n";

    return @matches;
}

sub create_fix_sites_string(\@ \@ \$)
{
    my ($user_fix_ary_ref, $match_fix_ary_ref, $seq_ref) = @_;
    my $fixed_sites = "";
    my $cur = 0;
    my $u_fix_pos;
    my $add_next;

    # check for sites in both sets: user and match
    foreach (@{$match_fix_ary_ref})
    {
        #print $_,substr($$seq_ref, $_, 1),"\n";

        $add_next = substr($$seq_ref, $_, 1).':'.$_;

        for (; $cur <= $#{$user_fix_ary_ref}; $cur++)
        {
            $u_fix_pos = substr(${$user_fix_ary_ref}[$cur], 2);

            if ($u_fix_pos > $_)
            {
                last;
            }
            elsif ($u_fix_pos == $_)
            {
                $add_next = ${$user_fix_ary_ref}[$cur]
            }
            else
            {
                $fixed_sites .= " -n ".${$user_fix_ary_ref}[$cur];
            }
        }
        
        $fixed_sites .= " -n ".$add_next;

        #print $fixed_sites."\n";
    }

    for (; $cur <= $#{$user_fix_ary_ref}; $cur++)
    {
        $fixed_sites .= " -n ".${$user_fix_ary_ref}[$cur];
        #print $fixed_sites."\n";
    }

    return $fixed_sites;
}
# FUNCTIONS        - END


# MAIN             - BEGIN
my %arg_hash;
my $brot_params = '';
my $fixed_sites = '';
my $seq;
my $struct_len;
my $pred_structure;
my $mfe;
my $ens_en;
my @matching_sites;
my $turns;

# parse commandline
parse_args(%arg_hash);

# assemble brot parameter string (wo fixed sites)
$brot_params = assemble_brot_param_string(%arg_hash);

# init fixed nucs
if ($#{$arg_hash{fix}} >=0)
{
    @{$arg_hash{fix}} = map  { $_->[0] }            
                        sort { $a->[1] cmp $b->[1] }
                        map  { [$_, substr($_, 2)] }
                        @{$arg_hash{fix}};

    $fixed_sites = ' -n '.join(' -n ', @{$arg_hash{fix}}).' ';
}

$struct_len = length($arg_hash{structure});

$turns = 3 * $struct_len;

while ($turns > 0)
{
    $seq = predict_sequence($brot_params,
                            $fixed_sites,
                            $arg_hash{structure},
                            $struct_len);
    
    ($pred_structure, $mfe, $ens_en)  = predict_structure ($seq, $struct_len);
    
    @matching_sites = compare_structures($arg_hash{structure}, $pred_structure);
    if ($struct_len - scalar(@matching_sites) > 0)
    {
        $fixed_sites = create_fix_sites_string(@{$arg_hash{fix}},
                                               @matching_sites,
                                               $seq);
    }
    else
    {
        last;
    }

    print "$turns: $seq ", $pred_structure, " ",$mfe, " ",$ens_en,"\n";
    $turns--;
}

if ($turns != 0)
{
    print "$seq ($mfe)\n";
}
else
{
    print "Could not find a sequence for \"$arg_hash{structure}\"!\n";
}

# check if result

# when do we verify structure? -> brot verifies, catch error
# MAIN             - END


__END__


=head1 NAME

brot_iterative - Run B<brot>, fold structure, compare in- and output, and
re-run.

=head1 SYNOPSIS

B<brot_iterative> [options] STRUCTURE

=head1 DESCRIPTION

B<brot_iterative> does the whole sequence design cycle for you: Design a
sequence from a structure (the input), predict the sturcture of the designed
sequence (the output), compare input and output and re-run B<brot> to fix
unwanted/ missing interactions. 

Re-running B<brot means>, that the "good sites" are set to a fixed state, while
unintended interactions are left open.

- strategy for wrong positions?
- basically same options as B<brot>

For structure prediction and comparison, the Vienna RNA package is used
(RNAfold, RNAdistance).

Please note: All options also work if invoked with an unique prefix. The
default values of the options are the same as for the B<brot> binary.

=head1 OPTIONS

=over 2

=head2 Options handed over to B<brot>

=over 6

=item B<-scoring <NAME>>

Use a certain energy model to score the RNA sequence. Possible NAMEs are `NN',
`nussinov' and `simpleNN'. `nussinov' uses a Nussinov model, only counting
Hbonds. `NN' invokes the Nearest Neighbour model. `simpleNN' uses a very coarse
grained approximation of the Nearest Neighbour model, scoring every base pair
as a stack, averaging mismatches in stacked pairs and without any loop
parameters.

=item B<-steps <INT>>

Number of iterations for the update of site probabilities.

=item B<-temp <FLOAT>>

Initial temperature of the system.

=item B<-seed <INT>>

Random seed used for adding thermal noise to the Nearest Neighbour model. This
is used for producing different answers for the same structure. Adding thermal
noise means adding small random numbers to all parameters. The range of those
numbers is [0.005, -0.005] and should be beyond the level of significance.
Since we use a correlated random number generator you can reproduce answers by
using the same seed. Only takes effect when using with `NN' as scoring scheme.

=item B<-negative-design-scaling <FLOAT>>

Scaling factor for the negative design term.

=item B<-heterogenity-term-scaling <FLOAT>>

Scaling factor for the sequence heterogenity term.

=item B<-entropy-output <FILENAME>>

Write the changes of the sequence matrix entropy, short and long term avg.'s
and temperature to given file. Only values from the simulation of interest are
written.

=item B<-window-size <INT>>

Size of the window to the left and right of a base to be considered for
calculating the heterogeneity term when using the `NN' scoring scheme. Please
note that this is always only one half of the window. Only takes effect when
using with `NN' as scoring scheme.

=item B<-sm-entropy <FLOAT>>

If the entropy of the sequence matrix drops below this value, the simulation
will stop.

=item B<-lambda <FLOAT>>

Describes the portion of an old probability to be mixed with the new one. The
old probability gets l a share of l, the new (1 - l). Used to avoid oscillation
in the system.

=item B<-beta-long <FLOAT>>

Share of the current long term avg. entropy to be used for the next step. New
value calculates from o * S_long + (1-o) * S_current.

=item B<-beta-short <FLOAT>>

Share of the current short term avg. entropy to be used for the next step. New
value calculates from o * S_short + (1-o) * S_current.

=item B<-speedup-threshold <FLOAT>>

If the ratio of current short- and long term entropy drops below this value, we
slow down cooling, above we speed up.

=item B<-min-cool <FLOAT>>

If the cooling factor drops below this value we do no further speedups.

=item B<-scale-cool <FLOAT>>

Speeding up cooling is done via (c * (c * q)).

=item B<-fixed-nuc <NUCLEOTIDE:INT>>

Use a fixed nucleotide in a position in the sequence during the simulation. As
NUCLEOTIDE the whole 15 letter RNA alphabet is allowed. The position INT has to
be in the range of the structure. May be called repetitive.

=back

=back

=over 2

=head2 Options controling the script itself

=over 6

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

=back

=back

=over 2

=head2 Options controling help texts

=over 6

=item B<-man>

Print a man page for the program. This is a more detailed description of the
script.

=item B<-help>

Print a help message and exit. This only contains the synopsis and the options
description.

=cut

# compare function using rnadistance
#    my ($ref_ref, $pred_ref) = @_;
#    my @err_run;
#
#    unless(open(FH,
#                "echo \"$$ref_ref\n$$pred_ref\" | $CmdRnaDistance -B -DP 2>&1 |"))
#    {
#        msg_error_and_die ("Could not start $CmdRnaDistance: $!\n");
#    }
#
#    foreach (<FH>)
#    {
#        print $_;
#        if ($_ =~ /P\:\s*(\d+)/)
#        {
#            close(FH);
#            return $1;
#        }
#        push(@err_run, $_);
#    }
#
#    close(FH);
#
#    msg_error_and_die("Running RNAdistance failed, output of\n"
#                 ."\`echo \"$$ref_ref\\n$$pred_ref\" | $CmdRnaDistance -DP\`: ",
#                      @err_run);

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
