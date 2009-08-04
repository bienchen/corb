#!@PERL@
# -*- perl -*-
# @configure_input@
# Last modified: 2009-07-06.16


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
use POSIX;
# PERL packages    - END


# PRIVATE packages - BEGIN
use lib '@CORB_PERL5LIB@';
#use PBar;
use CorbIO qw(:all);
use RNA qw(:all);
use Simplex qw(:all);
# PRIVATE packages - END

# CONSTANTS        - BEGIN
my %SMALL_TESTS = ('a' => '((....))',
                   'aa' => '((((((...(((......)))...(((...)))...))))))',
                   'aaa' => '((((....(((...)))...(((((((((....)))))))))......))))',
                   'aaaa' => '(((...(((.......((((.....))))))))))',
                   #'aaaaa' => '(((((...(((...)))...((((...))))...(((...(...))))...(((((....)))))..(...)....(((...)))...)))))',
                   #'aaaaaa' => '...(((...)))(((())))...(((...(...))))',
                   #'aaaaaaa' => '(((((....)))))..(...).((()))',
                   #'aaaaaaaa' => '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
    );
# CONSTANTS        - END

# GLOBALS          - BEGIN
# GLOBALS          - END


# FUNCTIONS        - BEGIN
# parse the arguments of the script and store them in a hash
#   parseargs(argument_hashref, error_msgref)
sub parseargs(\%)
{
    my ($argument_hashref) = @_;
    my  $optcatchresult = 0;
    my  $help;
    my  $man;
    my  $verbose;

    # wrap the signal handler for warnings to copy them to our message space
    local $SIG{__WARN__} = sub
                           {
                               msg_error_and_die("@_");
                           };

    # set defaults
    #$argument_hashref->{max} = 10.0;
    #$argument_hashref->{min} = 0.0;
    #$argument_hashref->{step} = 0.01;
    #$argument_hashref->{n} = 10000;

    # parse @ARGV
    $optcatchresult = GetOptions(
        'out=s'              => \$argument_hashref->{out},
        #'min=f'              => \$argument_hashref->{min},  
        #'step=f'             => \$argument_hashref->{step}, 
        #'n=i'                => \$argument_hashref->{n}, 
        'verbose!'           => \$verbose,
        'help'               => \$help,
        'man'                => \$man
                                );

    if ($optcatchresult == 0) { return 0 }

    if (defined($help)) { return 2 }

    if (defined($man)) { return 3 }

    if ($verbose)
    {
#        PBar::enable();
        enable_verbose;
    }

    return 1;
}

sub predict_sequence(\$ $ $ $ $ $ $ $ $ $ $ $)
{
    my ($struct_ref, $neg, $het, $T, $ed, $lambda, $b_l, $b_s, $spu, $mc, $sc, $steps)
        = @_;
    my $seq;
    my $struct_len = length($$struct_ref);

    unless(open(FH, "./corb \"brot -q $sc -j $mc -u $spu -i $b_s -o $b_l -l $lambda -e $ed -t $T -s $steps -d $neg -h $het $$struct_ref\" |"))
    {
        msg_error_and_die ("Could not start brot\n");
    }

    foreach (<FH>)
    {
        if ($_ =~ /([ACGUacgu]{$struct_len})/)
        {
            close(FH);
            return $1;
        }
    }
    close(FH);

    msg_error_and_die ("No sequence found for \"brot -s $steps -d $neg -h "
                      ."$het $$struct_ref\"\n");
}

sub eval_sequence(\$ \$)
{
    my ($struct_ref, $seq_ref) = @_;

    unless(open(FH, "./corb \"er2de $$seq_ref $$struct_ref\" |"))
    {
        msg_error_and_die ("Could not start er2de\n");
    }

    foreach (<FH>)
    {
        if ($_ =~ /G\s+\=\s+(\-?\d+\.\d+)/)
        {
            close (FH);
            return $1;
        }
    }
    close(FH);

    #msg_warning ("Could not evaluate energy of \"$$seq_ref\",".
    #                   " \"$$struct_ref\"\n");
    return 1000000;
}

sub predict_mfe_structure (\$)
{
    my ($seq_ref) = @_;
    my $seq_len = length($$seq_ref);
    my $result_count = 0;
    my $structure;
    my $mfe_en;
    my $ens_en;

    unless(open(FH, "echo \"$$seq_ref\" | ./RNAfold -d2 -p0 2>&1 |"))
    {
        msg_error_and_die ("Could not start RNAfold\n");
    }

    foreach (<FH>)
    {
        #print $_;
        if ($_ =~ /([\(\.\)]{$seq_len})\s+\(\s*(\-?\d+\.\d+)\)/)
        {
            $result_count++;
            $structure = $1;
            $mfe_en = $2;
        }
        elsif ($_ =~ /free\s+energy\s+of\s+ensemble\s+\=\s+(\-\d+\.\d+)/)
        {
            $result_count++;

            $ens_en = $1;
        }
    }

    close(FH);

    if (!defined ($structure))
    {
        msg_error_and_die ("No structure found found for sequence "
                          ."\"$$seq_ref\"\n");
    }
    elsif (!defined ($ens_en))
    {
        msg_error_and_die ("No ensemble energy found found for sequence "
                          ."\"$$seq_ref\"\n");        
    }

    #print "\n\n".$structure.$mfe_en.$ens_en."\n\n";

    return ($structure, $mfe_en, $ens_en);
}

sub compare_structures (\$ \$)
{
    my ($ref_ref, $pred_ref) = @_;

    unless(open(FH, "echo \"$$ref_ref\n$$pred_ref\" | ./RNAdistance  2>&1 |"))
    {
        msg_error_and_die ("Could not start RNAdistance\n");
    }

    foreach (<FH>)
    {
        #print $_;
        if ($_ =~ /f\:\s*(\d+)/)
        {
            close(FH);
            return $1;
        }
    }

    close(FH);

    msg_error_and_die ("Could not compare \"$$ref_ref\" and \"$$pred_ref\"");   
}

sub cost_func (\@ \@ \@)
{
    my ($param_aryref, $aryref_2, $aryref_3) = @_;
    my ($neg,
        $het,
        $T,
        $ed,
        $lambda,
        $b_l,
        $b_s,
        $spu,
        $mc, $sc, $steps)
        = @$param_aryref;
    my $predicted_seq;
    my $predicted_struct;
    my @prediction;
    my $dist = 0;
    my $costs = 0;
    my $target_p;
    my $target_en;
    my $fh;
    my $nerr = 0;

    #$steps = floor ($steps);

    printf "Testing %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ".
           "%5.2f %d: ",
    $neg, # 0.746094,
    $het,     # 1
    $T, # 13.6719,
    $ed,      # 0.2,
    $lambda,  # 0.625
    $b_l,     # 0.95,
    $b_s,     # 0.6
    $spu,      # 0.99,
    0.85,     #   $mc,
    0.99,     #   $sc;
    50000;

    $fh = open_or_die('artificial_small_tests', "<");

    #foreach (keys(%SMALL_TESTS))
    while (<$fh>)
    {
        #print $_;
        chomp $_;
        # design sequence
        #print "Testing $neg $het $steps on \"$SMALL_TESTS{$_}\"\n";
        $predicted_seq = predict_sequence(#$SMALL_TESTS{$_},
                                          $_,
                                          $neg,     # 0.746094,
                                          $het,     # 1
                                          $T,       # 13.6719,
                                          $ed,      # 0.2,
                                          $lambda,  # 0.625,
                                          $b_l,     # 0.95,
                                          $b_s, # 0.6,
                                          $spu, # 0.99,
                                          0.85,     #     $mc,
                                          0.99,     #     $sc,
                                          50000     #     $steps
            );

        $target_en = eval_sequence(#$SMALL_TESTS{$_},
                                   $_,
                                   $predicted_seq);

        #print "Seq: $predicted_seq\n";

        # predict structure
        @prediction = predict_mfe_structure ($predicted_seq);
        $predicted_struct = $prediction[0];
        #print "MFE: \"$predicted_struct\", $prediction[1], $prediction[2]\n";

        # evaluate
        $dist = compare_structures (#$SMALL_TESTS{$_},
                                    $_,
                                    $predicted_struct);

        # ("target_p         = %f\n", exp((ens_en - target_en) / RT));
        if ($target_en != 1000000)
        {
            $target_p = exp (($prediction[2] - $target_en) / 1);
        }
        else
        {
            $target_p = 0;
            #print "$target_en\n";
            #return $target_en;
        }
        #print "en: $target_en p: $target_p\n";

        # costs
        # dist * (1-target_p)
        $costs += ($dist * (1 - $target_p));

        if ($dist == 0)
        {
            $nerr++;
        }
        else
        {
            msg_warning ("Seq/ Struct w err's:\n$_\n$predicted_seq\n");
        }
    }
    print "$costs $nerr\n";

    close ($fh);

    #print @$aryref_1;
    #print "\n\n";
    #print $#$aryref_2;
    #print "\n\n";
    #print $#$aryref_3;
    #print "\n\n";

    return $costs;
}

# FUNCTIONS        - END


# MAIN             - BEGIN
my %arg_hash;
my $pod_verbose = 1;
my $ret_val;
my %result;
#                d   h   t    e         l    o            i    u     j     q s
#my @ini_guess = (5,  5,  10,  0.25,      0.6, 0.95,        0.5, 0.99, 0.85, 0.99, 500000);
my @ini_range = (
                 0.05,           # neg       
                 0.01,           # het     0.25  
                 0.5,               # temp      
                 0.00,        # s_drop    
                 0.01,         # lambda
                 0.005,        # b_long    
                 0.01,         # b_short   
                 0.1,        # speedup   
#                 0.05,        # min-cool  
#                 0.05,        # cool-scale
#                 500000       # steps     
    );

my @ini_guess = (
                 0.25,          # neg       
                 1.075,            # het 1     
                 30,           # temp      
                 0.05,          # s_drop   0.19375 0.2
                 0.83,        # lambda 0.85
                 0.95,        # b_long    
                 0.55,         # b_short   
                 1.0,        # speedup   0.75
#                 0.75,        # min-cool  
#                 0.75,        # cool-scale
#                 50000        # steps     
                 );

my @low_bound = (
                 0.1,           # neg
                 1,           # het       
                 29,           # temp      
                 0.05,       # s_drop    
                 0.81,           # lambda  0.7  
                 0.9,         # b_long    
                 0.5,         # b_short   
                 0.5,         # speedup   
#                 0.5,         # min-cool  
#                 0.5,         # cool-scale
#                 10000        # steps
    );

my @upp_bound = (
                 0.4,           # neg       
                 1.2,          # het 5     
                 31,         # temp      
                 0.06,         # s_drop    
                 0.85,      # lambda    0.9
                 0.99,        # b_long    0.99
                 0.6,        # b_short   
                 1.0,   # speedup   
#                 0.99,        # min-cool  
#                 0.9990000,   # cool-scale
#                 1000000    # steps
    );

my @names     = ('neg', 'het', 'temp', 's_drop', 'lambda', 'b_long', 'b_short',
                 'speedup', 'min-cool', 'cool-scale', 'steps');
my %s_arg     = (
#   These are obligatory entries
    func         => \&cost_func,
    ini_pt       => \@ini_guess,
    max_iter     => 10000,
    max_restart  => 3,
    f_tol        => 1e-5,
#       These are optional
    names        => \@names,
    ini_range    => \@ini_range,
    lower        => \@low_bound,
    upper        => \@upp_bound,
    scatter      => 0.20,
    );

# parse commandline
$ret_val = parseargs(%arg_hash);
if ($ret_val > 1)
{
    if ($ret_val == 3) { $pod_verbose = 2 }
    pod2usage(-exitval => 0, -verbose => $pod_verbose); 
}

if (defined ($arg_hash{out}))
{
    $s_arg{o_file} = $arg_hash{out};
}

# output params
for (my $i = 0; $i <= $#ini_guess; $i++)
{
    print "$names[$i] $ini_guess[$i] $ini_range[$i] $low_bound[$i] "
         ."$upp_bound[$i]\n";
}

# init random number generator
# srand ();

# run simplex
%result = simplex (%s_arg);

# evaluate result
if ($result {success} == $SPLX_SUCCESS)
{
    print "Simplex happy \n";
}
elsif ($result {success} == $SPLX_TOO_MANY)
{
    print "Simplex too many\n";
}
elsif ($result {success} == $SPLX_BROKEN)
{
    msg_warning ("Simplex broken\n");
}
else
{
    msg_error_and_die ("Unresolved programming bug\n");
}

# feedback
print $result{restart}, ' restarts and ', $result{ncycle}, " cycles\n";
print $result{value}, " best value of cost function\n";
#$DB::single = 1;
for (my $i = 0; $i < $s_arg{n_dim}; $i++)
{
    printf "%5s: %4g\n", $names[$i], $result{best}[$i];
}
print "\n";

# MAIN             - END


__END__

=head1 NAME

simplex_brot_params - Parameter optimisation using the simplex algorithm.

=head1 SYNOPSIS

B<simplex_brot_params> [options]

=head1 DESCRIPTION

bla

=head1 OPTIONS

=over 8

=item B<-out C<prefix>>

Prefix for the output of the simplex runs. If given, two files C<prefix>_hi.out
and C<prefix>_low.out are created. The files may not exist before the script is
started.

=item B<-exteriorsonly>

Only test exterior loops.

=item B<-bulgesonly>

Only test bulge loops.

=item B<-interiorsonly>

Only test interior loops.

=item B<-hairpinsonly>

Only test hairpin loops.

=item B<-multionly>

Only test multiloops.

=item B<-testsonly>  

Only process set of test sequences.

=item B<-man>

Print a man page for the program.

=item B<-help>

Print a help message and exit.

=cut

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
