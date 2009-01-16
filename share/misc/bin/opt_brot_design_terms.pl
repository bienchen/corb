#!@PERL@ -w
# -*- perl -*-
# @configure_input@
# Last modified: 2009-01-13.17


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
use PBar;
use CorbIO qw(:all);
use RNA qw(:all);
# PRIVATE packages - END

# CONSTANTS        - BEGIN
my %SMALL_TESTS = ('a' => '((....))',
                   'aa' => '((((((...(((......)))...(((...)))...))))))',
                   'aaa' => '((((....(((...)))...(((((((((....)))))))))......))))',
                   'aaaa' => '(((...(((.......((((.....))))))))))',
                   'aaaaa' => '(((((...(((...)))(((())))...(((...(...))))(((((....)))))..(...).((()))...)))))',
                   'aaaaaa' => '...(((...)))(((())))...(((...(...))))',
                   'aaaaaaa' => '(((((....)))))..(...).((()))',
                   'aaaaaaaa' => '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....');
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
    $argument_hashref->{max} = 10.0;
    $argument_hashref->{min} = 0.0;
    $argument_hashref->{step} = 0.01;
    $argument_hashref->{n} = 10000;

    # parse @ARGV
    $optcatchresult = GetOptions(
        'max=f'              => \$argument_hashref->{max},
        'min=f'              => \$argument_hashref->{min},  
        'step=f'             => \$argument_hashref->{step}, 
        'n=i'                => \$argument_hashref->{n}, 
        'verbose!'           => \$verbose,
        'help'               => \$help,
        'man'                => \$man
                                );

    if ($optcatchresult == 0) { return 0 }

    if (defined($help)) { return 2 }

    if (defined($man)) { return 3 }

    if ($verbose)
    {
        PBar::enable();
        enable_verbose;
    }

    return 1;
}

sub predict_sequence(\$ $ $ $)
{
    my ($struct_ref, $neg, $het, $steps) = @_;
    my $seq;
    my $struct_len = length($$struct_ref);

    unless(open(FH, "../../../src/corb \"brot -s $steps -d $neg -h $het $$struct_ref\" |"))
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

sub predict_mfe_structure (\$)
{
    my ($seq_ref) = @_;
    my $seq_len = length($$seq_ref);

    unless(open(FH, "echo \"$$seq_ref\" | ./RNAfold -d2  2>&1 |"))
    {
        msg_error_and_die ("Could not start RNAfold\n");
    }

    foreach (<FH>)
    {
        if ($_ =~ /([\(\.\)]{$seq_len})/)
        {
            close(FH);
            return $1;
        }
    }

    close(FH);

    msg_error_and_die ("No structure found found for sequence \"$$seq_ref\"\n");
}

sub compare_structures (\$ \$)
{
    my ($ref_ref, $pred_ref) = @_;

    unless(open(FH, "echo \"$$ref_ref\n$$pred_ref\" | ./RNAdistance  2>&1 |"))
    {
        msg_error_and_die ("Could not start RNAfold\n");
    }

    foreach (<FH>)
    {
        print $_;
        if ($_ =~ /f\:\s*(\d+)/)
        {
            close(FH);
            return $1;
        }
    }

    close(FH);

    msg_error_and_die ("Could not compare \"$$ref_ref\" and \"$$pred_ref\"");   
}

sub evaluate_params($ $ $)
{
    my ($neg, $het, $steps) = @_;
    my $result = 0;
    my $predicted_seq;
    my $predicted_struct;
    my $dist = 0;

    # for all test sequences
    foreach (keys(%SMALL_TESTS))
    {
        print "Testing $neg $het on \"$SMALL_TESTS{$_}\"\n";
        $predicted_seq = predict_sequence($SMALL_TESTS{$_}, $neg, $het, $steps);
        print "Seq: $predicted_seq\n";
        $predicted_struct = predict_mfe_structure ($predicted_seq);
        print "MFE: \"$predicted_struct\"\n";
        $dist += compare_structures ($SMALL_TESTS{$_}, $predicted_struct);
        #print "end\n";
    }

    $result = $dist * 0.5;

    return $result;
}
# FUNCTIONS        - END


# MAIN             - BEGIN
my %arg_hash;
my $pod_verbose = 1;
my %param_hash;
my $ret_val;
my ($p, $q, $r);
my $result;
my $best_result = 100000;
my $best_neg;
my $best_het;
my %pbar;

# parse commandline
$ret_val = parseargs(%arg_hash);
if ($ret_val > 1)
{
    if ($ret_val == 3) { $pod_verbose = 2 }
    pod2usage(-exitval => 0, -verbose => $pod_verbose); 
}

msg_verbose("Start of parameter range:  ${arg_hash{min}}\n");
msg_verbose("End of parameter range:    ${arg_hash{max}}\n");
msg_verbose("Step size:                 ${arg_hash{step}}\n");
msg_verbose("No. of steps per run:      ${arg_hash{n}}\n");

# single evaluation

msg_verbose("\nEvaluating negative design term, only:");
PBar::start(%pbar, $arg_hash{max});
for ($p = $arg_hash{min}; $p <= $arg_hash{max}; $p += $arg_hash{step})
{
    #print $p."\n";
    $result = evaluate_params($p, 0, $arg_hash{n});

    if ($result < $best_result)
    {
        $best_result = $result;
        $best_neg = $p;
    }

    PBar::update(%pbar, $p);
}
PBar::finish(%pbar);
print "Best parameter for negative design term, single run: $best_neg "
     ."($best_result)\n";

msg_verbose("\nEvaluating heterogenity term, only:");
$best_result = 100000;
PBar::start(%pbar, $arg_hash{max});
for ($p = $arg_hash{min}; $p <= $arg_hash{max}; $p += $arg_hash{step})
{
    $result = evaluate_params(0, $p, $arg_hash{n});

    if ($result < $best_result)
    {
        $best_result = $result;
        $best_het = $p;
    }

    PBar::update(%pbar, $p);
}
PBar::finish(%pbar);
print "Best parameter for heterogenity term, single run: $best_het "
     ."($best_result)\n";

msg_verbose("\nEvaluating both terms:");
PBar::start(%pbar, $arg_hash{max}*$arg_hash{max});
for ($p = $arg_hash{min}; $p <= $arg_hash{max}; $p += $arg_hash{step})
{
    for ($q = $arg_hash{min}; $q <= $arg_hash{max}; $q += $arg_hash{step})
    {
        $result = evaluate_params($p, $q, $arg_hash{n});
        
        if ($result < $best_result)
        {
            $best_result = $result;
            $best_neg = $p;
            $best_het = $q;
        }
        $r += $arg_hash{step};
        PBar::update(%pbar, $r);
    }
}
PBar::finish(%pbar);
print "Best parameters: neg($best_neg), het($best_het) ($best_result)\n";

# MAIN             - END


__END__

=head1 NAME

ex_cmp_er2de_rnaeval - Exhaustive comp. of er2de and RNAeval...

=head1 SYNOPSIS

B<ex_cmp_er2de_rnaeval> [options]

=head1 DESCRIPTION

B<ex_cmp_er2de_rnaeval> should be used to test the implementation of the Turner
energies in  B<corb> against the Vienna Package as gold standard. All possible
combinations of sequences/ structures are tested as are available in
Vienna 1.7.2. Additionally a few hundred test sequences are also checked.
Please note that the "only" options may actually be combined.

=head1 OPTIONS

=over 8

=item B<-stackssonly>

Only test stacking parameters.

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
