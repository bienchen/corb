#!@PERL@ -w
# -*- perl -*-
# @configure_input@
# Last modified: 2009-03-07.10


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
use sigtrap qw{handler my_sig_handler normal-signals error-signals};
# PERL packages    - END


# PRIVATE packages - BEGIN
use lib '@CORB_PERL5LIB@';
use PBar;
use CorbIO qw(:all);
use RNA qw(:all);
# PRIVATE packages - END

# CONSTANTS        - BEGIN
# CONSTANTS        - END

# GLOBALS          - BEGIN
my %sig_handler_params;
# GLOBALS          - END


# FUNCTIONS        - BEGIN
sub my_sig_handler
{
    my $sig = shift;
    my $fh;
    my $file = "brute_force_testset_small";

    $sig = defined $sig ? $sig : "????";

    msg_warning ("Signal caught: SIG$sig.\n");

    if (   (defined ($sig_handler_params{'current'}))
        && (defined ($sig_handler_params{'names'}))
        && (defined ($sig_handler_params{'winners'}))
        )
    {
        print "Current state: ";
        foreach (@{$sig_handler_params{'names'}})
        {
            print "$_: $sig_handler_params{'current'}->{$_} ";
        }
        print "\n";
    }

    if (   (defined ($sig_handler_params{'names'}))
        && (defined ($sig_handler_params{'winners'}))
        )
    {
        report_winners ($sig_handler_params{'winners'},
                        $sig_handler_params{'names'});
    }

    die ("Bye!\n");
}

sub increment_params (\@ \% \% \% \%)
{
    my ($names_aryref, $current_hashref, $stop_hashref, $start_hashref,
        $increment_hashref) = @_;
    my $n_names = $#{$names_aryref};
    my $i;
    my $param;

    for ($i = $n_names; $i >= 0; $i--)
    {
        $param = $$names_aryref[$i];
        #print "$param\n";

        $current_hashref->{$param} += $increment_hashref->{$param};

        if (   $current_hashref->{$param} > $stop_hashref->{$param})
        {
            $current_hashref->{$param} = $start_hashref->{$param};
        }
        else
        {
            return 0;
        }
    }

    return 1;
}

sub push_win_params (\% \% $)
{
    my ($target_hashref, $source_hashref, $costs) = @_;

    foreach (keys (%{$source_hashref}))
    {
        push (@{$target_hashref->{$_}}, $source_hashref->{$_});
    }
    push (@{$target_hashref->{'__min'}}, $costs);
}

sub report_winners (\% \@)
{
    my ($win_hashref, $names_aryref) = @_;
    my $results = -1;
    my $i;

    if (defined ($win_hashref->{'__min'}))
    {
        $results = $#{$win_hashref->{'__min'}};
    }

    print "Good param combinations: ".($results + 1)."\n";

    for ($i = 0; $i <= $results; $i++)
    {
        print "#".($i + 1), ": ${$win_hashref->{'__min'}}[$i]\n";
        foreach (@{$names_aryref})
        {
            print "$_: ${$win_hashref->{$_}}[$i]\n";
        }
    }

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

sub predict_mfe_structure (\$)
{
    my ($seq_ref) = @_;
    my $seq_len = length($$seq_ref);

    unless(open(FH, "echo \"$$seq_ref\" | ./RNAfold -d2 -p0 2>&1 |"))
    {
        msg_error_and_die ("Could not start RNAfold\n");
    }

    foreach (<FH>)
    {
        #print $_;
        if ($_ =~ /([\(\.\)]{$seq_len})\s+\(\s*(\-?\d+\.\d+)\)/)
        {
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

sub eval_params (\% \@)
{
    my ($param_hashref, $tests_aryref) = @_;
    my $pseq;
    my $pstruct;
    my $dist;
    my $i = 0;
    my $tmp;
    my $succount = 0;
    my $neg = $param_hashref->{'neg'};
    my $het = $param_hashref->{'het'};
    my $T   = $param_hashref->{'T'};
    my $ed  = 0.05;
    my $l   = $param_hashref->{'l'};
    my $b_l = $param_hashref->{'b_long'};
    my $b_s = $param_hashref->{'b_short'};
    my $spu = $param_hashref->{'sp_thresh'};
    my $mc  = $param_hashref->{'min_cool'};
    my $sc  = $param_hashref->{'scale_cool'};
    my $s   = 10000;

    #print "Testing $neg $het $ed $l $b_l $b_s $spu $mc $sc $s\n";

    foreach (<@{$tests_aryref}>)
    {
        #print "$i: $_ $#{$tests_aryref}\n";
        # predict sequence
        $pseq = predict_sequence($_,
                                 $neg,
                                 $het,
                                 $T,
                                 $ed,
                                 $l,
                                 $b_l,
                                 $b_s,
                                 $spu,
                                 $mc,
                                 $sc,
                                 $s);
        
        # predict mfe structure
        $pstruct = predict_mfe_structure ($pseq);

        # compare mfe structure
        $dist = compare_structures ($_, $pstruct);

        #print $_."\n$pseq\n$pstruct\n$dist\n";
        if ($dist > 0)
        {
            #print $_."\n$pseq\n$pstruct\n$dist\n";
            #print "$dist\n";
            # put the failed structure at the beginning of the list
            #$tmp = $_;
            if ($i != 0)
            { 
                $tmp = splice (@{$tests_aryref}, $i, 1);
                unshift (@{$tests_aryref}, $tmp);
            }

            return $succount;
        }
        if ($dist == 0)
        {
            $succount++;
            #print $_."\n$pseq\n$pstruct\n$dist\n"; 
        }

        $i++;
    }

    return $succount;
}

# FUNCTIONS        - END


# MAIN             - BEGIN
my $full_stop = 0;
my $max = 0.1;
my $cur = 0;
my %win_params;
my @tests;
my $fh;
my $timestamp;

# params
my @names = (
#        'ed',
        'l',
        'b_long',
        'b_short',
        'sp_thresh',
        'min_cool',
        'scale_cool',
        'T',
        'neg',
        'het',
    );
my %start = (
        'neg'        => 0,
        'het'        => 0,
        'T'          => 70, #0,
#        'ed'        => 0,
        'l'          => 0,
        'b_long'     => 0,
        'b_short'    => 0,
        'sp_thresh'  => 0.9,
        'min_cool'   => 0.5,
        'scale_cool' => 0.1,
    );
my %stop = (
        'neg'        => 10,
        'het'        => 10,
        'T'          => 75, #100,
#        'ed'        => 0.1,
        'l'          => 0.9,
        'b_long'     => 1,
        'b_short'    => 1,
        'sp_thresh'  => 1,
        'min_cool'   => 1,
        'scale_cool' => 1,
    );
my %increment = (
        'neg'        => 1,
        'het'        => 1,
        'T'          => 5,
#        'ed'        => 0.02,
        'l'          => 0.1,
        'b_long'     => 0.05,
        'b_short'    => 0.05,
        'sp_thresh'  => 0.01,
        'min_cool'   => 0.05,
        'scale_cool' => 0.05,
    );
my %current;

# 'min-cool', 'cool-scale', 'steps'

# init current param values
if ($#ARGV == 17)
{
    %current = (
        'neg'        => $ARGV[0],
        'het'        => $ARGV[1],
        'T'          => $ARGV[2],
        'l'          => $ARGV[3],
        'b_long'     => $ARGV[4],
        'b_short'    => $ARGV[5],
        'sp_thresh'  => $ARGV[6],
        'min_cool'   => $ARGV[7],
        'scale_cool' => $ARGV[8],
        );

    %stop = (
        'neg'        => $ARGV[9],
        'het'        => $ARGV[10],
        'T'          => $ARGV[11],
        'l'          => $ARGV[12],
        'b_long'     => $ARGV[13],
        'b_short'    => $ARGV[14],
        'sp_thresh'  => $ARGV[15],
        'min_cool'   => $ARGV[16],
        'scale_cool' => $ARGV[17],
    );
}
else
{
    foreach (@names)
    {
        $current{$_} = $start{$_};
    }
}

print "Start params for search:\n";
print "neg:        $current{'neg'}\n";
print "het:        $current{'het'}\n";
print "T:          $current{'T'}\n";
print "l:          $current{'l'}\n";
print "b_long:     $current{'b_long'}\n";
print "b_short:    $current{'b_short'}\n";
print "sp_thresh:  $current{'sp_thresh'}\n";
print "min_cool:   $current{'min_cool'}\n";
print "scale_cool: $current{'scale_cool'}\n\n";
print "Stop params for search:\n";
print "neg:        $stop{'neg'}\n";
print "het:        $stop{'het'}\n";
print "T:          $stop{'T'}\n";
print "l:          $stop{'l'}\n";
print "b_long:     $stop{'b_long'}\n";
print "b_short:    $stop{'b_short'}\n";
print "sp_thresh:  $stop{'sp_thresh'}\n";
print "min_cool:   $stop{'min_cool'}\n";
print "scale_cool: $stop{'scale_cool'}\n";

#%current = (
#        'neg'        => 7,
#        'het'        => 0,
#        'T'          => 40,
##        'ed'        => 0.02,
#        'l'          => 0,
#        'b_long'     => 0.8,
#        'b_short'    => 0.25,
#        'sp_thresh'  => 0.91,
#        'min_cool'   => 0.95,
#        'scale_cool' => 0.95,
#    );

$current{'l'} = 0.8;

$sig_handler_params{'names'} = \@names;
$sig_handler_params{'current'} = \%current;
$sig_handler_params{'winners'} = \%win_params;

# read list of test structures
$fh = open_or_die('artificial_small_tests', "<");
@tests = <$fh>;
close ($fh);

# fetch first timestamp
$timestamp = time;

# while not all params reached their upper boundary
while (! $full_stop)
{

    # evaluate
    $cur = eval_params (%current, @tests);

    if ($cur == $max)
    {
        push_win_params (%win_params, %current, $cur);
        $max = ($cur + $max)/2;
    }
    elsif ($cur > $max)
    {
        %win_params = ();
        push_win_params (%win_params, %current, $cur);

       $max = $cur;
    }

    # output after certain time
    if ((time - $timestamp) > 3600)
    {
        print "Progress report (". localtime() .")\n";
        print "------------------------------------------\n";
        report_winners (%win_params, @names);
        print "Tested until: ";
        foreach (@names)
        {
            print "$_: $current{$_} ";
        }
        print "\n";

        $timestamp = time;
    }
        
    # increment
    $full_stop = increment_params (@names, %current, %stop, %start, %increment);
}

# report
report_winners (%win_params, @names);

# MAIN             - END

#    if (($cur < ($min + $delta)) && ($cur > ($min - $delta)))
#    {
#        push_win_params (%win_params, %current, $cur);
#        $min = ($cur + $min)/2;
#    }
#    elsif ($cur < $min)
#    {
#        %win_params = ();
#        push_win_params (%win_params, %current, $cur);
#
#       $min = $cur;
#    }

__END__

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
