#!@PERL@ -w
# -*- perl -*-
# @configure_input@
# Last modified: 2009-01-05.10


# Copyright (C) 2008 Stefan Bienert
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

# Todo:
#      - make changing ALL params by default
#        Problem: If not, you will forget changing the stack_mm params which
#        are not checked by ex_... (er2de vs. RNAeval)!!!

# SETTINGS         - BEGIN
BEGIN
{
  # try to get bash compatible shell
  $ENV{'SHELL'} = '@SHELL@' if exists $ENV{'DJGPP'};
  $|=1; # dissable line buffer for print statements
}
# SETTINGS         - END


# PERL packages    - BEGIN
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use File::Basename;
#use File::Copy;
# PERL packages    - END


# PRIVATE packages - BEGIN
use CorbIO qw(:all);
# PRIVATE packages - END

# CONSTANTS        - BEGIN
my $A = 0;
my $C = 1;
my $G = 2;
my $U = 3;
my @ALPHA = ($A, $C, $G, $U);
my @BASEPAIRS = ('cg', 'gc', 'gu', 'ug', 'au', 'ua');
my @BASES = ('a', 'c', 'g', 'u');
my $UNDEF_ENERGY = 'FLOAT_UNDEF';
# CONSTANTS        - END

# GLOBALS          - BEGIN
# GLOBALS          - END


# FUNCTIONS        - BEGIN
# calculate the decadic logarithm
#   log10(number)
sub log10($)
{
    my ($n) = @_;
    return log($n)/log(10);
}

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

    # parse @ARGV
    $optcatchresult = GetOptions(
        'stacking_energies!' => \$argument_hashref->{stacking_energies},
        'int11_energies!'    => \$argument_hashref->{int11_energies},
        'int21_energies!'    => \$argument_hashref->{int21_energies},
        'int22_energies!'    => \$argument_hashref->{int22_energies},
        'mismatch_interior!' => \$argument_hashref->{mismatch_interior},
        'mismatch_hairpin!'  => \$argument_hashref->{mismatch_hairpin},
        'mismatch_stack!'    => \$argument_hashref->{mismatch_stack},
        'internals!'         => \$argument_hashref->{internals},
        'hairpins!'          => \$argument_hashref->{hairpins},
        'bulges!'            => \$argument_hashref->{bulges},      
        'dangle5!'           => \$argument_hashref->{dangle5},
        'dangle3!'           => \$argument_hashref->{dangle3},
        'tetraloops!'        => \$argument_hashref->{tetraloops},
        'sourcefile=s'       => \$argument_hashref->{sourcefile},
        'verbose!'           => \$verbose,
        'help'               => \$help,
        'man'                => \$man
                                );

    if ($optcatchresult == 0) { return 0 }

    if (defined($help)) { return 2 }

    if (defined($man)) { return 3 }

    if ($verbose)
    {
        enable_verbose;
    }

    # check that one param file was given
    if ($#ARGV != 0)
    {
        msg_error_and_die("Parameter file missing.\n"
                         ."Try \"-help\" or \"-man\" for more information.\n");
    }

    # fetch file name
    if ( ! -r $ARGV[0])
    {
        msg_error_and_die("Param file does not exist or is not readable: "
                         ."$ARGV[0]\n");
    }    
    $argument_hashref->{paramfile} = $ARGV[0];

    if (defined ($argument_hashref->{sourcefile}))
    {
        if ((! -r $argument_hashref->{sourcefile})
            || (! -w $argument_hashref->{sourcefile}))
        {
            msg_error_and_die("File \'$argument_hashref->{sourcefile}\'does "
                             ."not exist or is not read- and/ or "
                             ."writeable.\n");
        }
    }

    return 1;
}

sub slice_tags_or_die($ $ \@ \@)
{
    my ($open_tag, $close_tag, $file, $pslice) = @_;
    my $i;
    my $start;
    my $end;

    for ($i = 0; $i < $#{$file} + 1; $i++)
    {
        if (@{$file}[$i] =~ /${open_tag}/)
        {
            $start = $i + 1;
        }
        elsif (@{$file}[$i] =~ /${close_tag}/)
        {
            $end = $i;
        }
    }

    if ((defined ($start)) && (defined ($end)))
    {
        # splice everything between the tags from array
        splice(@{$file}, $start, $end - $start, @{$pslice});
    }
    else
    {
        msg_error_and_die("Tags \"$open_tag\" and \"$close_tag\" not found in "
                          ."source file.\n");
    }
}

sub read_en_matrix_offbone(\*)
{
    my ($FH) = @_;
    my $i = 0;
    my @table = ();

    while (<$FH>)
    {
        if ($i > 0)
        {
            if ($_ =~ /\s+-?\d+\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/)
            {
                $table[$i-1][$A] = $1;
                $table[$i-1][$C] = $2;
                $table[$i-1][$G] = $3;
                $table[$i-1][$U] = $4;
            }
            else
            {
                msg_error_and_die ("Unrecognised line found: $_");
            }
        }
        if ($i == 4)
        {
            return @table;
        }
        $i++;
    }
}

sub read_en_matrix_stream_offbone(\* $)
{
    my ($FH, $rows) = @_;
    my $i = 1;
    my @table = ();
    my ($a, $c, $g, $u);

    while (<$FH>)
    {
        if ($i > 0)
        {
            if ($_ =~ /\s+(?:-?\d+|INF)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+
                       (-?\d+)/x)
            {
                $table[$i-1][$A] = $1;
                $table[$i-1][$C] = $2;
                $table[$i-1][$G] = $3;
                $table[$i-1][$U] = $4;
            }
            else
            {
                msg_error_and_die ("Unrecognised line found: $_");
            }
        }
        if ($i == $rows)
        {
            return @table;
        }
        $i++;
    }
}

sub read_en_matrix(\*)
{
    my ($FH) = @_;
    my $i = 0;
    my @table = ();

    while (<$FH>)
    {
        if ($_ =~ /\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/)
        {
            $table[$i][$A] = $1;
            $table[$i][$C] = $2;
            $table[$i][$G] = $3;
            $table[$i][$U] = $4;
        }
        else
        {
            msg_error_and_die ("Unrecognised line found: $_");
        }

        if ($i == 3)
        {
            return @table;
        }
        $i++;
    }
}

sub fetch_stacking_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my $i = 0;

    while (<$FH>)
    {
        # start reading the matrix after the comment line
        if ($_ =~ /\s*\/\*\s*[AUGCaugc]{2}\s*[AUGCaugc]{2}\s*[AUGCaugc]{2}
                          \s*[AUGCaugc]{2}\s*[AUGCaugc]{2}\s*[AUGCaugc]{2}
                          \s*\@\s*\*\//x)
        {
            # just skip comment line
        }
        elsif ($_ =~/\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)
                     \s+(-?\d+)/x)
        {
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[0]}
            = $1;
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[1]}
            = $2;
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[2]}
            = $3;
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[3]}
            = $4;
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[4]}
            = $5;
            $param_hashref->{'stacking_energies'}{$BASEPAIRS[$i]}{$BASEPAIRS[5]}
            = $6;

            $i++;
            if ($i == 6)
            {
                return;
            }
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in stacking_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_int11_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;
    my ($bp1, $bp2);

    while (<$FH>)
    {
        # fetch bp + base
        #/* GC..GC */
        if ($_ =~ /\/\*\s*([ACGUacgu\@\s]{2})\.\.([ACGUacgu\@\s]{2})\s*\*\//x)
        {
            @table = read_en_matrix_offbone(*$FH);

            $bp1 = $1;
            $bp2 = $2;

            if (lc($bp1) eq "gc")
            {
                $bp1 = "cg"
            }
            elsif (lc($bp1) eq "cg")
            {
                $bp1 = "gc";
            }
            if (lc($bp2) eq "gc")
            {
                $bp2 = "cg"
            }
            elsif (lc($bp2) eq "cg")
            {
                $bp2 = "gc";
            }

            @{$param_hashref->{'int11_energies'}{lc($bp1)}{lc($bp2)}}
            = @table;
        }
        elsif ($_ =~ /^\s*$/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in int11_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_int21_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;

    while (<$FH>)
    {
        # fetch bp + base
        if ($_ =~ /\/\*\s+([ACGUacgu\@\s]{2})\.([ACGUacgu\@])\.+
                   ([ACGUacgu\@\s]{2})\s+/x)
        {
            @table = read_en_matrix_offbone(*$FH);


            @{$param_hashref->{'int21_energies'}{reverse(lc($1))}
                                                {reverse(lc($3))}{lc($2)}}
            = @table;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in int21_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_int22_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;
    my ($bp1, $bp2);
    my ($b1, $b2);

    while (<$FH>)
    {
        # fetch bp + base
        if ($_ =~ /\/\*\s+([\s\@ACGUacgu]{2})\.([ACGUacgu])([ACGUacgu])\.+
                   ([ACGUacgu\s\@]{2})\s+/x)
        {
            @table = read_en_matrix(*$FH);

            $bp1 = $1;
            $bp2 = $4;
            $b1 = $2;
            $b2 = $3;

            # possibilities:
            # $1 $4                      caacaaagacg 6.8 6.3    
            # reverse($1) $4                
            # $1 reverse($4)              
            # reverse($1) reverse($4)       
            # $4 $1                      
            # reverse($4) $1                    
            # $4 reverse($1)                 
            # reverse($4) reverse($1)       

            # transform bp due to error in Vienna1.72 int22 table
            if (lc($bp1) eq "gc")
            {
                $bp1 = "cg";
            }
            elsif (lc($bp1) eq "cg")
            {
                $bp1 = "gc";
            }
            if (lc($bp2) eq "gc")
            {
                $bp2 = "cg";
            }
            elsif (lc($bp2) eq "cg")
            {
                $bp2 = "gc";
            }
                                                 
            @{$param_hashref->{'int22_energies'}{lc($bp1)}
                                                {lc($bp2)}               
                                                {lc($b1)}{lc($b2)}}
            = @table;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in int22_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_mismatch_interior_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;
    my $bpcount = 0;

    while (<$FH>)
    {
        if ($_ =~ /\s+-?\d+\s+-?\d+\s+-?\d+\s+-?\d+\s+-?\d+/)
        {
            if ($bpcount > 5)
            {
                return;
            }

            @table = read_en_matrix_stream_offbone(*$FH, 4);

            @{$param_hashref->{'mismatch_interior'}{$BASEPAIRS[$bpcount]}}
            = @table;

            $bpcount++;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in mismatch_interior "
                             ."table: $_\n");
        }
    }
}

sub fetch_mismatch_hairpin_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;
    my $bpcount = 0;

    while (<$FH>)
    {
        if ($_ =~ /\s+-?\d+\s+-?\d+\s+-?\d+\s+-?\d+\s+-?\d+/)
        {
            if ($bpcount > 5)
            {
                return;
            }

            @table = read_en_matrix_stream_offbone(*$FH, 4);

            @{$param_hashref->{'mismatch_hairpin'}{$BASEPAIRS[$bpcount]}}
            = @table;

            $bpcount++;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in mismatch_hairpin "
                             ."table: $_\n");
        }
    }
}

sub fetch_dangle5_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;

    while (<$FH>)
    {
        # start reading the matrix after the comment line
        if ($_ =~ /\s*\/\*\s*\@\s*[AUGCaugc]\s*[AUGCaugc]\s*[AUGCaugc]
                          \s*[AUGCaugc]\s*\*\//x)
        {
            # just skip comment line
        }
        elsif ($_ =~/\s+(?:INF\s*){5}/x)
        {
            @table = read_en_matrix_stream_offbone(*$FH, 6);

            @{$param_hashref->{'dangle5'}} = @table;

            return;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in dangle5 table: "
                             ."$_\n");
        }
    }
}

sub fetch_dangle3_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @table;

    while (<$FH>)
    {
        # start reading the matrix after the comment line
        if ($_ =~ /\s*\/\*\s*\@\s*[AUGCaugc]\s*[AUGCaugc]\s*[AUGCaugc]
                          \s*[AUGCaugc]\s*\*\//x)
        {
            # just skip comment line
        }
        elsif ($_ =~/\s+(?:INF\s*){5}/x)
        {
            @table = read_en_matrix_stream_offbone(*$FH, 6);

            @{$param_hashref->{'dangle3'}} = @table;

            return;
        }
        elsif ($_ =~ /^\s*\n/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in dangle3 table: "
                             ."$_\n");
        }
    }
}

sub fetch_loop_energies($ \@)
{
    my ($FH, $param_aryref) = @_;
    my $pattern;
    my $plen;
    my $i;
    my $val;

    while (<$FH>)
    {
        # determine pattern size
        $pattern = '';
        $plen = 0;
        while ($_ =~ /${pattern}\s+(INF|-?\d+)/)
        {
            $pattern .= '\s+(INF|-?\d+)';
            $plen++;
        }
        if ($plen == 0)
        {
            seek($FH, - length($_), 1); # unread line
            return;
        }
        if ($_ =~ /${pattern}/x)
        {
            for ($i = 1; $i <= $plen; $i++)
            {
                $val = substr($_, $-[$i], $+[$i] - $-[$i]);

                if ($val == 'INF')
                {
                    $val = $UNDEF_ENERGY;
                }
                push @{$param_aryref}, $val;
            }
        }
        else
        {
            msg_error_and_die("Unregocnised line found in loop list: "
                             ."$_\n");
        }
    }
}

sub fetch_tetraloop_energies($ \%)
{
    my ($FH, $param_hashref) = @_;
    my @loop;

    while (<$FH>)
    {
        if ($_ =~ /\s*([aucgAUCG]{6})\s*(-?\d+\.?\d*)/)
        {
            @loop = split //, lc($1);

            ${$param_hashref->{tetraloops}}{$loop[0]}
                                           {$loop[1]}
                                           {$loop[2]} 
                                           {$loop[3]}
                                           {$loop[4]} 
                                           {$loop[5]} = $2;
        }
        elsif ($_ =~ /^\s*$/)
        {
            return;
        }
        else
        {
            msg_error_and_die("Unregocnised line found in tetra loop list: "
                              ."$_\n");
        }

    }

}

sub fetch_parameters($ \% \%)
{
    my ($filename, $param_hashref, $arg_hashref) = @_;

    # open file
    my $fh = open_or_die($filename);

    # search file for param tables
    while (<$fh>)
    {
        # read tables in separated functions
        if ($_ =~ /\#\s+stack_energies\s*$/)
        {
            msg_verbose("  Reading stack_energies table...");
            fetch_stacking_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+int11\_energies\s*$/)
        {
            msg_verbose("  Reading int11_energy table...");
            fetch_int11_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+int21\_energies\s*$/)
        {
            msg_verbose("  Reading int21_energy table...");
            fetch_int21_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+int22_energies/)
        {
            msg_verbose("  Reading int22_energy table...");
            fetch_int22_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+mismatch\_interior/)
        {
            msg_verbose("  Reading mismatch_interior table...");
            fetch_mismatch_interior_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+mismatch\_hairpin/)
        {
            msg_verbose("  Reading mismatch_hairpin table...");
            fetch_mismatch_hairpin_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+dangle5\s*$/)
        {
            msg_verbose("  Reading dangle5 table...");
            fetch_dangle5_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+dangle3\s*$/)
        {
            msg_verbose("  Reading dangle3 table...");
            fetch_dangle3_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+internal_loop\s*$/)
        {
            msg_verbose("  Reading internal loop table...");
            fetch_loop_energies($fh, @{$param_hashref->{'internals'}});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+hairpin\s*$/)
        {
            msg_verbose("  Reading hairpin table...");
            fetch_loop_energies($fh, @{$param_hashref->{'hairpins'}});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+bulge\s*$/)
        {
            msg_verbose("  Reading bulge table...");
            fetch_loop_energies($fh, @{$param_hashref->{'bulges'}});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s*Tetraloops\s*$/)
        {
            msg_verbose("  Reading tetra loop scores...");
            fetch_tetraloop_energies($fh, %{$param_hashref});
            msg_verbose("done\n");
        }
    }

    close($fh);
}

sub output_int11_energies(\%)
{
    my ($param_hashref) = @_;
    my ($bp1, $bp2);
    my ($i, $j);
    my $tmp = 1;
    my $cur;
    my $dig = 0;
    my $c1 = 0;
    my $c2 = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
                foreach (@ALPHA)
                {
                    $i = $_;
                    foreach (@ALPHA)
                    {
                        $j = $_;
                        $cur = ${$param_hashref->{$bp1}{$bp2}}[$i][$j];

                        if (($cur * $cur) > ($tmp * $tmp))
                        {
                            $tmp = $cur;
                        }
                    }
                }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        push @out, "   /* ".uc($bp1)." */\n";
        push @out, "   bp1 = this->bp_idx["
             .substr($bp1, 0, 1)."]["
             .substr($bp1, 1, 1)."];\n";
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            push @out, "   /*    ".uc($bp2)." */\n";
            push @out, "   bp2 = this->bp_idx["
                .substr($bp2, 0, 1)."]["
                .substr($bp2, 1, 1)."];\n";

            foreach (@ALPHA)
            {
                $i = $_;
                foreach (@ALPHA)
                {
                    $j = $_;
                    push @out, sprintf("   this->G_int11[bp1][bp2][".
                                       $BASES[$i]."][".$BASES[$j]."] = %*i - offset;",
                                $dig, ${$param_hashref->{$bp1}{$bp2}}[$i][$j]);
                    
                    push @out, " /* ";
#                    if ($c1 == 0)
#                    {
#                        push @out, uc($b);
#                    }
                    if ($c2 == 0)
                    {
                        push @out, uc($BASES[$i]);
                    }
                    else
                    {
                        push @out, " ";
                    }
                    
                    push @out, " ".uc($BASES[$j])." */\n";

                    $c1++;
                    if ($c1 == 16)
                    {
                        $c1 = 0;
                    }
                        $c2++;
                    if ($c2 == 4)
                    {
                        $c2 = 0;
                    }
                }
                #print "\n";                    
            }
        }
        push @out, "\n";
    }

    return @out;
}

sub output_int21_energies(\%)
{
    my ($param_hashref) = @_;
    my ($bp1, $bp2, $b);
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my $c1 = 0;
    my $c2 = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            foreach (@BASES)
            {
                $b = $_;
                foreach (@ALPHA)
                {
                    $i = $_;
                    foreach (@ALPHA)
                    {
                        $j = $_;
                        $cur = ${$param_hashref->{$bp1}{$bp2}{$b}}[$i][$j];

                        if (($cur * $cur) > ($tmp * $tmp))
                        {
                            $tmp = $cur;
                        }
                    }
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        push @out, "   /* ".uc($bp1)." */\n";
        push @out, "   bp1 = this->bp_idx["
             .substr($bp1, 0, 1)."]["
             .substr($bp1, 1, 1)."];\n";
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            push @out, "   /*   ".uc($bp2)." */\n";
            push @out, "   bp2 = this->bp_idx["
                .substr($bp2, 0, 1)."]["
                .substr($bp2, 1, 1)."];\n";

            foreach (@BASES)
            {
                $b = $_;
                foreach (@ALPHA)
                {
                    $i = $_;
                    foreach (@ALPHA)
                    {
                        $j = $_;
                        push @out, sprintf("   this->G_int21[bp1][bp2][$b][".
                               $BASES[$i]."][".$BASES[$j]."] = %*i - offset;",
                             $dig, ${$param_hashref->{$bp1}{$bp2}{$b}}[$i][$j]);

                        push @out, " /* ";
                        if ($c1 == 0)
                        {
                            push @out, uc($b);
                        }
                        else
                        {
                            push @out, " ";
                        }
                        push @out, " ";
                        if ($c2 == 0)
                        {
                            push @out, uc($BASES[$i]);
                        }
                        else
                        {
                            push @out, " ";
                        }

                        push @out, " ".uc($BASES[$j])." */\n";

                        $c1++;
                        if ($c1 == 16)
                        {
                            $c1 = 0;
                        }
                        $c2++;
                        if ($c2 == 4)
                        {
                            $c2 = 0;
                        }
                    }
                    #print "\n";                    
                }
            }
        }
        push @out, "\n";
    }

    return @out;
}

sub output_int22_energies(\%)
{
    my ($param_hashref) = @_;
    my ($bp1, $bp2, $b1, $b2);
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my $c1 = 0;
    my $c2 = 0;
    my $c3 = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            foreach (@BASES)
            {
                $b1 = $_;
                foreach (@BASES)
                {
                    $b2 = $_;
                    foreach (@ALPHA)
                    {
                        $i = $_;
                        foreach (@ALPHA)
                        {
                            $j = $_;
                       $cur = ${$param_hashref->{$bp1}{$bp2}{$b1}{$b2}}[$i][$j];

                           if (($cur * $cur) > ($tmp * $tmp))
                           {
                               $tmp = $cur;
                           }
                        }
                    }
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        push @out, "   /* ".uc($bp1)." */\n";
        push @out,  "   bp1 = this->bp_idx["
             .substr($bp1, 0, 1)."]["
             .substr($bp1, 1, 1)."];\n";
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            push @out, "   /*   ".uc($bp2)." */\n";
            push @out, "   bp2 = this->bp_idx["
                .substr($bp2, 0, 1)."]["
                .substr($bp2, 1, 1)."];\n";

            foreach (@BASES)
            {
                $b1 = $_;
                foreach (@BASES)
                {
                    $b2 = $_;
                    foreach (@ALPHA)
                    {
                        $i = $_;
                        foreach (@ALPHA)
                        {
                            $j = $_;
                  push @out, sprintf("   this->G_int22[bp1][bp2][$b1][$b2]["
                         .$BASES[$i]."][".$BASES[$j]."] = %*i - offset;",
                      $dig, ${$param_hashref->{$bp1}{$bp2}{$b1}{$b2}}[$i][$j]);

                            push @out, " /* ";
                            if ($c1 == 0)
                            {
                                push @out, uc($b1);
                            }
                            else
                            {
                                push @out, " ";
                            }
                            $c1++;
                            if ($c1 == 64)
                            {
                                $c1 = 0;
                            }
                            push @out, " ";
                            if ($c2 == 0)
                            {
                                push @out, uc($b2);
                            }
                            else
                            {
                                push @out, " ";
                            }
                            $c2++;
                            if ($c2 == 16)
                            {
                                $c2 = 0;
                            }
                            push @out, " ";
                            if ($c3 == 0)
                            {
                                push @out, uc($BASES[$i]);
                            }
                            else
                            {
                                push @out, " ";
                            }
                            $c3++;
                            if ($c3 == 4)
                            {
                                $c3 = 0;
                            }
                            push @out, " ".uc($BASES[$j])." */\n";
                        }                
                    }
                }
            }
        }
        push @out, "\n";
    }

    return @out;
}

sub output_mismatch_interior_energies(\%)
{
    my ($param_hashref) = @_;
    my $bp1;
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
                $cur = ${$param_hashref->{$bp1}}[$i][$j];
                
                if (($cur * $cur) > ($tmp * $tmp))
                {
                    $tmp = $cur;
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        push @out, "   /* ".uc($bp1)." */\n";
        push @out, "   bp1 = this->bp_idx["
            .substr($bp1, 0, 1)."]["
            .substr($bp1, 1, 1)."];\n";
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
#                if (${$param_hashref->{$bp1}}[$i][$j] != 0)
#                {
                    push @out, sprintf("   this->G_mismatch_interior[bp1][".
                                 $BASES[$i]."][".$BASES[$j]."] = %*i - offset;",
                                       $dig, ${$param_hashref->{$bp1}}[$i][$j]);
#                }
#                else
#                {
#                    push @out, sprintf("   this->G_mismatch_interior[bp1][".
#                                $BASES[$i]."][".$BASES[$j]."] = %*s - offset;",
#                                       $dig, "0.0f");
#                }

                push @out, " /* ".uc($BASES[$i]).uc($BASES[$j])." */\n";
            }
            push @out, "\n";
        }
    }

    return @out;
}

sub output_mismatch_hairpin_energies(\%)
{
    my ($param_hashref) = @_;
    my $bp1;
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
                $cur = ${$param_hashref->{$bp1}}[$i][$j];
                
                if (($cur * $cur) > ($tmp * $tmp))
                {
                    $tmp = $cur;
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        push @out, "   /* ".uc($bp1)." */\n";
        push @out, "   bp1 = this->bp_idx["
            .substr($bp1, 0, 1)."]["
            .substr($bp1, 1, 1)."];\n";
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
                push @out, sprintf("   this->G_mismatch_hairpin[bp1][".
                       $BASES[$i]."][".$BASES[$j]."] = %*i - offset;",
                       $dig, ${$param_hashref->{$bp1}}[$i][$j]);
                
                push @out, " /* ".uc($BASES[$i]).uc($BASES[$j])." */\n";
            }
            push @out, "\n";
        }
    }

    return @out;
}

sub output_mismatch_stack_energies(\% \%)
{
    my ($hairpin_hashref, $interior_hashref) = @_;
    my $bp1;
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
                $cur = int((${$hairpin_hashref->{$bp1}}[$i][$j]
                            + ${$interior_hashref->{$bp1}}[$i][$j])/2);
                
                if (($cur * $cur) > ($tmp * $tmp))
                {
                    $tmp = $cur;
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    push @out, "   /* stacks containing a mismatch */\n";
    push @out, "   /* mi: param from mismatch_interior table */\n";
    push @out, "   /* mh: param from mismatch_hairpin */\n\n";
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;

        push @out, "   bp1 = this->bp_idx["
            .substr($bp1, 0, 1)."]["
            .substr($bp1, 1, 1)."];\n";

        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;

                push @out, "   /* ".uc($bp1);
                push @out, " ".uc($BASES[$i].$BASES[$j])." */\n";
                push @out, "   /* mh: ".${$hairpin_hashref->{$bp1}}[$i][$j];
                push @out, " mi: ".${$interior_hashref->{$bp1}}[$i][$j];
                push @out, "*/\n";
                
                push @out, sprintf("   this->G_mm_stack[bp1][(int) "
                                   ."this->bp_idx["
                       .$BASES[$i]."][".$BASES[$j]."]] = %*i - offset;\n",
                       $dig, int((${$hairpin_hashref->{$bp1}}[$i][$j]
                                  + ${$interior_hashref->{$bp1}}[$i][$j])/2));
            }
            push @out, "\n";
        }
    }

    return @out;
}

sub output_stacking_energies(\%)
{
    my ($param_hashref) = @_;
    my ($bp1, $bp2);
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    foreach (@BASEPAIRS)
    {
        $bp1 = $_;
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;

            $cur = $param_hashref->{$bp1}{$bp2};

            if (($cur * $cur) > ($tmp * $tmp))
            {
                $tmp = $cur;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    foreach (@BASEPAIRS)
    {
        $bp1 = $_;

        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            push @out, "   /* ".uc($bp1)." ".uc($bp2)." */\n";           
            push @out, "   /* 5'- ".uc(substr($bp1, 0, 1).substr($bp2, 1, 1))."\n";
            push @out, "          ".uc(substr($bp1, 1, 1).substr($bp2, 0, 1))
                ." - 5' */\n";

            push @out, "   this->G_stack[(int) this->bp_idx[(int)"
                .substr($bp1, 0, 1)."][(int)"
                .substr($bp1, 1, 1)."]]\n";
            push @out, "                [(int) this->bp_idx[(int)"
                .substr($bp2, 0, 1)."][(int)"
                .substr($bp2, 1, 1)."]] = ";
            push @out, sprintf("%*i - offset;\n",
                               $dig,
                               $param_hashref->{$bp1}{$bp2});

            push @out, "\n";            
        }
    }
    return @out;
}

sub output_dangle5_energies(\@)
{
    my ($param_aryref) = @_;
    my ($bp, $b);
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    for ($i = 0; $i <= $#BASEPAIRS; $i++)
    {
        for ($j = 0; $j <= $#BASES; $j++)
        {
            $cur = ${$param_aryref}[$i][$j];

            if (($cur * $cur) > ($tmp * $tmp))
            {
                $tmp = $cur;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    for ($i = 0; $i <= $#BASEPAIRS; $i++)
    {
        $bp = $BASEPAIRS[$i];

        push @out, "   /* ".uc($bp)." */\n";

        for ($j = 0; $j <= $#BASES; $j++)
        {
            $b = $BASES[$j];

            push @out, "   this->G_dangle5[(int)this->bp_idx["
                .substr($bp, 0, 1)."]["
                .substr($bp, 1, 1)."]]["
                .$b."] = ";
            push @out, sprintf("%*i - offset;\n", $dig, ${$param_aryref}[$i][$j]);      
        }
        push @out, "\n";
    }
    return @out;
}

sub output_dangle3_energies(\@)
{
    my ($param_aryref) = @_;
    my ($bp, $b);
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;
    my @out;

    # fetch largest no.
    for ($i = 0; $i <= $#BASEPAIRS; $i++)
    {
        for ($j = 0; $j <= $#BASES; $j++)
        {
            $cur = ${$param_aryref}[$i][$j];

            if (($cur * $cur) > ($tmp * $tmp))
            {
                $tmp = $cur;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    for ($i = 0; $i <= $#BASEPAIRS; $i++)
    {
        $bp = $BASEPAIRS[$i];

        push @out, "   /* ".uc($bp)." */\n";

        for ($j = 0; $j <= $#BASES; $j++)
        {
            $b = $BASES[$j];

            push @out, "   this->G_dangle3[(int)this->bp_idx["
                .substr($bp, 0, 1)."]["
                .substr($bp, 1, 1)."]]["
                .$b."] = ";
            push @out, sprintf("%*i - offset;\n", $dig, ${$param_aryref}[$i][$j]);      
        }
        push @out, "\n";
    }
    return @out;
}

sub output_hairpins_energies(\@)
{
    my ($param_aryref) = @_;
    my $loop_len = 0;
    my $tmp = 0;
    my $dig = 0;
    my $dig_idx = 0;
    my @out;

    # fetch largest no.
    foreach (<@{$param_aryref}>)
    {
        if ($_ ne $UNDEF_ENERGY)
        {
            if (($_ * $_) > ($tmp * $tmp))
            {
                $tmp = $_;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    if ($#{$param_aryref} > 0)
    {
        $dig_idx = int(log10($#{$param_aryref}) + 1);
    }

    foreach (<@{$param_aryref}>)
    {
        push @out, "   this->G_hairpin_loop[".
            sprintf("%*i", $dig_idx, $loop_len)
            ."] = ";

        if ($_ eq $UNDEF_ENERGY)
        {
           push @out, $UNDEF_ENERGY;
        }
        else
        {
            push @out, sprintf("%*i - offset", $dig, $_);
        }

        push @out, ";\n";
        $loop_len++;
    }
    return @out;
}

sub output_internals_energies(\@)
{
    my ($param_aryref) = @_;
    my $loop_len = 0;
    my $tmp = 0;
    my $dig = 0;
    my $dig_idx = 0;
    my @out;

    # fetch largest no.
    foreach (<@{$param_aryref}>)
    {
        if ($_ ne $UNDEF_ENERGY)
        {
            if (($_ * $_) > ($tmp * $tmp))
            {
                $tmp = $_;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    if ($#{$param_aryref} > 0)
    {
        $dig_idx = int(log10($#{$param_aryref}) + 1);
    }

    foreach (<@{$param_aryref}>)
    {
        push @out, "   this->G_internal_loop[".
            sprintf("%*i", $dig_idx, $loop_len)
            ."] = ";

        if ($_ eq $UNDEF_ENERGY)
        {
           push @out, $UNDEF_ENERGY;
        }
        else
        {
            push @out, sprintf("%*i - offset", $dig, $_);
        }

        push @out, ";\n";
        $loop_len++;
    }
    return @out;
}

sub output_bulges_energies(\@)
{
    my ($param_aryref) = @_;
    my $loop_len = 0;
    my $tmp = 0;
    my $dig = 0;
    my $dig_idx = 0;
    my @out;

    # fetch largest no.
    foreach (<@{$param_aryref}>)
    {
        if ($_ ne $UNDEF_ENERGY)
        {
            if (($_ * $_) > ($tmp * $tmp))
            {
                $tmp = $_;
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    if ($#{$param_aryref} > 0)
    {
        $dig_idx = int(log10($#{$param_aryref}) + 1);
    }

    foreach (<@{$param_aryref}>)
    {
        push @out, "   this->G_bulge_loop[".
            sprintf("%*i", $dig_idx, $loop_len)
            ."] = ";

        if ($_ eq $UNDEF_ENERGY)
        {
           push @out, $UNDEF_ENERGY;
        }
        else
        {
            push @out, sprintf("%*i - offset", $dig, $_);
        }

        push @out, ";\n";
        $loop_len++;
    }
    return @out;
}

sub output_tetraloop_energies(\%)
{
    my ($tl_hashref) = @_;
    my @out;
    my ($b1, $b2, $b3, $b4, $b5, $b6);
    my $tmp = 1;
    my $dig = 0;
    my $dig_idx = 1;
    my $no = 0;

    # determine max. no. of digits
    foreach (keys %{$tl_hashref})
    {
        $b1 = $_;
        foreach (keys %{$tl_hashref->{$b1}})
        {
            $b2 = $_;
            foreach (keys %{$tl_hashref->{$b1}{$b2}})
            {
                $b3 = $_;
                foreach (keys %{$tl_hashref->{$b1}{$b2}{$b3}})
                {
                    $b4 = $_;
                    foreach (keys %{$tl_hashref->{$b1}{$b2}{$b3}{$b4}})
                    {
                        $b5 = $_;
                        foreach(keys %{$tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}})
                        {
                            $b6 = $_;

                            if (  ($tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}{$b6}
                                 * $tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}{$b6})
                                  > ($tmp * $tmp))
                            {
                                $tmp =
                                   $tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}{$b6};
                            }
                            $no++;
                        }
                    }
                }
            }
        }
    }
    if ($tmp < 0)
    {
        $dig++;
        $tmp = $tmp * (-1);
    }
    $dig += int(log10($tmp) + 1);

    if ($no > 0)
    {
        $dig_idx = int(log10($no) + 1);
    }

    # create output
    $no = 0;
    foreach (sort keys %{$tl_hashref})
    {
        $b1 = $_;
        foreach (sort keys %{$tl_hashref->{$b1}})
        {
            $b2 = $_;
            foreach (sort keys %{$tl_hashref->{$b1}{$b2}})
            {
                $b3 = $_;
                foreach (sort keys %{$tl_hashref->{$b1}{$b2}{$b3}})
                {
                    $b4 = $_;
                    foreach (sort keys %{$tl_hashref->{$b1}{$b2}{$b3}{$b4}})
                    {
                        $b5 = $_;
                        foreach(sort keys %{$tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}})
                        {
                            $b6 = $_;

                            push @out, "   /* ".uc("$b1$b2$b3$b4$b5$b6")
                                .sprintf(" %*i */\n", $dig,
                                  $tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}{$b6});
                            push @out, sprintf("   l = this->tetra_loop[%*i];"
                                               ."\n", $dig_idx, $no);
                            push @out, "   l[0] = $b1; l[1] = $b2; "
                                      ."l[2] = $b3; l[3] = $b4; l[4] = $b5; "
                                      ."l[5] = $b6; l[6] = \'\\0\';\n";
                            push @out, sprintf("   this->G_tetra_loop[%*i] = ",
                                               $dig_idx, $no)
                                .sprintf("%*i - offset;\n\n", $dig,
                                  $tl_hashref->{$b1}{$b2}{$b3}{$b4}{$b5}{$b6});

                            $no++;
                        }
                    }
                }
            }
        }
    }

    return @out;
}
# FUNCTIONS        - END

# insert params at certain labels (for substitution take opening and closing)
# add assertion action
# read/ output more params

# MAIN             - BEGIN
my %arg_hash;
my $pod_verbose = 1;
my %param_hash;
my $ret_val;
my $filehandle = *STDOUT;
my $i;
my @file;
my @output;

# parse commandline
$ret_val = parseargs(%arg_hash);
if ($ret_val > 1)
{
    if ($ret_val == 3) { $pod_verbose = 2 }
    pod2usage(-exitval => 0, -verbose => $pod_verbose); 
}

# fetch parameters
fetch_parameters($arg_hash{paramfile},
                 %param_hash,
                 %arg_hash);

if (defined ($arg_hash{sourcefile}))
{
    $filehandle = open_or_die("$arg_hash{sourcefile}");
    @file = <$filehandle>;
    close($filehandle);
}

# write parameters
if (! defined ($arg_hash{stacking_energies}))
{
    @output = output_stacking_energies(
        %{$param_hash{'stacking_energies'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_STACKING_ENERGIES", "END_STACKING_ENERGIES",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (!defined ($arg_hash{int11_energies}))
{
    @output = output_int11_energies(%{$param_hash{'int11_energies'}});

    if (defined ($arg_hash{sourcefile}))
    {
        slice_tags_or_die("BEGIN_INT11_ENERGIES", "END_INT11_ENERGIES",
                          @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{int21_energies}))
{
    @output = output_int21_energies(%{$param_hash{'int21_energies'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_INT21_ENERGIES", "END_INT21_ENERGIES",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{int22_energies}))
{
    @output = output_int22_energies(%{$param_hash{'int22_energies'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_INT22_ENERGIES", "END_INT22_ENERGIES",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (!defined ($arg_hash{mismatch_interior}))
{
    @output = output_mismatch_interior_energies(
        %{$param_hash{'mismatch_interior'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_MISMATCH_INTERIOR", "END_MISMATCH_INTERIOR",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{mismatch_hairpin}))
{
    @output = output_mismatch_hairpin_energies(
        %{$param_hash{'mismatch_hairpin'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_MISMATCH_HAIRPIN", "END_MISMATCH_HAIRPIN",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{mismatch_stack}))
{
    @output = output_mismatch_stack_energies(
       %{$param_hash{'mismatch_hairpin'}}, %{$param_hash{'mismatch_interior'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_MISMATCH_STACK", "END_MISMATCH_STACK",
                         @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{tetraloops}))
{
    @output = output_tetraloop_energies(%{$param_hash{'tetraloops'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_TETRA_LOOPS", "END_TETRA_LOOPS",
                          @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{dangle5}))
{
    @output = output_dangle5_energies(@{$param_hash{'dangle5'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_DANGLE_5", "END_DANGLE_5", @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{dangle3}))
{
    @output = output_dangle3_energies(@{$param_hash{'dangle3'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_DANGLE_3", "END_DANGLE_3", @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{internals}))
{
    @output = output_internals_energies(@{$param_hash{'internals'}});

    if (defined ($arg_hash{sourcefile}))
    {
        slice_tags_or_die("BEGIN_INTERNALS", "END_INTERNALS", @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{hairpins}))
{
    @output = output_hairpins_energies(@{$param_hash{'hairpins'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_HAIRPINS", "END_HAIRPINS", @file, @output);
    }
    else
    {
        print @output;
    }
}
if (! defined($arg_hash{bulges}))
{
    @output = output_bulges_energies(@{$param_hash{'bulges'}});

    if (defined ($arg_hash{sourcefile}))
    {
        # seek opening and closing tag
        slice_tags_or_die("BEGIN_BULGES", "END_BULGES", @file, @output);
    }
    else
    {
        print @output;
    }
}

if (defined ($arg_hash{sourcefile}))
{
    $filehandle = open_or_die(">$arg_hash{sourcefile}");
    print $filehandle @file;
    close($filehandle);
}
# MAIN             - END


__END__

=head1 NAME

B<read_viennaparams> - Read a parameter file of the Vienna package and output/
store the parameters.

=head1 SYNOPSIS

B<read_viennaparams> [options] <param file>

=head1 DESCRIPTION

B<read_viennaparams> reads a parameter file of the Vienna RNA package in the
first place. The parameters are then printed on C<STDOUT> or integrated
into a C source file. 

For applying the parameters to a source file, we use certain tags within the
code. This gives you the chance to develop your code without paying attention
to new parameter releases. When a new set arrives, B<read_viennaparams> is able
to replace parameters inbetween opening and closing tags automatically. The
tags have to be embraced by an ANSI C comment: C<\* TAG *\>.

Here is the list of tags:

=over

=item C<E<lt>BEGIN|ENDE<gt>_STACKING_ENERGIES> - Energies for stacked pairs.

=item C<E<lt>BEGIN|ENDE<gt>_INT11_ENERGIES> - Energies for internal 1x1 loops.

=item C<E<lt>BEGIN|ENDE<gt>_INT21_ENERGIES> - Energies for internal 2x1 loops.

=item C<E<lt>BEGIN|ENDE<gt>_INT22_ENERGIES> - Energies for internal 2x2 loops.

=item C<E<lt>BEGIN|ENDE<gt>_MISMATCH_INTERIOR - Closing bp en. for int.loops.

=item C<E<lt>BEGIN|ENDE<gt>_MISMATCH_HAIRPIN - Closing bp en. for hairpin loops.

=item C<E<lt>BEGIN|ENDE<gt>_MISMATCH_STACK - Closing bp en. for mismatch stacks.

=item C<E<lt>BEGIN|ENDE<gt>_INTERNALS - Internal loop parameters.

=item C<E<lt>BEGIN|ENDE<gt>_HAIRPINS - Hairpin loop parameters.

=item C<E<lt>BEGIN|ENDE<gt>_BULGES - Bulge loop parameters.

=item C<E<lt>BEGIN|ENDE<gt>_DANGLE_5 - 5' dangling end energies.

=item C<E<lt>BEGIN|ENDE<gt>_DANGLE_3 - 3' dangling end energies.

=item C<E<lt>BEGIN|ENDE<gt>_TETRA_LOOPS - Tetra loop energies.

=back

Parameters are only stored in file if option B<-sourcefile> is given.

=head1 OPTIONS

=over 8

=item B<-stacking_energies>

Do not output/ edit energies for stacked base pairs.

=item B<-int11_energies>

Do not output/ edit energies for internal 1x1 loops.

=item B<-int21_energies>

Do not output/ edit energies for internal 2x1 loops.

=item B<-int22_energies>

Do not output/ edit energies for internal 2x2 loops.

=item B<-mismatch_interior>

Do not output/ edit energies for the opening/ closing base pair of internal
loops with
size E<gt> 2.

=item B<-mismatch_hairpin>

Do not output/ edit energies for the opening/ closing base pair of hairpin
loops.

=item B<-mismatch_stack>

Do not output/ edit energies for the artificial mismatch stack energies.

=item B<-internals>

Do not output/ edit energies for internal loops.

=item B<-hairpins>

Do not output/ edit energies for hairpin loops.

=item B<-bulges>

Do not output/ edit energies for bulge loops.

=item B<-tetraloops>

Do not output/ edit energies for tetra loops.

=item B<-dangle5>

Do not output/ edit energies for the 5' dangling end.

=item B<-dangle3>

Do not output/ edit energies for the 3' dangling end.

=item B<-sourcefile C<file>>

Path to the C source file to store parameters in.

=item B<-verbose>

Be verbose.

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
