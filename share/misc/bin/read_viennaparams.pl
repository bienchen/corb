#!@PERL@ -w
# -*- perl -*-
# @configure_input@
# Last modified: 2008-12-11.09


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


# SETTINGS         - BEGIN
BEGIN
{
  # try to get bash compatible shell
  $ENV{'SHELL'} = '@SHELL@' if exists $ENV{'DJGPP'};
  # $|=1; # dissable line buffer for print statements
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
# CONSTANTS        - END

# GLOBALS          - BEGIN
# GLOBALS          - END


# FUNCTIONS        - BEGIN
# write error
#   error_msg(message)
sub error_msg_and_die($)
{
    my ($message) = @_;

    die("$0:ERROR:".$message);
}

# open a file or die tryin'
#   open_or_die(filename, filehandle)
sub open_or_die($ \*)
{
    my ($filename, $filehandle) = @_;

    open($filehandle, $filename) or
        error_msg_and_die("Unable to open file \"$filename\": $!\n");
}

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
                               error_msg_and_die("@_");
                           };

    # set defaults

    # parse @ARGV
    $optcatchresult = GetOptions(
        'int21_energies!'    => \$argument_hashref->{int21_energies},
        'int22_energies!'    => \$argument_hashref->{int22_energies},
        'mismatch_interior!' => \$argument_hashref->{mismatch_interior},
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
        error_msg_and_die("Parameter file missing.\n"
                         ."Try \"-help\" or \"-man\" for more information.\n"); 
    }

    # fetch file name
    if ( ! -r $ARGV[0])
    {
        error_msg_and_die("Param file does not exist or is not readable: "
                         ."$ARGV[0]\n");
    }    
    $argument_hashref->{paramfile} = $ARGV[0];

    return 1;
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
                error_msg_and_die ("Unrecognised line found: $_");
            }
        }
        if ($i == 4)
        {
            return @table;
        }
        $i++;
    }
}

sub read_en_matrix_stream_offbone(\*)
{
    my ($FH) = @_;
    my $i = 1;
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
                error_msg_and_die ("Unrecognised line found: $_");
            }
        }
        if ($i == 4)
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
            error_msg_and_die ("Unrecognised line found: $_");
        }

        if ($i == 3)
        {
            return @table;
        }
        $i++;
    }
}

sub fetch_int21_energies(\* \%)
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
            error_msg_and_die("Unregocnised line found in int21_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_int22_energies(\* \%)
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
            error_msg_and_die("Unregocnised line found in int22_energies "
                             ."table: $_\n");
        }
    }
}

sub fetch_mismatch_interior_energies(\* \%)
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

            @table = read_en_matrix_stream_offbone(*$FH);

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
            error_msg_and_die("Unregocnised line found in mismatch_interior "
                             ."table: $_\n");
        }
    }
}

sub fetch_parameters($ \%)
{
    my ($filename, $param_hashref) = @_;

    # open file
    open_or_die($filename, *FH);

    # search file for param tables
    while (<FH>)
    {
        # read tables in separated functions
        if ($_ =~ /\#\s+int21\_energies\s*$/)
        {
            msg_verbose("  Reading int21_energy table...");
            fetch_int21_energies(*FH, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+int22_energies/)
        {
            msg_verbose("  Reading int22_energy table...");
            fetch_int22_energies(*FH, %{$param_hashref});
            msg_verbose("done\n");
        }
        if ($_ =~ /\#\s+mismatch\_interior/)
        {
            msg_verbose("  Reading mismatch_interior table...");
            fetch_mismatch_interior_energies(*FH, %{$param_hashref});
            msg_verbose("done\n");
        }
    }

    close(FH);
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
        print "   /* ".uc($bp1)." */\n";
        print "   bp1 = this->bp_idx["
             .substr($bp1, 0, 1)."]["
             .substr($bp1, 1, 1)."];\n";
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            print "   /*   ".uc($bp2)." */\n";
            print "   bp2 = this->bp_idx["
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
                        printf("   this->G_int21[bp1][bp2][$b][".
                               $BASES[$i]."][".$BASES[$j]."] = %*i;",
                             $dig, ${$param_hashref->{$bp1}{$bp2}{$b}}[$i][$j]);

                        print " /* ";
                        if ($c1 == 0)
                        {
                            print uc($b);
                        }
                        else
                        {
                            print " ";
                        }
                        print " ";
                        if ($c2 == 0)
                        {
                            print uc($BASES[$i]);
                        }
                        else
                        {
                            print " ";
                        }

                        print " ".uc($BASES[$j])." */\n";

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
        print "\n";
    }
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
        print "   /* ".uc($bp1)." */\n";
        print "   bp1 = this->bp_idx["
             .substr($bp1, 0, 1)."]["
             .substr($bp1, 1, 1)."];\n";
        foreach (@BASEPAIRS)
        {
            $bp2 = $_;
            print "   /*   ".uc($bp2)." */\n";
            print "   bp2 = this->bp_idx["
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
                  printf("   this->G_int22[bp1][bp2][$b1][$b2]["
                         .$BASES[$i]."][".$BASES[$j]."] = %*i;",
                      $dig, ${$param_hashref->{$bp1}{$bp2}{$b1}{$b2}}[$i][$j]);

                            print " /* ";
                            if ($c1 == 0)
                            {
                                print uc($b1);
                            }
                            else
                            {
                                print " ";
                            }
                            $c1++;
                            if ($c1 == 64)
                            {
                                $c1 = 0;
                            }
                            print " ";
                            if ($c2 == 0)
                            {
                                print uc($b2);
                            }
                            else
                            {
                                print " ";
                            }
                            $c2++;
                            if ($c2 == 16)
                            {
                                $c2 = 0;
                            }
                            print " ";
                            if ($c3 == 0)
                            {
                                print uc($BASES[$i]);
                            }
                            else
                            {
                                print " ";
                            }
                            $c3++;
                            if ($c3 == 4)
                            {
                                $c3 = 0;
                            }
                            print " ".uc($BASES[$j])." */\n";
                        }                
                    }
                }
            }
        }
        print "\n";
    }
}

sub output_mismatch_interior_energies(\%)
{
    my ($param_hashref) = @_;
    my $bp1;
    my ($i, $j);
    my $tmp = 0;
    my $cur;
    my $dig = 0;

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
        print "   /* ".uc($bp1)." */\n";
        print "   bp1 = this->bp_idx["
            .substr($bp1, 0, 1)."]["
            .substr($bp1, 1, 1)."];\n";
        foreach (@ALPHA)
        {
            $i = $_;
            foreach (@ALPHA)
            {
                $j = $_;
                printf("   this->G_mismatch_interior[bp1][".
                       $BASES[$i]."][".$BASES[$j]."] = %*i;",
                       $dig, ${$param_hashref->{$bp1}}[$i][$j]);
                
                print " /* ".uc($BASES[$i]).uc($BASES[$j])." */\n";
            }
            print "\n";
        }
    }
}

# FUNCTIONS        - END


# MAIN             - BEGIN
my %arg_hash;
my $pod_verbose = 1;
my %param_hash;
my $ret_val;

# parse commandline
$ret_val = parseargs(%arg_hash);
if ($ret_val > 1)
{
    if ($ret_val == 3) { $pod_verbose = 2 }
    pod2usage(-exitval => 0, -verbose => $pod_verbose); 
}

# fetch parameters
fetch_parameters($arg_hash{paramfile}, %param_hash);

# write parameters
    if ($arg_hash{int21_energies})
    {
        output_int21_energies(%{$param_hash{'int21_energies'}});
    }
    if ($arg_hash{int22_energies})
    {
        output_int22_energies(%{$param_hash{'int22_energies'}});
    }
    if ($arg_hash{mismatch_interior})
    {
        output_mismatch_interior_energies(%{$param_hash{'mismatch_interior'}});
    }
# MAIN             - END


__END__

=head1 NAME

read_viennaparams - ...

=head1 SYNOPSIS

B<read_viennaparams> [options] <param file>

=head1 DESCRIPTION

B<reformat> ...  B<corb> project

=head1 OPTIONS

=over 8

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
