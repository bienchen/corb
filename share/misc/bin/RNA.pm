# Last modified: 2009-11-06.15
#
#
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

=pod

=head1 NAME

RNA - Everything RNA. Alphabets, base pairs, test sequences...

=head1 SYNOPSIS

    use RNA qw(:all);

=head1 DESCRIPTION

This module should provide everything you need to work with RNA in Perl. This
includes the standard alphabet, Watson-Crick base pairs, base to number codings
and many more. For testing we also include a set of sequence and structure
tupels, here.
 
To avoid the pollution of a scripts namespace nothing is exported by default.
Therefore everything has to be imported individually or via tags.
Following tags are defined:

=over 4

=item * sigma - everything related to the RNA alphabet
    (L<C<$A>|"$A, $C, $G, $U">,
    L<C<$G>|"$A, $C, $G, $U">,
    L<C<$C>|"$A, $C, $G, $U">,
    L<C<$U>|"$A, $C, $G, $U">,
    L<C<@ALPHA>|"@ALPHA">,
    L<C<@BASEPAIRS>|"@BASEPAIRS">,
    L<C<$N_BPAIRS>|"$N_BPAIRS">,
    L<C<@BASES>|"@BASES">,
    L<C<$N_BASES>|"$N_BASES">)

=item * ViennaRNA - all the wrappers for functions from the Vienna RNA Package
    (L<C<RNAfold()>|"RNAfold">)
    (L<C<%TESTS>|"%TESTS">)

=item * all - everything above

=back

Since direct calling of functions via package naming is also allowed, the
function names are not prefixed with the package name.
(does this hold for vars?)

=cut


package RNA;

use strict;
use warnings;

# no idea if this is the recommended way to do it, but since Perl < 5.10 does
# not come with Module::Load::Conditional, it seems to be the simplest
our $DB_ON = 0;
eval { require DBI };
if ( !$@ )
{
    eval { require DBD::SQLite };
    if ( !$@ ) 
    {
        $DB_ON = 1;
    }
}

BEGIN {
    use Exporter();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    $VERSION     = 0.01;

    @ISA         = qw(Exporter);
    @EXPORT      = ();
    %EXPORT_TAGS = (sigma => [qw($A $G $C $U
                                 @ALPHA
                                 @BASEPAIRS $N_BPAIRS
                                 @BASES $N_BASES)],
                    ViennaRNA => [qw(&RNAfold)]
                   );

    # add all the other tags to the ":all" tag, deleting duplicates
    {
        my %seen;
        
        push @{$EXPORT_TAGS{all}},
        grep {!$seen{$_}++} @{$EXPORT_TAGS{$_}} foreach keys %EXPORT_TAGS;
    }

    # automatically add all tagged functions to the EXPORT_OK list
    Exporter::export_ok_tags('all');
}
#our @EXPORT_OK; # do we need this while it is already defined above?


# EXPORTED GLOBALS     - BEGIN <- automatically exported
# EXPORTED GLOBALS     - END

# NON-EXPORTED GLOBALS - BEGIN <- not automatically but still exportABLE!

=pod

=head1 GLOBAL VARIABLES

=head2 $A, $C, $G, $U

These guys hold the translation of bases into numbers.

=cut

our ($A, $C, $G, $U);
*A = \0;
*C = \1;
*G = \2;
*U = \3;

=pod

=head2 @ALPHA

The RNA alphabet, translated to numbers. Is defined using the numbers of
L<C<$A>|"$A, $C, $G, $U">, L<C<$G>|"$A, $C, $G, $U">, L<C<$C>|"$A, $C, $G, $U">
and L<C<$U>|"$A, $C, $G, $U">.

=cut

our @ALPHA = ($A, $C, $G, $U);

=head2 @BASEPAIRS

The list of canonical base pairs including whobble C<GU>. Each entry is a string
of 2 characters. All pairs are written in lower case.

=cut

our @BASEPAIRS = ('cg', 'gc', 'gu', 'ug', 'au', 'ua');

=head2 $N_BPAIRS

Size of list L<C<@BASEPAIRS>|"@BASEPAIRS">. That is the number of canonical
base pairs including whobble C<GU>.

=cut

#my $N_BPAIRS;

our $N_BPAIRS;
*N_BPAIRS = \($#BASEPAIRS + 1);

=head2 @BASES

List of RNA standard nucleotides. All characters are written in lower case.

=cut

our @BASES = ('a', 'c', 'g', 'u');

=head2 $N_BASES

Size of list L<C<@BASES>|"@BASES">. That is the number of standard RNA bases.

=cut

our $N_BASES;
*N_BASES = \($#BASES + 1);

# NON-EXPORTED GLOBALS - END

# PRIVATE GLOBALS      - BEGIN <- no access from outside
# PRIVATE GLOBALS      - END

=pod

=head1 FUNCTIONS

=head2 RNAfold

Folds an RNA sequence into a 2D structure using RNAfold of the Vienna RNA
Package. The return value is a hash with C<structure> and C<mfe> as keys,
pointing to the structure and the "minimum free energy" of it.

Since folding can take some time, we provide a databse mechanism, which is able
to store all your results for reuse. The idea is to store each call to RNAfold
with its parameters and on recurring calls just send back the answer from the
database instead of calculating it again. An entry in the database contains the
sequence, a unified parameterstring, ...path... and the results of the call. 

The database will be stored in a single file.... The database system we use is
SQLite, therefore it should be shareable among different operating systems.

Recording may be turned of by calling ... but cannot be re-enabled.

=over 4

=item ARGUMENTS

=over 4

=item sequence

The sequence to be folded. Has to be a scalar.

=item parameters

Parameters to be passed to RNAfold. All options have to be passed in a single
string in the same way they would be used in a command line call to RNAfold.
May be omitted if not used.

=item command

As default, this functions just sends "RNAfold" as command to a shell. This
option can be used to change the call to whatever you want. Just omit it to get
the default behaviour.

=back

=item EXAMPLE

    my $verbose
    my $optcatchresult = GetOptions('verbose!' => \$verbose);

    if ($optcatchresult == 0) { return 0 }

    enable_verbose() if $verbose;

=back

=cut

sub RNAfold
{
    my ($seq, $params, $cmd) = @_;
    my $len = length($seq);
    my $i;
    my @err_run;
    my %result;

    if (!defined($cmd))
    {
        $cmd = 'RNAfold';
    }

    if (!defined($params))
    {
        $params = '';
    }

    # start folding
    unless(open(FH, "echo \"$seq\" | $cmd $params 2>&1 |"))
    {
        msg_error_and_die ("Could not start $cmd $params: $!\n");
    }

    @err_run = ();
    
    foreach (<FH>)
    {
        if ($_ =~ /([\(\.\)]{$len})\s+\(\s*(\-?\d+\.\d+)\)/)
        {
            $result{structure} = $1;
            $result{mfe} = $2;
        }
        
        push(@err_run, $_);
    }

    close(FH);

    if (! defined ($result{structure}))
    {
        msg_error_and_die("Running RNAfold failed, output of\n"
                          ."\`echo \"$seq\" | $cmd $params\"\`: ",
                          @err_run);
    }
        

    return %result;
}

1;

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:

