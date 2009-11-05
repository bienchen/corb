# Last modified: 2009-11-05.11
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

=item * tests - everything related to testing
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
                    tests => [qw(%TESTS)]
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


1;

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:

