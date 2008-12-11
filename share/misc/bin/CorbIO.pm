# Last modified: 2008-12-10.17
#
#
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

=pod

=head1 NAME

CorbIO - Simple input and output routines.

=head1 SYNOPSIS

    use CorbIO qw(:verbose);

    enable_verbose();

    msg_verbose("Only seen in verbose mode");

    if (is_verbose())
    {
      print("Also only seen in verbose mode");      
    }

    example for open or die, warn and error
    use exporter in this example

=head1 DESCRIPTION

This module provides simple message and file handling functions. The idea is to
provide functionality which is needed everywhere, like error messages, warnings
and file opening, as one-line-function. Of course this means, that on problems
the calling script is stopped (via C<die()>). As a rule, functions which might
exit a caller, tell so in their names.

To avoid the pollution of a scripts namespace nothing is exported by default.
Therefore each function to be used has to be imported individually or via tags.
Following tags are defined:

=over 4

=item * verbose - all functions related to the verbose mode
    (L<C<enable_verbose()>|"enable_verbose">,
    L<C<disable_verbose()>|"disable_verbose">,
    L<C<is_verbose()>|"is_verbose">,
    L<C<msg_verbose()>|"msg_verbose">)

=item * msg - Additional messaging functions (L<C<msg_warning()>|"msg_warning">)

=item * all - all functions from the tags above

=back

Since direct calling of functions via package naming is also allowed, the
function names are not prefixed with the package name.

=cut


package CorbIO;

use strict;
use warnings;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    $VERSION     = 0.01;
    
    @ISA         = qw(Exporter);
    @EXPORT      = ();
    %EXPORT_TAGS = (verbose => [qw(enable_verbose
                                   disable_verbose
                                   is_verbose
                                   msg_verbose)],
                    msg => [qw(msg_warning)]
                   );

    # add all the other tags to the ":all" tag,
    # deleting duplicates
    {
        my %seen;
        
        push @{$EXPORT_TAGS{all}},
        grep {!$seen{$_}++} @{$EXPORT_TAGS{$_}} foreach keys %EXPORT_TAGS;
    }

    # automatically add all tagged functions to the EXPORT_OK list
    Exporter::export_ok_tags('all');
}
our @EXPORT_OK;


# EXPORTED GLOBALS     - BEGIN
# EXPORTED GLOBALS     - END

# NON-EXPORTED GLOBALS - BEGIN
# NON-EXPORTED GLOBALS - END

# PRIVATE GLOBALS      - BEGIN
my $private_func_msg_verbose  = sub {};
my $private_func_is_verbose  = sub { return 0; };
# PRIVATE GLOBALS      - END

=pod

=head1 FUNCTIONS

=head2 enable_verbose

This function enables the verbose messaging mode. As default this is disabled,
meaning that calls to L<C<msg_verbose()>|"msg_verbose">, produce no output at
all. Additionally, L<C<is_verbose()>|"is_verbose"> returns 0 if disabled, 1
otherwise. The idea is to enable the verbose mode if a script is called to be
verbose by the user.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

    my $verbose
    my $optcatchresult = GetOptions('verbose!' => \$verbose);

    if ($optcatchresult == 0) { return 0 }

    enable_verbose() if $verbose;

=back

=cut

sub enable_verbose
{
    $private_func_msg_verbose = sub {
        my ($msg) = @_;

        print($msg);
    };

    $private_func_is_verbose = sub { return 1; };
}

=head2 disable_verbose

This function disables the verbose mode. After calling, calls to
L<C<msg_verbose()>|"msg_verbose"> produce no output and
L<C<is_verbose()>|"is_verbose"> returns 0.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

    disable_verbose();

=back

=cut

sub disable_verbose
{
    $private_func_msg_verbose = sub {};

    $private_func_is_verbose = sub { return 0; };
}

=pod

=head2 is_verbose

Returns 1 if verbose mode is enabled (former call to
L<C<enable_verbose()>|"enable_verbose">), 0 otherwise.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

    if (is_verbose())
    {
      PBar::enable;
    }

=back

=cut

sub is_verbose
{
    return &$private_func_is_verbose;
}

=head2 msg_verbose

Writes a message to C<STDOUT> depending on the value of
L<C<is_verbose()>|"is_verbose">. Utilises the standard C<print> function and
hence output might be buffered.

=over 4

=item ARGUMENTS

=over 4

=item message

The message to be written if in verbose mode. Has to be provided as scalar.

=back

=item EXAMPLE

     my $formatted_float = sprintf("%3.3f", 1/7);

     msg_verbose("1/7 = ${formatted_float}\n");

=back

=cut

sub msg_verbose
{
    return &$private_func_msg_verbose (@_);
}

=head2 msg_warning

Just writes a warning message to C<STDERR>. There is nothing special about this
but the message is preceeded by "WARNING:". This goes into the direction of the
GNU coding standards and should give as a unified look for all messages of our
scripts.

=over 4

=item ARGUMENTS

=over 4

=item message

The message to be written. Has to be provided as scalar.

=back

=item EXAMPLE

     msg_warning("This is a warning.");

=back

=cut

sub msg_warning($)
{
    my ($msg) = @_;

    print STDERR "WARNING:".$msg;
}

=pod

=head1 AUTHOR

Stefan Bienert (bienert@zbh.uni-hamburg.de)

=head1 COPYRIGHT

Copyright 2008 (C) Stefan Bienert

This module is part of CoRB.

CoRB is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CoRB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoRB.  If not, see <http://www.gnu.org/licenses/>.

=cut

1;

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:

