# Last modified: 2008-12-09.13
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
                                   msg_verbose)]
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

The message to be written if in verbose mode.

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

=head2 finish / pbar_finish

This function finalises the progress bar. Any former printing of the bar to
STDOUT is overwritten by a new bar with 100% values. Should be called after a
loop. Calling without issuing L<C<start()>|"start / pbar_start"> before,
will lead to warnings/ errors. This is just because we are accessing some
values which have to be initialised, first. Calling without the use of the 
L<C<update()>|"update / pbar_update"> function in the loop does not make much
sense. In this case you will just see a jump in the progress bar from 0%
(L<C<start()>|"start / pbar_start">) to 100%
(L<C<finish()>|"finish / pbar_finish">).  Only produces some output after
L<C<enable()>|"enable / pbar_enable"> was issued.

=over 4

=item ARGUMENTS

=over 4

=item progressbar

A hash to store all data of a particular progress bar. Is modified by this
function. Has to be initialised by L<C<start()>|"start / pbar_start"> before.
Has to be provided as reference while using
L<C<pbar_finish()>|"finish / pbar_finish">.

=back

=item EXAMPLE

     my %pbar;
     my $count = 100;

     PBar::start(%pbar, $count);

     for (my $i = 0; $i < $count; $i++)
     {
        print $i;
        PBar::update(%pbar, $i);
     }

     PBar::finish(%pbar);

=back

=cut

sub finish(\%)
{
#    return &$private_func_finish (@_);
}

sub pbar_finish
{
#    return &$private_func_finish (@_);
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

