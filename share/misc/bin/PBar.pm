# Last modified: 2008-11-10.10
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


=head1 NAME

PBar - A very simple progress bar.

=head1 SYNOPSIS

    use PBar;

    my %progress;
    my $count = 100;

    PBar::enable();
    PBar::start(%progress, $count);

    for (my $i = 0; $i < $count; $i++)
    {
      PBar::update(%progress, $i);
    }

    PBar::finish(%progress);

=head1 DESCRIPTION

This module provides a very simple progress bar for use in your loops. Where
"simple" means that we have only one layout for the bar (% finished, bar
itself, ETA) and that it works for a U***x shell. We do not care about thread
safety, here, or any other operating system. The idea is to have something
which just works in our environment without any DEPENDENCIES to other
non-standard packages. We use it mostly for generating parameters which go into
our software. So this is intended to work on the developer side and not on all
sides in general. For a really cool progress bar please refer to the
Term::ProgressBar package from CPAN. This is definitively what you should use in
productive scripts.

=cut


package PBar;

use strict;
use warnings;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    $VERSION     = 0.01;
    
    @ISA         = qw(Exporter);
    @EXPORT      = qw(&pbar_enable);
    %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],
    
    # your exported package globals go here,
    # as well as any optionally exported functions
    @EXPORT_OK   = qw();
}
our @EXPORT_OK;


# EXPORTED GLOBALS     - BEGIN
# EXPORTED GLOBALS     - END

# NON-EXPORTED GLOBALS - BEGIN
# NON-EXPORTED GLOBALS - END

# PRIVATE GLOBALS      - BEGIN
my $private_start = sub {};
# PRIVATE GLOBALS      - END

=head1 FUNCTIONS

=over 4

=item enable / pbar_enable

This function enables the visual progress bar. As default the progress bar is
disabled, meaning that calls to L<B<start>|/start>, L<B<update>|/update> and
L<B<finish>|/finish> produce no output at all. The idea is to only enable the
progressbar if a script is called in verbose mode.

B<EXAMPLE>

    if ($Verbose)
    {
       PBar::enable;
    }

=cut

sub enable
{
    $private_start = sub {
        my ($pbar_hashref, $count) = @_;
        $pbar_hashref->{count} = $count;
        $pbar_hashref->{bar} = "";
        $pbar_hashref->{port} = 0;
        $pbar_hashref->{lw} = 47;

        # set up bar
        $pbar_hashref->{barrierbar} = ($count / ($pbar_hashref->{lw} - 1));
        $pbar_hashref->{nextbar} = $pbar_hashref->{barrierbar};
        $pbar_hashref->{stepbar} = (1 / $pbar_hashref->{barrierbar});

        # set up percentage counter
        $pbar_hashref->{barrierport} = ($count / 99);
        $pbar_hashref->{nextport} = $pbar_hashref->{barrierport};
        $pbar_hashref->{stepport} = (1 / $pbar_hashref->{barrierport});

        # set up time
        $pbar_hashref->{time} = time;
        $pbar_hashref->{tdiff} = 0;

        printf("\n %3.0f%% [%s%*s ETA --:--:-- YY-MM-DD", $pbar_hashref->{port},
               $pbar_hashref->{bar},
               ($pbar_hashref->{lw} - length($pbar_hashref->{bar}) + 1), "]");
    };
}

sub pbar_enable
{
    return enable();
}

=item start

This function enables the visual progress bar. As default the progress bar is

B<EXAMPLE>

  my $FH = open_file_for_read("arbitrary.file");

  while(<$FH>)
  {
    print $_;
  }

  close($FH);

=cut

sub start(\% $)
{
    return &$private_start (@_);
}

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

