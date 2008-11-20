# Last modified: 2008-11-20.17
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
itself, ETA) and that it works for a shell. We do not care about thread
safety, here, or any other operating system. The idea is to have something
which just works in our environment without any DEPENDENCIES to other
non-standard packages. We use it mostly for generating parameters which go into
our software. So this is intended to work on the developer side and not on all
sides in general. For a really cool progress bar please refer to the
B<C<Term::ProgressBar>> package from CPAN. This is definitively what you should
use in productive scripts.

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
    @EXPORT      = ();
    %EXPORT_TAGS = (PBarControl => [qw(&pbar_enable
                                       &pbar_disable
                                       &pbar_start
                                       &pbar_update
                                       &pbar_finish)]);
    @EXPORT_OK   = qw(&pbar_enable
                      &pbar_disable
                      &pbar_start &pbar_update
                      &pbar_finish);
}
our @EXPORT_OK;


# EXPORTED GLOBALS     - BEGIN
# EXPORTED GLOBALS     - END

# NON-EXPORTED GLOBALS - BEGIN
# NON-EXPORTED GLOBALS - END

# PRIVATE GLOBALS      - BEGIN
my $private_func_start  = sub {};
my $private_func_update = sub {};
my $private_func_finish = sub {};
# PRIVATE GLOBALS      - END

=pod

=head1 FUNCTIONS

=head2 enable / pbar_enable

This function enables the visual progress bar. As default the progress bar is
disabled, meaning that calls to L<C<start()>|"start / pbar_start">,
L<C<update()>|"update / pbar_update"> and L<C<finish()>|"finish / pbar_finish">
produce no output at all. The idea is to only enable the progress bar if a
script is called in verbose mode.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

     if ($Verbose)
     {
       PBar::enable;
     }

=back

=cut

sub enable
{
    $private_func_start = sub {
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

    $private_func_update = sub {
        my ($pbar_hashref, $current, ) = @_;
        my $need_rewrite = 0;

        if ($current > $pbar_hashref->{nextport})
        {
            $need_rewrite = 1;

            $pbar_hashref->{nextport} += $pbar_hashref->{barrierport};  
            $pbar_hashref->{port} = ($current * $pbar_hashref->{stepport});
        }
        if ($current > $pbar_hashref->{nextbar})
        {
            $need_rewrite = 1;

            $pbar_hashref->{bar} = sprintf("% *s",
                                           ($current * $pbar_hashref->{stepbar})
                                           , "");
            $pbar_hashref->{bar} =~ s/\s/=/g;
            $pbar_hashref->{nextbar} += $pbar_hashref->{barrierbar};
        }

        if ($need_rewrite)
        {
            $pbar_hashref->{tdiff} += (time - $pbar_hashref->{time});
            my $eta = (  (99 - $pbar_hashref->{port})
                       * $pbar_hashref->{tdiff})     /$pbar_hashref->{port};
            my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
                = localtime(time + $eta);

            #print "$current $pbar_hashref->{nextbar}\n";
            printf("\r %3.0f%% [%s%*s ETA %02d:%02d:%02d %02d-%02d-%02d",
                   $pbar_hashref->{port},
                   $pbar_hashref->{bar},
                   ($pbar_hashref->{lw} - length($pbar_hashref->{bar}) + 1),
                   "]", $hour, $min, $sec, $year % 100,$mon + 1, $mday);

            $pbar_hashref->{time} = time;
        }
    };

    $private_func_finish = sub {
        my ($pbar_hashref) = @_;
        
        $pbar_hashref->{tdiff} += (time - $pbar_hashref->{time});
        
        my $sec = $pbar_hashref->{tdiff} % 60;
        my $min = int(($pbar_hashref->{tdiff} % 3600) / 60);
        my $hou = int($pbar_hashref->{tdiff} / 3600);

        if ($hou > 100)
        {
            $hou = 99;
            $sec = 99;
            $min = 99;
        }

        $pbar_hashref->{bar} = sprintf("% *s", $pbar_hashref->{lw}, "");
        $pbar_hashref->{bar} =~ s/\s/=/g;
        printf("\r 100%% [%s] Elapsed time %02d:%02d:%02d\n",
               $pbar_hashref->{bar},
               $hou, $min, $sec);
    };
}

sub pbar_enable
{
    return enable();
}

=head2 disable / pbar_disable

This function disables the visual progress bar. After calling this function,
calls to L<C<start()>|"start / pbar_start">, L<C<update()>|"update / pbar_update">
and L<C<finish()>|"finish / pbar_finish"> produce no output.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

       PBar::disable;

=back

=cut

sub disable
{
    $private_func_start = sub {};

    $private_func_update = sub {};

    $private_func_finish = sub {};
}

sub pbar_disable
{
    return disable();
}

=pod

=head2 start / pbar_start

This function starts a new progress bar. This means that it is initially
written to the terminal and all variables are set. Only produces some output
after L<C<enable()>|"enable / pbar_enable"> was issued.

=over 4

=item ARGUMENTS

=over 4

=item progressbar

A hash to store all data of a particular progress bar. Is initialised by this
function. Has to be provided as reference while using
L<C<pbar_start()>|"start / pbar_start">.

=item count

Number of loop iterations. Used to calculate step size of the progressbar.

=back

=item EXAMPLE

     my %pbar;
     my $count = 100;

     PBar::start(%pbar, $count);

     for (my $i = 0; $i < $count; $i++)
     {
        print $i;
     }

=back

=cut

sub start(\% $)
{
    return &$private_func_start (@_);
}

sub pbar_start
{
    return &$private_func_start (@_);
}

=head2 update / pbar_update

This function updates the progress bar. Any former printing of the bar to
STDOUT is overwritten by a new bar with updated values. Should be called within
a loop. Calling without issuing L<C<start()>|"start / pbar_start"> before,
will lead to warnings/ errors. This is just because we are accessing some
values which have to be initialised, first. Only produces some output
after L<C<enable()>|"enable / pbar_enable"> was issued.

=over 4

=item ARGUMENTS

=over 4

=item progressbar

A hash to store all data of a particular progress bar. Is modified by this
function. Has to be initialised by L<C<start()>|"start / pbar_start"> before.
Has to be provided as reference while using
L<C<pbar_update()>|"update / pbar_update">.

=item current

No. of steps already processed in the loop.

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

=back

=cut

sub update(\% $)
{
    return &$private_func_update (@_);
}

sub pbar_update
{
    return &$private_func_update (@_);
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
    return &$private_func_finish (@_);
}

sub pbar_finish
{
    return &$private_func_finish (@_);
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

