# Last modified: 2009-02-13.15
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

Simplex - multidimensional optimiser.

=head1 SYNOPSIS

  use Simplex;

=head1 DESCRIPTION

This module provides functions for a multidimensional simplex.

It needs
=over 4

=item *
reference to a function

=item *
array of starting points

=item *
maximum number of steps

=item *
maximumn number of (re)starts

=item *
tolerance for recognising convergence

=back

It can also use

=over 4

=item *
an array of lower bounds

=item *
an array of upper bounds

=item *
output files for recording progress

=item *
amount of spread/scatter on initial simplex

=back

=head1 PARAMETERS

Parameters are passed by a reference to a single hash like this

    my %simplex_arg;
    $simplex_arg {max_iter} = 1000;

    my %result = simplex (\%simplex_arg);

In detail:

=head2 Obligatory Parameters

=over

=item *
function

  %simplex_arg {func} = \&myfunc;

This is a reference to a function to be minimised.

   sub myfunc () {
       my $a = shift;
       my $b = shift;
       return ($a * $a + ($b - 2) * ($b - 2);
   }
   ...
       %simplex_arg{func} = \&myfunc;

=item *
Initial guess

    my @initial_guess = (1, 3, 5, 4);
    %simplex_arg {ini_pt} = \@initial_guess

This must be a reference to an array.
The number of dimensions should be correct for the function being optimised.

=item *
Maximum iterations

This is the maximum number of iterations per restart. The
optimisation may restart many times.

    %simplex_arg {max_iter} = 1000;

This is single integer;

=item *
Maximum restarts

This is the maximum number of times the simplex may be
started.
It is not really the number of restarts.
If you set it to 1, you will get an initial minimisation only
and no restarts.

=item *
Tolerance

Convergence criteria are always fun.  Currently, we stop when

=over 8

=item

The difference between best and worst corners of the simplex
is less than C<f_tol> and

=item

The difference between the worst value and previous worst value
is less than C<f_tol>.

=back

To set this

    %simplex_args {f_tol} = 10e-5;

=back

=head2 Optional Parameters

=over

=item *
Fixed parameter array

If the hash of simplex arguements contains a member called
B<fix_param>, it will be taken as a reference to an array of
fixed parameters. These are things that will not change from step
to step of minimising.  The calling mechanism is explained below.

=item *
Lower bounds

You can specify lower bounds for the search. If the simplex
tries to go below this in any dimension, it will reject the
point.

    my @lower = (-2, -3, -4, -5);
    %simplex_arg {lower} = \@lower;

This must be a reference to an array. The dimensions must
agree with those of C<@ini_pt>.

=item *
Upper bounds.

This has corresponding behaviour to the lower bounds array.

    my @upper = (10, 20, 10, 3000);
    %simplex_arg {upper} = \@upper;

=item *
Output file

This is not essential. If you want to record the progress,
then set this like

    %simplex_arg {o_file} = 'splx_out';

This example would end up creating two files

=over

=item splx_out_hi.out

=item splx_out_low.out

=back

This first lists the cycle number, function value and test
point for the worst corner of the simplex. The second gives
the same information for the best point on the simplex

When reading this, the worst (highest) point should change from
step to step. The best point will only change every so often.

=item *
Initial Scatter

The simplex is initialised by taking your guess and spreading
corners of the simplex around it, perfectly evenly. In two dimensions,
This would correspond to surrounding your initial guess by a
triangle.

This parameter controls the width of the spread.

    %simplex_args{scatter} = 0.2;

Would result in a spread of 20 % of the size in each
dimension, going S<10 %> up and S<10 %> down. If your initial
guess is C<10>, then the simplex would span a range of C<9> to
C<11> in this dimension. If you have a two dimensional
problem, the initial values would be S<9, 10 and 11>.
If you have a four dimensional problem, the values would be
S<9, 9.5, 1, 10.5, and 11>.

If you do not specify a value, some default like 0.2 will be used.

=back

=head1 RETURN

The routine returns a reference to a hash with three elements

=over

=item *
success / failure

=item *
array containing best point

=item *
value at best point

=back

Access them like this:

    my $result = simplex (\%simplex_arg);
    if ($$result {success} == $SPLX_SUCCESS) {
        print "Simplex happy \n"; }

    my $success = $$result{success};
    my @best = @{${$s_arg{result}}{best}};
    my $best_value = $$result {value};
    print "best value of $best_value at \n@best\n";

The element, C<$$result{success}> can have one of three values
which are exported to the caller:

=over

=item $SPLX_SUCCESS

No problems.

=item $SPLX_TOO_MANY

The routine did not converge within the maximum number of iterations.

=item $SPLX_BROKEN

A programming bug.

=back


=head1 NOTES and OPERATION

=over

=item *

The code is taken from Numerical recipes, but much changed.

=item *

On each restart, the simplex is centred at the best value
previously found.

=item *

If the simplex hits a plateau, nothing terrible should happen.
If the best and worst points are on the plateau, it will just
return, which is a bummer. If some of the points are on a
plateau, the whole simplex will just contract about the best
point.

=item *

Do not look for a parameter with the number of dimensions.
Since perl arrays know how big they are, there is no point
in adding another parameter. We get the dimensionality by
looking at the size of the C<@ini_pt> array.

=head1 EXAMPLE

We have a small function like this:

    sub toyfunc ( \@)
    {
        my ($a, $b) = @_;
        return (($a + 8) * ($a + 8) + ($a - 40) * ($a - 40) + 30 * $a * sin($a)
                + $b * $b * $b * $b_);

    }

Then try
    my $fref = \&toyfunc;
    my @guess = (1, 14);
    my @lower = (-10, -3000);
    my @upper = (20, 230);
    my $max_iter = 1000;
    my $max_restart = 5;
    my $f_tol = 10e-7;

    my %result;
    my %s_arg = (
        func        => \&toyfunc,
        ini_pt      => \@guess,
        lower       => \@lower,
        upper       => \@upper,
        max_iter    => $max_iter,
        max_restart => $max_restart,
        o_file      => 'splx_out',
        scatter     => 0.20,
        f_tol       => $f_tol,
        result      => \%result
                 
    );
    
    my $result = simplex (\%s_arg);

    if ($result {success} == $SPLX_SUCCESS) {
        print "Simplex happy \n";
    } elsif ($result {success} == $SPLX_TOO_MANY) {
        print "Simplex too many\n"; }
    for (my $i = 0; $i < $s_arg{n_dim}; $i++) {
        printf '%4g ', "${${$s_arg{result}}{best}}[$i]"; }
    print "\n";

=head1 ANOTHER EXAMPLE

Imagine we have some parameters which should be passed to the
cost function, but do not vary.  In the caller, we have

    my %s_arg
    my %fix_param;
    $fix_param { num_days } = 5;
    $fix_param { colour } = 'red';

    my %s_arg = (
        func        => \&align_cost,
        ini_pt      => \@guess,
        lower       => \@lower,
        upper       => \@upper,
        max_iter    => $max_iter,
        max_restart => $max_restart,
        o_file      => 'splx_out',
        scatter     => 0.20,
        f_tol       => $f_tol,
        result      => \%result,
        fix_param   => \%fix_param
    );
   [... lots of code ..]
    use lib 'path_to_Simplex';
    use Simplex2;
    my $ret = simplex2 (\%s_arg);

In the cost function, we have something like

   sub align_cost (\@ \%)
    {
        my $var_param = shift;
        my $fix_param = shift;

        my $num_days = $$fix_param { num_days };
        my $colour   = $$fix_param { colour };
        
        my $first_variable = $$var_param[0];
            [.... calculate and return a cost ..]


=head1 BUGS



=cut


package Simplex;

use strict;
use warnings;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use CorbIO qw(:msg :file);

BEGIN {
    use Exporter();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

    # set the version for version checking
    $VERSION     = 0.01; # after everything is adopted to CoRB AND Andrews doc
                         # fixes a here, set to 1.00

    @ISA         = qw(Exporter);
    @EXPORT      = ();
    %EXPORT_TAGS = (all => [qw(&simplex
                               $SPLX_SUCCESS
                               $SPLX_TOO_MANY
                               $SPLX_BROKEN)]
                   );

    # automatically add all tagged functions to the EXPORT_OK list
    Exporter::export_ok_tags('all');

    disable_msg_caller('Simplex');
}
#our @EXPORT_OK; # do we need this while it is already defined above?

# ----------------------- Constants     -----------------------------
# INI_SCATTER is the width of the range over which values will be spread.
# If our initial value is 10, and INI_SCATTER = 0.2, then initial
# values will run from 9 .. 11.
our ($INI_SCATTER, $ALPHA, $BETA, $GAMMA);
*INI_SCATTER = \0.20;

# Some values, using names from numerical recipes version
*ALPHA = \-1.0;
*BETA  = \0.5;
*GAMMA = \2.0;

# When checking bounds, we can either fix them or just report them.
our ($BOUND_OK, $BOUND_BROKEN, $BOUND_REPORT, $BOUND_FIX);
*BOUND_OK     = \0;
*BOUND_BROKEN = \1;

*BOUND_REPORT = \0;
*BOUND_FIX    = \1;

our ($ACCEPT, $REJECT, $SPLX_TOO_MANY, $SPLX_SUCCESS, $SPLX_BROKEN);
*ACCEPT = \0;
*REJECT = \1;

*SPLX_SUCCESS  = \0;
*SPLX_TOO_MANY = \1;
*SPLX_BROKEN   = \2;

# ----------------------- print_simplex  ----------------------------
# For debugging
sub
print_simplex (\@)
{
    my ($simplex_aryref) = @_;

    my $n_dim = $#$simplex_aryref;
    for (my $row = 0; $row <= $n_dim; $row++)
    {
        for (my $i = 0; $i < $n_dim; $i++)
        {
            printf ('%7g ', $$simplex_aryref[$row][$i]);
        }
        print "\n";
    }
}

# ----------------------- rand_array  -------------------------------
# Cute way to take elements of an array and return the same array
# with elements in scrambled order
sub
rand_array (\@)
{
    my ($in_aryref) = @_;

    my @tmp = @$in_aryref;
    my $todo = $#tmp;
    my $end  = $#tmp;
    for ( my $i = 0; $i < $todo; $i++, $end--)
    {
        my $ndx = int (rand ($end + 1));
        my $t = splice (@tmp, $ndx, 1);
        $$in_aryref[$i] = $t;
    }
    $$in_aryref[$todo] = $tmp[0];
}

# ----------------------- spread_simplex ----------------------------
# Given a simplex and central point, spread values evenly around in
# all dimensions.
# Sometimes, we want to spread the simplex by scattering points
# over a range like +- 20 %.
# On other days, we want to use a range specified for each point.
# If the array, ini_range is defined, use it via the spread_range().
# function. If it is not defined, use spread_product().
# As a last frill, we avoid some correlations between parameters by
# taking the range for each and randomising the order.

sub
spread_simplex (\@ \% $)
{
    my ($simplex, $s_arg, $scatter) = @_;
#    my $simplex   = shift;
#    my $s_arg     = shift;
#    my $scatter   = shift;

    my $n_dim  = $$s_arg {n_dim};
    my $delta  = 1.0 / $n_dim;
    my $ini_pt = $$s_arg {ini_pt};
    my $ini_range = undef;
    if (exists ($$s_arg {ini_range}))
    {
        $ini_range = $$s_arg {ini_range};
    }

    my @tmp_row;

    if ( ! defined ($ini_range))
    {          # This approach does not
        for (my $i = 0; $i < $n_dim; $i++)
        { # like non-zero values
            if ($$ini_pt[$i] == 0.0)
            {
                $$ini_pt[$i] = 0.0001;
            } 
        }
    }
    
    for (my $i = 0; $i <= $n_dim; $i++)
    {
        $tmp_row[$i] = $i * $delta ;          # Spread evenly across range
        $tmp_row[$i] -= 0.5;                  # Shift to range -0.5 .. 0 .. 0.5
    }

#   Shift is the full range variation for each parameter.
    my @shift;                     # This could be a property of each parameter
    if (defined ($ini_range))
    {
        for (my $i = 0; $i < $n_dim; $i++)
        {
            if (defined ($$ini_range[$i] ))
            {  # Use specific value if defined
                $shift[$i] = $$ini_range [$i];
            } # else the default
            else
            {
                $shift[$i] = $scatter;
            }
        }
    }
    else
    {
        for (my $i = 0; $i < $n_dim; $i++)
        {
            $shift [$i] = $scatter * $$ini_pt[$i];
        }
    }
    rand_array (@tmp_row);
    for (my $row = 0; $row <= $n_dim; $row++)
    {
        for ( my $i = 0; $i < $n_dim; $i++)
        {
            my $tmp = $$ini_pt [$i] + $tmp_row[$i] * $shift[$i];
            $$simplex[$row][$i] = $tmp;
        }
        push (@tmp_row, shift (@tmp_row));
    }

    if (bound_check ($s_arg, $simplex, $BOUND_FIX) == $BOUND_BROKEN)
    {
        msg_warning ("Bound violation setting up initial simplex\n");
    }
   #print_simplex (@$simplex);
}

# ----------------------- init_y   ----------------------------------
# At start of the day, you have to calculate the function value for
# each of the $ndim+1 points on the simplex.  Return them in the
# $y array.
# This is our first call to the cost function which is quite dangerous.
# The return value is pretty important.
# undef means broken.
sub
init_y ($ \@ $ $ $)
{
    my ($fix_param, $simplex, $y, $func, $names) = @_;
    my $n_dim = $#$simplex;
    my $error = "init_y: cost function seems undefined. i and simplex ";
    for (my $i = 0; $i < $n_dim + 1; $i++)
    {  # @{$$.. nothing easier ?
        my $a  = $$y[$i] = &$func (\@{$$simplex[$i]}, $fix_param, $names);
        if ( ! defined ($a))
        {
            msg_warning ("$error $i, ", "@{$$simplex[$i]}\n"); return undef;
        }
    }
    return 1;
}

# ----------------------- get_three     -----------------------------
# Given a list of numbers, return the indices of the highest,
# second highest and lowest - in that order.
sub
get_three (@)
{
    my $list = \@_;
    my $high = 0;
    my $next = 1;
    my $low  = 1;
    if ($$list[$next] > $$list[$high])
    {
        $high = 1;
        $next = 0;
        $low  = 0;
    }

    for (my $i = 2; $i <= $#$list; $i++)
    {
        if ($$list[$i] > $$list[$high])
        {
            $next = $high;
            $high = $i;
        }
        elsif ($$list[$i] >= $$list[$next] )
        {
            $next = $i;
        }
        elsif ($$list[$i] < $$list[$low])
        {
            $low = $i;
        }
    }
    return ($high, $next, $low);
}

# ----------------------- point_check -------------------------------
# Check if a trial point is within bounds. We may not even have an
# array of bounds, in which case, we skip the loop. For real
# generality, we even allow for an array with just some of the
# members defined and only check bounds when they seem to exist.
sub
point_check (\% \@)
{
    my ($s_arg, $point) = @_;
    #my $s_arg = shift;
    #my $point = shift;

    if (defined ($$s_arg{lower}))
    {
        my $lower = $$s_arg {lower};
        for (my $i = 0; $i <= $#$point; $i++)
        {
            if (defined ($$lower [$i]))
            {
                if ($$point [$i] < $$lower [$i])
                {
                    return $BOUND_BROKEN;
                }
            }
        }
    }
    if (defined ($$s_arg{upper}))
    {
        my $upper = $$s_arg {upper};
        for (my $i = 0; $i <=$#$point; $i++)
        {
            if (defined ($$upper [$i]))
            {
                if ($$point [$i] > $$upper [$i])
                {
                    return $BOUND_BROKEN;
                }
            }
        }
    }
}

# ----------------------- bound_check -------------------------------
# We are given a simplex as well as upper and lower bounds.
# The last argument tells us if we should fix bounds, or just report
# them.
sub
bound_check (\% \@ $)
{
    my ($s_arg, $simplex, $todo) = @_;

    #my $s_arg   = shift;
    #my $simplex = shift;
    #my $todo    = shift;
    my $result  = $BOUND_OK;

    if (defined ($$s_arg{lower}))
    {
        my $lower = $$s_arg {lower};
        for (my $i = 0; $i < $#$simplex + 1; $i++)
        {
            for (my $j = 0; $j < $#$simplex; $j++)
            {
                if (defined ($$lower [$j]))
                {
                    if ($$simplex[$i][$j] < $$lower[$j])
                    {
                        $result = $BOUND_BROKEN;
                        if ($todo == $BOUND_FIX)
                        {
                            $$simplex[$i][$j] = $$lower[$j];
                        }
                    }
                }
            }
        }
    }
    if (defined ($$s_arg{upper}))
    {
        my $upper = $$s_arg{upper};
        for (my $i = 0; $i < $#$simplex + 1; $i++)
        {
            for (my $j = 0; $j < $#$simplex; $j++)
            {
                if (defined ($$upper [$j]))
                {
                    if ($$simplex[$i][$j] > $$upper[$j])
                    {
                        $result = $BOUND_BROKEN;
                        if ($todo == $BOUND_FIX)
                        {
                            $$simplex[$i][$j] = $$upper[$j];
                        }
                    }
                }
            }
        }
    }

    return $result;
}

# ----------------------- get_psum   --------------------------------
sub
get_psum (\@ $)
{
    my ($simplex, $n_dim) = @_;
    my @psum;
    for (my $j = 0; $j < $n_dim; $j++)
    {
        my $sum = 0.0;
        for (my $i = 0; $i < $n_dim + 1; $i++)
        {
            $sum += $$simplex [$i][$j];
            $psum [$j] = $sum;
        }
    }
    return @psum;
}

# ----------------------- splx_out  ---------------------------------
# Output success or path of simplex.
# We want to be able to look at accepted points and best points at
# each step. These are not the same.
# At any cycle, we may replace the highest point. The new point is
# part of the simplex, but it is not necessarily the best point.
# We also need to keep track of the best point at each step.
# Lets write them to separate files.
sub
splx_out ($ $ $ \@ \@ $ $ $)
{
    my ($step,
        $low,
        $high,
        $simplex,
        $y,
        $bestfile,     # For best current point
        $highfile,     # Worst (highest) current point
        $move) = @_;

    #my $step     = shift;
    #my $low      = shift;
    #my $high     = shift;
    #my $simplex  = shift;
    #my $y        = shift;
    #my $bestfile = shift;       # For best current point
    #my $highfile = shift;       # Worst (highest) current point
    #my $move     = shift;

    my $gfmt = ' %.3f';
    my $scorefmt = ' %.4f';
    printf $bestfile "%d $scorefmt",$step, $$y[$low];
    for (my $i = 0; $i < $#$simplex; $i++)
    {
        printf $bestfile $gfmt, $$simplex[$low][$i];
    }
    printf $bestfile "\n";

    printf $highfile "%d $scorefmt",$step, $$y[$high];
    for (my $i = 0; $i < $#$simplex; $i++)
    {
        printf $highfile $gfmt, $$simplex[$high][$i];
    }

    print $highfile  " $move\n";
}


# ----------------------- amotry      -------------------------------
sub
amotry (\% \@ \@ \@ $ $)
{
    my ($s_arg, $simplex, $y, $psum, $high, $fac) = @_;

    #my $s_arg   = shift;
    #my $simplex = shift;
    #my $y       = shift;
    #my $psum    = shift;
    #my $high    = shift;
    #my $fac     = shift;
    my $func    = $$s_arg {func};
    my $n_dim = $#$y;
    my @ptry;  # This will be our trial point
    my $ytry;  # Function value at trial point
    my ($fac1, $fac2);
    my $accept = $REJECT;
    $fac1 = (1.0 - $fac) / $n_dim;
    $fac2 = $fac1 - $fac;
    for ( my $j = 0; $j < $n_dim; $j++)
    {
        $ptry [$j] = $$psum [$j] * $fac1 - $$simplex [$high][$j] * $fac2;
    }
    if (point_check (%$s_arg, @ptry) == $BOUND_BROKEN)
    {
        msg_warning ("Point out of bounds\n", return ($$y[$high] + 1, $accept));
    }

    $ytry = &$func ( \@ptry, $$s_arg { fix_param }, \@{$$s_arg{names}});

    if ( $ytry < $$y[$high])
    {   # This is good, so replace highest point
        $accept = $ACCEPT;
        $$y[$high] = $ytry;
        for (my $j = 0; $j < $n_dim; $j++)
        {
            $$psum [$j] += $ptry [$j] - $$simplex[$high] [$j];
        }
        @{$$simplex[$high]} = @ptry;
    }
    return ($ytry, $accept);
}

# ----------------------- simplex_once  -----------------------------
sub
simplex_once (\%)
{
    my ($s_arg) = @_;
    my $func        = $$s_arg {func};
    my $n_dim       = $$s_arg {n_dim};
    my $ini_pt      = $$s_arg {ini_pt};
    my $f_tol       = $$s_arg {f_tol};
    my $max_iter    = $$s_arg {max_iter};
    my $o_file_low  = $$s_arg {o_file_low};
    my $o_file_hi   = $$s_arg {o_file_hi};
    my %result;#      = $$s_arg {result};
    my $fix_param     = undef;
    if (defined ( $$s_arg { fix_param }))
    {
        $fix_param = $$s_arg { fix_param}; 
    }

    my @psum;         # Geometric hack stores vertex sums (saves a nanosecond)
    my $ncycle = 0;   # count of number of cycles through main loop

    my @simplex;
    $simplex[$n_dim][$n_dim - 1] = 0.0;

    my $scatter;
    if (defined ($$s_arg {scatter}))
    {
        $scatter = $$s_arg {scatter};
    }
    else
    {
        $scatter = $$s_arg {scatter} = $INI_SCATTER / 2.0;
    }
#   Now, lets build the simplex.
    spread_simplex (@simplex, %$s_arg, $scatter);
    my @y;                 # Array of function values

    if (! init_y ($fix_param, @simplex, \@y, $func, \@{$$s_arg{names}}))
    {
        msg_warning ("bad return from init_y\n");
        return undef;
    }
    @psum = get_psum (@simplex, $n_dim);
    my $prev_worst = $y[0];
    my $move = '';
    my $prev_best;
    while (1)
    {
        my $ytry;
        my ($high, $next, $ilo) = get_three (@y);
        splx_out( $ncycle, $ilo, $high, @simplex, @y, $o_file_low,
                  $o_file_hi, $move);
        if (defined ($$s_arg {test_fix_param} )) {
            if (($y[$ilo] < $prev_best) || ( ! defined ($prev_best)))
            {
                my $n = \@{$$s_arg{names}};
                my $tp = \%{$$s_arg{test_fix_param}};
                my $t = &$func (\@{$simplex[$ilo]}, $tp, $n);
                print O_FILE_TEST "$ncycle $t\n";
            }
        }
        $prev_best = $y[$ilo];
        $move = '';
        my @foo = @{$simplex[$ilo]};  # Save best point so far
        $result{best} =   \@foo;
        $result{value} =  $y[$ilo];
        $result{ncycle} = $ncycle;
        {
            my $diff = $y[$high] - $y[$ilo];
            if ($diff < $f_tol)
            {
                $diff = abs ($prev_worst - $y[$high]);
                if ($diff < $f_tol)
                {
                    $result { success } = $SPLX_SUCCESS;
                    return (%result) ;
                }
            }
        }

        if ($ncycle >= $max_iter)
        {
            $result {success} = $SPLX_TOO_MANY;
            return (%result);
        }
        $prev_worst = $y[$high];
        my $a;
        ($ytry, $a) = amotry (%$s_arg, @simplex, @y, @psum, $high, $ALPHA);
        if ($a == $ACCEPT)
        {
            $move = 'r';
        }
        if ($ytry <= $y[$ilo])
        {       # Better than best point, try extension
            ($ytry, $a) =amotry (%$s_arg,@simplex, @y, @psum, $high, $GAMMA);
            if ($a == $ACCEPT)
            {
                $move .= 'e'
            }
        } elsif ($ytry >= $y[$next])
        { # Little progress, so try 1d contract
            my $ysave = $y[$high];
            ($ytry,$a) = amotry (%$s_arg, @simplex, @y, @psum, $high, $BETA);
            if ($a == $ACCEPT)
            {
                $move .= 'c'
            }
            if ($ytry >= $ysave)
            {     # No improvement, do complete contract
                $move .= 'a';
                for (my $i = 0; $i < $n_dim + 1; $i++)
                {
                    my @tmp;
                    if ( $i != $ilo)
                    {
                        for (my $j = 0; $j < $n_dim; $j++)
                        {
                            $tmp[$j] = 0.5 * ($simplex[$i][$j] +
                                              $simplex[$ilo][$j] );
                        }
                        @{$simplex[$i]} = @tmp;
                        $y[$i] = &$func (\@tmp, $fix_param,\@{$$s_arg{names}});
                    }
                }
                @psum = get_psum (@simplex, $n_dim);
            }
        }
        $ncycle++;
    }
}

# ----------------------- check_param   -----------------------------
sub
check_param (\% $)
{
    my ($s_arg, $name) = @_;
    if ( ! defined ($$s_arg { $name }))
    {
        msg_error_and_die ("Simplex parameter \'$name\' not defined\n");
    }

}

# ----------------------- input_sane    -----------------------------
# Check for obligatory arguments for simplex
sub
input_sane (\%)
{
    my ($s_arg) = @_;
    foreach my $n ('func', 'ini_pt', 'max_iter', 'max_restart', 'f_tol')
    {
        check_param (%$s_arg, $n); 
    }
}

# ----------------------- simplex       -----------------------------
# This is the callable simplex interface.
# It will call simplex_once() for up to $max_restart times.
sub
simplex (\%)
{
    use IO::Handle;
    my ($s_arg) = @_;
    my %r;

    input_sane (%$s_arg);

    my $max_restart = $$s_arg {max_restart};

    if (defined ($$s_arg {o_file}) )
    {
        my $o_file = $$s_arg {o_file };
        my $o_file_hi  = $o_file . '_hi.out';
        my $o_file_low = $o_file . '_low.out';

        *O_FILE_HI = open_or_die ("$o_file_hi", ">");
        *O_FILE_LOW = open_or_die ("$o_file_low", ">");

        O_FILE_HI->autoflush(1);
        O_FILE_LOW->autoflush(1);
    }

    if (defined ($$s_arg {test_fix_param} ))
    {
        print "Doing test file\n";
        my $o_file;
        if (defined ($$s_arg {o_file}) )
        {
            $o_file = $$s_arg {o_file } . '_test.out';
        }
        else
        {
            $o_file = 'default_test.out';
        }
        *O_FILE_TEST = open_or_die ("$o_file", ">");
        O_FILE_TEST->autoflush (1);
    }

    if ( ! defined ($$s_arg {n_dim}))
    {
        $$s_arg {n_dim} = $#{$$s_arg{ini_pt}} + 1;
    }

    $$s_arg{o_file_low}  = \*O_FILE_LOW;
    $$s_arg{o_file_hi}   = \*O_FILE_HI;
    my $prev_best = undef;
    my $i;
    for ($i = 0; $i < $max_restart; $i++)
    {
        %r = simplex_once (%$s_arg);

        if ( ! %r )
        {
            return ({ 'success' => $SPLX_BROKEN });     # an anonymous hash
        }
        @{$$s_arg {ini_pt}} = @{$r{best}};
        my $improve = abs( $prev_best - $r{value});
        if ($improve < $$s_arg{f_tol}) 
        {
            last;
        }
        $prev_best = $r{value};
        $$s_arg {scatter} /= 2.0;
        if (exists ($$s_arg {ini_range}))
        {
            my $ini_range = $$s_arg {ini_range};
            for ( my $j = 0; $j < $#$ini_range; $j++)
            {
                $$ini_range [$j] /= 2.0; 
            }
        }

        print O_FILE_HI "\n";
        print O_FILE_LOW "\n";
    }
    close (O_FILE_HI);
    close (O_FILE_LOW);
    if (defined ($$s_arg {test_fix_param} ))
    {
        close (O_FILE_TEST);
    }
    $r{restart} = $i;

    return %r;
}

1;

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
