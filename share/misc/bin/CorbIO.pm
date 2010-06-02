# Last modified: 2010-05-21.15
#
#
# Copyright (C) 2010 Stefan Bienert
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

    use CorbIO qw(:all);

    enable_verbose();

    msg_verbose("Only seen in verbose mode");

    if (is_verbose())
    {
      print("Also only seen in verbose mode");      
    }
    else
    {
      msg_warning("Verbose mode is disabled!");
    }

    my $fh = open_or_die("/file/to.read", '<');

    foreach (<$fh>)    
    {
      print $_;
    }

    close($fh);

    my @lines = file_2_array_or_die("/file/to.read");

    msg_warning("This prints a warning to STDERR");

    my ($header, $sequence) = read_single_fasta_or_die("sequence.fasta");

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

=item * msg - Additional messaging functions 
    (L<C<msg_warning()>|"msg_warning">,
    L<C<msg_error_and_die()>|"msg_error_and_die">,
    L<C<disable_msg_caller()>|"disable_msg_caller">)

=item * file - File handling
    (L<C<open_or_die()>|"open_or_die">,
    L<C<file_2_array_or_die()>|"file_2_array_or_die">,
    L<C<read_single_fasta_or_die()>|"read_single_fasta_or_die">)

=item * misc - Everything only weakly related to IO
    (L<C<gnuplot_histogram()>|"gnuplot_histogram">)

=item * all - all functions from the tags above

=back

Since direct calling of functions via package naming is also allowed, the
function names are not prefixed with the package name.

=cut


package CorbIO;

use strict;
use warnings;

BEGIN {
    use Exporter();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    $VERSION     = 0.02;
    
    @ISA         = qw(Exporter);
    @EXPORT      = ();
    # we try to remove the '&' in front of the functions: Once they are in use,
    # try to run without '&', if everything is O.K.: remove
    # Reason: Documentation states, loading is faster wo '&'
    %EXPORT_TAGS = (verbose => [qw(enable_verbose
                                   &disable_verbose
                                   &is_verbose
                                   msg_verbose)],
                    msg     => [qw(msg_warning
                                   msg_error_and_die
                                   &disable_msg_caller)],
                    file    => [qw(open_or_die
                                   &file_2_array_or_die
                                   &read_single_fasta_or_die)],
                    misc    => [qw(gnuplot_histogram)]
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


# EXPORTED GLOBALS     - BEGIN <- automatically exported
# EXPORTED GLOBALS     - END


# NON-EXPORTED GLOBALS - BEGIN <- not automatically but still exportABLE!
# NON-EXPORTED GLOBALS - END


# PRIVATE GLOBALS      - BEGIN <- no access from outside
my $Private_func_msg_verbose  = sub {};
my $Private_func_is_verbose   = sub { return 0; };
my %Disable_msg_for = ('CorbIO' => 1);
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
    $Private_func_msg_verbose = sub {
        my (@msg) = @_;

        print @msg;
    };

    $Private_func_is_verbose = sub { return 1; };
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
    $Private_func_msg_verbose = sub {};

    $Private_func_is_verbose = sub { return 0; };
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
    return &$Private_func_is_verbose;
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

     # boring example
     my $formatted_float = sprintf("%3.3f", 1/7);

     msg_verbose("1/7 = ${formatted_float}\n");

     # instead of working with is_verbose() for verbose printing, we can 
     # compute everything we need just in the msg_verbose() statement. This
     # saves us from checking if we are in verbose state and then using
     # msg_verbose(), checking again, to print a message to the right stream
     # and of correct format.
     my %pets = (
         'birds'     => [ 'Canary', 'Great tit', 'Stork' ],
         'goggies'   => [ 'Beagle', 'Dachshund', 'Hot', 'Rottweiler' ],
         'dinosaurs' => [ 'Tyrannosaurus', 'Grandma', 'Apatosaurus' ]
         );
     msg_verbose("Pets in shop:\n", &{ sub
                                       {
                                           my @lop;

                                           foreach (keys(%pets))
                                           {
                                               push(@lop, "   ${_}:\n");

                                               foreach (@{$pets{$_}})
                                               {
                                                   push(@lop, "      ${_}\n");
                                               }
                                           }

                                           return @lop;
                                       }}, "\n");

=back

=cut

sub msg_verbose
{
    return &$Private_func_msg_verbose (@_);
}

=head2 disable_msg_caller

This function disables a given package/ module/ script from being referenced by
L<C<msg_warning()>|"msg_warning"> and
L<C<msg_error_and_die()>|"msg_error_and_die">. If a disabled item is found as a
caller of a message function, we go one up in the stack trace.

=over 4

=item ARGUMENTS

=over 4

=item caller

Caller of a message function to be excluded from being named as a caller. Has
to be provided in scalar context.

=back

=item EXAMPLE

    disable_msg_caller('CorbIO');

=back

=cut

sub disable_msg_caller($)
{
    my ($caller) = @_;

    $Disable_msg_for{$caller}++;
}

# internal routine
# get the line and file of a subroutine call.
# takes no arguments, just returns file and lineno..
sub s_corbio_get_call_info()
{
    my $level = 0;
    my %call_info;

    {
        @call_info{ qw(pack file line) } = caller($level);
        
        #print "   caller(${level}):".caller($level)."\n";
        $level++;
        redo if $Disable_msg_for{$call_info{pack}};
    }
    unless (defined $call_info{pack})
    {
        die("A somewhat impossible error occured: Could not find calling "
            ."package.");
    }
    
    #print "   pack:       ".$call_info{pack}."\n";
    #print "   file:       ".$call_info{file}."\n";
    #print "   line:       ".$call_info{line}."\n";

    return ($call_info{file}, $call_info{line});
}

=head2 msg_warning

Just writes a warning message to C<STDERR>. There is nothing special about this
but the message is preceded by "file:line:WARNING:". This goes into the
direction of the GNU coding standards and should give as a unified look for all
messages of our scripts.

=over 4

=item ARGUMENTS

=over 4

=item message list

The list of strings to be written.

=back

=item EXAMPLE

     msg_warning("This is a warning.");

=back

=cut

sub msg_warning(@)
{
    my (@msg) = @_;
    my ($file, $line) = s_corbio_get_call_info();

    print STDERR "${file}:${line}:WARNING: ", @msg;
}


=head2 msg_error_and_die

Writes an error message to C<STDERR> and exits the script. There is nothing
special about this but the message is preceded by "file:line:ERROR:". This goes
into the direction of the GNU coding standards and should give as a unified
look for all messages of our scripts. 

=over 4

=item ARGUMENTS

=over 4

=item message

The list of strings to be written.

=back

=item EXAMPLE

     msg_warning("This is a warning.");

=back

=cut

sub msg_error_and_die(@)
{
    my (@msg) = @_;
    my ($file, $line) = s_corbio_get_call_info();
    my $status = $!;

    print STDERR "${file}:${line}:ERROR: ", @msg;

    # simulate die() behaviour
    if ($status == 0)
    {
        $status = $? >> 8;
    }
    if ($status == 0)
    {
        $status = 255;
    }
    exit $status;
}


=head2 open_or_die

Returns a file handle to a file opened in a given mode. If opening fails, the
script dies. If mode is omitted, we assume that the file should be opened for
reading.

=over 4

=item ARGUMENTS

=over 4

=item filename

The file to be opened.

=item mode

Valid modes are:

=over

=item * ">"

    Open a file for writing. If the file does not exist it will be created. If
    the file already exists, the script will C<die>. This is not usual B<Perl>
    behaviour and could be a little bit annoying while developing new scripts.
    On the other hand you will name your first-born after us, the first time
    this prevents loss of data from an overnight run. Just believe us, we have
    learned our lessons...

    ...several times.

=item * "<"

    Open a file for reading.

=item * ">>"

    Open a file for appending. If file exists, data written to it occurs at
    the end. If file does not exists, it will be created.

=item * "+>" "+<"

    Open a file for reading and writing. Using "+>" would lead to clobbering
    of the file first using B<Perl> C<open>. Since we do not overwrite nothing
    here, trying to do so will exit the script.

=item * "|-"

    The file to be opened is a shell command. The command is executed and a
    pipe towards it is established.

=item * "-|"

    The file to be opened is a shell command. The command is executed and an
    outgoing pipe is established.

=back

=back

=item EXAMPLE

     my $file_2_read = open_or_die("/file/to/be.read", "<");
     my $file_2_write = open_or_die("file/to/copy.to", ">");

     foreach (<$file_2_read>)
     {
       print $file_2_write $_;
     }

     close($file_2_read);
     close($file_2_write);

=back

=cut

sub open_or_die
{
    my ($filename, $mode) = @_;
    local *FH;

    if (! defined($filename))
    {
        msg_error_and_die ("Filename missing in attempt to open file.\n");        
    }

    if (! defined($mode))
    {
        $mode = '<';
    }

    # check if we are trying to open an existing file in write/ write-read mode
    if ((($mode eq '>') || ($mode eq '+>')) && (-e $filename))
    {
        msg_error_and_die ("Unable to create file \"$filename\": File "
                          ."already exists\n");
    }

# SB: 18.08.09 changed from open(FH, $mode, $filename) to be able to open
#     STDOUT as a filehandle
    open(FH, $mode.$filename) or
        msg_error_and_die("Unable to open file \"$filename\": $!\n");

    return *FH;
}

=head2 file_2_array_or_die

Returns an array filled with the content of a file. If the file can not be
opened, the script dies.

=over 4

=item ARGUMENTS

=over 4

=item filename

The file to be opened.

=back

=item EXAMPLE

     my @lines;

     @lines = file_2_array_or_die("testfile.sth");

=back

=cut

sub file_2_array_or_die
{
    my $filehandle = open_or_die("@_", '<');

    my @lines = <$filehandle>;

    close($filehandle);

    return @lines;
}

=head2 read_single_fasta_or_die

Read a FASTA file carrying a single sequence and return the sequence with its
header.

=over 4

=item ARGUMENTS

=over 4

=item filename

The FASTA file. The file given is assumed to be in FASTA format, hence its
extension is arbitrary.

=back

=item EXAMPLE

     my ($header, $sequence) = read_single_fasta_or_die("sequence.fas");

=back

=cut

sub read_single_fasta_or_die
{
    my $filehandle = open_or_die("@_", '<');
    my $header;
    my $sequence = '';

    foreach(<$filehandle>)
    {
	if($_ =~ />(.*)/)
	{   
            $header = $1;
	}
	elsif($_ =~ /(.*)/)
	{
            $sequence .= $1;
	}
    }

    close($filehandle);

    if (!defined($header))
    {
        msg_error_and_die ("FASTA file \"@_\" is missing a header ",
                           "line.\n");
    }

    return ($header, $sequence);
}

=head2 gnuplot_histogram

Simple interface to C<gnuplot> for creating histograms.

Histograms can be produced just by providing the data, but can also be tuned
using C<%settings> parameter. Everything is written to a file with the terminal
(file type) determined from the extension of a filename provided.

Since plotting with C<gnuplot> is not an easy business, this function does not
cover all possible errors that might occur. Hence, we do not recommend using
this function in a productive script but for data evaluation. If you stumble
upon an uncovered problem, feel free to add a solution to this function.

=over 4

=item ARGUMENTS

=over 4

=item data

A hash containing the data where the keys are the labels for the x-axis and the
values are the data points. Since the hash is used as a reference, no anonymous
hashes allowed, here. The keys will be used in sorted order.

As value each entry has to carry an array of values belonging to its key. This
will be presented as a group of bars in the plot.

The members of an array may be named using the C<legend> setting.

=item filename

Name where to store the histogram. The extension determines the terminal used.
Currently supported are:

=over 4

=item * pdf, triggered by '*.pdf'

=item * postscript, triggered by '*.ps'

=back

=item settings

A hash to fine tune the plot. Keys describe the item to tune. Currently
supported are:

=over 4

=item * title, Heading of the plot

=item * xlabel, Naming the x-axis

=item * ylabel, Naming the y-axis

=item * legend, Array of names for the data lines; Note that you have to name all lines or none.

=item * termsettings, A string carrying options for the terminal, use as in C<gnuplot>

=item * xtics, Optionstring for xtics, use as in C<gnuplot>

=item * yupperrange, Defines the upper y range in the C<plot> command

=back

=back

=item EXAMPLE

    my %data = (0.1 =>[10], 0.2 =>[20], 0.3 =>[30], 0.4 =>[40], 0.5 =>[50]);
    gnuplot_histogram(%data, "data1.ps",(xlabel => 'Portion',
                                         ylabel => 'Occurrences',
                                         termsettings => 'color'));

    %data = (0.1 => [10, 50, 5],
             0.2 => [20, 40, 30],
             0.3 => [30, 30, 55],
             0.4 => [40, 20, 80],
             0.5 => [50, 10, 105]);
    gnuplot_histogram(%data, "data2.ps",(xlabel => 'Portion',
                                         ylabel => 'Occurrences',
                                         legend => ['rise', 'fall', 'misc'],
                                         termsettings => 'color'));

=back

=cut

sub gnuplot_histogram(\% $; %)
{
    my ($data_hashref, $filename, %settings) = @_;
    my $run = 1;
    my $mode = '';
    my $msg = '';
    my $trm = 'not set';
    my $i;
    my @xtic;
    my $using;
    my @legend;
    my ($data1, $data2);
    my $termopts = '';
    my $add_xtics = '';
    my $gplt_script = '';
    my $yup = '';

    if (keys(%{$data_hashref}) == 0)
    {
        return;
    }

    # determine gnuplot terminal
    if ($filename =~ /\.pdf$/)
    {
        $trm = 'pdf';
    }
    elsif ($filename =~ /\.ps$/)
    {
        $trm = 'postscript';
    }

    #    eval %settings
    if (%settings)
    {
        if (defined($settings{title}))
        {
            $gplt_script .= "set title \\\"$settings{title}\\\"\n";
        }
        if (defined($settings{xlabel}))
        {
            $gplt_script .= "set xlabel \\\"$settings{xlabel}\\\"\n";
        }
        if (defined($settings{ylabel}))
        {
            $gplt_script .= "set ylabel \\\"$settings{ylabel}\\\"\n";
        }
        if (defined($settings{legend}))
        {
            @legend = @{$settings{legend}};
        }

        # terminal settings
        if (defined($settings{termsettings}))
        {
            $termopts = $settings{termsettings};
        }

        # xtics addons
        if (defined($settings{xtics}))
        {
            $add_xtics = $settings{xtics};
        }

        # y-ranges
        if (defined($settings{yupperrange}))
        {
            $yup = $settings{yupperrange};
        }
    }

    # create script
    #    set terminal & output file
    $gplt_script .= "set term $trm $termopts\nset out '$filename'\n";

    #    standard settings for histograms
    $gplt_script .= "set style histogram clustered gap 1\n"
                   ."set style data histograms\n"
                   ."set style fill solid 1.00 border -1\n";

    #    add xtics
    @xtic = sort(keys(%{$data_hashref}));
    $gplt_script .= "set xtics $add_xtics (\\\"$xtic[0]\\\" 0";
    for($i = 1; $i <= $#xtic; $i++)
    {
        $gplt_script .= ", \\\"$xtic[$i]\\\" $i";
    }
    $gplt_script .= ")\n";

    # create data
    $data1 = '';
    foreach(sort(keys(%{$data_hashref})))
    {
        unless (uc(ref($data_hashref->{$_})) eq 'ARRAY')
        {
            msg_error_and_die("Value for key \"$_\" in data hash is not an ",
                              "array in ", (caller(0))[3], ".\n");
        }

        $data1 .= join(' ', @{$data_hashref->{$_}})."\n";
    }
    $data1 .= "e\n";

    # create 'using' part of plot command
    if (!@legend)
    {
        foreach (@{$data_hashref->{$xtic[0]}})
        {
            push(@legend, '');
        }
    }
    else
    {
        $gplt_script .= "set key reverse Left left\n";
    }

    $using = '';
    $data2 .= $data1;
    for ($i = 1; $i <= $#legend; $i++)
    {
        $using .= ', \'-\' u '.($i + 1).' t \''.$legend[$i].'\'';
        $data2 .= $data1;
    }

    $gplt_script .= "plot [:][0:$yup] '-' using 1 t '".$legend[0]."'".$using
        ."\n";

    $gplt_script .= $data2;

    # execute script
    while ($run)
    {
        $run = 0;

        unless(open(FH, "echo \"$gplt_script\" | gnuplot 2>&1 |"))
        {
            msg_error_and_die ('Failure at running gnuplot with script ',
                               "\"$gplt_script\": $!\n");
        }

        # catch error messages
        foreach (<FH>)
        {
            #print $_;
            if ($_ =~ /line\s+\d+:\s+unknown\s+or\s+ambiguous\s+terminal\s+
                      type;\s+type\s+just\s+'set\s+terminal'\s+for\s+a\s+list/x)
            {
                $gplt_script = 'set terminal';
                $mode = 'termtypes';
                $run = 1;
                last;
            }
            elsif (    (($mode eq '') or ($mode eq 'uncatched'))
                   and ($_ !~ /^\s*\^\s*$/)
                   and ($_ !~ /^$/)
                   and ($_ !~ /gnuplot\>\s+set\s+term\s+/))
            {
                $mode = 'uncatched';
                $msg .= $_;
            }
            elsif ($mode eq 'termtypes')
            {
                $msg .= $_;
            }
        }
        close(FH);
    }

    if ($mode eq 'termtypes')
    {
        msg_error_and_die ("Unsupported format of output file \"$filename\".",
                           $msg);
    }
    elsif ($mode eq 'uncatched')
    {
        msg_error_and_die ("Uncatched error occured for script\n\n",
                           $gplt_script,
                           "\ngnuplot output:\n",
                           $msg);
    }
    #print $gplt_script;
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

