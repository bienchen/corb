#!/usr/bin/perl
#!@PERL@
# -*- perl -*-
# @configure_input@
# Last modified: 2009-09-30.23


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


# SETTINGS         - BEGIN
BEGIN
{
  # try to get bash compatible shell
  $ENV{'SHELL'} = '@SHELL@' if exists $ENV{'DJGPP'};
  $|++;
}
# SETTINGS         - END


# PERL packages    - BEGIN
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
# PERL packages    - END


# PRIVATE packages - BEGIN
use lib '@CORB_PERL5LIB@';
use CorbIO qw(:all);
# PRIVATE packages - END


# CONSTANTS        - BEGIN
#our %Targets = ( 'txi' => 1 );
#our %Modes = ( 'man' => 1 );
# CONSTANTS        - END


# GLOBALS          - BEGIN
# GLOBALS          - END


# FUNCTIONS        - BEGIN
# Message wrapper to be used for warning messages so we are able to simulate
# -Werror behaviour: Instead of msg_warning msg_error_and_die will be called.
my $Private_func_msg_warning  = \&msg_warning;
sub w_msg_warning
{
    return &$Private_func_msg_warning (@_);
}

# parse the arguments of the script and store them in a hash
#   parseargs(argument_hashref)
sub parse_args(\%)
{
    my ($argument_hashref) = @_;
    my $werror;
    my $i;

    # wrap the signal handler for warnings to copy them to our message space
    local $SIG{__WARN__} = sub { msg_error_and_die("@_") };

    # set defaults

    # parse @ARGV
    GetOptions(
        # data
        'preamble=s' => \$argument_hashref->{preamble},
        'macros=s'   => \$argument_hashref->{macros},
        # strip settings
#        'mode=s'          => \$argument_hashref->{mode},
        # output settings
        'out=s'      => \$argument_hashref->{out},
#        'target-format=s' => \$argument_hashref->{target},
        # misc
        'Werror'     => \$werror,
        # info
        'verbose!'   => sub { enable_verbose; },
        'help'       => sub { pod2usage(-exitval => 0, -verbose => 1) },
        'man'        => sub { pod2usage(-exitval => 0, -verbose => 2) }
              );

    # first message of what we are doing, here
    msg_verbose('CoRB Documentation tool', "\n", '-----------------------',
                "\n", ' Stripping user manuals from source...', "\n"); 

    # Warning options
    if (defined($werror))
    {
        $Private_func_msg_warning  = \&msg_error_and_die;
    }

    # verify option arguments
    if (defined($argument_hashref->{preamble}))
    {
        if (! -r $argument_hashref->{preamble})
        {
            msg_error_and_die("Preamble file (option preamble)",
                              "\"$argument_hashref->{preamble}\" does not ",
                              "exist or is not readable.\n");
        }
    }
    else
    {
        msg_error_and_die("Option preamble is mandatory but missing\n");
    }

    if (defined($argument_hashref->{macros}))
    {
        if (! -r $argument_hashref->{macros})
        {
            msg_error_and_die("Macro definition file (option macros)",
                              "\"$argument_hashref->{macros}\" does not ",
                              "exist or is not readable.\n");
        }
    }

    if (defined($argument_hashref->{out}))
    {
        if (-r $argument_hashref->{out})
        {
            w_msg_warning("Output file (option out) ",
                          "\"$argument_hashref->{out}\" already exists.\n");
        }
    }
    else
    {
        msg_error_and_die("Option out is mandatory but missing\n");
    }

    # get command line argument (an argument given without option)
    # fetch file name
    if ($#ARGV < 0)
    {
        msg_error_and_die("Input file(s) missing. Try \"--help\" or \"--man\" ",
                          "for more information.\n");
    }

    for ($i = 0; $i <= $#ARGV; $i++)
    {
        if ( ! -r $ARGV[$i])
        {
            msg_error_and_die("Input file does not exist or is not readable: ",
                              "$ARGV[$i]\n");
        }
    }

    $argument_hashref->{input} = \@ARGV;

    # target
#    if (defined($argument_hashref->{target}))
#    {
#        if (! defined($Targets{$argument_hashref->{target}}))
#        {
#            msg_error_and_die('Unknown target "', $argument_hashref->{target},
#                              '" to option target-format, allowed: "',
#                              join('", "', keys(%Targets)),'"', "\n");
#        }
#    }

    # doc mode
#    if (defined($argument_hashref->{mode}))
#    {
#        if (! defined($Targets{$argument_hashref->{mode}}))
#        {
#            msg_error_and_die('Unknown strip mode "', $argument_hashref->{mode},
#                              '" to option mode, allowed: "',
#                              join('", "', keys(%Modes)),'"', "\n");
#        }
#    }
#    else
#    {
#        msg_error_and_die("Mandatory option mode missing.\n");
#    }
}

# validate_texinfo_file_or_die($file)
# Verify a file to be a valid texinfo chapter.
# Has no return value, exits on error.
#   validate_texinfo_file_or_die($file)
sub validate_texinfo_file_or_die($)
{
    my ($file) = @_;
    my @lines = file_2_array_or_die($file);
    my $i;

    for ($i = 0; $i <= $#lines; $i++)
    {
    #  print $lines[$i];
# search files for certain items
# verify texinfo code -> iterative function, which gets lines and one object for states; e.g. if in "example" env, wait-status for "end exapmple" is set in object
    }
}
# FUNCTIONS        - END


### MAIN           - BEGIN
my %arg_hash;
my %priv_macros;

# fetch options/ arguments
parse_args(%arg_hash);

# Fetch additional macros
%priv_macros = read_texinfo_macros_or_die($arg_hash{macros});
# NEXT: Load & understand additional macro file
#     - load macros for macro check
#     - substitution
#       - existing macros can still be used
#       - changes to stuff written by developer
#       - changed macros: might be followed by special tag added to modified
#         header
#     - check for known/ unknown macros

# NEXTNEXT: if option "out" is omitted, write to stdout

# NEXTNEXTNEXT: verify preamble file
# Verify preamble to be valid texinfo chapter
msg_verbose('  Preamble file: ', $arg_hash{preamble}, "\n");
validate_texinfo_file_or_die($arg_hash{preamble});

# strip documentation

### MAIN           - END


__END__


=head1 NAME

doctool - ...

=head1 SYNOPSIS

B<doctool> [options] ...

=head1 DESCRIPTION

Describe first task.

First Task of the doctool: Help producing the CORB User Manual
- Assemble a file pulled in by the master texinfo file (has to be chapter)
- requires a preamble from file
- strips from all files given
- strips in alphabetical order
- creates menues/ updates menue if already section node/ menue exists
- does not link chapter node into master
- checks file
- doctool -preamble <file> -out <file> <input files>

Later:
- ARGV missing file message: Determined from expected input

=head1 OPTIONS

You do not need to provide the full option string, an unambiguous prefix is
enough.

=over 8

=item B<--preamble <file>>

Give a text to be placed at the beginning of an output file. Format has to fit
the output.

For stripping a tool-description for a user manual, this has to be the beginning
of a chapter which can be pulled in by the master Texinfo file. Please note
that the preamble is I<NOT> the actual file to be pulled, but the file named by
B<--out>.

=item B<--macros <file>>

Here you may name a file which carries additional macro definitions for
Texinfo. Since all Texinfo input is checked for known and unknown macros, this
is mandatory if you use your own macros.

Furthermore, macros from this file are used for substitution in chapter and
section headers as well as in node names. Since Texinfo has problems
recursively resolving own macros, we do this in headings. For node names we
stay with your name but give it a new menu entry name.

=item B<--out <file>>

Name a file to write to. Existing files will be overwritten.

=item B<--target-format <txi>>

Define the output format for stripped content. If omitted, the target is
guessed by B<doctool>. Allowed values are

=over

=item C<txi>

Produce Texinfo output.

=back

=item B<--mode <man>>

Define what should be stripped: Api documentation for developers or an user
manual. This option is mandatory. Allowed values are

=over

=item C<man>

Strip user manual content.

=back

=item B<--verbose>

Be verbose.

=item B<--man>

Print a man page for the program. This is a more detailed description of the
script.

=item B<--help>

Print a help message and exit. This only contains the synopsis and the options
description.

=cut

# old considerations: 
#- strip info from c/ pl/ m4/ ... source files
#  - check data: linewidth, right use of item, ...
#  - modes: tool, api -> check completeness, form -> mode
#  - strip files
#- assemble into texinfo files to be included
#  - allow for prefix, suffix
#  - assembly modes: alphabetically or as handed in
#  - define output file name
#  - guess format from file name -> target-format
#- update menues/ nodes of texinfo stuff
#  - by providing single files or toplevel parent
#  - create node label from section name if formatted
#- produce man pages from striped content
#  - no. has to fit mode api or tool
#- pay attention to 'xxx' comments (via Perl-module to be reusable by reformat)
#- all checks via module to be reusable


# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
