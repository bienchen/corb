#!@PERL@ -w
# -*- perl -*-
# @configure_input@
# Last modified: 2008-04-30.00

# Copyright (c) 2007 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.


# SETTINGS         - BEGIN
BEGIN
{
  # try to get bash compatible shell
  $ENV{'SHELL'} = '@SHELL@' if exists $ENV{'DJGPP'};
  $|=1; # dissable line buffer for print statements
}
# SETTINGS         - END


# PERL packages    - BEGIN
use strict;
use warnings;
use File::Basename;
use File::Copy;
use Getopt::Long;
use Pod::Usage;
# use sigtrap; on aborting signals restore backup file?
# PERL packages    - END


# PRIVATE packages - BEGIN
# PRIVATE packages - END


# GLOBALS          - BEGIN
my %Known_Extensions   = (
                           'm4' => 'M4 macro language',
                         );
my %Format_Preamble    = (
                           'm4' => \&format_preamble_m4,
                         );
my %Define_Preamble    = (
                           'm4' => \&define_preamble_m4,
                         );
my %Define_Codesection = (
                           'm4' => \&define_codesection_m4,
                         );
my %Define_End         = (
                           'm4' => \&define_end_m4,
                         );
my $VERBOSE          = 0;
my $BCKEXT           = "bkp";
my $DWIDTH           = 1;
# GLOBALS          - END


# FUNCTIONS        - BEGIN
##########################
###   M4 formatting    ###
##########################
# define the start and end line of a M4 source file preamble
#   define_preamble_m4($start_ref, $end_ref, @lines_ref)
sub define_preamble_m4(\$ \$ \@)
{
    my ($start_ref, $end_ref, $lines_ref) = @_;
    my $started = 0;

    $$start_ref = 1;

    for (my $i = 0; $i <= $#{$lines_ref}; $i++)
    {
        if (@{$lines_ref}[$i] =~ /(^\#\s*$|
                                   Last\s+modified\:|
                                   Copyright\s+\(C\)\s+\d+|
                                   See\s+COPYING\s+file\s+in\s+the\s+top\s+
                                   level\s+directory\s+of\s+this\s+tree\s+for
                                   \s+licence|
                                   ^\s*$)/xi)
        {
            # print @{$lines_ref}[$i];
            $started = 1;
        }
        else
        {
            if ($started)
            {   
                $$end_ref = $i;
                return;
            }
        }
    }

    $$end_ref = $#{$lines_ref};
}

# define the start and end line of a M4 source file code section
#   define_codesection_m4($start_ref, $end_ref, @lines_ref)
sub define_codesection_m4(\$ \$ \@)
{
    my ($start_ref, $end_ref, $lines_ref) = @_;
    my $started = 0;

    for (my $i = $$start_ref; $i <= $#{$lines_ref}; $i++)
    {
        if (@{$lines_ref}[$i] !~ /^\s*$/)
        {
            if (@{$lines_ref}[$i] =~ /^\s*dnl\#\s*Local\s+variables\:/i)
            {
               $$end_ref = $i;
               return;
            }
            else
            {
                if (! $started)
                {
                    $$start_ref = $i + 1;
                    $started = 1;
                }
            }
        }
    }

    $$end_ref = $#{$lines_ref};
}

# define the start and end line of a M4 source file end section
#   define_end_m4($start_ref, $end_ref, @lines_ref)
sub define_end_m4(\$ \$ \@)
{
    my ($start_ref, $end_ref, $lines_ref) = @_;
    my $started = 0;

    for (my $i = $$start_ref; $i <= $#{$lines_ref}; $i++)
    {
        if (@{$lines_ref}[$i] =~ /(^\s*dnl\#\s*Local\s+variables\:|
                                   ^\s*$)/xi)
        {
          if (! $started)
          {
              $$start_ref = $i + 1;
              $started = 1;
          }
        }
    }

    $$end_ref = $#{$lines_ref} + 1;
}

# format the preamble of a M4 source file, returns the no. of lines added to
# the preamble
#   format_preamble_m4($start, $end, $file, @lines_aryref)
sub format_preamble_m4($ $ $ \@)
{
    my ($start, $end, $file, $lines_aryref) = @_;
    my $i;
    my $line;

    if ($end > ($#{$lines_aryref} + 1))
    {
        $end = $#{$lines_aryref} + 1;
    }

#    a
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-26.23
#

    $i = $start;
    if (@{$lines_aryref}[$i - 1] !~ /^\#$/)
    {
        # leading spaces
        # trailing spaces
        # empty line
        # text
        form_violation_msg("BLA\n", $file, $i, 0);
    }

    for ($i = $start; $i <= $end; $i++)
    {
        $line = @{$lines_aryref}[$i - 1];
        warning_msg(sprintf ("line %*d: $line", $DWIDTH, $i));
    }

    return 0;
}


##########################
###   Arb. functions   ###
##########################
# Sets the width of the max. lineno. as default width for writting to stdout.
#   set_digits(n)
sub set_dwidth($)
{
    my ($n) = @_;

    $n = log($n)/log(10);
    $DWIDTH =  sprintf("%.0f",($n + 0.5));    
}

# copy a file with a unique file name or find a file with the same content
#   cpy_file(filename, error_msgref)
sub cpy_file($ \$)
{
    my ($file, $error_msgref) = @_;
    my $cpy_name = "";
    my $dir;
    my $filename;
    my $bck = 1;

    # strip directory from file
    ($filename, $dir, undef) = fileparse($file, qr//);

    # search for filename.\d+.bck in dir list
    unless (opendir (DIR, $dir))
    {
        $$error_msgref = "Could not open directory \"${dir}\"\n";
        return "";
    }
    foreach (readdir(DIR))
    {
        if ($_ =~ /$filename\.(\d+)\.$BCKEXT/)
        {
            verbose_msg("    backup found: \"$_\"\n");
            if ($bck <= $1)
            {
                $bck = $1 + 1;
            }
        }
    }
    closedir(DIR);

    # create copy $filename.$bck.$BCKEXT
    $cpy_name = $file."\.".$bck."\.".$BCKEXT;
    unless (copy($file, $cpy_name))
    {
        $$error_msgref = "Could not copy \"$filename\" to \"$cpy_name\"\n";
    }

    return $cpy_name;
}

# write message if in verbose mode
#   verbose_msg(message)
sub verbose_msg($)
{
    my ($message) = @_;

    print $message if $VERBOSE;
}

# write warning
#   warning_msg(message)
sub warning_msg($)
{
    my ($message) = @_;

    print STDERR "WARNING:".$message;
}

# write message for format violation
#   form_violation_msg($message, $file, $line, $col)
sub form_violation_msg($ $ $ $)
{
    my ($message, $file, $line, $col) = @_;

    printf(STDERR "${file}:%*d:%2d: $message\n", $DWIDTH, $line, $col);
}

# check whether a known file format is given
#   validate_file_format(file-extension, error_msg_ref)
sub validate_file_format($ \$)
{
    my ($extension, $error_msgref) = @_;

    if (! defined($Known_Extensions{$extension}))
    {
        $$error_msgref = "Unknwon file format: ${extension}\n";
        return 0;
    }

    return 1;
}

# strip extension from file and return a format identifier. Returns "" on error.
#   get_file_format(filename, error_msg_ref)
sub get_file_format($ \$)
{
    my ($filename, $error_msgref) = @_;
    my $format = "";
    my(undef, undef, $extension) = fileparse($filename, qr{\..*});

    if ($extension ne "")
    {
        $extension = join('', (split(/^\./, $extension)));

        # get extension
        if (validate_file_format($extension, $$error_msgref))
        {
            $format = $extension;
        }
    }
    else
    {
        $$error_msgref = "Could not get extension from: \"${filename}\"\n";
    }


    return $format;
}

# parse the arguments of the script and store them in a hash
#   parseargs(argument_hashref, error_msgref)
sub parseargs(\% \$)
{
    my ($argument_hashref, $error_msgref) = @_;
    my  $optcatchresult = 0;
    my  $help;
    my  $man;

    # wrap the signal handler for warnings to copy them to our message space
    local $SIG{__WARN__} = sub
                           {
                               $$error_msgref = "@_";
                               return 1;
                           };

    # set defaults

    # parse @ARGV
    $optcatchresult = GetOptions(
                                 'format=s' => \$argument_hashref->{format},
                                 'auto-fix!'  => \$argument_hashref->{fix},
                                 'verbose!' => \$VERBOSE,
                                 'help'     => \$help,
                                 'man'      => \$man
                                );

    if ($optcatchresult == 0) { return 0 }

    if (defined($help)) { return 2 }

    if (defined($man)) { return 3 }

    # verify options/ arguments
    if (defined $argument_hashref->{format})
    {
        unless (validate_file_format($argument_hashref->{format},
                                     $$error_msgref))
        {
            return 0;
        }
    }

    # check that at least the file is given
    if ($#ARGV < 0)
    {
         $$error_msgref = "At least one name of a source file has to be ".
                          "given. Try \"-help\" or \"-man\" for more ".
                          "information.\n";
        return 0;       
    }

    # fetch file names
    foreach (@ARGV)
    {
        if ( ! -r $_)
        {
            $$error_msgref = "Source file does not exist or is not readable: ".
                             "$_\n";
            return 0;
        }
        push(@{$argument_hashref->{files}}, $_);
    }

    return 1;
}

# FUNCTIONS        - END


# MAIN             - BEGIN
my $ret_val;
my %arg_hash;
my $pod_verbose = 1;
my $current_file;
my $format;
my $filecpy = "";
my @lines;
my $extended = 0;
my %file_structure;
my $error_msg;

# parse commandline
$ret_val = parseargs(%arg_hash, $error_msg);

if ($ret_val == 0)
{
    die("$0: $error_msg");
}
elsif ($ret_val > 1)
{
    if ($ret_val == 3) { $pod_verbose = 2 }
    pod2usage(-exitval => 0, -verbose => $pod_verbose); 
}

# for each file
foreach $current_file (@{$arg_hash{files}})
{
    verbose_msg("Processing file \"${current_file}\"...\n");

    # determin mime type
    if (! defined($arg_hash{format}))
    {
        verbose_msg("  determining mime type... ");
        $format = get_file_format($current_file, $error_msg);

        if ($format eq "") { die("$0: $error_msg") }        
    }
    else
    {
        $format = $arg_hash{format};
        verbose_msg("  forced mime type... ");
    }
    verbose_msg("${Known_Extensions{$format}}\n");    

    # secure copy of file
    if ($arg_hash{fix})
    {
        verbose_msg("  creating backup of file...\n");
        $filecpy = cpy_file($current_file, $error_msg);

        if ($filecpy eq "") { die("$0: $error_msg") }
        verbose_msg("  new backup: \"$filecpy\"\n");
    }
    
    # load file
    verbose_msg("  loading \"${current_file}\"...");
    open(FILE, "<", $current_file) or die("\n$0: Could not open ".
                                          "\"${current_file}\"\n");
    
    @lines = <FILE>;
    close(FILE);
    verbose_msg(" done\n");

    verbose_msg("  checking format...\n");

    set_dwidth($#lines + 1);

    # first determine file structure
    $file_structure{'preamble_start'} = 0;
    $file_structure{'preamble_end'} = 0;
    &{$Define_Preamble{$format}}(\$file_structure{'preamble_start'},
                                 \$file_structure{'preamble_end'},
                                 \@lines);
    verbose_msg("    found preamble from line ".
                "$file_structure{'preamble_start'} to line ".
                "$file_structure{'preamble_end'}\n");

    $file_structure{'code_start'} = $file_structure{'preamble_end'};
    $file_structure{'code_end'} = 0;
    &{$Define_Codesection{$format}}(\$file_structure{'code_start'},
                                    \$file_structure{'code_end'},
                                    \@lines);
    verbose_msg("    found code section from line ".
                "$file_structure{'code_start'} to line ".
                "$file_structure{'code_end'}\n");

    $file_structure{'end_start'} = $file_structure{'code_end'};
    $file_structure{'end_end'} = 0;
    &{$Define_End{$format}}(\$file_structure{'end_start'},
                            \$file_structure{'end_end'},
                            \@lines);
    verbose_msg("    found end section from line $file_structure{'end_start'} ".
                "to line $file_structure{'end_end'}\n");

    # check preamble
    if ($file_structure{'preamble_end'} != 0)
    {
        $extended = &{$Format_Preamble{$format}}(
                                              $file_structure{'preamble_start'},
                                              $file_structure{'preamble_end'},
                                              ${current_file},
                                              \@lines
                                                );
    }

    # secure storing reformated code as file

    verbose_msg("finished processing file \"${current_file}\"\n");
}

# MAIN             - END


__END__

=head1 NAME

reformat - Fix the format (indentation, spacing, ...) of a source file.

=head1 SYNOPSIS

B<reformat> [options] <source file 1> <source file 2> ...

=head1 DESCRIPTION

B<reformat> checks the formatting of a source file to fit the common style of
the B<corb> project. With option B<auto-fix>, all issues are corrected in a
given file. B<THIS MEANS YOUR SOURCE CODE WILL BE AUTOMATICALLY CHANGED>!
Before any changes are applied, a backup copy of the original file is created.
This copy is extended by ".bkp" and a number. If anything goes wrong during
reformatting, just copy the backup to the orignal file name.

- spell check: check spelling, excluded from auto-fix, unknown words have to be
  put into defined dictionary

=head1 OPTIONS

=over 8

=item B<-f FORMAT, --format>

Format has to be one of "m4". Define the format of given source files. Without
this option, the source code format is determined on the file extension.
Supported so far: "m4" for the M4 macro language.

=item B<-a, --auto-fix>

Fix format problems in a source file automatically. With this option, beside
making considerations on the formatting of code, the changes are written to
the file. 

=item B<-v, --verbose>

Enables the verbose mode. In this mode, a little bit of information is
provided as output on each step during reformating.

=item B<-m, --man>

Print a help message plus a description for the program.

=item B<-m, --help>

Print a help message and exit.

=back

=cut


# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
