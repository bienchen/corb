# Last modified: 2010-03-31.14
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

RNA - Everything RNA. Alphabets, base pairs...

=head1 SYNOPSIS

    use RNA qw(:all);

=head1 DESCRIPTION

This module should provide everything you need to work with RNA in Perl. This
includes the standard alphabet, Watson-Crick base pairs, base to number codings
and many more.
 
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
    (L<C<RNAfold()>|"RNAfold">,
     L<C<RNAdistance()>|"RNAdistance">,
     L<C<set_vienna_path()>|"set_vienna_path">)

=item * DB - all functions dealing with DB functionality
    (L<C<set_db_name()>|"set_db_name">,
     L<C<disable_db_chache()>|"disable_db_cache">)

=item * all - everything above

=back

Since direct calling of functions via package naming is also allowed, the
function names are not prefixed with the package name.
(does this hold for vars?)

=cut


package RNA;

use strict;
use warnings;
use Fcntl ':flock';

# PRIVATE packages - BEGIN
use CorbIO qw(:all);
# PRIVATE packages - END

# no idea if this is the recommended way to do it, but since Perl < 5.10 does
# not come with Module::Load::Conditional, it seems to be the simplest
our $DB_ON = 0;
eval { require DBI };
if ( !$@ )
{
    eval { require DBD::SQLite };
    if ( !$@ ) 
    {
        eval { require Digest::SHA1 };        
        if ( !$@ )
        {
            $DB_ON = 1;
        }
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
                    ViennaRNA => [qw(&RNAfold
                                     &RNAdistance
                                     &set_vienna_path)],
                    DB => [qw(&set_db_name
                              &disable_db_cache)]
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

# NON-EXPORTED GLOBALS - BEGIN <- not automatically but still exportable!

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

The list of canonical base pairs including wobble C<GU>. Each entry is a string
of 2 characters. All pairs are written in lower case.

=cut

our @BASEPAIRS = ('cg', 'gc', 'gu', 'ug', 'au', 'ua');

=head2 $N_BPAIRS

Size of list L<C<@BASEPAIRS>|"@BASEPAIRS">. That is the number of canonical
base pairs including wobble C<GU>.

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
my $Private_db_file          = 'corb-calls.sdb';
my $Private_dbh;
my $Private_func_RNAfold     = sub {};
my $Private_func_RNAdistance = sub {};
my %Private_vienna = (
    'RNAfold' => 'RNAfold',
    'RNAdistance' => 'RNAdistance',
    );
# PRIVATE GLOBALS      - END/home/bienert/projects/wetlabdesign.git/bin

# INTERNAL ROUTINES    - BEGIN
# internal routine for starting the database and creating ALL tables
sub s_start_db
{
    if (!defined($Private_dbh))
    {
        my $sth;
        my $ref;
        my %tables;
        my $found = 4;
        my $lock_handle;

        # locking the database file
        $lock_handle = open_or_die($Private_db_file);
        flock($lock_handle, LOCK_EX) or
            msg_error_and_die('Unable to acquire exclusive lock on database ',
                              "file \"$Private_db_file\": $!");

        $Private_dbh = DBI->connect("dbi:SQLite:dbname=$Private_db_file",
                                    '', '' , { PrintError => 0 ,
                                               AutoCommit => 1,
             # we open the DB with an error handler which is inherited by
             # everything derived from $Private_dbh
             HandleError => sub {msg_error_and_die ('Database problem ',
                                                    "(\"$Private_db_file\"): ",
                                                            $DBI::errstr, "\n")}
                                    });

        # check for seq. and struct. tables
        $sth = $Private_dbh->table_info(undef, '%', '%', 'TABLE');
        
        while ($ref = $sth->fetchrow_hashref())
        {
            if (defined($ref->{'TABLE_NAME'}))
            {
                #print "\n|", $ref->{'TABLE_NAME'}, "|";
                if ($ref->{'TABLE_NAME'} eq 'sequence')
                {
                    $tables{'sequence'} = 1;
                    if ($found == 0)
                    {
                        last;
                    }
                    $found--;
                }
                elsif ($ref->{'TABLE_NAME'} eq 'structure')
                {
                    $tables{'structure'} = 1;
                    if ($found == 0)
                    {
                        last;
                    }
                    $found--;
                }
                elsif ($ref->{'TABLE_NAME'} eq 'RNAfold')
                {
                    $tables{'RNAfold'} = 1;
                    if ($found == 0)
                    {
                        last;
                    }
                    $found--;
                }
                elsif ($ref->{'TABLE_NAME'} eq 'RNAdistance')
                {
                    $tables{'RNAdistance'} = 1;
                    if ($found == 0)
                    {
                        last;
                    }
                    $found--;
                }
            }
        }

        # create tables
        if (!defined($tables{'sequence'}))
        {
          $Private_dbh->do('CREATE TABLE sequence ('.
                           'id INTEGER PRIMARY KEY NOT NULL, '.
                           'sha1 CHAR(40) UNIQUE NOT NULL, '.
                           'string BLOB UNIQUE NOT NULL, '.
                           'length INTEGER NOT NULL)');
        }
        if (!defined($tables{'structure'}))
        {
            $Private_dbh->do('CREATE TABLE structure ('.
                             'id INTEGER PRIMARY KEY NOT NULL, '.
                             'sha1 CHAR(40) UNIQUE NOT NULL, '.
                             'string BLOB UNIQUE NOT NULL, '.
                             'length INTEGER NOT NULL)');
        }
        if (!defined($tables{'RNAfold'}))
        {
            $Private_dbh->do('CREATE TABLE RNAfold ('.
                            'id INTEGER PRIMARY KEY NOT NULL, '.
                            'seq_id INTEGER NOT NULL, '.
                            'struct_id INTEGER NOT NULL, '.
                            'parameters TEXT, '.
                            'mfe REAL, '.
                            'eoe REAL, '.
                            'pomfe REAL, '.
                            'FOREIGN KEY(seq_id) REFERENCES sequence(id), '.
                            'FOREIGN KEY(struct_id) REFERENCES structure(id), '.
                             'UNIQUE(seq_id, struct_id, parameters)'.
                            ')');

            $Private_dbh->do('CREATE INDEX RNAfoldindex ON RNAfold(seq_id)');
        }
        if (!defined($tables{'RNAdistance'}))
        {
            $Private_dbh->do('CREATE TABLE RNAdistance ('.
                            'id INTEGER PRIMARY KEY NOT NULL, '.
                            'struct1_id INTEGER NOT NULL, '.
                            'struct2_id INTEGER NOT NULL, '.
                            'parameters TEXT, '.
                            'distance REAL, '.
                           'FOREIGN KEY(struct1_id) REFERENCES structure(id), '.
                           'FOREIGN KEY(struct2_id) REFERENCES structure(id), '.
                            'UNIQUE(struct1_id, struct2_id, parameters)'.
                            ')');

   $Private_dbh->do('CREATE INDEX RNAdistanceindex ON RNAdistance(struct1_id)');
        }

        close($lock_handle);
    }
}
# INTERNAL ROUTINES    - END

=pod

=head1 FUNCTIONS

=head2 set_db_name

Some wrappers in this package use a SQLite database as a permanent storage.
This function can be used to swap to another database file. 

=over 4

=item ARGUMENTS

=over 4

=item file

New filename to be used for the database.

=back

=item EXAMPLE

set_db_name('new_db.sdb');

=back

=cut

sub set_db_name($)
{
    my ($file) = @_;

    if (defined($Private_dbh))
    {
        $Private_dbh->disconnect  or 
            msg_error_and_die ("Could not disconnect database ",
                               "\"$Private_db_file\": $DBI::errstr\n");
    }

    $Private_db_file = $file;

    if ($DB_ON)
    {
        $Private_func_RNAfold     = \&first_db_RNAfold;
        $Private_func_RNAdistance = \&first_db_RNAdistance;
    }
}

=head2 disable_db_cache

Turns off ALL DB functionality for ALL functions of this module. Please note,
that the DB connection can not be re-enabled, once it was closed.

=over 4

=item ARGUMENTS

none.

=item EXAMPLE

disable_db_cache();

=back

=cut

sub disable_db_cache
{
    if ($DB_ON == 1)
    {
        $Private_func_RNAfold     = \&default_RNAfold;
        $Private_func_RNAdistance = \&default_RNAdistance;
    }

    $DB_ON = 0;
}

=head2 set_vienna_path

Set a path to the binaries of the Vienna package.

If all binaries of the Vienna package are in the same place, can be used to
omit individual commands in calls of the wrapper functions.

=over 4

=item ARGUMENTS

=over 4

=item path

The path to the binaries. Trailing backslash may be omitted.

=back

=item EXAMPLE

set_vienna_path(/home/user/vienna);

=back

=cut

sub set_vienna_path($)
{
    my ($path) = @_;

    if (substr($path, -1, 1) ne '/') { $path .='/' }

    foreach (keys(%Private_vienna))
    {
        $Private_vienna{$_} = $path.$Private_vienna{$_};
    }
}

=head2 RNAfold

Folds an RNA sequence into a 2D structure using RNAfold of the Vienna RNA
Package. The return value is a hash with C<structure>, C<mfe>, C<eoe> and 
C<pomfe> as keys, pointing to the structure, its "minimum free energy", the
"free energy of ensemble" and the "frequency of mfe structure in ensemble".

Here is a list of abbreviations, used as keys:

=over 4

=item * mfe = "Minimum free energy"

=item * eoe = "Energy of ensemble"

=item * pomfe = "Probability of mfe structure"

=back

Since folding can take some time, we provide a databse mechanism, which is able
to store all your results for reuse. The idea is to store each call to RNAfold
with its parameters and on recurring calls just send back the answer from the
database instead of calculating it again. An entry in the database contains the
sequence, a unified parameterstring, and the results of the call. 

The database will be stored in a single file defaulting to "corb-calls.sdb".
The database system we use is SQLite, therefore it should be shareable among
different operating systems.

Recording may be turned of by calling
L<C<disable_db_cache()>|"disable_db_cache"> but cannot be re-enabled.

Please note that we do not check for the binary used to produce results.

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

my %prediction = RNAfold('AUGCAUGC', '-noPS -d2');

=back

=cut

sub default_RNAfold
{
    my ($seq, $params, $cmd) = @_;
    my $len = length($seq);
    my @err_run;
    my %result;
    
    if (!defined($cmd))
    {
        $cmd = $Private_vienna{'RNAfold'};
    }
    
    if (!defined($params))
    {
        $params = '';
    }
    elsif ($params =~ /(\;)/)
    {
        msg_error_and_die('Parameterstring for RNAfold contains ',
                          "malicious characters: $1");
    }
    
    # start folding
    unless(open(FH, "echo \"$seq\" | $cmd $params 2>&1 |"))
    {
        msg_error_and_die ("Could not start $cmd $params: $!\n");
    }

    foreach (<FH>)
    {
        #print $_;
        if ($_ =~ /([\(\.\)]{$len})\s+\(\s*(\-?\d+\.\d+)\)/)
        {
            $result{structure} = $1;
            $result{mfe} = $2;
        }
        elsif ($_ =~ /free\s+energy\s+of\s+ensemble\s+=\s*(\-?\d+\.\d+)\s*kcal/)
        {
            $result{eoe} = $1;
        }
        elsif ($_ =~ /frequency\s+of\s+mfe\s+structure\s+in\s+ensemble\s+(\-?\d+\.\d+(?:e\-?\d+)?);/)
        {
            $result{pomfe} = $1;
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

sub db_RNAfold
{
    my ($seq, $params, $cmd) = @_;
    my $seq_len = length($seq);                  # sequence length
    my %result;                                  # container for folding result
    my $seq_id;                                  # id from table sequence
    my $struct_id;                               # id from table structure
    my $cparams = '';                            # condensed parameter string
    my $seq_sha1 = Digest::SHA1::sha1_hex($seq); # sha1 hash of sequence
    my $ref;                                     # fetch all kinds of refs

    if (defined($params))
    {
        $cparams = $params;
        # we have a non-empty param.string, sort it
        if ($cparams =~ /^(.*\-d)\s*(\d.*)$/) #                 -d[0|1|2|3]
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-p)\s*(\d.*)$/) #                 -p[0|2]
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-T)\s*(-?\d+(?:\.\d+)?.*)$/) #    -T temp 
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-e)\s*(\d.*)$/) #                 -e 1|2
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-P)\s*(.*)$/) #                   -P paramfile
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-nsp)\s*(.*)$/) #                 -nsp pairs
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }
        if ($cparams =~ /^(.*\-S)\s*(-?\d+(?:\.\d+)?.*)$/) #    -S scale 
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }

        # split and sort string
        my @parry = sort(split(/ /, $cparams));
        $cparams = join('', @parry);
    }

    # db is running, check if experiment is already stored
    my $sth = $Private_dbh->prepare(q{
                                 SELECT id FROM sequence WHERE sequence.sha1=?;
                                     });
    $sth->execute($seq_sha1);

    if ($ref = $sth->fetchrow_hashref())
    {
        $seq_id = $ref->{'id'};

        $sth = $Private_dbh->prepare(q{
        SELECT struct_id, mfe, eoe, pomfe FROM RNAfold WHERE seq_id=? AND parameters=?;
                                      });
        $sth->execute($seq_id, $cparams);

        if ($ref = $sth->fetchrow_hashref())
        {
            $struct_id = $ref->{'struct_id'};
            $result{mfe} = $ref->{'mfe'};
            $result{eoe} = $ref->{'eoe'};
            $result{pomfe} = $ref->{'pomfe'};
        }
    }
    else
    {
      $sth = $Private_dbh->prepare(
           "INSERT INTO sequence(sha1,string,length) VALUES (?,?,?)");
      $sth->execute($seq_sha1, $seq, $seq_len);
    }

    if (defined($struct_id))
    {
        $sth = $Private_dbh->prepare(q{
                  SELECT string FROM structure WHERE id=?;
                                      });
        $sth->execute($struct_id);

        if ($ref = $sth->fetchrow_hashref())
        {
            $result{structure} = $ref->{'string'};
        }
        else
        {
            msg_error_and_die ('Database problem ', "(\"$Private_db_file\"): ",
                              "Entry for struct with id $struct_id not found!");
        }

        return %result; #default_RNAfold($seq, $params, $cmd);
    }
    else # store experiment
    {
        # get seq.id: We fetch it from the table, again, because of possible
        # changes of the db since the first SELECT statement.
        $sth = $Private_dbh->prepare(q{
                                 SELECT id FROM sequence WHERE sequence.sha1=?;
                                     });
        $sth->execute($seq_sha1);

        if ($ref = $sth->fetchrow_hashref())
        {
            $seq_id = $ref->{'id'};
        }
        else
        {
            msg_error_and_die ('Database problem ', "(\"$Private_db_file\"): ",
                               "Entry for sequence \"$seq\" not found!");
        }

        # check whether we have to store structure
        %result = default_RNAfold($seq, $params, $cmd);
        my $struct_sha1 = Digest::SHA1::sha1_hex($result{structure});
        $sth = $Private_dbh->prepare(q{
                                SELECT id FROM structure WHERE structure.sha1=?;
                                     });
        $sth->execute($struct_sha1);
        if ($ref = $sth->fetchrow_hashref())
        {
            $struct_id = $ref->{'id'};
        }
        else
        {
            $sth = $Private_dbh->prepare(
                "INSERT INTO structure(sha1,string,length) VALUES (?,?,?)");
            $sth->execute($struct_sha1,
                          $result{structure},
                          length($result{structure}));

            $sth = $Private_dbh->prepare(q{
                                SELECT id FROM structure WHERE structure.sha1=?;
                                     });
            $sth->execute($struct_sha1);
            if ($ref = $sth->fetchrow_hashref())
            {
                $struct_id = $ref->{'id'};
            }
            else
            {
                msg_error_and_die ('Database problem ',
                                   "(\"$Private_db_file\"): Entry for ",
                                   "structure \"$result{structure}\" not ",
                                   'found!');
            }
        }
        #print "ENTRY: $seq_id, $struct_id\n";

        # store RNAfold experiment
        $sth = $Private_dbh->prepare(
            "INSERT INTO RNAfold(seq_id, struct_id, parameters, mfe, eoe, pomfe) "
            ."VALUES (?,?,?,?,?,?)");
        $sth->execute($seq_id, $struct_id, $cparams, $result{mfe},
                      $result{eoe},$result{pomfe});
    }

    return %result;
}

sub first_db_RNAfold
{
    $Private_func_RNAfold = \&db_RNAfold;
    s_start_db();
    return db_RNAfold(@_);
}

if (!$DB_ON)
{
    msg_warning ("Database cache offline.\n");
    $Private_func_RNAfold = \&default_RNAfold;
}
else
{
    $Private_func_RNAfold = \&first_db_RNAfold;
}

sub RNAfold
{
    return &$Private_func_RNAfold(@_);
}

=head2 RNAdistance

Compare two RNA secondary structures using RNAdistance of the Vienna RNA
Package. The return value is a hash with the base pair distance, pointed to by
C<DP> as key. If you need anything else than the base pair distance, feel free
to extend this function.

Since comparing structures can take some time, we provide a databse mechanism,
which is able to store all your results for reuse. The idea is to store each
call to RNAdistance with its parameters and on recurring calls just send back
the answer from the database instead of calculating it again. An entry in the
database contains the two structures, a unified parameterstring, and the
results of the call. 

The database will be stored in a single file defaulting to "corb-calls.sdb".
The database system we use is SQLite, therefore it should be shareable among
different operating systems.

Recording may be turned of by calling
L<C<disable_db_cache()>|"disable_db_cache"> but cannot be re-enabled.

Please note that we do not check for the binary used to produce results.

=over 4

=item ARGUMENTS

=over 4

=item structure1

First of the two structures.

=item structure2

Second of the two structures.

=item parameters

Parameters to be passed to RNAdistance. All options have to be passed in a
single string in the same way they would be used in a command line call to
RNAdistance. May be omitted if not used.

=item command

As default, this functions just sends "RNAdistance" as command to a shell. This
option can be used to change the call to whatever you want. Just omit it to get
the default behaviour.

=back

=item EXAMPLE

my %comparison = RNAdistance('(((...)))', '(((......)))', '-DP');

=back

=cut

sub default_RNAdistance
{
    my ($struct1, $struct2, $params, $cmd) = @_;
    #my $len1 = length($struct1);
    #my $len2 = length($struct2);
    my @err_run;
    my %result;
    
    if (!defined($cmd))
    {
        $cmd = $Private_vienna{'RNAdistance'};
    }
    
    if (!defined($params))
    {
        $params = '';
    }
    elsif ($params =~ /(\;)/)
    {
        msg_error_and_die('Parameterstring for RNAdistance contains ',
                          "malicious characters: $1");
    }

    # start comparison
    unless(open(FH, "echo \"$struct1\n$struct2\" | $cmd $params 2>&1 |"))
    {
        msg_error_and_die ("Could not start $cmd $params: $!\n");
    }

    foreach (<FH>)
    {
        if ($_ =~ /P\:\s*(\d+)/)
        {
            $result{DP} = $1;
            next;
        }
        
        push(@err_run, $_);
    }
    
    close(FH);
    
    if (! defined ($result{DP}))
    {
        msg_error_and_die("Running RNAdistance failed, output of\n"
                 ."\`echo \"$struct1\n$struct2\" | $cmd $params\`: ",
                          @err_run);
    }

    return %result;
}

sub db_RNAdistance
{
    my ($struct1, $struct2, $params, $cmd) = @_;
    my $len1 = length($struct1);                         # structure size
    my $len2 = length($struct2);
    my %result;                                          # result container
    my $struct1_id;                                      # id of table structure
    my $struct2_id;
    my $cparams = '';                                    # condensed parameters
    my $struct1_sha1 = Digest::SHA1::sha1_hex($struct1); # sha1 key of structure
    my $struct2_sha1 = Digest::SHA1::sha1_hex($struct2);
    my $ref;                                             # all kinds of refs
    my $sth_s;                                           # handle for SELECTs
    my $sth_i;                                           # handle for INSERTs

    if (defined($params))
    {
        $cparams = $params;
        # we have a non-empty param.string, sort it
        if ($cparams =~ /^(.*\-D)\s*([fhwcFHWCP\s]*)(.*)$/) #      -D[fhwcFHWCP]
        {
            $cparams = substr($cparams, 0, length($1));
            if (defined($2))
            {
                $cparams .= join('', sort(split(//,$2)));
            }
            $cparams .= $3;
        }
        if ($cparams =~ /^(.*\-X)\s*([pmfc\s]*)(.*)$/)      #      -X[p|m|f|c]
        {
            $cparams = substr($cparams, 0, length($1));
            if (defined($2))
            {
                $cparams .= join('', sort(split(//,$2)));
            }
            $cparams .= $3;
        }
        if ($cparams =~ /^(.*\-B)\s*([^\s]+.*)$/)         #      -B <file>
        {
            $cparams = substr($cparams, 0, length($1));
            $cparams .= $2;
        }

        # split and sort string
        my @parry = sort(split(/ /, $cparams));
        $cparams = join('', @parry);
    }

    # check if experiment is already stored
    # fetch id of first struct
    $sth_s = $Private_dbh->prepare(q{
                                SELECT id FROM structure WHERE structure.sha1=?;
                                    });
    $sth_s->execute($struct1_sha1);

    if ($ref = $sth_s->fetchrow_hashref())
    {
        $struct1_id = $ref->{'id'};
    }
    else
    {
        $sth_i = $Private_dbh->prepare(q{
                       INSERT INTO structure(sha1,string,length) VALUES (?,?,?);
                                        });
        $sth_i->execute($struct1_sha1, $struct1, $len1);

        $sth_s->execute($struct1_sha1);

        if ($ref = $sth_s->fetchrow_hashref())
        {
            $struct1_id = $ref->{'id'};
        }
        else
        {
            msg_error_and_die ('Database problem ', "(\"$Private_db_file\"): ",
                               "Stored structure \"$struct1\" not found!");
        }
    }

    # fetch id of 2nd struct
    $sth_s->execute($struct2_sha1);

    if ($ref = $sth_s->fetchrow_hashref())
    {
        $struct2_id = $ref->{'id'};
    }
    else
    {
        if (!defined($sth_i))
        {
            $sth_i = $Private_dbh->prepare(q{
                       INSERT INTO structure(sha1,string,length) VALUES (?,?,?);
                                            });
        }

        $sth_i->execute($struct2_sha1, $struct2, $len2);

        $sth_s->execute($struct2_sha1);

        if ($ref = $sth_s->fetchrow_hashref())
        {
            $struct2_id = $ref->{'id'};
        }
        else
        {
            msg_error_and_die ('Database problem ', "(\"$Private_db_file\"): ",
                               "Stored structure \"$struct2\" not found!");
        }
    }

    # now we have ids for both structures, try to fetch experiment
    $sth_s = $Private_dbh->prepare(
                  'SELECT distance FROM RNAdistance WHERE struct1_id=? '.
                                                     'AND struct2_id=? '.
                                                     'AND parameters=?;');

    $sth_s->execute($struct1_id, $struct2_id, $cparams);

    if ($ref = $sth_s->fetchrow_hashref())
    {
        $result{DP} = $ref->{'distance'};
        return %result; #default_RNAdistance(@_);
    }

    #print "STRUCT1: $struct1_id STRUCT2: $struct2_id\n";

    %result = default_RNAdistance(@_);

    $sth_i = $Private_dbh->prepare('INSERT INTO RNAdistance(struct1_id, '.
                                                           'struct2_id, '.
                                                           'parameters, '.
                                                           'distance) '.
                                   'VALUES(?,?,?,?)');
    $sth_i->execute($struct1_id, $struct2_id, $cparams, $result{DP});

    return %result;
}

sub first_db_RNAdistance
{
    $Private_func_RNAdistance = \&db_RNAdistance;
    s_start_db();
    return db_RNAdistance(@_);
}

if (!$DB_ON)
{
    $Private_func_RNAdistance = \&default_RNAdistance;
}
else
{
    $Private_func_RNAdistance = \&first_db_RNAdistance;
}

sub RNAdistance
{
    return &$Private_func_RNAdistance(@_);
}

1;

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
