package Chemistry::File::SMILES;

$VERSION = "0.33";
# $Id: SMILES.pm,v 1.7 2004/06/11 00:10:45 itubert Exp $

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol;
use Carp;


=head1 NAME

Chemistry::File::SMILES - SMILES linear notation parser/writer

=head1 SYNOPSYS

    #!/usr/bin/perl
    use Chemistry::File::SMILES;

    # parse a SMILES string
    my $s = 'C1CC1(=O)[O-]';
    my $mol = Chemistry::Mol->parse($s, format => 'smiles');

    # print a SMILES string
    print $mol->print(format => 'smiles');

    # parse a SMILES file
    my @mols = Chemistry::Mol->read("file.smi", format => 'smiles');

    # write a multiline SMILES file
    Chemistry::Mol->write("file.smi", mols => [@mols]);


=head1 DESCRIPTION

This module parses a SMILES (Simplified Molecular Input Line Entry
Specification) string. This is a File I/O driver for the PerlMol project.
L<http://www.perlmol.org/>. It registers the 'smiles' format with
Chemistry::Mol.

This parser interprets anything after whitespace as the molecule's name;
for example, when the following SMILES string is parsed, $mol->name will be
set to "Methyl chloride":

    CCl	 Methyl chloride

The name is not included by default on output. However, if the C<name> option
is defined, the name will be included after the SMILES string, separated by a
tab.

    print $mol->print(format => 'smiles', name => 1);

=head2 Multiline SMILES and SMILES files

A file or string can contain multiple molecules, one per line.

    CCl	 Methyl chloride
    CO	 Methanol

Files with the extension '.smi' are assumed to have this format.

=cut

# INITIALIZATION
Chemistry::Mol->register_format('smiles');
my $Smiles_parser = __PACKAGE__->new;

#=begin comment
#
#=over
#
#=cut

sub file_is {
    my $self = shift;
    $self->name_is(@_);
}

sub name_is {
    my ($self, $name) = @_;
    $name =~ /\.smi/;
}

sub parse_string {
    my ($self, $string, %opts) = @_;
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || "Chemistry::Atom";
    my $bond_class = $opts{bond_class} || "Chemistry::Bond";

    my (@lines) = split /(?:\n|\r\n?)/, $string;
    my @mols;
    for my $line (@lines) {
        my ($smiles, $name) = split " ", $line, 2;
        my $mol = $mol_class->new;
        $Smiles_parser->parse($smiles, $mol);
        $mol->name($name);
        return $mol unless wantarray;
        push @mols, $mol;
    }
    return @mols;
}

### The contents of the original Chemistry::Smiles module start below

my $Symbol = qr/
    s|p|o|n|c|b|Zr|Zn|Yb|Y|Xe|W|V|U|Tm|Tl|Ti|Th|
    Te|Tc|Tb|Ta|Sr|Sn|Sm|Si|Sg|Se|Sc|Sb|S|Ru|Rn|Rh|Rf|Re|Rb|Ra|
    Pu|Pt|Pr|Po|Pm|Pd|Pb|Pa|P|Os|O|Np|No|Ni|Ne|NdNb|Na|N|Mt|Mt|
    Mo|Mn|Mg|Md|Lu|Lr|Li|La|Kr|K|Ir|In|I|Hs|Hs|Ho|Hg|Hf|He|H|Ge|
    Gd|Ga|Fr|Fm|Fe|F|Eu|Es|Er|Dy|Ds|Db|Cu|Cs|Cr|Co|Cm|Cl|Cf|Ce|
    Cd|Ca|C|Br|Bk|BiBh|Be|Ba|B|Au|At|As|Ar|Am|Al|Ag|Ac|\*
/x; # Order is reverse alphabetical to ensure longest match

my $Simple_symbol = qr/Br|Cl|B|C|N|O|P|S|F|I|s|p|o|n|c|b/;

my $Bond = qr/(?:[-=#:.\/\\])?/; 
my $Simple_atom = qr/($Simple_symbol)/;   #3
my $Complex_atom = qr/
    (?:
        \[                          #begin atom
        (\d*)                       #4 isotope
        ($Symbol)                   #5 symbol
        (\@{0,2})                   #6 chirality
        (?:(H\d*))?                 #7 H-count
        (\+{2,}|-{2,}|\+\d*|-\d*)?  #8 charge
        \]                          #end atom 
    )
/x;

my $Digits = qr/(?:($Bond)(?:\d|%\d\d))*/; 
my $Chain = qr/
    \G(                                     #1
        (?: 
            ($Bond)                         #2
            (?:$Simple_atom|$Complex_atom)  #3-8
            ($Digits)                       #9
        ) 
        |\( 
        |\)
        |.+
    )
/x;

my $digits_re = qr/($Bond)(\%\d\d|\d)/;

my %type_to_order = (
    '-' => 1,
    '=' => 2,
    '#' => 3,
    '/' => 1,
    '\\' => 1,
    '' => 1, # not strictly true
    '.' => 0,
);

#=item Chemistry::Smiles->new([add_atom => \&sub1, add_bond => \&sub2])
#
#Create a SMILES parser. If the add_atom and add_bond subroutine references
#are given, they will be called whenever an atom or a bond needs to be added
#to the molecule. If they are not specified, default methods, which
#create a Chemistry::Mol object, will be used.
#
#=cut

sub new {
    my $class = shift;
    my %opts = @_;
    require Chemistry::Mol unless $opts{add_atom} && $opts{add_bond};
    my $self = bless {
        add_atom => $opts{add_atom} || \&add_atom,
        add_bond => $opts{add_bond} || \&add_bond,
    }, $class;
}

#=item $obj->parse($string, $mol)
#
#Parse a Smiles $string. $mol is a "molecule state object". It can be anything;
#the parser doesn't do anything with it except sending it as the first parameter
#to the callback functions. If callback functions were not provided when
#constructing the parser object, $mol must be a Chemistry::Mol object, because
#that's what the default callback functions require.
#
#=cut

sub parse {
    my $self = shift;
    my ($s, $mol) = @_;
    $self->{stack} = [ undef ];
    $self->{digits} = {};

    while ($s =~ /$Chain/g) {
        #my @a = ($1, $2, $3, $4, $5, $6, $7, $8);
        #print Dumper(\@a);
        my ($all, $bnd, $sym, $iso, $sym2, $chir, $hcnt, $chg, $dig) 
            = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
        if ($all eq '(') {
            $self->start_branch();
        } elsif ($all eq ')') {
            $self->end_branch();
        } elsif ($sym) { # Simple atom
            no warnings;
            my @digs = parse_digits($dig);
            $self->atom($mol, $bnd, '', $sym, '', '', '', \@digs);
        } elsif ($sym2) { # Complex atom
            no warnings;
            my @digs = parse_digits($dig);
            if ($hcnt eq 'H') { 
                $hcnt = 1;
            } else {
                $hcnt =~ s/H//;
            }
            unless ($chg =~ /\d/) {
                $chg = ($chg =~ /-/) ? -length($chg) : length($chg);
            }
            $self->atom($mol, $bnd, $iso, $sym2, $chir, $hcnt || 0, $chg, \@digs);
        } else {
            croak "SMILES ERROR: '$all'\n";
        }
    }
    $mol;
}

sub parse_digits {
    my ($dig) = @_;
    my @digs;
    while ($dig && $dig =~ /$digits_re/g) {
        push @digs, {bnd=>$1, dig=>$2};
    }
    @digs;
}

sub atom {
    my $self = shift;
    my ($mol,$bnd,$iso,$sym,$chir,$hcount,$chg,$digs) = @_;
    #{no warnings; local $" = ','; print "atom(@_)\n"}
    my $a = $self->{add_atom}($mol,$iso,$sym,$chir,$hcount,$chg);
    if($self->{stack}[-1]) {
        $self->{add_bond}($mol, $bnd, $self->{stack}[-1], $a);
    }
    for my $dig (@$digs) {
        if ($self->{digits}{$dig->{dig}}) {
            if ($dig->{bnd} && $self->{digits}{$dig->{dig}}{bnd}
                &&  $dig->{bnd} ne $self->{digits}{$dig->{dig}}{bnd}){
                die "SMILES: Inconsistent ring closure\n";
            }
            $self->{add_bond}($mol, 
                $dig->{bnd} || $self->{digits}{$dig->{dig}}{bnd}, 
                $self->{digits}{$dig->{dig}}{atom}, $a);
            delete $self->{digits}{$dig->{dig}};
        } else {
            $self->{digits}{$dig->{dig}} = {atom=>$a, bnd=>$dig->{bnd}};
        }
    }
    $self->{stack}[-1] = $a;
}

#=back
#
#=head1 CALLBACK FUNCTIONS
#
#=over
#
#=item $atom = add_atom($mol, $iso, $sym, $chir, $hcount, $chg)
#
#Called by the parser whenever an atom is found. The first parameter is the
#state object given to $obj->parse(). The other parameters are the isotope,
#symbol, chirality, hydrogen count, and charge of the atom. Only the symbol is
#guaranteed to be defined. Mnemonic: the parameters are given in the same order
#that is used in a SMILES string (such as [18OH-]). This callback is expected to
#return something that uniquely identifies the atom that was created (it might
#be a number, a string, or an object).
#
#=cut

# Default add_atom callback 
sub add_atom {
    my ($mol, $iso, $sym, $chir, $hcount, $chg) = @_;
    my $atom = $mol->new_atom(symbol=>$sym);
    $iso && $atom->attr('smiles/isotope' => $iso);
    $chir && $atom->attr('smiles/chirality' => $chir);
    length $hcount && $atom->attr('smiles/h_count' => $hcount);
    $chg && $atom->attr('smiles/charge' => $chg);
    $atom;
}

#=item add_bond($mol, $type, $a1, $a2)
#
#Called by the parser whenever an bond needs to be created. The first parameter
#is the state object given to $obj->parse(). The other parameters are the bond
#type and the two atoms that need to be bonded. The atoms are identified using
#the return values from the add_atom() callback.
#
#=back
#
#=end comment
#
#=cut

# Default add_bond callback 
sub add_bond {
    my ($mol, $type, $a1, $a2) = @_;
    my $order = $type_to_order{$type} or return; # don't add bonds of order 0
    my $bond = $mol->new_bond(type=>$type, atoms=>[$a1, $a2], order=>$order);
    $bond->attr("smiles/type" => $type);
    $bond;
}

sub start_branch {
    my $self = shift;
    #print "start_branch\n";
    push @{$self->{stack}}, $self->{stack}[-1];
}

sub end_branch {
    my $self = shift;
    #print "end_branch\n";
    pop @{$self->{stack}};
}

##### SMILES WRITER ########

my %ORDER_TO_TYPE = (
    2 => '=', 1 => '', 3 => '#',
);

my %ORGANIC_ELEMS = (
    Br => 1, Cl => 1, B => 1, C => 1, N => 1, O => 1, P => 1, S => 1, 
    F => 1, I => 1, s => 1, p => 1, o => 1, n => 1, c => 1, b => 1,
);

sub write_string {
    my ($self, $mol_ref, %opts) = @_;

    my $eol;
    my @mols;
    if ($opts{mols}) {
        @mols = @{$opts{mols}};
        $eol = "\n";
    } else {
        @mols = $mol_ref; 
        $eol = "";
    }

    my $smiles;
    for my $mol (@mols) {
        $mol = $mol->clone; 
        collapse_hydrogens($mol);

        my $atom = $mol->atoms(1);
        my $ring_atoms = {};
        find_ring_bonds($mol, $atom, undef, {}, $ring_atoms);
        $smiles .= branch($mol, $atom, undef, {}, $ring_atoms);
        if ($opts{name}) {
            $smiles .= "\t" . $mol->name;
        }
        $smiles .= $eol;
    }
    return $smiles;
}

sub find_ring_bonds {
    my ($mol, $atom, $from_bond, $visited, $ring_atoms) = @_;

    $visited->{$atom}  = 1;
    for my $bn ($atom->bonds_neighbors) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond}  = 1;
        if ($visited->{$nei}) { # closed ring
            #print "closing ring\n";
            $ring_atoms->{$nei}++;
        } else {
            find_ring_bonds($mol, $nei, $bond, $visited, $ring_atoms);
        }
    }
}

sub branch {
    my ($mol, $atom, $from_bond, $visited, $digits) = @_;

    my $prev_branch = "";
    my $smiles;
    $smiles .= $ORDER_TO_TYPE{$from_bond->order} if $from_bond;
    $digits->{count}++;
    if ($ORGANIC_ELEMS{$atom->symbol}) {
        $smiles .= $atom->symbol;
    } else {
        my $h_count = $atom->attr("smiles/h_count");
        $h_count = $h_count ? ($h_count > 1 ? "H$h_count" : 'H') : '';
        $smiles .= "[" . $atom->symbol . "$h_count]";
    }
    if ($digits->{$atom}) {  # opening a ring
        my @d;
        for (1 .. $digits->{$atom}) {
            push @d, next_digit($digits);
        }
        $digits->{$atom} = \@d;
        $smiles .= join "", map { $_ < 10 ? $_ : "%$_"} @d;
    }

    $visited->{$atom}  = 1;
    for my $bn ($atom->bonds_neighbors) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond} = 1;
        if ($visited->{$nei}) { # closed a ring
            my $digit = shift @{$digits->{$nei}};
            $smiles .= $ORDER_TO_TYPE{$bond->order};
            $smiles .= $digit < 10 ? $digit : "%$digit";
            $digits->{used_digits}[$digit] = 0; # free for future use
            $visited->{$bond} = 1;
        } else {
            my $branch = branch($mol, $nei, $bond, $visited, $digits);
            if ($prev_branch) {
                $smiles .= "($prev_branch)";
            }
            $prev_branch = $branch;
        }
    }
    $smiles .= "$prev_branch";
    $smiles;
}

sub next_digit {
    my ($digits) = @_;
    for (my $i = 1; $i < 100; $i++) {
        unless ($digits->{used_digits}[$i]) {
            $digits->{used_digits}[$i] = 1;  # mark as used
            return $i;
        }
    }
    die "no more available smiles digits!";  # shouldn't happen
}

sub collapse_hydrogens {
    my ($mol) = @_;

    for my $atom (grep {$_->symbol eq 'H'} $mol->atoms) {
        my ($neighbor) = $atom->neighbors;
        $atom->delete;
        my $h_count = $neighbor->attr("smiles/h_count");
        $h_count++;
        $neighbor->attr("smiles/h_count", $h_count);
    }
}


1;

=head1 CAVEATS

Branches that start before an atom, such as (OC)C, which should be equivalent
to C(CO) and COC, according to some variants of the SMILES specification. Many
other tools don't implement this rule either.

Not yet available: Proper handling of aromatic atoms;  unique (canonical)
SMILES output.

=head1 VERSION

0.33

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>

The SMILES Home Page at http://www.daylight.com/dayhtml/smiles/
The Daylight Theory Manual at 
http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

