package Chemistry::File::SMILES;

$VERSION = "0.40";
# $Id: SMILES.pm,v 1.8 2004/06/16 19:43:04 itubert Exp $

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol;
use Chemistry::Bond::Find 'assign_bond_orders';
use List::Util 'first';
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

=head1 OPTIONS

=over

=item aromatic

On output, detect aromatic atoms and bonds by means of the Chemistry::Ring
module, and represent the organic aromatic atoms with lowercase symbols.

=item unique

When used on output, canonicalize the structure if it hasn't been canonicalized
already and generate a unique SMILES string. This option implies "aromatic".

=back

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
    %opts = (kekulize => 1, %opts);
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";

    # these two are not used yet...
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;
    my $bond_class = $opts{bond_class} || $mol_class->bond_class;

    my (@lines) = split /(?:\n|\r\n?)/, $string;
    my @mols;
    for my $line (@lines) {
        my ($smiles, $name) = split " ", $line, 2;
        my $mol = $mol_class->new;
        $Smiles_parser->parse($smiles, $mol);
        $mol->name($name);
        add_implicit_hydrogens($mol);
        if ($opts{kekulize}) {
            assign_bond_orders($mol, method => "itub", use_coords => 0, 
                scratch => 0, charges => 0);
        }
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

my %ORGANIC_ELEMS = (
    Br => 1, Cl => 1, B => 3, C => 4, N => 3, O => 2, P => 3, S => 2, 
    F => 1, I => 1, s => 1, p => 1, o => 1, n => 1, c => 1, b => 1,
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
    my ($self, $s, $mol) = @_;
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
            $self->atom($mol, $bnd, '', $sym, '', undef, '', \@digs);
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
    my $atom = $mol->new_atom(symbol => ucfirst $sym);
    $iso && $atom->attr('smiles/isotope' => $iso);
    $iso && $atom->mass($iso);
    $chir && $atom->attr('smiles/chirality' => $chir);
    defined $hcount && $atom->hydrogens($hcount);
    $chg && $atom->formal_charge($chg);
    if ($sym =~ /^[a-z]/) {
        $atom->attr("smiles/aromatic", 1);
    }
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
 
sub add_implicit_hydrogens {
    my ($mol) = @_;
    for my $atom ($mol->atoms) {
        #print "H=".$atom->hydrogens."\n";
        unless (defined $atom->hydrogens) {
            my $h_count = $ORGANIC_ELEMS{$atom->symbol} - $atom->valence;
            $h_count = 0 if $h_count < 0;
            if ($atom->attr("smiles/aromatic") and $atom->symbol =~ /^[CN]$/) {
                $h_count--;
            }
            $atom->hydrogens($h_count);
        }
    }
}

##### SMILES WRITER ########

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
        my @atoms = $mol->atoms; 

        if ($opts{unique}) {
            unless ($atoms[0]->attr("canon/class")) {
                require Chemistry::Canonicalize;
                Chemistry::Canonicalize::canonicalize($mol);
            }
            $opts{aromatic} = 1; # all unique smiles have to be aromatic
            @atoms = sort {
                $a->attr("canon/class") <=> $b->attr("canon/class")
            } @atoms;
        }

        aromatize($mol) if $opts{aromatic};

        my $visited = {};
        my @s;
        for my $atom (@atoms) {
            next if $visited->{$atom};
            my $ring_atoms = {};

            # first pass to find and number the ring bonds
            find_ring_bonds($mol, \%opts, $atom, undef, {}, $ring_atoms);

            # second pass to actually generate the SMILES string
            push @s, branch($mol, \%opts, $atom, undef, $visited, $ring_atoms);
        }
        $smiles .= join '.', @s;

        if ($opts{name}) {
            $smiles .= "\t" . $mol->name;
        }
        $smiles .= $eol;
    }
    return $smiles;
}

sub find_ring_bonds {
    my ($mol, $opts, $atom, $from_bond, $visited, $ring_atoms) = @_;

    $visited->{$atom}  = 1;
    for my $bn (sorted_bonds_neighbors($atom, $opts)) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond}  = 1;
        if ($visited->{$nei}) { # closed ring
            #print "closing ring\n";
            $ring_atoms->{$nei}++;
        } else {
            find_ring_bonds($mol, $opts, $nei, $bond, $visited, $ring_atoms);
        }
    }
}

sub branch {
    my ($mol, $opts, $atom, $from_bond, $visited, $digits) = @_;

    my $prev_branch = "";
    my $smiles;
    $smiles .= bond_symbol($from_bond, $opts);
    $digits->{count}++;
    $smiles .= format_atom($atom, $opts);
    if ($digits->{$atom}) {  # opening a ring
        my @d;
        for (1 .. $digits->{$atom}) {
            push @d, next_digit($digits);
        }
        $digits->{$atom} = \@d;
        $smiles .= join "", map { $_ < 10 ? $_ : "%$_"} @d;
    }

    $visited->{$atom}  = 1;
    for my $bn (sorted_bonds_neighbors($atom, $opts)) {
        my $nei  = $bn->{to};
        my $bond = $bn->{bond};
        next if $visited->{$bond};
        $visited->{$bond} = 1;
        if ($visited->{$nei}) { # closed a ring
            my $digit = shift @{$digits->{$nei}};
            $smiles .= bond_symbol($bond, $opts);
            $smiles .= $digit < 10 ? $digit : "%$digit";
            $digits->{used_digits}[$digit] = 0; # free for future use
            $visited->{$bond} = 1;
        } else {
            my $branch = branch($mol, $opts, $nei, $bond, $visited, $digits);
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

sub aromatize {
    my ($mol) = @_;

    require Chemistry::Ring::Find;

    for my $atom ($mol->atoms) {
        my @rings = Chemistry::Ring::Find::find_ring(
            $atom, all => 1, min => 5, max => 7);
        for my $ring (@rings) {
            if ($ring->is_aromatic) {
                for my $ring_elem ($ring->atoms, $ring->bonds) {
                    $ring_elem->aromatic(1);
                }
            }
        }
    }
}

sub sorted_bonds_neighbors {
    my ($atom, $opts) = @_;
    my @bn = $atom->bonds_neighbors;
    if ($opts->{unique}) {
        @bn = sort { 
            $a->{to}->attr("canon/class") <=> $b->{to}->attr("canon/class") 
        } @bn;
    }
    @bn;
}

my %ORDER_TO_TYPE = (
    2 => '=', 1 => '', 3 => '#',
);

sub bond_symbol {
    my ($bond, $opts) = @_;
    return '' unless $bond;
    return '' if $opts->{aromatic} && $bond->aromatic;
    return $ORDER_TO_TYPE{$bond->order};
}

sub format_atom {
    my ($atom, $opts) = @_;
    my $s;
    if ($ORGANIC_ELEMS{$atom->symbol} && !$atom->formal_charge) {
        $s = $atom->symbol;
        if ($opts->{aromatic} && $atom->aromatic) {
            $s = lc $s;
        }
    } else {
        my $h_count = $atom->hydrogens;
        my $charge  = $atom->formal_charge;
        my $symbol  = $atom->symbol;
        my $iso     = $atom->attr("smiles/isotope") || '';

        $charge  = $charge ? sprintf("%+d", $charge): '';
        $h_count = $h_count ? ($h_count > 1 ? "H$h_count" : 'H') : '';

        $s = "[$iso$symbol$h_count$charge]";
    }
    $s;
}


1;

=head1 CAVEATS

Reading branches that start before an atom, such as (OC)C, which should be
equivalent to C(OC) and COC, according to some variants of the SMILES
specification. Many other tools don't implement this rule either.

=head1 VERSION

0.40

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

