package Chemistry::File::SMILES;

$VERSION = "0.21";
use 5.006001;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol;
use Carp;


=head1 NAME

Chemistry::File::SMILES - SMILES parser

=head1 SYNOPSYS

    #!/usr/bin/perl
    use Chemistry::File::SMILES;

    my $s = 'C1CC1(=O)[O-]';
    my $mol = Chemistry::Mol->parse($s, format => 'smiles');

=head1 DESCRIPTION

This module parses a SMILES (Simplified Molecular Input Line
Entry Specification) string. It registers the 'smiles' format with 
Chemistry::Mol.

=cut

# INITIALIZATION
Chemistry::Mol->register_format('smiles');
my $Smiles_parser = __PACKAGE__->new;

=begin comment

=over

=cut


sub parse_string {
    my ($self, $string, %opts) = @_;
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || "Chemistry::Atom";
    my $bond_class = $opts{bond_class} || "Chemistry::Bond";

    my $mol = $mol_class->new;
    $Smiles_parser->parse($string, $mol);
    return $mol;
}

### The contents of the original Chemistry::Smiles module start below

my $Symbol = qr/
    s|p|o|n|c|b|Zr|Zn|Yb|Y|Xe|W|V|U|Tm|Tl|Ti|Th|
    Te|Tc|Tb|Ta|Sr|Sn|Sm|Si|Sg|SeSc|Sb|S|Ru|Rn|Rh|Rf|Re|Rb|Ra|
    Pu|Pt|Pr|Po|Pm|Pd|Pb|Pa|P|Os|O|Np|No|Ni|Ne|NdNb|Na|N|Mt|Mt|
    Mo|Mn|Mg|Md|Lu|Lr|Li|La|Kr|K|Ir|In|I|Hs|Hs|Ho|Hg|Hf|He|H|Ge
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

=item Chemistry::Smiles->new([add_atom => \&sub1, add_bond => \&sub2])

Create a SMILES parser. If the add_atom and add_bond subroutine references
are given, they will be called whenever an atom or a bond needs to be added
to the molecule. If they are not specified, default methods, which
create a Chemistry::Mol object, will be used.

=cut

sub new {
    my $class = shift;
    my %opts = @_;
    require Chemistry::Mol unless $opts{add_atom} && $opts{add_bond};
    my $self = bless {
        add_atom => $opts{add_atom} || \&add_atom,
        add_bond => $opts{add_bond} || \&add_bond,
    }, $class;
}

=item $obj->parse($string, $mol)

Parse a Smiles $string. $mol is a "molecule state object". It can be anything;
the parser doesn't do anything with it except sending it as the first parameter
to the callback functions. If callback functions were not provided when
constructing the parser object, $mol must be a Chemistry::Mol object, because
that's what the default callback functions require.

=cut

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

=back

=head1 CALLBACK FUNCTIONS

=over

=item $atom = add_atom($mol, $iso, $sym, $chir, $hcount, $chg)

Called by the parser whenever an atom is found. The first parameter is the
state object given to $obj->parse(). The other parameters are the isotope,
symbol, chirality, hydrogen count, and charge of the atom. Only the symbol is
guaranteed to be defined. Mnemonic: the parameters are given in the same order
that is used in a SMILES string (such as [18OH-]). This callback is expected to
return something that uniquely identifies the atom that was created (it might
be a number, a string, or an object).

=cut

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

=item add_bond($mol, $type, $a1, $a2)

Called by the parser whenever an bond needs to be created. The first parameter
is the state object given to $obj->parse(). The other parameters are the bond
type and the two atoms that need to be bonded. The atoms are identified using
the return values from the add_atom() callback.

=back

=end comment

=cut

# Default add_bond callback 
sub add_bond {
    my ($mol, $type, $a1, $a2) = @_;
    my $order = $type_to_order{$type};
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

1;

=head1 BUGS

The SMILES specification is not fully implemented yet. For example, branches
that start before an atom (such as (OC)C, which should be equivalent to C(CO)
and COC).

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>

The SMILES Home Page at http://www.daylight.com/dayhtml/smiles/
The Daylight Theory Manual at 
http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

