#!/home/ivan/bin/perl

use blib;
use Chemistry::File::SMILES;
use Chemistry::File::MDLMol;

my $mol = Chemistry::Mol->parse($ARGV[0] || 'CCC', format => 'smiles');
#my $mol = Chemistry::Mol->read($ARGV[0] || die "pls give a .mol\n");
printf "%s (%s)\n", $mol->name, $mol->formula;
my $smiles = $mol->print(format => 'smiles');
print "$smiles\n";
