use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SMILES',
    'VERSION_FROM'	=> 'SMILES.pm', # finds $VERSION
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol              0.24
        Chemistry::Canonicalize     0.10
        Chemistry::Bond::Find       0.21
        Chemistry::Ring::Find       0.10
        List::Util                  0
    )},
);
