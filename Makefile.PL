use 5.006001;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SMILES',
    'VERSION_FROM'	=> 'SMILES.pm', # finds $VERSION
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => "0.10",
    },
);
