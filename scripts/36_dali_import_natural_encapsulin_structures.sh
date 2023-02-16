#!/bin/bash

mkdir structures/experimental

wget "https://files.rcsb.org/download/7MU1.pdb";
wget "https://files.rcsb.org/download/6NJ8.pdb";
wget "https://files.rcsb.org/download/6X8M.pdb";
wget "https://files.rcsb.org/download/7S2T.pdb";

mv *.pdb structures/experimental/;

perl ../DALI/DaliLite.v5/bin/import.pl --pdbfile structures/experimental/7MU1.pdb --pdbid 7MU1 --dat DALI/data/ --clean;
perl ../DALI/DaliLite.v5/bin/import.pl --pdbfile structures/experimental/6NJ8.pdb --pdbid 6NJ8 --dat DALI/data/ --clean;
perl ../DALI/DaliLite.v5/bin/import.pl --pdbfile structures/experimental/6X8M.pdb --pdbid 6X8M --dat DALI/data/ --clean;
perl ../DALI/DaliLite.v5/bin/import.pl --pdbfile structures/experimental/7S2T.pdb --pdbid 7S2T --dat DALI/data/ --clean;