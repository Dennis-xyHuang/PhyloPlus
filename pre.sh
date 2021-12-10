#!/bin/bash

dir_resolve()
{
    cd "$1" 2>/dev/null || return $?
    echo "`pwd -P`"
}
mypwd=`pwd -P`
pattern=" |'"
scripts="./scripts"
NCBI_FTP_SERVER="ftp://ftp.ncbi.nlm.nih.gov"
NCBI_TAXONOMY="${NCBI_FTP_SERVER}/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
NCBI_ASSEMBLY_GB="${NCBI_FTP_SERVER}/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
NCBI_ASSEMBLY_RS="${NCBI_FTP_SERVER}/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"

if [[ $mypwd =~ $pattern ]]; then
    echo "PATH ERROR: Check there are no spaces in the absolute path." && exit
fi;

cd $mypwd
if [ ! -d "./NCBI_dmp_files" ]; then
    echo "mkdir -p ./NCBI_dmp_files" && mkdir -p ./NCBI_dmp_files;
fi;
NCBI_dmp_files_rel="./NCBI_dmp_files"
NCBI_dmp_files_abs="`dir_resolve "$NCBI_dmp_files_rel"`"

cd ./NCBI_dmp_files
echo "Downloading taxonomy dump files from NCBI ftp server..." && wget -q $NCBI_TAXONOMY
gunzip -c new_taxdump.tar.gz | tar xf -
echo "Downloading genome assembly summary files from NCBI ftp server..." && wget -q $NCBI_ASSEMBLY_GB $NCBI_ASSEMBLY_RS

cd ..
echo ">>> python $scripts/summarize_dmp_files.py $NCBI_dmp_files_abs" && python $scripts/summarize_dmp_files.py $NCBI_dmp_files_abs

cd ./NCBI_dmp_files
find . -type f -not -name '*.pkl' -delete
