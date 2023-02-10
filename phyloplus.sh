#!/bin/bash

dir_resolve()
{
    cd "$1" 2>/dev/null || return $?
    pwd -P
}

help()
{
  echo    ""
  echo -e "\033[1mSyntax:\033[0m"
  echo    "phyloplus.sh [options] <arguments>"
  echo    ""
  echo -e "\033[1mOptions:\033[0m"
  echo    "-m, --mode       Set the mode to determine which step(s) to run. Choose from"
  echo -e "                 \033[3m\"download\"\033[0m, \033[3m\"build\"\033[0m or \033[3m\"test\"\033[0m."
  echo    "-r, --ref        Location of the reference directory. Choose one of the child"
  echo    "                 directories in ./reference."
  echo    "-i, --input      Location of input taxonomy ID text file."
  echo    "-e, --email      Set the email address per NCBI Entrez requirements. Visit"
  echo -e "                 \033[3mhttps://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_\033[0m"
  echo -e "                 \033[3mGuidelines_and_Requirement\033[0m for more details."
  echo    "-o, --output     Directory to place all generated outputs. Will create this"
  echo    "                 directory if it does not exist."
  echo    "-t, --taxrank    Taxonomic rank to display query taxa in the final output."
  echo -e "                 Choose from \033[3m\"species\"\033[0m, \033[3m\"genus\"\033[0m or \
\033[3m\"family\"\033[0m (default: species)."
  echo    "-t1              Threshold used to determine outlier tips if a query taxon is"
  echo    "                 mapped at the species level (t1 ≥ 0, default: 1)."
  echo    "-t2              Threshold used to determine outlier tips if a query taxon is"
  echo    "                 mapped at the species group level (t2 ≥ 0, default: 2)."
  echo    "-t3              Threshold used to determine outlier tips if a query taxon is"
  echo    "                 mapped at the genus level or above (t3 ≥ 0, default: 2)."
  echo    "-t4              Threshold used to determine the outlier tip if only two"
  echo    "                 reference tips exist to locate a query taxon (0 < t4 < 1,"
  echo    "                 default: 0.75)."
  echo    "-h, --help       Print this help message."
  echo    ""
  echo -e "\033[1mUsages:\033[0m"
  echo -e "\033[4mDownload and summarize NCBI dump files:\033[0m"
  echo    "phyloplus.sh -m download"
  echo    ""
  echo -e "\033[4mGenerate phylogeny with the input taxonomy ID file:\033[0m"
  echo    "phyloplus.sh -m build -r <reference directory> -i <input file> -o <output directory> -e <email address> \
[-t <taxonomic rank> -t1 <threshold 1> -t2 <threshold 2> -t3 <threshold 3> -t4 <threshold 4>]"
  echo    ""
  echo -e "\033[4mRun a quick test:\033[0m"
  echo    "phyloplus.sh -m test -e <email address>"
  echo    "This is equivalent to: phyloplus.sh -m build -r ./reference/bacterial/207 -i ./sample_data/sample_taxIDs.txt\
 -o ./sample_data/sample_output -e <email address> -t species -t1 1 -t2 2 -t3 2 -t4 0.75"
  echo    ""
  echo -e "\033[7mEND\033[0m"
}

positional=()
boolean_mode=false
boolean_download=false
boolean_build=false
boolean_test=false
boolean_taxrank=false
boolean_t1=false
boolean_t2=false
boolean_t3=false
boolean_t4=false

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -m|--mode)
      mode="$2"
      case "$mode" in
        *download*)
          boolean_mode=true
          boolean_download=true
          ;;
        *build*)
          boolean_mode=true
          boolean_build=true
          ;;
        *test*)
          boolean_mode=true
          boolean_test=true
          ;;
        *)
          boolean_mode=false
          ;;
      esac
      shift 2
      ;;
    -i|--input)
      input_rel="$2"
      shift
      ;;
    -r|--ref)
      ref_dir_rel="$2"
      shift
      ;;
    -e|--email)
      email="$2"
      shift
      ;;
    -o|--output)
      output_dir_rel="$2"
      shift
      ;;
    -t|--taxrank)
      boolean_taxrank=true
      taxrank="$2"
      shift
      ;;
    -t1)
      boolean_t1=true
      t1="$2"
      shift
      ;;
    -t2)
      boolean_t2=true
      t2="$2"
      shift
      ;;
    -t3)
      boolean_t3=true
      t3="$2"
      shift
      ;;
    -t4)
      boolean_t4=true
      t4="$2"
      shift
      ;;
    -h|--help)
      help
      exit
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

if [ $boolean_mode = false ]; then
  exit
fi

pattern=" |'"
working_dir=$(dir_resolve "$(dirname "${BASH_SOURCE:-$0}")")
if [[ $working_dir =~ $pattern ]]; then
  echo "PathError: check there are no spaces in the absolute path: $working_dir" && exit
else
  cd "$working_dir" && echo "cd $working_dir"
fi

scripts_rel="./scripts/"
scripts_abs="$(dir_resolve "$scripts_rel")"
ncbi_dmp_dir_rel="./NCBI_dmp_file"
file1="assembly_summary.csv"
file2="taxonomy_summary.csv"

if [ $boolean_download = true ]; then

  ncbi_ftp_server="https://ftp.ncbi.nlm.nih.gov"
  assembly_gb_current="${ncbi_ftp_server}/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
  assembly_gb_deprecated="${ncbi_ftp_server}/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
  assembly_rs_current="${ncbi_ftp_server}/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
  assembly_rs_deprecated="${ncbi_ftp_server}/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
  taxonomy_dump="${ncbi_ftp_server}/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"

  if [ ! -d "$ncbi_dmp_dir_rel" ]
  then
      echo "mkdir -p $ncbi_dmp_dir_rel" && mkdir -p "$ncbi_dmp_dir_rel"
  fi

  ncbi_dmp_dir_abs="$(dir_resolve "$ncbi_dmp_dir_rel")"

  cd "$ncbi_dmp_dir_abs" || exit
  echo "Downloading taxonomy dump files from NCBI ftp server..."
  wget -q $taxonomy_dump
  gunzip -c new_taxdump.tar.gz | tar xf -
  echo "Downloading genome assembly summary files from NCBI ftp server..."
  wget -q $assembly_gb_current $assembly_gb_deprecated $assembly_rs_current $assembly_rs_deprecated

  cd "$working_dir" || exit
  echo ">>> python ./scripts/summarize_dmp_files.py $ncbi_dmp_dir_abs" &&
  python "$scripts_abs"/summarize_dmp_files.py "$ncbi_dmp_dir_abs"

  cd "$ncbi_dmp_dir_abs" || exit
  echo "Removes intermediate files..."
  rm -f ./*genbank*.txt ./*refseq*.txt ./*.dmp ./new_taxdump.tar.gz

else

  if [ ! -d "$ncbi_dmp_dir_rel" ]; then
    echo "FileNotFoundError: dump summary files not found, please run 'download' mode first." && exit
  else
    ncbi_dmp_dir_abs="$(dir_resolve "$ncbi_dmp_dir_rel")"
    if [ ! -f "$ncbi_dmp_dir_abs/$file1" ] || [ ! -f "$ncbi_dmp_dir_abs/$file2" ]; then
      echo "FileNotFoundError: dump summary files not found, please run 'download' mode first." && exit
    fi
  fi

  if [ $boolean_taxrank = false ]; then
    taxrank="species"
  fi
  if [ $boolean_t1 = false ]; then
    t1=1
  fi
  if [ $boolean_t2 = false ]; then
    t2=2
  fi
  if [ $boolean_t3 = false ]; then
    t3=2
  fi
  if [ $boolean_t4 = false ]; then
    t4=0.75
  fi

  if [ $boolean_test = true ]; then
    output_dir_rel="./sample_data/sample_output"
    ref_dir_rel="./reference/bacterial/207"
    input_rel="./sample_data/sample_taxIDs.txt"
  fi

  if [ $boolean_test = true ] || [ $boolean_build = true ]; then
    if [ ! -d "$output_dir_rel" ]; then
      echo "mkdir -p $output_dir_rel" && mkdir -p "$output_dir_rel"
    fi
    output_dir_abs="$(dir_resolve "$output_dir_rel")"
    if [[ $output_dir_abs =~ $pattern ]]; then
      echo "PathError: Check there are no spaces in the absolute path: $output_dir_abs" && exit
    fi
    if [ ! -d "$ref_dir_rel" ]; then
      echo "PathError: Invalid reference directory path: $ref_dir_rel" && exit
    fi
    ref_dir_abs="$(dir_resolve "$ref_dir_rel")"
    if [[ $ref_dir_abs =~ $pattern ]]; then
      echo "PathError: Check there are no spaces in the absolute path: $ref_dir_abs" && exit
    fi
    if [ ! -f "$input_rel" ]; then
      echo "PathError: Invalid input file path: $input_rel" && exit
    fi
    input_path="$(dirname "$input_rel")"
    input_file="$(basename "$input_rel")"
    input_res_path="$(dir_resolve "$input_path")"
    input_path_abs="$input_res_path/$input_file"
    if [[ $input_path_abs =~ $pattern ]]; then
      echo "PathError: Check there are no spaces in the absolute path: $input_path_abs" && exit
    fi
    printf ">>> python ./scripts/extract_lineage.py %s %s %s %s %s %s\n" \
    "$ref_dir_abs" "$input_path_abs" "$output_dir_abs" "$email" "$ncbi_dmp_dir_abs" "$taxrank" &&
    python "$scripts_abs"/extract_lineage.py "$ref_dir_abs" "$input_path_abs" "$output_dir_abs" "$email" \
    "$ncbi_dmp_dir_abs" "$taxrank" &&
    printf ">>> Rscript ./scripts/generate_phylogeny.R %s %s %s %s %s %s\n" \
    "$ref_dir_abs" "$output_dir_abs" "$t1" "$t2" "$t3" "$t4" &&
    Rscript "$scripts_abs"/generate_phylogeny.R "$ref_dir_abs" "$output_dir_abs" "$t1" "$t2" "$t3" "$t4"
    rm -f "$output_dir_abs"/*lineage.csv
  else
    exit
  fi
fi

exit