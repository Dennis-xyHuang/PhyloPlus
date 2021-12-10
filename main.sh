#!/bin/bash

columnize () {
    indent=$1;
    collen=$(($(tput cols)-indent));
    keyname="$2";
    value=$3;
    while [ -n "$value" ] ; do
        printf "%-20s %-${indent}s\n" "$keyname" "${value:0:$collen}";
        keyname="";
        value=${value:$collen};
    done
}

Help()
{
    # Display Help
    echo ""
    echo "Syntax: ./main.sh [-m <mode>][-p <phylogeny path>][-o <output directory>] [-i <input taxID path> | -e <e-mail address> | -t1 <threshold 1> | -t2 <threshold 2> | -t3 <threshold 3> | -t4 <threshold 4>]"
    echo ""
    echo "Flags:"
    columnize 25 "-m, --mode" "Set the mode to determine which step(s) to run. Choose <mode> from \"extract\", \"expand\" or \"both\"."
    columnize 25 "-p, --phylo_path" "Location of the reference phylogeny (default: ./sample_data/sample_phylo.tree)."
    columnize 25 "-o, --output_dir" "Directory to place all generated outputs. Will create adirectory if it does not exist (default: ./sample_data/sample_output)."
    columnize 25 "-i, --taxid_path" "Location of input taxID text file (default: ./sample_data/sample_taxIDs.txt)."
    columnize 25 "-e, --email" "Set the email address per NCBI Entrez batch processing requirements. See https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requirement for  more details."
    columnize 25 "-t1" "Set the threshold for Z-score above which is used to   determine outlier tips within a group of reference tipsto locate species that can be mapped at the species le-vel (default: 1)."
    columnize 25 "-t2" "Set the threshold for Z-score above which is used to   determine outlier tips within a group of reference tipsto locate species that can be mapped at the species gr-oup level (default: 2)."
    columnize 25 "-t3" "Set the threshold for Z-score above which is used to   determine outlier tips within a group of reference tipsto locate species that can be mapped at the genus level(default: 3)."
    columnize 25 "-t4" "Set the threshold for the fraction above which is used to detect presence of an outlier tip when only two ref-erence tips are available to locate a species. The fra-ction used is defined as (the distance to the MRCA node) / (the distance to the base root) for the more dista-nt tip (default: 0.75)."
    echo ""
    echo "Options:"
    echo "Retrieve complete lineage information:"
    echo "Usage: ./main.sh -m extract -p <phylogeny path> -o <output directory> -i <input taxID path> - e <e-mail address>"
    echo ""
    echo "Locate species within the reference phylogeny:"
    echo "Usage: ./main.sh -m expand -p <phylogeny path> -o <output directory> -t1 <threshold 1> -t2 <threshold 2> -t3 <threshold 3> -t4 <threshold 4>"
    echo ""
    echo "Perform both steps consecutively:"
    echo "Usage: ./main.sh -m both -p <phylogeny path> -o <output directory> -i <input taxID path> - e <e-mail address> -t1 <threshold 1> -t2 <threshold 2> -t3 <threshold 3> -t4 <threshold 4>"
    echo "";
}

dir_resolve()
{
    cd "$1" 2>/dev/null || return $?
    echo "`pwd -P`"
}
pattern=" |'"
mypwd=`pwd -P`
POSITIONAL=()
M=false
I=false
P=false
R=false
Phylo=false
O=false
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -m|--mode)
      mode="$2"
      case "$mode" in
        *extract*)
        P=true
	;;
      esac
      case "$mode" in
        *expand*)
        R=true
        ;;
      esac
      case "$mode" in
        *both*)
        R=true
	    P=true
	;;
      esac
      shift
      shift
      ;;
    -i|--taxid_path)
      I=true
      input_rel="$2"
      shift
      ;;
    -p|--phylo_path)
      Phylo=true
      phylop_file="$2"
      shift
      ;;
    -e|--email)
      email="$2"
      shift
      ;;
    -o|--output_dir)
      O=true
      output_rel="$2"
      shift
      ;;
    -t1)
      t_1=true
      t1="$2"
      shift
      ;;
    -t2)
      t_2=true
      t2="$2"
      shift
      ;;
    -t3)
      t_3=true
      t3="$2"
      shift
      ;;
    -t4)
      t_4=true
      t4="$2"
      shift
      ;;
    -h|--help)
      Help
      exit;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

# Default directory values
def_sample_data_rel="./sample_data/"
def_output_rel="./sample_data/sample_output/"
def_phylo_rel="sample_phylo.tree"
def_input_rel="sample_taxIDs.txt"
scripts="./scripts"
def_sample_data_abs="`dir_resolve "$def_sample_data_rel"`"
def_phylop_abs_path="$def_sample_data_abs/`basename \"$def_phylo_rel\"`"
def_i_abs_path="$def_sample_data_abs/`basename \"$def_input_rel\"`"
def_output_abs_path="`dir_resolve \"$def_output_rel\"`"
f1="tip_lineage.tsv"
f2="phylo_spp_lineage.tsv"
f3="added_spp_lineage.tsv"
dmp_abs="${mypwd}/NCBI_dmp_files"

if $O ; then
    cd $mypwd;
    if [ ! -d $output_rel ]; then
        echo "mkdir -p $output_rel" && mkdir -p $output_rel;
    fi;
    output="`dir_resolve \"$output_rel\"`"
else
    cd $mypwd;
    if [ ! -d $def_output_abs_path ]; then
        echo "mkdir -p $def_output_abs_path" && mkdir -p $def_output_abs_path;
    fi;
    output="$def_output_abs_path"
fi

if [[ $output =~ $pattern ]]; then
    echo "OUTPUT PATH ERROR: Check there are no spaces in the absolute path." && exit
fi;

if $Phylo ; then
    phylop_abs_file="`dirname \"$phylop_file\"`"
    phylop_name="`basename \"$phylop_file\"`"
    phylop_path="`dir_resolve \"$phylop_abs_file\"`"
    phylop="$phylop_path/$phylop_name"
    # echo "phylop is $phylop"
    # for this check [[ ]] is required, [ ] will not work
else
    phylop="$def_phylop_abs_path";
fi

if [[ $phylop =~ $pattern ]]; then
    echo "PHYLO PATH ERROR: Check there are no spaces in the absolute path." && exit
fi;

if $I ; then
    i_path="`dirname \"$input_rel\"`"
    i_file="`basename \"$input_rel\"`"
    i_res_path="`dir_resolve \"$i_path\"`"
    i_path="$i_res_path/$i_file"
else
    i_path="$def_i_abs_path";
fi

if [[ $i_path =~ $pattern ]]; then
    echo "Input Taxid PATH ERROR: Check there are no spaces in the absolute path." && exit
fi;

#echo "i_path: `dirname \"$input_rel\"`"
#echo "phylop: $phylop is used: $Phylo"
#echo "input: $i_path is used: $I"
#echo "output: $output is used: $O"
# Default values for thresholds
if [ ! $t_1 ]; then
    t1=1
fi
if [ ! $t_2 ]; then
    t2=2
fi
if [ ! $t_3 ]; then
    t3=3
fi
if [ ! $t_4 ]; then
    t4=0.75
fi

# Check if user files for output in R only are valid
if ! $P; then
    if $R; then
        if [ ! -f "$output/$f1" ] || [ ! -f "$output/$f2" ] || [ ! -f "$output/$f3" ]; then
            echo "Lineage information not found, please run 'extract' mode first." && exit;
        fi;
    fi;
fi

#echo "phylop: $phylop is used: $Phylo"
#echo "input: $i_path is used: $I"
#echo "output: $output is used: $O"
# values have been funneled into these variables and tested for whitespaces
if $P && $R; then
    printf ">>> python $scripts/extract_lineage.py $phylop $i_path $output $email $dmp_abs \n&&\n>>> Rscript $scripts/expand_phylogeny.R $phylop $output $t1 $t2 $t3 $t4\n" && python $scripts/extract_lineage.py $phylop $i_path $output $email $dmp_abs && Rscript $scripts/expand_phylogeny.R $phylop $output $t1 $t2 $t3 $t4
else
    if $P; then
        echo ">>> python $scripts/extract_lineage.py $phylop $i_path $output $email $dmp_abs" && python $scripts/extract_lineage.py $phylop $i_path $output $email $dmp_abs
    fi;

    if $R; then
        echo ">>> Rscript $scripts/expand_phylogeny.R $phylop $output $t1 $t2 $t3 $t4" && Rscript $scripts/expand_phylogeny.R $phylop $output $t1 $t2 $t3 $t4;
    fi;
fi
