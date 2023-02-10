#!/usr/bin/env python
import os
import sys
import pandas as pd
import phyloplus
from Bio import Phylo

input_reference_dir_path = sys.argv[1]
input_user_file_path = sys.argv[2]
output_dir_path = sys.argv[3]
email = sys.argv[4]
NCBI_dmp_dir_path = sys.argv[5]
taxonomic_rank = sys.argv[6]

assert email != "", "Please provide an email address."
assert os.path.exists(NCBI_dmp_dir_path), "NCBI dump files have not been downloaded and processed."
assert taxonomic_rank in ["species", "genus", "family"], \
    "Specified taxonomic rank must be one of 'species', 'genus', or 'family'"

if not os.path.exists(output_dir_path):
    os.makedirs(output_dir_path)

reference_files = os.listdir(input_reference_dir_path)
input_reference_file = [name for name in reference_files if name.startswith("bac") or name.startswith("arc")][0]
input_phylo_file = [name for name in reference_files if name.startswith("GTDB")][0]
input_reference_path = os.path.join(input_reference_dir_path, input_reference_file)
input_phylo_path = os.path.join(input_reference_dir_path, input_phylo_file)
output_tip_lineage_file_path = os.path.join(output_dir_path, "tip_lineage.csv")
output_user_lineage_file_path = os.path.join(output_dir_path, "user_lineage.csv")
output_note_file_path = os.path.join(output_dir_path, "note.txt")

print("Loading and processing input files...")
# reference phylogeny
phylo_tree = Phylo.read(input_phylo_path, "newick")
phylo_tip_labels = set([clade.name for clade in phylo_tree.find_clades() if clade.name])
tip_labels = pd.read_csv(input_reference_path, header=None, sep="\t", names=["TipLabel", "Lineage"])
tip_labels.loc[:, "Accession"] = tip_labels.loc[:, "TipLabel"].apply(lambda x: x.strip()[3:])
tip_label_groups = dict((j, i) for (i, j) in enumerate(tip_labels["Lineage"].sort_values().unique()))
tip_labels.loc[:, "TipGroup"] = tip_labels.loc[:, "Lineage"].apply(lambda x: tip_label_groups[x])
tip_group_members = tip_labels.groupby(by=['TipGroup'])["TipLabel"].apply(list)
group_member_series = tip_labels.loc[:, "TipGroup"].apply(lambda x: tip_group_members.loc[x])
representative_tip_series = group_member_series.apply(lambda x: list(set(x) & phylo_tip_labels)[0])
tip_labels.loc[:, "TipOnPhylo"] = representative_tip_series
tip_labels = tip_labels.loc[:, ["TipLabel", "Accession", "TipGroup", "TipOnPhylo"]]
# user input
user_file = open(input_user_file_path)
unique_user_taxids = []
for line in user_file:
    line = line.strip()
    taxid = line.split(sep="\t")[-1].strip()
    if taxid != '':
        unique_user_taxids.append(taxid)
unique_user_taxids = [int(taxid) for taxid in set(unique_user_taxids)]
# dump summary files
dump_assembly = pd.read_csv(os.path.join(NCBI_dmp_dir_path, "assembly_summary.csv"))
dump_taxonomy = pd.read_csv(os.path.join(NCBI_dmp_dir_path, "taxonomy_summary.csv"))

print("*** Step 1-1 ***")
print("Extracting lineage information for the reference phylogeny...")
print("Converting genome assembly accession numbers extracted from phylogenetic tip labels into NCBI taxonomy IDs...")
print("Retrieving genome assembly information from downloaded dump files...")
output_assembly_filter = dump_assembly.loc[:, "Accession"].isin(tip_labels.loc[:, "Accession"])
output_assembly_1 = dump_assembly.loc[output_assembly_filter, :]
remaining_accessions = set(tip_labels.loc[:, "Accession"]) - set(output_assembly_1.loc[:, "Accession"])
if len(remaining_accessions):
    print("Retrieving remaining genome assembly information from the NCBI server...")
    output_assembly_2, abnormal_accessions = phyloplus.convert_acc_to_taxid_main(remaining_accessions, email)
    if len(abnormal_accessions):
        print("The following genome assembly accession numbers can not be found by Entrez:\n" +
              ", ".join(abnormal_accessions))
    output_assembly = pd.concat([output_assembly_1, output_assembly_2], axis=0, ignore_index=True)
else:
    output_assembly = output_assembly_1
output_assembly = tip_labels.merge(output_assembly, on="Accession", how="left")
print("Summarization of genome assembly information is done!")

print("*** Step 1-2 ***")
print("Extracting full lineage information for taxonomy IDs found in the reference phylogeny...")
unique_ref_taxids = list(output_assembly.loc[:, "TaxID"].dropna().unique())
print("Retrieving taxonomic information from downloaded dump files...")
output_taxonomy_ref_filter = dump_taxonomy.loc[:, "TaxID"].isin(unique_ref_taxids)
output_taxonomy_ref_1 = dump_taxonomy.loc[output_taxonomy_ref_filter, :]
remaining_taxids_ref = set(unique_ref_taxids) - set(output_taxonomy_ref_1.loc[:, "TaxID"])
if len(remaining_taxids_ref):
    print("Retrieving remaining taxonomic information from the NCBI server...")
    output_taxonomy_ref_2, abnormal_taxids_ref = \
        phyloplus.convert_taxid_to_lineage_main(remaining_taxids_ref, email, phyloplus.taxonomic_ranks_cols)
    if len(abnormal_taxids_ref):
        print("The following taxonomy IDs can not be found by Entrez:\n" +
              ", ".join([str(taxid) for taxid in abnormal_taxids_ref]))
    output_taxonomy_ref = pd.concat([output_taxonomy_ref_1, output_taxonomy_ref_2], axis=0, ignore_index=True)
else:
    output_taxonomy_ref = output_taxonomy_ref_1
output_taxonomy_ref = output_taxonomy_ref.reset_index(drop=True)
print("Summarization of taxonomic information is done!")

print("*** Step 1-3 ***")
print("Generating the summary csv file for reference phylogeny...")
final_output_ref = output_assembly.merge(output_taxonomy_ref, on="TaxID", how="left")
final_output_ref = final_output_ref.drop(columns=["TaxID"]).rename(columns={"CurrentTaxID": "TaxID"})
final_output_ref.to_csv(output_tip_lineage_file_path, sep=",", index=False)
print("File output completed!")

print("*** Step 2-1 ***")
print("Extracting full lineage information for taxonomy IDs found in the user input file...")
print("Retrieving taxonomic information from downloaded dump files...")
output_taxonomy_user_filter = dump_taxonomy.loc[:, "TaxID"].isin(unique_user_taxids)
output_taxonomy_user_1 = dump_taxonomy.loc[output_taxonomy_user_filter, :]
remaining_taxids_user = set(unique_user_taxids) - set(output_taxonomy_user_1.loc[:, "TaxID"])
abnormal_taxids_user = []
if len(remaining_taxids_user):
    print("Retrieving remaining taxonomic information from the NCBI server...")
    output_taxonomy_user_2, abnormal_taxids_user = \
        phyloplus.convert_taxid_to_lineage_main(remaining_taxids_user, email, phyloplus.taxonomic_ranks_cols)
    if len(abnormal_taxids_user):
        print("The following taxonomy IDs can not be found by Entrez:" +
              ", ".join([str(taxid) for taxid in abnormal_taxids_user]))
    output_taxonomy_user = pd.concat([output_taxonomy_user_1, output_taxonomy_user_2], axis=0, ignore_index=True)
else:
    output_taxonomy_user = output_taxonomy_user_1
output_taxonomy_user = output_taxonomy_user.reset_index(drop=True)
print("Summarization of taxonomic information is done!")

print("*** Step 2-2 ***")
print("Updating taxonomy IDs found in the user input that are outdated or merged (if any) ...")
outdated_filter = output_taxonomy_user.loc[:, "TaxID"] != output_taxonomy_user.loc[:, "CurrentTaxID"]
outdated_taxonomy = output_taxonomy_user.loc[outdated_filter, ["TaxID", "CurrentTaxID"]]
output_taxonomy_user = output_taxonomy_user.drop_duplicates(subset=["CurrentTaxID"], ignore_index=True)

with open(output_note_file_path, "w") as text_file:
    if len(abnormal_taxids_user) == 0:
        text_file.write("All taxonomy IDs from user input can be found in NCBI taxonomy database.\n")
    else:
        text_file.write("The following taxonomy IDs from user input cannot be found in NCBI taxonomy database and "
                        "are removed before further processing:\n")
        for taxid in abnormal_taxids_user:
            text_file.write("{0}\n".format(taxid))
    if len(outdated_taxonomy) == 0:
        text_file.write("\nAll taxonomy IDs from user input are up-to-date, no changes are made.\n")
    else:
        text_file.write("\nThe following taxonomy IDs from user input are outdated. They are changed or merged into "
                        "their corresponding current IDs. Only current taxonomy IDs are used in further processing:\n"
                        "Input taxonomy ID\t--->\tCurrent taxonomy ID\n")
        for outdated_id, current_id in zip(outdated_taxonomy.loc[:, "TaxID"], outdated_taxonomy.loc[:, "CurrentTaxID"]):
            text_file.write("{0}\t---->\t{1}\n".format(int(outdated_id), int(current_id)))

print("*** Step 2-3 ***")
print("Generating the csv file summarizing lineage information for unique taxa at the user-defined taxonomic rank that "
      "are identified in the user input.")
target_taxids = output_taxonomy_user.loc[:, phyloplus.taxonomic_ranks_cols[taxonomic_rank]].dropna().unique().tolist()
print("Retrieving taxonomic information from downloaded dump files...")
output_taxonomy_user_filter_updated = dump_taxonomy.loc[:, "TaxID"].isin(target_taxids)
output_taxonomy_user_1_updated = dump_taxonomy.loc[output_taxonomy_user_filter_updated, :]
remaining_taxids_user_updated = set(target_taxids) - set(output_taxonomy_user_1_updated.loc[:, "TaxID"])
abnormal_taxids_user_updated = []
if len(remaining_taxids_user_updated):
    print("Retrieving remaining taxonomic information from the NCBI server...")
    output_taxonomy_user_2_updated, abnormal_taxids_user_updated = \
        phyloplus.convert_taxid_to_lineage_main(remaining_taxids_user_updated, email, phyloplus.taxonomic_ranks_cols)
    if len(abnormal_taxids_user_updated):
        print("The following taxonomy IDs can not be found by Entrez:" +
              ", ".join([str(taxid) for taxid in abnormal_taxids_user_updated]))
    output_taxonomy_user_updated = pd.concat([output_taxonomy_user_1_updated, output_taxonomy_user_2_updated],
                                             axis=0, ignore_index=True)
else:
    output_taxonomy_user_updated = output_taxonomy_user_1_updated
output_taxonomy_user_updated = output_taxonomy_user_updated.reset_index(drop=True)
print("Summarization of taxonomic information is done!")
final_output_user = output_taxonomy_user_updated.drop(columns="TaxID").rename(columns={"CurrentTaxID": "TaxID"})
final_output_user.to_csv(output_user_lineage_file_path, sep=",", index=False)
print("File output completed!")
