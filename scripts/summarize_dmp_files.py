#!/usr/bin/env python
import sys
import os
import pandas as pd
import phyloplus

NCBI_dmp_dir_path = sys.argv[1]

# Processing assembly summary reports. Note that Entrez query of an outdated assembly accession number will return
# document summary of the outdated accession number, rather than the one with updated versions.

print("Processing genome assembly summary reports from both GenBank and RefSeq sources...")
NCBI_assembly_reports = ["assembly_summary_genbank.txt", "assembly_summary_genbank_historical.txt",
                         "assembly_summary_refseq.txt", "assembly_summary_refseq_historical.txt"]
NCBI_assembly_reports_path = {name: os.path.join(NCBI_dmp_dir_path, name) for name in NCBI_assembly_reports}
assembly_output_path = os.path.join(NCBI_dmp_dir_path, "assembly_summary.csv")
NCBI_assembly_dataframes = []
for assembly_report_name, assembly_report_path in NCBI_assembly_reports_path.items():
    temp_assembly_df = pd.read_csv(assembly_report_path, sep="\t", skiprows=[0], low_memory=False)
    subset_colnames = []
    subset_colnames.append(next(filter(lambda x: 'assembly_accession' in x, temp_assembly_df.columns)))
    subset_colnames.append(next(filter(lambda x: 'species_taxid' in x, temp_assembly_df.columns)))
    temp_assembly_df = temp_assembly_df[subset_colnames]
    NCBI_assembly_dataframes.append(temp_assembly_df)
final_assembly_df = pd.concat(NCBI_assembly_dataframes)
final_assembly_df = final_assembly_df.dropna(axis=0, how="any").reset_index(drop=True)
final_assembly_df = final_assembly_df.rename(columns={subset_colnames[0]: "Accession", subset_colnames[1]: "TaxID"})
final_assembly_df.to_csv(assembly_output_path, header=True, index=False)
print("Assembly summary file output completed!")

# Processing taxonomic summary reports. Note that Entrez query of a merged taxonomy ID will automatically return
# document summary of the updated taxonomy ID.

print("Processing taxonomy summary reports...")
dump_names_path = os.path.join(NCBI_dmp_dir_path, "names.dmp")
dump_nodes_path = os.path.join(NCBI_dmp_dir_path, "nodes.dmp")
dump_lineage_path = os.path.join(NCBI_dmp_dir_path, "taxidlineage.dmp")
dump_merged_path = os.path.join(NCBI_dmp_dir_path, "merged.dmp")
taxonomy_output_path = os.path.join(NCBI_dmp_dir_path, "taxonomy_summary.csv")

# Read taxonomic dump files and preprocess them.
# Retrieve scientific names: names.dmp has three fields: tax_id, name_txt, unique name, name class.
# Retrieve taxonomic ranks: nodes.dmp has the following fields: tax_id, parent tax_id, rank, embl code, etc.
# Retrieve lineage information: taxidlineage.dmp has two fields: tax_id, lineage.
# Retrieve merged TaxIDs: merged.dmp has two fields: old_tax_id, new_tax_id.

dump_names = phyloplus.read_dump(dump_names_path, {3: "scientific name"}, {0: "TaxID", 1: "ScientificName"})
dump_nodes = phyloplus.read_dump(dump_nodes_path, None, {0: "TaxID", 2: "TaxonomicRank"})
dump_lineage = phyloplus.read_dump(dump_lineage_path)
# Convert sequence of TaxIDs to LineageEx-like form.
LineageEx_map = dict(zip(dump_nodes.loc[:, "TaxID"],
                         dump_nodes.apply(lambda x: {"TaxId": x.loc["TaxID"], "Rank": x.loc["TaxonomicRank"]}, axis=1)))
# Split LineageEx-like list of dictionaries into different columns.
dump_lineage.loc[:, "Lineage"] = dump_lineage.iloc[:, 1].apply(phyloplus.build_lineage, args=(LineageEx_map,))
dump_lineage = phyloplus.convert_lineage_to_ranks(dump_lineage, "Lineage", phyloplus.taxonomic_ranks_cols)
taxo_rank_cols = phyloplus.taxonomic_ranks_cols.values()
dump_lineage = dump_lineage.rename(columns={0: "TaxID"}).loc[:, ["TaxID", *taxo_rank_cols]]
dump_merged = phyloplus.read_dump(dump_merged_path, None, {0: "TaxID", 1: "CurrentTaxID"})

# Combine all taxonomy-related information
final_taxid_df = pd.DataFrame({"TaxID": dump_names.loc[:, "TaxID"], "CurrentTaxID": dump_names.loc[:, "TaxID"]})
final_taxid_df = pd.concat([final_taxid_df, dump_merged], axis=0, ignore_index=True)
final_taxid_df = final_taxid_df.merge(dump_names, left_on="CurrentTaxID", right_on="TaxID", suffixes=(None, "_1"))
final_taxid_df = final_taxid_df.merge(dump_nodes, left_on="CurrentTaxID", right_on="TaxID", suffixes=(None, "_2"))
final_taxid_df = final_taxid_df.merge(dump_lineage, left_on="CurrentTaxID", right_on="TaxID", suffixes=(None, "_3"))
final_taxid_df = final_taxid_df.drop(columns=["TaxID_1", "TaxID_2", "TaxID_3"])
# LineageEx does not include the query TaxID itself in the list of dictionaries.
final_taxid_df = phyloplus.finalize_taxonomy_summary(final_taxid_df, "CurrentTaxID", "TaxonomicRank",
                                                     phyloplus.taxonomic_ranks_cols)
final_taxid_df.to_csv(taxonomy_output_path, header=True, index=False)
print("Taxonomy summary file output completed!")
