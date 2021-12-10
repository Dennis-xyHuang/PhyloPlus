#!/usr/bin/env python
import pandas as pd
import sys
import os
import gc

def build_lineage(taxid_list, mapping):
    """Given a list of taxIDs, generate a NCBI LineageEx-like list."""
    temp_list = []
    for taxid in taxid_list:
        temp_dict = {}
        if taxid == "":
            pass
        else:
            temp_dict["TaxId"] = taxid
            if taxid in mapping.keys():
                temp_dict["Rank"] = mapping[taxid]
            else:
                temp_dict["Rank"] = "NA"
        temp_list.append(temp_dict)
    return temp_list

NCBI_dmp_dir_path = sys.argv[1]
assembly_gb_path = os.path.join(NCBI_dmp_dir_path, "assembly_summary_genbank.txt")
assembly_rs_path = os.path.join(NCBI_dmp_dir_path, "assembly_summary_refseq.txt")
dump_names_path = os.path.join(NCBI_dmp_dir_path, "names.dmp")
dump_nodes_path = os.path.join(NCBI_dmp_dir_path, "nodes.dmp")
dump_taxidlineage_path = os.path.join(NCBI_dmp_dir_path, "taxidlineage.dmp")
assembly_to_taxID_out_path = os.path.join(NCBI_dmp_dir_path, "assembly_to_taxID.pkl")
taxid_to_lineage_out_path = os.path.join(NCBI_dmp_dir_path, "taxid_to_lineage.pkl")

print("Combining genome assembly records from both GenBank and RefSeq sources...")
assembly_gb = pd.read_csv(assembly_gb_path, sep = "\t", skiprows = [0], low_memory = False)
assembly_rs = pd.read_csv(assembly_rs_path, sep = "\t", skiprows = [0], low_memory = False)
assembly_gb_filtered = assembly_gb[["# assembly_accession", "taxid"]]
assembly_rs_filtered = assembly_rs[["# assembly_accession", "taxid"]]
combined_df = pd.concat([assembly_gb_filtered, assembly_rs_filtered], ignore_index=True)
combined_df = combined_df.rename(columns = {"# assembly_accession": "Accession", "taxid": "taxID"})
combined_df =  combined_df.astype({"Accession": "str", "taxID": "str"})
del assembly_gb, assembly_rs, assembly_gb_filtered, assembly_rs_filtered
gc.collect()
combined_df.to_pickle(assembly_to_taxID_out_path, protocol = 4)
print("File output completed!")


print("Retrieving scientific names for different taxIDs...")
dump_names = pd.read_csv(dump_names_path, sep = "\t\|\t", header = None, engine = "python")
dump_names.iloc[:, 0] = dump_names.iloc[:, 0].astype("str")
dump_names.iloc[:, 3] = dump_names.iloc[:, 3].str.rstrip("\t\|")
dump_names_sci = dump_names[dump_names.iloc[:, 3]=="scientific name"].iloc[:, [0, 1]]
dump_names_final = dump_names_sci.rename(columns = {0: "taxID", 1: "Scientific_Name"}).reset_index(drop = True)
del dump_names, dump_names_sci
gc.collect()

print("Retrieving taxonomic ranks for different taxIDs...")
dump_nodes = pd.read_csv(dump_nodes_path, sep = "\t\|\t", header = None, engine = "python")
dump_nodes.iloc[:, 0] = dump_nodes.iloc[:, 0].astype("str")
dump_nodes_final = dump_nodes.iloc[:, [0, 2]].rename(columns = {0: "taxID", 2: "Rank"})
taxid_to_rank_dict = dict(zip(dump_nodes_final.loc[:, "taxID"], dump_nodes_final.loc[:, "Rank"]))
del dump_nodes
gc.collect()

print("Retrieving complete taxonomic lineage information for different taxIDs...")
dump_taxidlineage = pd.read_csv(dump_taxidlineage_path, sep = "\t\|\t", header = None, engine = "python")
dump_taxidlineage.iloc[:, 0] = dump_taxidlineage.iloc[:, 0].astype("str")
dump_taxidlineage.iloc[:, 1] = dump_taxidlineage.iloc[:, 1].str.rstrip(" \t\|")
dump_taxidlineage["lineage_list"] = dump_taxidlineage.iloc[:, 1].str.split(" ")
dump_taxidlineage = dump_taxidlineage.rename(columns = {0: "taxID"})
dump_taxidlineage["Lineage"] = dump_taxidlineage["lineage_list"].apply(build_lineage, args = (taxid_to_rank_dict,))
dump_taxidlineage_final = dump_taxidlineage.loc[:, ["taxID", "Lineage"]]
del dump_taxidlineage
gc.collect()

print("Retrieving information for merged taxIDs...")
dump_merged_path = os.path.join(NCBI_dmp_dir_path, "merged.dmp")
dump_merged = pd.read_csv(dump_merged_path, sep = "\t\|\t", header = None, engine = "python")
dump_merged.iloc[:, 0] = dump_merged.iloc[:, 0].astype("str")
dump_merged.iloc[:, 1] = dump_merged.iloc[:, 1].str.rstrip("\t\|")
dump_merged = dump_merged.rename(columns = {0: "previous", 1: "current"})

print("Combining all gathered infomation...")
comp_summary = pd.DataFrame({"taxID": dump_names_final["taxID"],
                             "Scientific_Name": dump_names_final["Scientific_Name"]})
comp_summary = comp_summary.merge(dump_nodes_final, how = "inner", left_on = "taxID", right_on = "taxID")
comp_summary = comp_summary.merge(dump_taxidlineage_final, how = "inner", left_on = "taxID", right_on = "taxID")
comp_summary = comp_summary.astype({"taxID": "str", "Scientific_Name": "str", "Rank": "str"})
comp_summary["CurrentID"] = comp_summary["taxID"]
temp_cond_main = comp_summary["taxID"].isin(dump_merged["previous"])
temp_cond_merged = dump_merged["previous"].isin(comp_summary["taxID"])
temp_ind_main = comp_summary.index[temp_cond_main]
temp_ind_merged = dump_merged.index[temp_cond_merged]
comp_summary.loc[temp_ind_main, "CurrentID"] = dump_merged.loc[temp_ind_merged, "current"]
comp_summary.to_pickle(taxid_to_lineage_out_path, protocol = 4)
print("File output completed!")
