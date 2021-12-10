#!/usr/bin/env python
import pandas as pd
import numpy as np
from Bio import Entrez, Phylo
from math import log10, ceil
import time
import sys
import os

def keep_trying(func, max_tries, sleep = None, *args):
    counter = []
    while True:
        counter.append(1)
        try:
            return func(*args)
        except Exception:
        # Entrez batch search may encounter IndexError and HttpError.
            if len(counter) == max_tries:
                print("Maximum limit reached.")
                break
            if sleep is not None:
                time.sleep(sleep)
            else:
                break

def convert_acc_to_taxid_main(acc_list):
    """
    Given a list of genome assembly accession numbers, returns a dataframe summarizing corresponding taxIDs.
    """
    def from_acc_to_taxid(acc):
        Entrez.email = email
        handle = Entrez.esearch(db = "assembly", term = acc)
        record = Entrez.read(handle)
        id_key = record['IdList'][0]
        summ_handle = Entrez.esummary(db = "assembly", id = id_key, report = "full")
        summ_record = Entrez.read(summ_handle, validate = False)
        taxID = summ_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        return taxID
    acc_taxid = {}
    start_time = time.time()
    ave_proc_time = 0
    n = len(acc_list)
    if n == 0:
        print("", end = "\r")
    else:
        for i in range(n):
            print((">>> The {0:^{width}d}th genome being processed. Progress: {1:^6.2f}%; ETR: {2:^6.2F} min. <<<").format(
                   i + 1, (i + 1) * 100 / n, (n - i - 1) * ave_proc_time, width = ceil(log10(n + 1))), end = "\r")
            acc_key = acc_list[i]
            taxid = keep_trying(from_acc_to_taxid, 10, 5, acc_key)
            acc_taxid[acc_key] = taxid
            end_time = time.time()
            proc_time = (end_time - start_time) / 60
            ave_proc_time = proc_time / (i + 1)
            time.sleep(0.3)
        print("\n", end = "\r")
    print("Data retrieval completed!")
    df_acc = []
    df_taxid = []
    for i, j in acc_taxid.items():
        df_acc.append(i)
        df_taxid.append(j)
    tip_taxid = pd.DataFrame({"Accession": df_acc, "taxID": df_taxid})
    return tip_taxid

def convert_taxid_to_lineage_main(uniq_id_list):
    """
    Given a list of non-redudant NCBI taxIDs, returns a lineage summary dataframe containing the taxID itself, the
    scientific name, the taxonomic rank, the full lineage, and the current taxID assigned to the organism.
    """
    def from_taxid_to_lineage(taxid):
        info_dict = {}
        Entrez.email = email
        handle = Entrez.efetch(db = "Taxonomy", id = taxid, retmode = "xml")
        record = Entrez.read(handle)
        info_dict["Lineage"] = record[0]["LineageEx"]
        info_dict["ScientificName"] = record[0]["ScientificName"]
        info_dict["Rank"] = record[0]["Rank"]
        info_dict["CurrentID"] = record[0]["TaxId"]
        return info_dict
    def dict_slct(keyname, dictionary):
        return dictionary[keyname]
    uniq_sci_name = {}
    uniq_rank = {}
    uniq_lineage = {}
    uniq_currID = {}
    start_time = time.time()
    ave_proc_time = 0
    n = len(uniq_id_list)
    if n == 0:
        print("", end = "\r")
    else:
        for i in range(n):
            print((">>> The {0:^{width}d}th taxID being processed. Progress: {1:^6.2f}%; ETR: {2:^6.2F} min. <<<").format(
                   i + 1, (i + 1) * 100 / n, (n - i - 1) * ave_proc_time, width = ceil(log10(n + 1))), end = "\r")
            taxid = uniq_id_list[i]
            temp_dict = keep_trying(from_taxid_to_lineage, 10, 5, taxid)
            uniq_sci_name[taxid] = temp_dict["ScientificName"]
            uniq_rank[taxid] = temp_dict["Rank"]
            uniq_lineage[taxid] = temp_dict["Lineage"]
            uniq_currID[taxid] = temp_dict["CurrentID"]
            end_time = time.time()
            proc_time = (end_time - start_time) / 60
            ave_proc_time = proc_time / (i + 1)
            time.sleep(0.25)
        print("\n", end = "\r")
    print("Data retrieval completed!")
    sum_df = pd.DataFrame({"taxID": list(uniq_sci_name.keys()), "Scientific_Name": np.NaN, "Rank": np.NaN,
                           "Lineage": np.NaN, "CurrentID": np.NaN})
    sum_df["Scientific_Name"] = sum_df["taxID"].apply(dict_slct, args = (uniq_sci_name,))
    sum_df["Rank"] = sum_df["taxID"].apply(dict_slct, args = (uniq_rank,))
    sum_df["Lineage"] = sum_df["taxID"].apply(dict_slct, args = (uniq_lineage,))
    sum_df["CurrentID"] = sum_df["taxID"].apply(dict_slct, args = (uniq_currID,))
    return sum_df

def check_outdated_taxids(df, text):
    """
    Given a lineage summary dataframe, check if all extracted taxIDs are up to date. If not, change or merge them into
    current ones and delete duplicated records.
    """
    out_dict = {}
    if all(df["CurrentID"] == df["taxID"]) == True:
        print("All taxIDs identified in the " + text + " are up to date.")
    else:
        print("Not all taxIDs identified in the " + text + " are up to date:")
        print("Identified taxID   -->   Current taxID")
        condition = df["CurrentID"] != df["taxID"]
        temp_df = df[condition]
        for i, j in zip(temp_df["taxID"], temp_df["CurrentID"]):
            out_dict[i] = j
            print("{0:^16s}   -->   {1:^13s}".format(i, j))
            temp_cond = df["taxID"] == i
            temp_ind = df.index[temp_cond]
            df.loc[temp_ind, "taxID"] = j
        print("Records with outdated taxIDs are changed/merged into current ones.")
        df = df.drop_duplicates(subset = ["taxID"]).reset_index(drop = True)
    return (df, out_dict)

def convert_lineage_to_ranks(df):
    """
    Given a lineage summary dataframe, extract the taxID at different taxonomic ranks and store them in separate
    columns. Check if all organisms appearing in the dataframe have their species-level taxID.
    """
    def ext_rank(lineage_dict, rank):
        for i in range(len(lineage_dict)):
            # full lineage information is stored as a list of dictionaries.
            if lineage_dict[i]["Rank"] == rank:
                return lineage_dict[i]["TaxId"]
            else:
                continue
    df = df.reset_index(drop = True)
    df["superkingdom"] = df["Lineage"].apply(ext_rank, args = ("superkingdom",))
    df["phylum"] = df["Lineage"].apply(ext_rank, args = ("phylum",))
    df["class"] = df["Lineage"].apply(ext_rank, args = ("class",))
    df["order"] = df["Lineage"].apply(ext_rank, args = ("order",))
    df["family"] = df["Lineage"].apply(ext_rank, args = ("family",))
    df["genus"] = df["Lineage"].apply(ext_rank, args = ("genus",))
    df["species_group"] = df["Lineage"].apply(ext_rank, args = ("species group",))
    df["species"] = df["Lineage"].apply(ext_rank, args = ("species",))
    contradict_indices = []
    for i in range(len(df)):
        if df.loc[i,"Rank"] == "species":
            if df.loc[i,"species"] is None:
                df.loc[i,"species"] = df.loc[i,"taxID"]
            else:
                contradict_indices.append(i)
    if len(contradict_indices) >= 1:
        print("Potential contradictory record(s) found while converting the original lineage dataframe into a new " +
              "dataframe where taxIDs for different taxonomic rank are stored in separate columns:")
        for item in contradict_indices:
            print("Manual verification is required for ", df.loc[item,"Scientific_Name"], " (taxID: ",
                  df.loc[item,"taxID"], ") in terms of its species-level taxID.", sep = "")
    assert len(contradict_indices) == 0, ("Incorrect full lineage information retrieved for some taxIDs. See above " +
                                          "printed information for details.")
    return df

def uniq_spp_lineage(raw_lineage_df, uniq_spp_list, dmp_df):
    """
    Given a lineage summary dataframe, generate a dataframe where full lineage information for each of the unique
    species present in the original lineage dataframe is recorded.
    """
    temp_df = (raw_lineage_df[raw_lineage_df["Rank"] == "species"].drop_duplicates(subset = ["taxID"]).
               reset_index(drop = True)[["taxID", "Rank", "Scientific_Name", "Lineage"]])
    int_ids = list(set(uniq_spp_list).intersection(set(temp_df["taxID"])))
    diff_ids = list(set(uniq_spp_list).difference(set(temp_df["taxID"])))
    print("Inspecting all species-level taxIDs identified in the raw lineage summary dataframe:")
    if len(diff_ids) > 0:
        # Some genome assemblies / taxIDs may have been assigned taxonomic ranks lower than species (e.g., serovar,
        # subspecies, etc.). Their species-level taxIDs may or may not have already been directly extracted from
        # the same source. In the latter case, additional searches are needed for their corresponding species-level
        # taxIDs.
        print(("For species-level taxIDs present in the raw lineage summary dataframe, {0} out of {1} are " +
               "derived from lower level taxIDs rather than directly provided. Full lineage information for the " +
               "corresponding species need to be collected.").format(len(diff_ids),len(uniq_spp_list)))
        spp_df1 = temp_df[temp_df["taxID"].isin(int_ids)].reset_index(drop = True)
        print("Retrieving information from downloaded dump files...")
        spp_df2_1 = dmp_df[dmp_df["taxID"].isin(diff_ids)].reset_index(drop = True)
        spp_df2_1 = spp_df2_1[["taxID", "Rank", "Scientific_Name", "Lineage"]]
        diff_ids_rmng = list(set(diff_ids) - set(spp_df2_1["taxID"]))
        print("Data retrieval completed!")
        print("Retrieving remaining information from the NCBI server...")
        spp_df2_2 = convert_taxid_to_lineage_main(diff_ids_rmng)[["taxID", "Rank", "Scientific_Name", "Lineage"]]
        spp_df2 = spp_df2_1.append(spp_df2_2, ignore_index = True)
        spp_df = spp_df1.append(spp_df2, ignore_index = True)
        spp_df = spp_df.drop_duplicates(subset = "taxID").reset_index(drop = True)
    else:
        print("All species-level taxIDs identified in the raw lineage summary dataframe can be found in the provided " +
              "source of taxIDs directly.")
        spp_df = temp_df
    assert all(spp_df["Rank"] == "species") == True, ("Non-species taxIDs are included in the summary dataframe for " +
                                                      "non-redundant species.")
    return spp_df

input_phylo_path = sys.argv[1]
input_added_path = sys.argv[2]
output_dir_path = sys.argv[3]
email = sys.argv[4]
NCBI_dmp_path = sys.argv[5]

if not os.path.exists(output_dir_path):
    os.makedirs(output_dir_path)
output_tip_lineage_path = os.path.join(output_dir_path, "tip_lineage.tsv")
output_phylo_uniq_spp_path = os.path.join(output_dir_path, "phylo_spp_lineage.tsv")
output_added_uniq_spp_path = os.path.join(output_dir_path, "added_spp_lineage.tsv")

phylo_tree = Phylo.read(input_phylo_path, "newick")
tip_labels = pd.DataFrame([clade.name for clade in phylo_tree.find_clades() if clade.name], columns = ["Tip_Label"])
temp = tip_labels["Tip_Label"].str.split("_", expand = True)
tip_labels["Accession"] = temp[1].str.strip() + "_" + temp[2].str.strip()
acc_list = list(filter(None, tip_labels["Accession"]))

added_file = open(input_added_path)
added_taxid_list = []
for line in added_file:
    line = line.strip()
    taxid = line.split(sep = "\t")[-1].strip()
    if taxid != '':
        added_taxid_list.append(taxid)
added_taxid_list = list(set(added_taxid_list))

if os.path.exists(NCBI_dmp_path):
    print("NCBI dump files have already been downloaded.")
    dump_acc_taxid = pd.read_pickle(os.path.join(NCBI_dmp_path, "assembly_to_taxID.pkl"))
    dump_taxid_lineage = pd.read_pickle(os.path.join(NCBI_dmp_path, "taxid_to_lineage.pkl"))
else:
    print("NCBI dump files not downloaded.")
    dump_acc_taxid = pd.DataFrame({"Accession": pd.Series(dtype = "str"), "taxID": pd.Series(dtype = "str")})
    dump_taxid_lineage = pd.DataFrame({"taxID": pd.Series(dtype = "str"), "Scientific_Name": pd.Series(dtype = "str"),
                                       "Rank": pd.Series(dtype = "str"), "Lineage": np.NaN,
                                       "CurrentID": pd.Series(dtype = "str")})

print("*** Step 1/6 ***")
print("Converting genome assembly accession numbers extracted from phylogenetic tip labels into taxIDs under NCBI " +
      "taxonomic classification system.")
print("1a: Retrieving information from downloaded dump files...")
phylo_tip_taxids_1 = dump_acc_taxid[dump_acc_taxid["Accession"].isin(acc_list)].reset_index(drop = True)
print("Data retrieval completed!")
print("1b: Retrieving remaining information from the NCBI server...")
acc_list_rmng = list(set(acc_list) - set(phylo_tip_taxids_1["Accession"]))
phylo_tip_taxids_2 = convert_acc_to_taxid_main(acc_list_rmng)
phylo_tip_taxids = phylo_tip_taxids_1.append(phylo_tip_taxids_2, ignore_index = True)
phylo_tip_taxids = phylo_tip_taxids.drop_duplicates(subset = ["Accession"]).reset_index(drop = True)

print("\n*** Step 2/6 ***")
print("Extracting full lineage for taxIDs found in the phylogeny. Generating the raw lineage summary dataframe.")
uniq_tip_taxids = list(filter(None, phylo_tip_taxids["taxID"].unique()))
print("2a: Retrieving information from downloaded dump files...")
phylo_summary_1 = dump_taxid_lineage[dump_taxid_lineage["taxID"].isin(uniq_tip_taxids)].reset_index(drop = True)
print("Data retrieval completed!")
print("2b: Retrieving remaining information from the NCBI server...")
uniq_tip_taxids_rmng = list(set(uniq_tip_taxids) - set(phylo_summary_1["taxID"]))
phylo_summary_2 = convert_taxid_to_lineage_main(uniq_tip_taxids_rmng)
phylo_summary = phylo_summary_1.append(phylo_summary_2, ignore_index = True)
print("Checking if all taxIDs are up to date:")
(phylo_summary, phylo_updated_taxids) = check_outdated_taxids(phylo_summary, "phylogeny")
phylo_summary = phylo_summary.drop_duplicates(subset = ["taxID"]).reset_index(drop = True)

print("\n*** Step 3/6 ***")
added_list_overlap = list(set(added_taxid_list) & set(phylo_summary["taxID"]))
added_list_rmng = list(set(added_taxid_list) - set(phylo_summary["taxID"]))
print(("Extracting full lineage for taxIDs found in the user-provided file. An overlap of {0} taxIDs are present in " +
       "the phylogeny summary table, lineage information of the remaining {1} taxIDs need to be collected to " +
       "generate the raw lineage summary dataframe:").format(len(added_list_overlap), len(added_list_rmng)))
added_summary_1 = phylo_summary[phylo_summary["taxID"].isin(added_list_overlap)]
print("3a: Retrieving information from downloaded dump files...")
added_summary_2_1 = dump_taxid_lineage[dump_taxid_lineage["taxID"].isin(added_list_rmng)].reset_index(drop = True)
print("Data retrieval completed!")
print("3b: Retrieving remaining information from the NCBI server...")
added_list_rmng_rmng = list(set(added_list_rmng) - set(added_summary_2_1["taxID"]))
added_summary_2_2 = convert_taxid_to_lineage_main(added_list_rmng_rmng)
added_summary_2 = added_summary_2_1.append(added_summary_2_2, ignore_index = True)
print("Checking if all taxIDs are up to date:")
(added_summary, added_updated_taxids) = check_outdated_taxids(added_summary_2, "user-provided file")
added_summary_2 = added_summary_2.drop_duplicates(subset = ["taxID"]).reset_index(drop = True)
added_summary = added_summary_1.append(added_summary_2, ignore_index = True)
added_summary = added_summary.drop_duplicates(subset = ["taxID"]).reset_index(drop = True)

print("\n*** Step 4/6 ***")
print("Generating the TSV file summarizing lineage information for all tips present in the phylogeny.")
# Need to update taxIDs in dataframe phylo_tip_taxids as well before merging them
for i, j in phylo_updated_taxids.items():
    temp_cond = phylo_tip_taxids["taxID"] == i
    temp_ind = phylo_tip_taxids.index[temp_cond]
    phylo_tip_taxids.loc[temp_ind, "taxID"] = j
tip_lineage = tip_labels.merge(phylo_tip_taxids, left_on = "Accession", right_on = "Accession", how = "inner")
tip_lineage = tip_lineage.merge(phylo_summary, left_on = "taxID", right_on = "taxID", how = "inner")
tip_rank_df = convert_lineage_to_ranks(tip_lineage)[["Tip_Label", "Accession", "taxID", "Scientific_Name",
                                                     "superkingdom", "phylum", "class", "order", "family",
                                                     "genus", "species_group", "species"]]
tip_rank_df["Scientific_Name"] = tip_rank_df["Scientific_Name"].str.replace("'", "")
tip_rank_df["Scientific_Name"] = tip_rank_df["Scientific_Name"].str.replace("\"", "")
tip_rank_df.to_csv(output_tip_lineage_path, sep = "\t", index = False)
print("File output completed!")

print("\n*** Step 5/6 ***")
print("Generating the TSV file summarizing lineage information for each unique species identified in the phylogeny.")
print("< \"pseudo\" species which indeed refer to higher levels (e.g. Enterobacteriaceae bacterium, taxID: 1849603, " +
      "assigned rank: species) are defined here as those lacking genus-level descriptions and are removed from " +
      "the summary file. >")
phylo_uniq_spp = list(filter(None, tip_rank_df["species"].unique()))
phylo_uniq_spp_lineage = uniq_spp_lineage(phylo_summary, phylo_uniq_spp, dump_taxid_lineage)
print("TaxIDs at different ranks are now stored in different columns.")
phylo_uniq_spp_rank_df = convert_lineage_to_ranks(phylo_uniq_spp_lineage)
phylo_uniq_spp_rank_df = phylo_uniq_spp_rank_df[["taxID", "Scientific_Name", "superkingdom", "phylum", "class",
                                                 "order", "family", "genus", "species_group", "species"]]
temp_filter = ~phylo_uniq_spp_rank_df["genus"].isna()
phylo_uniq_spp_rank_df = phylo_uniq_spp_rank_df[temp_filter].reset_index(drop = True)
phylo_uniq_spp_rank_df["Scientific_Name"] = phylo_uniq_spp_rank_df["Scientific_Name"].str.replace("'", "")
phylo_uniq_spp_rank_df["Scientific_Name"] = phylo_uniq_spp_rank_df["Scientific_Name"].str.replace("\"", "")
phylo_uniq_spp_rank_df.to_csv(output_phylo_uniq_spp_path, sep = "\t", index = False)
print("File output completed!")

print("\n*** Step 6/6 ***")
print("Generating the TSV file summarizing lineage information for each unique species identified in the " +
      "user-provided file.")
added_rank_df = convert_lineage_to_ranks(added_summary)
added_rank_df = added_rank_df[["taxID", "Scientific_Name", "superkingdom", "phylum", "class", "order", "family",
                               "genus", "species_group", "species"]]
added_uniq_spp = list(filter(None, added_rank_df["species"].unique()))
added_uniq_spp_lineage = uniq_spp_lineage(added_summary, added_uniq_spp, dump_taxid_lineage)
print("TaxIDs at different ranks are now stored in different columns.")
added_uniq_spp_rank_df = convert_lineage_to_ranks(added_uniq_spp_lineage)
added_uniq_spp_rank_df = added_uniq_spp_rank_df[["taxID", "Scientific_Name", "superkingdom", "phylum", "class",
                                                 "order", "family", "genus", "species_group", "species"]]
added_uniq_spp_rank_df["Scientific_Name"] = added_uniq_spp_rank_df["Scientific_Name"].str.replace("'", "")
added_uniq_spp_rank_df["Scientific_Name"] = added_uniq_spp_rank_df["Scientific_Name"].str.replace("\"", "")
added_uniq_spp_rank_df.to_csv(output_added_uniq_spp_path, sep = "\t", index = False)
print("File output completed!")
