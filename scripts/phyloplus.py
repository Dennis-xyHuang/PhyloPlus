import time
import pandas as pd
import numpy as np
from Bio import Entrez
from math import log10, ceil
from urllib.error import URLError, HTTPError
from http.client import IncompleteRead


taxonomic_ranks_cols = {"superkingdom": "TaxSuperKingdom", "phylum": "TaxPhylum", "class": "TaxClass",
                        "order": "TaxOrder", "family": "TaxFamily", "genus": "TaxGenus",
                        "species group": "TaxSpeciesGroup", "species": "TaxSpecies"}


def read_dump(file_path, modify_rows=None, modify_cols=None):
    """
    Read dump files downloaded from NCBI FTP server and perform basic pre-processing.

    :param str file_path: path for NCBI dump file to be read.
    :param modify_rows: dictionary {column_index: keyword} indicating the filter to apply to the raw dataframe
        row-wise, default None keeps all rows.
    :type modify_rows: dict | None
    :param modify_cols: dictionary (column_index: new_name) providing rename information for columns, default None
        keeps all columns and makes no changes to column names.
    :type modify_cols: dict | None
    :return: processed dataframe of input dump file.
    :rtype: pd.DataFrame
    """
    dump_df = pd.read_csv(file_path, sep="|", header=None, engine="python")
    dump_df = dump_df.dropna(axis=1, how="all")
    dump_df_str_cols = dump_df.columns[dump_df.dtypes == "object"]
    dump_df.iloc[:, dump_df_str_cols] = dump_df.iloc[:, dump_df_str_cols].apply(lambda x: x.str.strip())
    if modify_rows:
        for index, keyword in modify_rows.items():
            dump_df = dump_df.loc[dump_df.loc[:, index] == keyword, :]
    if modify_cols:
        dump_df = dump_df.iloc[:, list(modify_cols.keys())]
        dump_df = dump_df.rename(columns=modify_cols)
    dump_df = dump_df.reset_index(drop=True)
    return dump_df


def build_lineage(taxids_str, mapping):
    """
    Takes in a string of a series of taxonomy IDs separated by ' ', generate a NCBI LineageEx-like list of dictionaries.
    Compared to LineageEx information returned by Entrez query, the value for key `TaxId` is an integer rather than a
    string, and the `ScientificName` item is not included.

    :param str taxids_str: a string of taxonomy IDs, separated by ' '.
    :param dict mapping: a dictionary of which keys are taxonomy IDs and values are corresponding LineageEx-like
        dictionaries; items for each LineageEx-like dictionary include {"TaxId": int, "Rank": str}.
    :return: a LineageEx-like list of dictionaries.
    :rtype: list of dict
    """
    if taxids_str == "":
        taxid_list = []
    else:
        taxid_list = [int(taxid) for taxid in taxids_str.split(" ")]
    temp_list = []
    for taxid in taxid_list:
        temp_dict = mapping.get(taxid, {"TaxId": taxid, "Rank": "NotFound"})
        temp_list.append(temp_dict)
    return temp_list


def convert_lineage_to_ranks(lineage_df, lineage_col, taxrank_cols_map):
    """
    Takes in a dataframe in which a column contains Entrez LineageEx or LineageEx-like list of dictionaries.
    Split lineage information and store them into individual taxonomic rank columns.

    :param pd.DataFrame lineage_df: a taxonomic summary dataframe that contains the 'Lineage' column which stores
        Entrez LineageEx (retrieved via Entrez search) or LineageEx-like (created from dump file) list of dictionaries
        for each taxonomy ID.
    :param str lineage_col: column name of which stores lineage information.
    :param dict taxrank_cols_map: dictionary where keys are taxonomic ranks reported by Entrez taxonomic query, and
        values are corresponding column names to be appended to the dataframe.
    :return: a modified summary dataframe of the input, information stored in column 'Lineage' is split and stored in
        individual taxonomic rank columns described in `taxrank_cols_map`.
    :rtype: pd.DataFrame
    """
    def extract_rank(lineage_dict, taxo_rank):  # Return TaxID at a given rank for a complete lineage record
        val = np.nan
        counter = 0
        while counter <= len(lineage_dict) - 1:
            if lineage_dict[counter]["Rank"] == taxo_rank:
                val = int(lineage_dict[counter]["TaxId"])  # TaxIDs are strings if lineage_dict come from Entrez search.
                break
            counter += 1
        return val

    for rank, colname in taxrank_cols_map.items():
        lineage_df.loc[:, colname] = lineage_df.loc[:, lineage_col].apply(extract_rank, args=(rank,))
    return lineage_df


def finalize_taxonomy_summary(lineage_df_with_rank, query_col, rank_col, taxrank_cols_map):
    """
    Takes in a dataframe in which a column describes taxonomic rank of the query TaxID, and individual taxonomic
    rank columns described in `taxrank_cols_map` contain TaxID of ancestral nodes for the query taxonomy ID. Finalize
    the dataframe so that the taxonomy ID itself is also included in the corresponding taxonomic column, if applicable.

    :param pd.DataFrame lineage_df_with_rank: a taxonomic summary dataframe that contains a column describing taxonomic
        ranks of query taxonomy IDs, and individual taxonomic rank columns.
    :param str query_col: column name of which describes query taxonomy IDs.
    :param str rank_col: column name of which describes taxonomic ranks of query taxonomy IDs.
    :param dict taxrank_cols_map: dictionary where keys are taxonomic ranks reported by Entrez taxonomic query, and
        values are corresponding column names in the dataframe.
    :return: a modified summary dataframe where the query taxonomy ID itself is also included in the corresponding
        taxonomic column, if applicable.
    :rtype: pd.DataFrame
    """
    for taxo_rank, colname in taxrank_cols_map.items():
        temp_condition = lineage_df_with_rank.loc[:, rank_col] == taxo_rank
        lineage_df_with_rank.loc[temp_condition, colname] = lineage_df_with_rank.loc[temp_condition, query_col]
    return lineage_df_with_rank


def keep_trying(func, max_tries, sleep=None, *args):
    """
    Keep trying calling a function until success or if maximum limit has been reached.

    :param function func: the function being called.
    :param int max_tries: the maximum limit of calling `func`.
    :param sleep: pause time between each call, default None does not sleep.
    :type sleep: float | None
    :param tuple args: positional arguments passed into `func`.
    :return: the return type of `fun` or None if reaches maximum limit.
    """
    count = 0
    while True:
        count += 1
        try:
            return func(*args)
        except (IndexError, URLError, HTTPError, IncompleteRead):
            if count >= max_tries:
                print("Maximum limit reached when calling the following function: " + " " * 15)
                print("--> {0}({1}).".format(func.__name__, ",".join([str(i) for i in args])))
                break
            if sleep is not None:
                time.sleep(sleep)
            else:
                break


def convert_acc_to_taxid_main(accessions, email, assembly_col="Accession", taxid_col="TaxID"):
    """
    Takes in an iterable of NCBI genome assembly accession numbers, returns a summary dataframe and a list of
    abnormal input.

    :param accessions: an iterable of NCBI assembly accession numbers.
    :type accessions: list | set
    :param str email: email address for Entrez search.
    :param str assembly_col: column name of which to store genome assembly accession numbers, default 'Accession'.
    :param str taxid_col: column name of which to store taxonomy IDs, default 'TaxID'.
    :return: [1] summary dataframe with columns for genome assembly accession numbers and their corresponding NCBI
        taxonomy IDs; [2] a list of abnormal input assembly accession numbers.
    :rtype: (pd.DataFrame, list)
    """
    def from_acc_to_taxid(acc):
        Entrez.email = email
        handle = Entrez.esearch(db="assembly", term=acc)
        record = Entrez.read(handle)
        id_key = record['IdList'][0]
        summ_handle = Entrez.esummary(db="assembly", id=id_key, report="full")
        summ_record = Entrez.read(summ_handle, validate=False)
        taxonomy_id = summ_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        taxonomy_id = np.nan if taxonomy_id == "" else int(taxonomy_id)
        return taxonomy_id
    abnormal_input = []
    acc_taxid = {}
    start_time = time.time()
    ave_proc_time = 0
    if len(accessions) == 0:
        print("", end="\r")
    else:
        for count, accession in enumerate(accessions):
            print("> The {0:^{width}d}th assembly being processed. Progress: {1:^6.2f}%; ETR: {2:^6.2f} min. <".format(
                count + 1, (count + 1) * 100 / len(accessions), (len(accessions) - count - 1) * ave_proc_time,
                width=ceil(log10(len(accessions) + 1))), end="\r")
            temp_taxid = keep_trying(from_acc_to_taxid, 10, 5, accession)
            if temp_taxid is None:
                abnormal_input.append(accession)
            else:
                acc_taxid[accession] = temp_taxid
            end_time = time.time()
            proc_time = (end_time - start_time) / 60
            ave_proc_time = proc_time / (count + 1)
            time.sleep(0.35)
        print("\n", end="\r")
    print("Data retrieval completed!")
    df_acc = []
    df_taxid = []
    for entrez_acc, entrez_taxid in acc_taxid.items():
        df_acc.append(entrez_acc)
        df_taxid.append(entrez_taxid)
    summary_df = pd.DataFrame({assembly_col: df_acc, taxid_col: df_taxid})
    return summary_df, abnormal_input


def convert_taxid_to_lineage_main(unique_taxids, email, taxrank_cols_map, taxid_col="TaxID",
                                  current_taxid_col="CurrentTaxID", sciname_col="ScientificName",
                                  rank_col="TaxonomicRank"):
    """
    Takes in an iterable of non-redundant NCBI taxonomy IDs, returns a summary dataframe and a list of abnormal input.

    :param unique_taxids: an iterable of non-redundant NCBI taxonomy IDs.
    :type unique_taxids: list | set
    :param str email: email address for Entrez search.
    :param str taxid_col: column name of which to store genome taxonomy IDs, default 'TaxID'.
    :param str current_taxid_col: column name of which to store current taxonomy IDs, default 'CurrentTaxID'.
    :param str sciname_col: column name of which to store scientific names, default 'ScientificName'.
    :param str rank_col: column name of which to store taxonomic ranks, default 'TaxonomicRank'.
    :param dict taxrank_cols_map: dictionary where keys are taxonomic ranks reported by Entrez taxonomic query, and
        values are corresponding column names in the dataframe.
    :return: [1] summary dataframe with columns for NCBI taxonomy IDs, the corresponding current IDs, scientific names,
        taxonomic ranks, and taxonomy IDs for ancestral nodes at specified rank levels indicated by `taxrank_cols_map`;
        [2] a list of abnormal input taxonomy IDs.
    :rtype: (pd.DataFrame, list)
    """
    def from_taxid_to_lineage(taxonomy_id):
        info_dict = {}
        Entrez.email = email
        handle = Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml")
        record = Entrez.read(handle)
        info_dict["Lineage"] = record[0].get("LineageEx", [])  # Taxonomy ID "1" does not have LineageEx data
        info_dict["ScientificName"] = record[0]["ScientificName"]
        info_dict["Rank"] = record[0]["Rank"]
        info_dict["CurrentTaxID"] = int(record[0]["TaxId"])
        return info_dict

    abnormal_input = []  # In case some user input contain taxonomy IDs outside NCBI database
    uniq_sci_name = {}
    uniq_rank = {}
    uniq_lineage = {}
    uniq_curr_id = {}
    start_time = time.time()
    ave_proc_time = 0
    if len(unique_taxids) == 0:
        print("", end="\r")
    else:
        for count, taxid in enumerate(unique_taxids):
            print("> The {0:^{width}d}th TaxID being processed. Progress: {1:^6.2f}%; ETR: {2:^6.2F} min. <".format(
                count + 1, (count + 1) * 100 / len(unique_taxids), (len(unique_taxids) - count - 1) * ave_proc_time,
                width=ceil(log10(len(unique_taxids) + 1))), end="\r")
            temp_dict = keep_trying(from_taxid_to_lineage, 10, 5, taxid)
            if temp_dict is None:
                abnormal_input.append(taxid)
            else:
                uniq_curr_id[taxid] = temp_dict["CurrentTaxID"]
                uniq_sci_name[taxid] = temp_dict["ScientificName"]
                uniq_rank[taxid] = temp_dict["Rank"]
                uniq_lineage[taxid] = temp_dict["Lineage"]
            end_time = time.time()
            proc_time = (end_time - start_time) / 60
            ave_proc_time = proc_time / (count + 1)
            time.sleep(0.35)
        print("\n", end="\r")
    print("Data retrieval completed!")
    summary_df = pd.DataFrame({taxid_col: list(uniq_sci_name.keys()), current_taxid_col: np.nan,
                               sciname_col: np.nan, rank_col: np.nan, "Lineage": np.nan})
    summary_df.loc[:, current_taxid_col] = summary_df.loc[:, taxid_col].apply(lambda x: uniq_curr_id.get(x))
    summary_df.loc[:, sciname_col] = summary_df.loc[:, taxid_col].apply(lambda x: uniq_sci_name.get(x))
    summary_df.loc[:, rank_col] = summary_df.loc[:, taxid_col].apply(lambda x: uniq_rank.get(x))
    summary_df.loc[:, "Lineage"] = summary_df.loc[:, taxid_col].apply(lambda x: uniq_lineage.get(x))
    summary_df = convert_lineage_to_ranks(summary_df, "Lineage", taxrank_cols_map)
    summary_df = summary_df.loc[:, [taxid_col, current_taxid_col, sciname_col, rank_col, *taxrank_cols_map.values()]]
    summary_df = finalize_taxonomy_summary(summary_df, "TaxID", "TaxonomicRank", taxrank_cols_map)
    summary_df = summary_df.reset_index(drop=True)
    return summary_df, abnormal_input
