'''
Contains all classes creating the CRISPR CAS dataset used for further processing from the database.

author: U.B.
'''

import pandas as pd
from crispr_cas_db.processing.JsonTables import CrisprLoci, CrisprLociRegions, CrisprRegions, CasCluster, Sequence

class CrisprCasDataset:
    '''Class representing the final CRISPR CAS dataset'''
    def __init__(self) -> None:
        '''INitializes a CrisprCasDAtaset object from the CRISPR dataset and the Cas dataset'''
        crispr_arrays, cas_clusters = CrisprDataset(), CasDataset()
        self._arrays = crispr_arrays.dataset
        self._cas_clusters = cas_clusters.dataset
        self._dataset = None

    @property
    def dataset(self) -> pd.DataFrame:
        '''Gets and returns the final dataset'''
        if self._dataset is None:
            self._dataset = self._arrays.copy()
            seq_class_dict  = dict(zip(self._cas_clusters["clustercas_sequence"], self._cas_clusters["clustercas_class"]))
            self._dataset["clustercas_class"] = self._dataset["crisprlocus_sequence"].map(seq_class_dict)
        return self._dataset
    
class CasDataset:
    '''Class representing the CAS dataset'''
    def __init__(self) -> None:
        '''Initializes a CasDataset object'''
        cas_clusters, sequence = CasCluster(), Sequence()
        self._cas_clusters = cas_clusters.data
        self._sequence = sequence.data
        self._dataset = None
    
    @property
    def dataset(self) -> pd.DataFrame:
        '''Gets and returns the final CAS dataset'''
        if self._dataset is None:
            self._dataset = self._merge_and_process()
        return self._dataset

    def _merge_sequences_cas_clusters(self) -> pd.DataFrame:
        '''Merges the CAS cluster and sequence dataframes and returns the merged dataframe'''
        df_merged =  self._cas_clusters.merge(
            self._sequence[["sequence_id", "sequence_strain"]], 
            left_on="clustercas_sequence", 
            right_on="sequence_id", 
            how="inner"
        ).drop(columns=["sequence_id"])
        return df_merged

    def _merge_and_process(self) -> pd.DataFrame:
        '''Filters the sequences out with one unique subtype from a single strain and returns the resulting dataframe'''
        df_merged = self._merge_sequences_cas_clusters()
        df_processed = df_merged.groupby("clustercas_sequence").filter(
            lambda x: len(x["clustercas_class"].unique()) == 1
            and x["clustercas_class"].iloc[0] != "CAS" 
            and pd.notna(x["clustercas_class"].iloc[0])
        ).groupby("sequence_strain").filter(
            lambda x: len(x["clustercas_class"].unique()) == 1
        )
        return df_processed
    
class CrisprDataset:
    '''Class representing the CRISPR dataset'''
    def __init__(self) -> None:
        '''Initializes a CrisprDataset'''
        regions, loci_regions, loci = CrisprRegions(), CrisprLociRegions(), CrisprLoci()
        self._regions = regions.data
        self._loci_regions = loci_regions.data
        self._loci = loci.data
        self._dataset = None

    @property
    def dataset(self) -> pd.DataFrame:
        '''Gets and returns the final CRISPR dataset'''
        if self._dataset is None:
            self._dataset = self._merge_and_process()
        return self._dataset 

    def _raw_crispr_array(self) -> pd.DataFrame:
        '''Gets the CRISPR arrays from the region, crisprlocus region and crisprlocus datatframes by merging them and return the resulting datatframe'''
        df_first_merge = self._loci_regions.merge(self._regions, left_on="crisprlocus_region_region", right_on="region_id", how="inner")
        df_second_merged = self._loci.merge(df_first_merge, left_on="crisprlocus_id", right_on="crisprlocus_region_crisprlocus", how="inner")
        return df_second_merged
    
    def _filter_evlvl4(self, df: pd.DataFrame) -> pd.DataFrame:
        '''Filters all evidence level 4 entries of the dataframe and returns the resulting dataframe'''
        if "crisprlocus_evidencelevel" in df.columns:
            df_evl4 = df[df["crisprlocus_evidencelevel"].isin([4])].reset_index(drop=True)
            return df_evl4
        else:
            print("Column 'crisprlocus_evidencelevel' does not exist!")
            return df
        
    def _filter_non_ambiguous_nt(self, df: pd.DataFrame) -> pd.DataFrame:
        '''Gets all sequences with non ambiguous nucleotides and returns the resulting dataframe'''
        valid_nucleotides = {'A', 'C', 'G', 'T'}
        df_valid_nucleotides = df[df['region_sequence'].apply(lambda seq: all(base in valid_nucleotides for base in seq))]
        return df_valid_nucleotides
    
    def _filter_correct_direction(self, df: pd.DataFrame) -> pd.DataFrame:
        '''Gets all entries of a dataframe which have an unambiguous sequence direction and returns the resulting dataframe'''
        df_direction = df[
            ((df["crisprlocus_orientation"] == 1.0) & (df["crisprlocus_potentialorientation"] == 1.0)) | 
            ((df["crisprlocus_potentialorientation"] == 2.0) & (df["crisprlocus_orientation"] == 2.0))
        ]
        return df_direction
    
    def _make_reverse_complementary(self, sequence: str) -> str:
        '''Creates and returns the reverse complementary of the entered sequence'''
        rev_sequence = sequence[::-1]
        complementary_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        rev_comp_seq = "".join([complementary_dict.get(nucleotide) for nucleotide in rev_sequence])
        return rev_comp_seq

    def _apply_reverse_complement(self, row: pd.Series) -> str:
        '''Creates the reverse complement of a row if the sequence is in reverse orientation'''
        if row["crisprlocus_potentialorientation"] == 2.0 and row["crisprlocus_orientation"] == 2.0:
            return self._make_reverse_complementary(row["region_sequence"])
        else:
            return row["region_sequence"]
        
    def _convert_reverse_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        '''Changes all sequences which are in reverse to the reverse complement in the specified dataframe'''
        df = df.copy()
        df["region_sequence"] = df.apply(self._apply_reverse_complement, axis=1)
        return df
    
    def _remove_extra_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        '''Removes redundant columns in the dataframe and returns the resulting dataframe'''
        columns_to_remove = [
            "crisprlocus_blastscore",
            "crisprlocus_evidencelevelreeval",
            "crisprlocus_spacerconservation",
            "crisprlocus_drconservation",
            "crisprlocus_trusted"
        ]
        cols_to_drop = [col for col in columns_to_remove if col in df.columns]
        return df.drop(columns=cols_to_drop).copy()
    
    def _merge_and_process(self) -> pd.DataFrame:
        '''Gets and processes the CRISPR array'''
        df = self._raw_crispr_array()
        df_lv4 = self._filter_evlvl4(df)
        df_non_amb_nt = self._filter_non_ambiguous_nt(df_lv4)
        df_cor_direc = self._filter_correct_direction(df_non_amb_nt)
        df_rev_comp = self._convert_reverse_sequences(df_cor_direc)
        df_rem_col = self._remove_extra_columns(df_rev_comp)
        return df_rem_col