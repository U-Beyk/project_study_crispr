'''
Classes for all the necessary tables. Parses the json-files and processes the content.

author: U.B.
'''

import pandas as pd
from abc import ABC, abstractmethod
from typing import override

FILEPATHS = {
    "crispr_regions": "crispr_cas_db/database_tables/region.json",
    "crispr_loci_regions": "crispr_cas_db/database_tables/crisprlocus_region.json",
    "crispr_loci": "crispr_cas_db/database_tables/crisprlocus.json",
    "cas_cluster": "crispr_cas_db/database_tables/clustercas.json",
    "sequence": "crispr_cas_db/database_tables/sequence.json"
}

class JsonTable(ABC):
    '''Abstract class for the Tables.'''
    def __init__(self, filepath: str) -> None:
        '''Initializes a table by reading the coresponding the json-file'''
        self._dataframe = pd.read_json(filepath)

    @property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        '''Returns the content of the json-file.'''
        pass


class CrisprRegions(JsonTable):
    '''Class for the region.json containing all the sequences of the CRISPR CAS database'''
    def __init__(self, filepath: str = FILEPATHS["crispr_regions"]) -> None:
        '''Initializes a CrisprRegions object'''
        super().__init__(filepath)

    @property
    @override
    def data(self) -> pd.DataFrame:
        '''Processes the content of the json-file and returns it as a dataframe'''
        df_reg_repeat_spacer = self._dataframe[(self._dataframe["region_category"] == 1) | (self._dataframe["region_category"] == 3)].reset_index(drop=True)
        df_reg_repeat_spacer["region_category"] = df_reg_repeat_spacer["region_category"].astype(str)
        df_reg_repeat_spacer["region_category"] = df_reg_repeat_spacer["region_category"].replace({"1.0":"Repeat", "3.0":"Spacer"})
        return df_reg_repeat_spacer

class CrisprLociRegions(JsonTable):
    '''Class for the crisprlocus_region.json containing all the sequences of the CRISPR CAS database'''
    def __init__(self, filepath: str = FILEPATHS["crispr_loci_regions"]) -> None:
        '''Initializes a CrisprLociRegions object'''
        super().__init__(filepath)

    @property
    @override
    def data(self) -> pd.DataFrame:
        '''Processes the content of the json-file and returns it as a dataframe'''
        df = self._dataframe.copy()
        df["crisprlocus_region_end"] = df["crisprlocus_region_start"] + df["crisprlocus_region_length"]
        return df
    
class CrisprLoci(JsonTable):
    '''Class for the crisprlocus.json containing all the sequences of the CRISPR CAS database'''
    def __init__(self, filepath: str = FILEPATHS["crispr_loci"]) -> None:
        '''Initializes a CrisprLoci object'''
        super().__init__(filepath)

    @property
    @override
    def data(self) -> pd.DataFrame:
        '''Processes the content of the json-file and returns it as a dataframe'''
        df = self._dataframe.copy()
        df["crisprlocus_end"] = df["crisprlocus_start"] + df["crisprlocus_length"]
        return df
    
class CasCluster(JsonTable):
    '''Class for the clustercas.json containing all the sequences of the CRISPR CAS database'''
    def __init__(self, filepath: str = FILEPATHS["cas_cluster"]) -> None:
        '''Initializes a CasCluster object'''
        super().__init__(filepath)

    @property
    @override
    def data(self) -> pd.DataFrame:
        '''Returns the content of the json-file as a dataframe'''
        return self._dataframe

class Sequence(JsonTable):
    '''Class for the sequence.json containing all the sequences of the CRISPR CAS database'''
    def __init__(self, filepath: str = FILEPATHS["sequence"]) -> None:
        '''Initializes a Sequence object'''
        super().__init__(filepath)

    @property
    @override
    def data(self) -> pd.DataFrame:
        '''Returns the content of the json-file as a dataframe'''
        return self._dataframe