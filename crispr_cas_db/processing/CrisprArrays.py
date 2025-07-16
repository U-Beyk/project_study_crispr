'''
Processes the CRISPR arrays from the CRISPR CAS dataset and stores the arrays with their corresponding subtypes in a fasta-file.

author: U.B.
'''

import os
import pandas as pd
from crispr_cas_db.processing.CrisprCasDataset import CrisprCasDataset

class CrisprArrays:
    '''A class processing and representing the CRISPR arrays'''
    def __init__(self) -> None:
        '''Initializes a CrisprArraySeq object with their arrays, repeats and spacers'''
        crispr_cas_dataset = CrisprCasDataset()
        self._arrays = crispr_cas_dataset.dataset
        self._repeats = None
        self._spacers = None
        self._subtyped_repeats = None

    @property
    def arrays(self) -> pd.DataFrame:
        '''Gets and returns the arrays as a dataframe'''
        return self._arrays

    @property
    def repeats(self) -> pd.DataFrame:
        '''Gets the repeats from the specified dataframe'''
        if self._repeats is None:
            self._repeats = self._arrays[self._arrays["region_category"] == "Repeat"].reset_index(drop=True)
        return self._repeats
    
    @property
    def subtyped_repeats(self) -> pd.DataFrame:
        '''Gets a list of repeats which can be assigned to a subtype'''
        if self._subtyped_repeats is None:
            valid_repeats = self.repeats.dropna(subset=["clustercas_class"])
            self._subtyped_repeats = valid_repeats.groupby("region_sequence")["clustercas_class"].apply(lambda x: list(set(x))).to_dict()
        return self._subtyped_repeats

    @property
    def spacers(self) -> pd.DataFrame:
        '''Gets the spacers from the specified dataframe'''
        if self._spacers is None:
            self._spacers = self._arrays[self._arrays["region_category"] == "Spacer"].reset_index(drop=True)
        return self._spacers
    
    def save_repeats_fasta(self, filename: str) -> None:
        '''Saves all the unique repeats in a fasta'''
        os.makedirs("crispr_cas_db/fasta_files", exist_ok=True)
        unique_sequences = self.repeats["region_sequence"].unique()
        fasta_filename = f"crispr_cas_db/fasta_files/{filename}.fasta"
        with open(fasta_filename, "w") as fasta_file:
            for index, sequence in enumerate(unique_sequences, start=1):
                subtypes = ",".join(self.subtyped_repeats.get(sequence, ["Undetermined"]))
                fasta_file.write(f">sequence_{index}|subtype:{subtypes}\n{sequence}\n")
        print(f"Fasta-file: {filename} was created")
    
    def check_number_unique_seq(self) -> None:
        '''Checks and prints the total number of all unique repeat and spacer sequences'''
        print(f"The total number of unique repeats is {len(self.repeats["region_sequence"].unique())}")
        print(f"The total number of unique spacers is {len(self.spacers["region_sequence"].unique())}", "\n")