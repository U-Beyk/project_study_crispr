'''
Processes the sequences further by building the mature CRISPR RNA and storing them with their correspponding subtypes as a fasta-file.

author: U.B.
'''

from crispr_cas_db.processing.CrisprRNA import MatureCrRNAs, CrRNA, CrisprArray
from crispr_cas_db.processing.CrisprArrays import CrisprArrays

class CrisprRNAs:
    '''Class representing all mature CRISPR RNAs in a list with their corresponding subtyps.'''
    def __init__(self, crispr_data: CrisprArrays) -> None:
        '''Initializes a CrisprRNAs object representing all mature CRISPR RNAs'''
        self._crispr_data = crispr_data
        self._crRNAs = None

    @property
    def crRNAs(self) -> list[CrRNA]:
        '''Property representing a list of CRIPR RNAs'''
        if self._crRNAs is None:
            self._crRNAs = {
                CrRNA(sequence=seq, subtype=",".join(sorted(subtypes)))
                for seq, subtypes in self._group_crRNAs().items()
                }
        return self._crRNAs
    
    def _assemble_crRNAs(self) -> list[CrRNA]:
        '''Assembles and builds all mature CRIPSR RNAs from the CRISPR arrays'''
        assembled_crRNAs, crispr_arrays = [], self._crispr_data.arrays
        for _, locus in crispr_arrays.groupby("crisprlocus_id"):
            array = CrisprArray(locus)
            if array.subtype != None:
                mature_crRNAs = MatureCrRNAs(array)
                assembled_crRNAs.extend(mature_crRNAs.crRNAs)
        return assembled_crRNAs
    
    def _group_crRNAs(self) -> dict[str, set[str]]:
        '''Removes duplicates and groupes the sequences with the same sequence, but different subtypes'''
        unique_sequences: dict[str, set[str]]  = dict()
        for crRNA in self._assemble_crRNAs():
            if crRNA.sequence in unique_sequences:
                unique_sequences[crRNA.sequence].add(crRNA.subtype)
            else:
                unique_sequences[crRNA.sequence] = {crRNA.subtype}
        return unique_sequences
    
    def save_crRNA_fasta(self, filename: str) -> None:
        '''Saves the processed and final CRISPR RNAs in a fasta-file with their corresponding subtypes'''
        fasta_filename = f"crispr_cas_db/fasta_files/{filename}.fasta"
        with open(fasta_filename, "w") as fasta_file:
            for index, crRNA in enumerate(self.crRNAs):
                fasta_file.write(f">sequence_{index}|subtype:{crRNA.subtype}\n{crRNA.sequence}\n")
        print(f"Fasta-file: {filename} was created")