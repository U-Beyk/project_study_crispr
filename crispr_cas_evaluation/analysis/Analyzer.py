'''
Contains all classes with the dataframes to be analyzed and its methods to filter necessary information.

author: U.B.
'''

from typing import Generic, TypeVar
from dataclasses import dataclass
from crispr_cas_evaluation.predictions.RNADataFrameAssembler import RNADataFrameAssembler
from crispr_cas_evaluation.predictions.RNARecordsAssembler import RNAmotiFoldRecordsAssembler, RNAmotiCesRecordsAssembler, CrRNAmotiFoldRecordsAssembler, CrRNAmotiCesRecordsAssembler

T = TypeVar("T", bound=RNAmotiFoldRecordsAssembler)
U = TypeVar("U", bound=RNAmotiCesRecordsAssembler)

@dataclass
class AnalyzerConfig:
    '''Dataclass containing the information of the fasta file and the RNAmotiFold and RNAHeliCes/RNAmotiCes prediction files.'''
    fasta_path: str
    motifold_csv_path: str
    motices_csv_path: str

class Analyzer(Generic[T, U]):
    '''Generic class for the different RNA assemblies.'''
    motifoldassembler_class: type[T] = RNAmotiFoldRecordsAssembler
    moticesassembler_class: type[U] = RNAmotiCesRecordsAssembler

    def __init__(self, config: AnalyzerConfig) -> None:
        '''Initializes an Analyzer object with an RNAmotiFold assembly and an RNAHeliCes/RNAmotiCes assembly.'''
        self._rna_motifold = self.motifoldassembler_class(
            RNADataFrameAssembler(config.fasta_path, config.motifold_csv_path)
        )
        self._rna_motices = self.moticesassembler_class(
            RNADataFrameAssembler(config.fasta_path, config.motices_csv_path)
        )

    @property
    def rna_motifold_assembly(self) -> T:
        '''Returns the RNAmotiFold assembly.'''
        return self._rna_motifold

    @property
    def rna_motices_assembly(self) -> U:
        '''Returns the RNAHeliCes/RNAmotiCes assembly.'''
        return self._rna_motices

    def filter_motifold_mfe_range(self) -> T:
        '''Filters the predictions with a mfe range in the RNAmotiFold assembly.'''
        return self._rna_motifold.filter_records(mfe_range=True)

    def filter_motices_mfe_range(self) -> U:
        '''Filters the predictions with a mfe range in the RNAHeliCes/RNAmotiCes assembly.'''
        return self._rna_motices.filter_records(mfe_range=True)

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (f"{self.__class__.__name__}("
                f"rna_motifold_assembly={self._rna_motifold}, "
                f"rna_motices_assembly={self._rna_motices})")


class CRISPRAnalyzer(Analyzer[CrRNAmotiFoldRecordsAssembler, CrRNAmotiCesRecordsAssembler]):
    '''Class for the CRISPR RNA assemblies.'''
    motifoldassembler_class = CrRNAmotiFoldRecordsAssembler
    moticesassembler_class = CrRNAmotiCesRecordsAssembler

    def __init__(self, config: AnalyzerConfig):
        '''Initializes a CRISPRAnalyzer object containing the RNA record assemblies.'''
        super().__init__(config)

    def filter_motifold_by_subtype(self, subtype: str, mfe_range: bool = False) -> CrRNAmotiFoldRecordsAssembler:
        '''Filters the RNAmotiFold records corresponding to a specified CRISPR subtype.
        Also allows to filter predictions within a mfe range.'''
        return self._rna_motifold.filter_records(
            condition=lambda record: subtype in record.subtypes,
            mfe_range=mfe_range
        )

    def filter_motices_by_subtype(self, subtype: str, mfe_range: bool = False) -> CrRNAmotiCesRecordsAssembler:
        '''Filters the RNAHeliCes/RNAmotiCes records corresponding to a specified CRISPR subtype.
        Also allows to filter predictions within a mfe range.'''
        return self._rna_motices.filter_records(
            condition=lambda record: subtype in record.subtypes,
            mfe_range=mfe_range
        )