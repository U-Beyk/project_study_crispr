'''
Contains all classes that assemble different types of RNA records into a dictionary.

author: U.B.
'''

import pandas as pd
from abc import ABC, abstractmethod
from typing import TypeVar, Generic, Callable, Self
from collections import Counter
from crispr_cas_evaluation.predictions.Prediction import RNAmotiFoldPrediction, RNAmotiCesPrediction
from crispr_cas_evaluation.predictions.RNARecord import RNARecord, RNAmotiFoldRecord, RNAmotiCesRecord, CrRNAmotiFoldRecord, CrRNAmotiCesRecord
from crispr_cas_evaluation.predictions.Visualization import BarChart, Histogram, ViolinPlot 

T = TypeVar("T", bound=RNARecord)

#TO-DO: Implement the visualization classes taking no dataframe types and instead a dict or count of motifs/motices for consistency in types.
class RNARecordsAssembler(ABC, Generic[T]):
    '''Generic class for an RNA record assembler.'''
    def __init__(self, rna_dataframe: pd.DataFrame):
        '''Initializes the assembler via specified dataframe.'''
        self._rna_dataframe = rna_dataframe
        self._rna_sequences: dict[str, T] = {}

    @property
    def rna_sequences(self) -> dict[str, T]:
        '''Assembles the RNA sequences if there are none in the dictionary, otherwise returns them.'''
        if not self._rna_sequences:
            self._rna_sequences = self._assemble_rna_sequences()
        return self._rna_sequences
    
    @property
    def lowest_mfe_value(self) -> float:
        '''Returns the lowest mfe value as a float.'''
        value = 0
        for record in self.rna_sequences.values():
            if record.mfe_value < value:
                value = record.mfe_value
        return value
    
    @property
    def motifs_set(self) -> set[str]:
        '''Returns a set of unique motifs of all RNA sequences.'''
        motifs = set()
        for prediction in self.rna_sequences.values():
            motifs.update(prediction.motifs_set)
        return motifs
    
    @property
    def motifs_count(self) -> dict[str, int]:
        '''Counts all unique motifs for each sequence. Also counts all sequences without any motifs.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.motifs_set)
        return dict(motifs_counter)
    
    @property
    def sequence_number(self) -> int:
        '''Gets and returns the number of seuquences.'''
        return len(self.rna_sequences)
    
    @property
    def prediction_number(self) -> int:
        '''Counts and returns the total number of predictions.'''
        number = 0
        for sequence in self.rna_sequences.values():
            number += sequence.prediction_number
        return number
    
    @property
    def distance_to_lowest_all_motifs(self) -> pd.DataFrame:
        '''Creates a dtaframe with the distance to the lowest mfe for each motif.'''
        distance_motifs = []
        for rna_seq in self.rna_sequences.values():
            for motif_data in rna_seq.distance_to_lowest_all_motifs:
                for motif, value in motif_data.items():
                    row = {"Motifs": motif, "Distance to mfe": value}
                    distance_motifs.append(row)
        return pd.DataFrame(distance_motifs)
    
    @property
    def median_distance_to_lowest_all_motifs(self) -> pd.Series:
        '''Calculates a series with the median of the distance to the lowest mfe for each motif.'''
        return self.distance_to_lowest_all_motifs.groupby("Motifs")["Distance to mfe"].median()
    
    def filter_records(self, condition: Callable[[T], bool] | None = None, mfe_range: bool = False) -> Self:
        '''Filters the records by a specified condition and optionally gets the predictions in a specified mfe range.
        Returns a new instance of the corresponding class.'''
        mfe_threshold = self.lowest_mfe_value * 0.1 if mfe_range else None
        new_sequences = {}
        for key, record in self.rna_sequences.items():
            filtered_record = record.filter_predictions(mfe_threshold) if mfe_threshold else record
            if condition is None or condition(filtered_record):
                new_sequences[key] = filtered_record
        return self._set_sequences(new_sequences)

    def _set_sequences(self, sequences: dict[str, T]) -> Self:
        '''Creates a new instance of the corresponding class, with a specified set of RNA sequences.'''
        new_instance = self.__class__(self._rna_dataframe)
        new_instance._rna_sequences = sequences
        return new_instance

    @abstractmethod
    def _assemble_rna_sequences(self):
        '''Assembles the RNA sequences.'''
        pass

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        lines = [f"{key}: {value}" for key, value in self.rna_sequences.items()]
        return "\n".join(lines)
    
U = TypeVar("U", bound=RNAmotiFoldRecord)

class RNAmotiFoldRecordsAssembler(RNARecordsAssembler[U], Generic[U]):
    '''Generic class for RNAmotiFold record assemblies.'''
    record_class: type[U] = RNAmotiFoldRecord

    def _assemble_rna_sequences(self):
        '''Assembles the RNA sequences with RNAmotiFold records.'''
        rna_sequences: dict[str, U] = {}
        for row in self._rna_dataframe:
            rna_sequence = rna_sequences.setdefault(row.ID, self.record_class(row.ID, row.sequence))
            rna_sequence.add_prediction(RNAmotiFoldPrediction(row.mfe, row.motBracket, row.Class))
        return rna_sequences
    
    def visualize_as_barchart(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a barchart.'''
        barchart = BarChart(
            description,
            self.motifs_count,
            self.sequence_number
        )
        barchart.save_plot(filepath)

    def visualize_as_violinplot(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a violinplot.'''
        plot_data = self.distance_to_lowest_all_motifs
        violinplot = ViolinPlot(
            description,
            plot_data,
            self.prediction_number
        )
        violinplot.save_plot(filepath)


class CrRNAmotiFoldRecordsAssembler(RNAmotiFoldRecordsAssembler[CrRNAmotiFoldRecord]):
    '''Class for RNAmotiFold assemblies from CRISPR RNA.'''
    record_class = CrRNAmotiFoldRecord 

    @property
    def unique_subtypes(self) -> set:
        '''Gets and returns all unique subtypes of the CRISPR RNA.'''
        subtypes = set()
        for rna_sequence in self.rna_sequences.values():
            subtypes.update(rna_sequence.subtypes)
        return subtypes

V = TypeVar("V", bound=RNAmotiCesRecord)

class RNAmotiCesRecordsAssembler(RNARecordsAssembler[V], Generic[V]):
    '''Generic class for RNAHeliCes/RNAmotiCes record assemblies.'''
    record_class: type[V] = RNAmotiCesRecord

    @property
    def potential_motifs_count(self) -> dict:
        '''Gets and returns all potential motifs computed via positional abstraction.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.potential_motifs_set)
        return dict(motifs_counter.most_common(10))
    
    def _assemble_rna_sequences(self):
        '''Assembles the RNA sequences with RNAHeliCes/RNAmotiCes records.'''
        rna_sequences: dict[str, V] = {}
        for row in self._rna_dataframe:
            rna_sequence = rna_sequences.setdefault(row.ID, self.record_class(row.ID, row.sequence))
            rna_sequence.add_prediction(RNAmotiCesPrediction(row.mfe, row.motBracket, row.Class))
        return rna_sequences
    
    def visualize_as_histogram(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a histogram.'''
        histogram = Histogram(
            description,
            self.potential_motifs_count,
            self.sequence_number
        )
        histogram.save_plot(filepath)

class CrRNAmotiCesRecordsAssembler(RNAmotiCesRecordsAssembler[CrRNAmotiCesRecord]):
    '''Class for RNAHeliCes/RNAmotiCes assemblies from CRISPR RNA.'''
    record_class = CrRNAmotiCesRecord

    @property
    def unique_subtypes(self) -> set:
        '''Gets and returns all unique subtypes of the CRISPR RNA.'''
        subtypes = set()
        for rna_sequence in self.rna_sequences.values():
            subtypes.update(rna_sequence.subtypes)
        return subtypes