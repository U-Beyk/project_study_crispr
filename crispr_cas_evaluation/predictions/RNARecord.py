'''
Contains all classes with different RNA record types to separate RNAmotiFold from RNAHeliCes/RNAmotiCes prediction.

author: U.B.
'''

from dataclasses import dataclass, field
from typing import Generic, TypeVar, Self
from crispr_cas_evaluation.predictions.Prediction import Prediction, RNAmotiCesPrediction, RNAmotiFoldPrediction

T = TypeVar("T", bound=Prediction)

@dataclass
class RNARecord(Generic[T]):
    '''Generic class for an RNA record, containing the correpsonding RNA sequence and all its predictions.'''
    sequence_header: str
    sequence: str
    _predictions: list[T] = field(default_factory=list)

    def __post_init__(self) -> None:
        '''After initializing an object the sequence ID is parsed from the header.'''
        self._sequence_id = self.sequence_header.split("|")[0]

    @property
    def sequence_id(self) -> str:
        '''Gets and returns the sequence ID.'''
        return self._sequence_id

    @property
    def predictions(self) -> list[T]:
        '''Gets and returns a list of the prediction.'''
        return self._predictions

    @property
    def mfe_prediction(self) -> T:
        '''Calculates and returns the prediction with the mfe value.'''
        return min(self.predictions, key=lambda pred: pred.free_energy, default=None)

    @property
    def mfe_value(self) -> float:
        '''Returns the mfe value as a float.'''
        return self.mfe_prediction.free_energy if self.mfe_prediction else float("inf")

    @property
    def motifs_set(self) -> set[str]:
        '''Returns a set of unique motifs of all predictions.'''
        return {motif for pred in self.predictions for motif in pred.motifs_set} or {"No motif"}

    @property
    def prediction_number(self) -> int:
        '''Returns the total number of all predictions.'''
        return len(self.predictions)
    
    @property
    def distance_to_lowest_all_motifs(self) -> list[dict[str, float]]:
        '''Returns a list of all motifs '''
        return [prediction.distance_to_mfe_by_motif for prediction in self.predictions]

    def filter_predictions(self, mfe_range: float) -> Self:
        '''Filters predictions within a specified mfe range from the mfe and returns a new RNARecord instance.'''
        mfe = self.mfe_value
        new_record = self.__class__(self.sequence_header, self.sequence)
        new_record._predictions = list(
            pred for pred in self._predictions if pred.free_energy <= mfe + abs(mfe_range)
        )
        new_record.update_mfe_distance()
        return new_record

    def add_prediction(self, prediction: T) -> None:
        '''Adds a prediction and updates its distances to the lowest mfe value.'''
        self._predictions.append(prediction)
        self.update_mfe_distance()

    def update_mfe_distance(self) -> None:
        '''Updates and calculates the distance to the mfe value for all predictions.'''
        if not self.predictions:
            return
        min_mfe = self.mfe_value
        for pred in self.predictions:
            pred.distance_to_mfe = min_mfe

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (
            f"{self.__class__.__name__}("
            f"sequence_id={self.sequence_id!r}, "
            f"sequence={self.sequence!r}, "
            f"predictions={self._predictions!r})"
        )

@dataclass
class RNAmotiFoldRecord(RNARecord[RNAmotiFoldPrediction]):
    '''RNA record for RNAmotiFold predictions.'''
    pass

@dataclass
class RNAmotiCesRecord(RNARecord[RNAmotiCesPrediction]):
    '''RNA record for RNAHeliCes/RNAmotiCes predictions.'''
    @property
    def potential_motifs_set(self) -> set[str]:
        '''Returns all the potential motifs that were computed via positional abstraction.'''
        potential_motifs = set()
        for prediction in self.predictions:
            prediction.compute_potential_motifs(self.sequence)
            potential_motifs.update(prediction.potential_motif_sequences)
        return potential_motifs

@dataclass
class CrRNAmotiFoldRecord(RNAmotiFoldRecord):
    '''RNA record of RNAmotiFold predictions for CRISPR RNA.'''
    subtypes: list[str] = field(init=False)

    def __post_init__(self):
        '''Additionally to the parent method, determines the subtype of the RNA record.'''
        super().__post_init__()
        subtype_section = self.sequence_header.split("|")[1]
        self.subtypes = subtype_section.replace("subtype:", "").split(",")

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return super().__repr__() + f", subtypes={self.subtypes!r}"

@dataclass
class CrRNAmotiCesRecord(RNAmotiCesRecord):
    '''RNA record of RNAHeliCes/RNAmotiCes predictions for CRISPR RNA.'''
    subtypes: list[str] = field(init=False)

    def __post_init__(self):
        '''Additionally to the parent method, determines the subtype of the RNA record.'''
        super().__post_init__()
        subtype_section = self.sequence_header.split("|")[1]
        self.subtypes = subtype_section.replace("subtype:", "").split(",")

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return super().__repr__() + f", subtypes={self.subtypes!r}"