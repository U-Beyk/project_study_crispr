'''
Contains all Prediction classes corresponding to either an RNAmotiFold prediction or an RNAmotiCes prediction.

author: U.B.
'''

import re
import pandas as pd
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

#TO-DO: Implement recognition of ambiguous motifs instead of hardcoding
AMBIGUOUS_MOTIFS = {"u": "GU", "g": "GT", "t": "GT"}

@dataclass
class Prediction(ABC):
    '''Abstract class of an RNA prediction.'''
    free_energy: float
    mot_bracket: str
    _distance_to_mfe: float = field(init=False, default=None)

    def __post_init__(self) -> None:
        '''Converts the free energy value into kcal/mol by dividing with 100.'''
        self.free_energy /= 100

    @property
    def distance_to_mfe(self) -> float | None:
        '''Returns the distance to the mfe value if present.'''
        return self._distance_to_mfe

    @distance_to_mfe.setter
    def distance_to_mfe(self, min_mfe: float) -> None:
        '''Sets the distance to the lowest mfe.'''
        self._distance_to_mfe = self.free_energy - min_mfe

    @property
    def distance_to_mfe_by_motif(self) -> dict[str, float]:
        '''Assigns to each motif the value for the distance to the mfe.
        If there are no motifs, it returns a corresponding dictionary.'''
        if not self.motifs_set:
            return {"No motif": self.distance_to_mfe}
        return {motif: self.distance_to_mfe for motif in self.motifs_set}

    @property
    @abstractmethod
    def motifs_set(self) -> set[str]:
        '''Returns a set of unique motifs in the prediction.'''
        pass
    
    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (f"{self.__class__.__name__}(mfe={self.free_energy:.2f}, "
                f"mot_bracket='{self.mot_bracket}', "
                f"distance_to_lowest_mfe={self.distance_to_mfe:.2f})")

@dataclass
class RNAmotiFoldPrediction(Prediction):
    '''Represents an RNAmotiFold prediction.'''
    motifs: str

    def __post_init__(self) -> None:
        '''Additionally to the parent method, turns the motifs attribute into an empty string if it has a "nan" value.'''
        super().__post_init__()
        self.motifs = "" if pd.isna(self.motifs) else self.motifs

    @property
    def motifs_set(self) -> set[str]:
        '''Gets and returns a set of unique motifs in the prediction, replacing ambiguous lowercase motifs.'''
        motifs = set()
        for motif in self.motifs:
            if motif.islower() and motif in AMBIGUOUS_MOTIFS:
                motifs.update(AMBIGUOUS_MOTIFS[motif])
            else:
                motifs.update(motif)
        return motifs

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (f"{super().__repr__()}, motifs='{self.motifs}')")

@dataclass
class RNAmotiCesPrediction(Prediction):
    '''Represents an RNAmotiCes prediction'''
    motices: str | float
    _potential_motif_sequences: list[str] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        '''Additionally to the parent method, turns the motices attribute into an empty string if it has a "nan" value.'''
        super().__post_init__()
        self.motices = "" if pd.isna(self.motices) else self.motices

    @property
    def potential_motif_sequences(self) -> list[str]:
        '''Gets and returns a list of potential motifs not corresponding to any known motif,
        computed via positional abstraction.'''
        return self._potential_motif_sequences

    def compute_potential_motifs(self, sequence) -> None:
        '''Computes the potential motif for the corresponding sequence of that prediction.'''
        for position in self.positions_without_motif:
            left_nt_pos = self._find_left_bracket_position(round(float(position)))
            right_nt_pos = self._find_right_bracket_position(round(float(position)))
            self._potential_motif_sequences.append(sequence[left_nt_pos: right_nt_pos])

    def _find_left_bracket_position(self, start_pos: int) -> int | None:
        '''Finds the nearest '(' bracket position to the left.'''
        for pos in range(start_pos, 0, -1):
            if self.mot_bracket[pos] == "(":
                return pos + 1
        return None

    def _find_right_bracket_position(self, start_pos: int) -> int | None:
        '''Finds the nearest ')' bracket position to the right.'''
        for pos in range(start_pos, len(self.mot_bracket)):
            if self.mot_bracket[pos] == ")":
                return pos
        return None

    @property
    def motifs_set(self) -> set[str]:
        '''Gets and returns a set of unique motifs, replacing ambiguous ones.'''
        motifs = set()
        for motif in str(self.motices):
            if motif.isalpha():
                if motif.islower() and motif in AMBIGUOUS_MOTIFS:
                    motifs.update(AMBIGUOUS_MOTIFS[motif])
                else:
                    motifs.update(motif)
        return motifs
    
    @property
    def positions_without_motif(self) -> list[float]:
        '''finds and returns all the secondary structures without a motif.'''
        return re.findall(r'\d+\.?\d+(?=_)', self.motices)

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (f"{super().__repr__()}, motices='{self.motices}')")