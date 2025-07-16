'''
Contains all necessary dataclasses and classes for the processing of the final mature CRISPR RNAs.

author: U.B.
'''

import pandas as pd
from dataclasses import dataclass
from crispr_cas_db.processing.Subtype import *

SUBTYPE_CLASS_DICT = {"CAS-TypeI-A": SubtypeI_A, "CAS-TypeI-B": SubtypeI_B, "CAS-TypeI-C": SubtypeI_C,
    "CAS-TypeI-D": SubtypeI_D, "CAS-TypeI-E": SubtypeI_E, "CAS-TypeI-F": SubtypeI_F, "CAS-TypeI-G": SubtypeI_G,
    "CAS-TypeII-A": SubtypeII_A, "CAS-TypeII-B": SubtypeII_B, "CAS-TypeII-C": SubtypeII_C,
    "CAS-TypeIII-A": SubtypeIII_A, "CAS-TypeIII-B": SubtypeIII_B, "CAS-TypeIII-D": SubtypeIII_D, 
    "CAS-TypeV-A": SubtypeV_A, "CAS-TypeV-F4": SubtypeV_F4}

@dataclass(frozen=True)
class ArrayElement:
    '''Represent an element in a CRISPR array.'''
    sequence: str

@dataclass(frozen=True)
class Repeat(ArrayElement):
    '''Represent a Repeat in a CRISPR array'''
    pass

@dataclass(frozen=True)
class Spacer(ArrayElement):
    '''Represent a Spacer in a CRISPR array'''
    pass

@dataclass(frozen=True)
class CrRNA:
    '''Represents a mature CRISPR RNA with its subtype'''
    sequence: str
    subtype: str


class CrisprArray:
    '''Class representing a single array'''
    def __init__(self, locus: pd.DataFrame) -> None:
        '''Initializes a CrisprArray object'''
        self._array = locus.sort_values("crisprlocus_region_start").reset_index(drop=True)
        self._subtype = None
        self._orientation = None
        self._elements = None

    @property
    def subtype(self) -> str:
        '''Gets and returns the subtype of the array, if a subtype can be clearly assigned'''
        if self._subtype is None:
            if self._array["clustercas_class"].nunique() == 1 and self._array["clustercas_class"].iloc[0] in SUBTYPE_CLASS_DICT:
                self._subtype = self._array["clustercas_class"].iloc[0]
        return self._subtype
    
    @property
    def orientation(self) -> str:
        '''Gets and returns the orientation of the array, either in forward or reverse direction'''
        if self._orientation is None:
            if self._array["crisprlocus_orientation"].iloc[0] == 1:
                self._orientation = "Forward"
            elif self._array["crisprlocus_orientation"].iloc[0] == 2:
                self._orientation = "Reverse"
        return self._orientation
    
    @property
    def elements(self) -> dict[int, Spacer | Repeat]:
        '''Represent all elements in an array (SPacer and Repeats)'''
        if self._elements is None:
            self._elements = self._construct_array_elements()
        return self._elements

    def _construct_array_elements(self) -> dict[int, Spacer | Repeat]:
        '''Constructs an indexed dictionary of all elements in the array with their corresponding type (Spacer/Repeat)'''
        elements = {}
        for index, row in self._array.iterrows():
            if row["region_category"] == "Spacer":
                elements[index] = Spacer(row["region_sequence"])
            elif row["region_category"] == "Repeat":
                elements[index] = Repeat(row["region_sequence"])
        return elements

class MatureCrRNAs:
    '''Class representing a the mature CRISPR RNAs of a single array after processing them'''
    def __init__(self, array: CrisprArray) -> None:
        '''Initializes a MatureCrRNAs object'''
        self._array_elements = array.elements
        self._subtype = array.subtype
        self._orientation = array.orientation
        self._crRNAs = None

    @property
    def crRNAs(self) -> list[CrRNA]:
        '''Gets and returns a list of all CRISPR RNAs in the array'''
        if self._crRNAs is None:
            self._crRNAs = self._assemble_crRNAs()
        return self._crRNAs
    
    def _assemble_crRNAs(self) -> list[CrRNA]:
        '''Assembles the final CRISPR RNAs by building the pre-crRNA and processing it depending on the subtype'''
        crRNAs = []
        for index, element in self._array_elements.items():
            if isinstance(element, Spacer):
                prev_element, next_element = self._extract_flanking_repeats(index)
                if prev_element is not None and next_element is not None:
                    CRISPRSubtype = SUBTYPE_CLASS_DICT.get(self._subtype)
                    crRNA: Subtype = CRISPRSubtype(prev_element.sequence, element.sequence, next_element.sequence)
                    crRNAs.append(CrRNA(crRNA.sequence, self._subtype))
        return crRNAs
    
    def _extract_flanking_repeats(self, index: int) -> tuple[Repeat, Repeat]:
        '''Extracts the flanking repeats of a Spacer in the array'''
        if self._orientation == "Forward":
            prev_element = self._array_elements.get(index - 1)
            next_element = self._array_elements.get(index + 1)
        elif self._orientation == "Reverse":
            prev_element = self._array_elements.get(index + 1)
            next_element = self._array_elements.get(index - 1)
        return prev_element, next_element