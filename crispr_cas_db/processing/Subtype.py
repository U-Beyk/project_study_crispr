'''
Processes the mature CRISPR RNA by its subtype. CRISPR CAS classes and subtypes are used synonymous.

author: U.B.
'''

from abc import ABC, abstractmethod
from typing import override

class Subtype(ABC):
    '''General abstract class processing mature CRISPR RNAs'''
    def __init__(self, repeat5: str, spacer: str, repeat3: str) -> None:
        '''Initializes CrisprSubtype object'''
        self._repeat5 = repeat5
        self._spacer = spacer
        self._repeat3 = repeat3

    @property
    @abstractmethod
    def sequence(self) -> str:
        '''Gets and returns the mature crRNA sequence'''
        pass

'''All subtypes process the mature crRNA depending on its type, detials why the subtype
processes it the way it does can be found in the text file "subtype_processing_sources.txt"'''
class SubtypeI_A(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"

class SubtypeI_B(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeI_C(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 11:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 11]}"
    
class SubtypeI_D(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeI_E(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeI_F(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeI_G(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeII_A(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._spacer[-20:] if len(self._spacer) > 20 else self._spacer}{self._repeat3[:19] if len(self._repeat3) >= 19 else self._repeat3}"
    
class SubtypeII_B(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._spacer[-20:] if len(self._spacer) > 20 else self._spacer}{self._repeat3[:19] if len(self._repeat3) >= 19 else self._repeat3}"
    
class SubtypeII_C(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._spacer[-20:] if len(self._spacer) > 20 else self._spacer}{self._repeat3[:19] if len(self._repeat3) >= 19 else self._repeat3}"
    
class SubtypeIII_A(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeIII_B(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 8:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 8]}"
    
class SubtypeIII_D(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat5[len(self._repeat5) - 11:]}{self._spacer}{self._repeat3[:len(self._repeat3) - 11]}"

class SubtypeV_A(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat3[-19:] if len(self._repeat3) >= 19 else self._repeat3}{self._spacer[:20] if len(self._spacer) > 20 else self._spacer}"

class SubtypeV_F4(Subtype):
    @property
    @override
    def sequence(self) -> str:
        return f"{self._repeat3[-17:] if len(self._repeat3) >= 17 else self._repeat3}{self._spacer[:20] if len(self._spacer) > 20 else self._spacer}"