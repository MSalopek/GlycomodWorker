from typing import NamedTuple, List


class GlycomodComposition(NamedTuple):
    theoretical_MH: float
    delta: float  
    long_notation: str
    short_notation: str
    theoretical_MTagH: float
    theoretical_MTagNa: float
    theoretical_MTagK: float
    theoretical_MTagH2: float
    theoretical_MTagHNa: float
    theoretical_MTagHK: float
    theoretical_MTagNa2: float
    theoretical_MTagNH4: float
    def __repr__(self):
        return f"Theoretical [MH]+: {self.theoretical_MH}, Delta: {self.delta}, Comp: {self.short_notation}"

 
class SubmittedMass(NamedTuple):
    # peak_number: int
    experimental_mass: float
    adduct: str  
    adduct_mass: float
    red_end_tag: str 
    red_end_tag_mass: float
    glycomod_structures: List[GlycomodComposition]

    def __repr__(self):
        return f"Experimental mass: {self.experimental_mass}\n" + \
            f"Adduct: {self.adduct} ({self.adduct_mass})\n" + \
            f"Reducing end: {self.red_end_tag} ({self.red_end_tag_mass})\n" + \
            f"Compositions: {len(self.glycomod_structures)}\n"

    def prep_csv_out(self):
        """Prepare tuples used for outputing as csv"""
        structure_list = []
        if len(self.glycomod_structures) > 0:
            for i in self.glycomod_structures:
                structure_list.append(
                    (
                        # self.peak_number,
                        self.experimental_mass,
                        self.red_end_tag,
                        self.red_end_tag_mass,
                        # self.adduct,
                        # self.adduct_mass,
                        i.short_notation,
                        i.long_notation,
                        i.theoretical_MH,
                        i.delta,
                        i.theoretical_MTagH,
                        i.theoretical_MTagNa,
                        i.theoretical_MTagK,
                        i.theoretical_MTagNH4,
                        i.theoretical_MTagH2,
                        i.theoretical_MTagHNa,
                        i.theoretical_MTagHK,
                        i.theoretical_MTagNa2,
                    )
                )
        return structure_list
