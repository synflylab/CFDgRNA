#!/usr/bin/env python3

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import math


class uDNA(Seq):
    def __init__(self, data, alphabet=IUPAC.unambiguous_dna):
        super().__init__(data, alphabet=IUPAC.unambiguous_dna)


class gRNA(uDNA):
    def __init__(self, data, alphabet=IUPAC.unambiguous_dna):
        if len(data) == 20:
            super().__init__(data)
        else:
            raise ValueError("gRNA targets must be 20 bp long")


class DistalFwdPrimer:
    def __init__(self, grna: gRNA):
        self.prefix = uDNA("GCGGCCCGGGTTCGATTCCCGGCCGATGCA")
        self.suffix = uDNA("GTTTTAGAGCTAGAAATAGCAAG")
        self.target = grna

    def seq(self):
        return str(self.prefix) + str(self.target) + str(self.suffix)


class DistalRevPrimer:
    def __init__(self, grna: gRNA):
        self.prefix = uDNA("ATTTTAACTTGCTATTTCTAGCTCTAAAAC")
        self.suffix = uDNA("TGCACCAGCCGGGAATCGAACCC")
        self.target = grna

    def seq(self):
        return str(self.prefix) + str(self.target.reverse_complement()) + str(self.suffix)


class InternalFwdPrimer:
    def __init__(self, grna: gRNA):
        self.suffix = uDNA("GTTTTAGAGCTAGAAATAGCAAG")
        self.target = grna

    def seq(self):
        return str(self.target) + str(self.suffix)


class InternalRevPrimer:
    def __init__(self, grna: gRNA):
        self.suffix = uDNA("TGCACCAGCCGGGAATCGAACCC")
        self.target = grna

    def seq(self):
        return str(self.target.reverse_complement()) + str(self.suffix)


class AssDistalFwdPrimer:
    def __init__(self, grna: gRNA):
        self.prefix = uDNA("GCGGCCCGGGTTCGATTCCCGGCCGATGCA")
        self.target = grna

    def seq(self):
        return str(self.prefix) + str(self.target)


class AssDistalRevPrimer:
    def __init__(self, grna: gRNA):
        self.prefix = uDNA("ATTTTAACTTGCTATTTCTAGCTCTAAAAC")
        self.target = grna

    def seq(self):
        return str(self.prefix) + str(self.target.reverse_complement())


class AssInternalFwdPrimer:
    def __init__(self, grna: gRNA):
        self.target = grna

    def seq(self):
        return str(self.target)


class AssInternalRevPrimer:
    def __init__(self, grna: gRNA):
        self.target = grna

    def seq(self):
        return str(self.target.reverse_complement())


class Assembler:
    max = 6

    def __init__(self, guides: list, prefix=None):
        self.prefix = prefix
        self.primers = {}
        self.guides = []
        for guide in guides:
            if isinstance(guide, gRNA):
                self.guides.append(guide)
            else:
                self.guides.append(gRNA(guide))

    def assemble(self):
        n_guides = len(self.guides)
        if n_guides <= self.max:
            for i, guide in enumerate(self.guides):
                if i == 0:
                    self.primers["PCR" + str(i + 1) + "fwd"] = DistalFwdPrimer(guide)
                elif i == len(self.guides) - 1:
                    self.primers["PCR" + str(i) + "rev"] = DistalRevPrimer(guide)
                else:
                    self.primers["PCR" + str(i) + "rev"] = InternalRevPrimer(guide)
                    self.primers["PCR" + str(i + 1) + "fwd"] = InternalFwdPrimer(guide)
        else:
            assemblies = []
            n_ass = math.ceil((n_guides - self.max) / (self.max - 1)) + 1
            n_full = n_ass - 2
            flanks = n_guides - n_full * (self.max - 1)
            par_div = math.ceil(flanks / 2)

            stop = par_div
            assemblies.append(Assembler(self.guides[0:par_div]))
            for i in range(0, n_full):
                start = par_div + i * (self.max - 1) - 1
                stop = par_div + (i + 1) * (self.max - 1)
                assemblies.append(Assembler(self.guides[start:stop]))
            assemblies.append(Assembler(self.guides[stop - 1:]))

            for i, ass in enumerate(assemblies):
                ass.assemble()
                if i == 0:
                    self.primers["ASS" + str(i + 1) + "fwd"] = AssDistalFwdPrimer(ass.guides[0]).seq()
                else:
                    self.primers["ASS" + str(i + 1) + "fwd"] = AssInternalFwdPrimer(ass.guides[0]).seq()
                for key, primer in ass.primers.items():
                    self.primers["ASS" + str(i + 1) + key] = primer.seq()
                if i == len(assemblies) - 1:
                    self.primers["ASS" + str(i + 1) + "rev"] = AssDistalRevPrimer(ass.guides[-1]).seq()
                else:
                    self.primers["ASS" + str(i + 1) + "rev"] = AssInternalRevPrimer(ass.guides[-1]).seq()


if __name__ == "__main__":
    ass = Assembler([
        "TCCTTTCGATCCGCGCCTTT",
        "GCAACATCAATTCATTAATT",
        "GTTAAGCTGAACCATTTCAT",
        "GAGTGAAAGTGACAGCTTGG",
        "AAGTTTTTAGCCTAAGATCT",
        "TCGACTGTGAATCTCTCACG",
        "AGGGGTCCCTATAATCAACA",
        "CAGCTTAGTCCGCCGCCGAG",
        "GACGAGGGTCATCATCCCGA",
        "AATCGCTGTAACAGGTGGAA"
    ], prefix="VM_")
    ass.assemble()
    for key, primer in ass.primers.items():
        print(ass.prefix + key,"\t", primer)
