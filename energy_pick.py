#!/usr/bin/env python3
"""energy_pick.py <gmx energy menu file> <md|inter>

Reads the term menu that `gmx energy` prints and returns the menu numbers of the
terms this pipeline needs, so terms are always chosen by name, never by a fixed number.
"""
import sys

menu = open(sys.argv[1], errors="ignore").read().splitlines()
mode = sys.argv[2]
avail = {}
for line in menu:
    t = line.split()
    if len(t) >= 2 and len(t) % 2 == 0 and all(t[i].isdigit() for i in range(0, len(t), 2)):
        for i in range(0, len(t), 2):
            avail.setdefault(t[i + 1], t[i])

if mode == "md":
    want = ["Bond", "Angle", "U-B", "Proper-Dih.", "Ryckaert-Bell.", "Fourier-Dih.", "Restricted-Dih.",
            "Improper-Dih.", "Per.-Imp.-Dih.", "CMAP-Dih.", "LJ-14", "Coulomb-14", "LJ-(SR)", "Disper.-corr.",
            "Coulomb-(SR)", "Coul.-recip.", "Potential", "Kinetic-En.", "Total-Energy", "Temperature",
            "Pressure", "Density", "Volume", "Box-X", "Box-Y", "Box-Z", "Box-XX", "Box-YY", "Box-ZZ",
            "Box-YX", "Box-ZX", "Box-ZY"]
else:
    want = ["Coul-SR:Protein-Ligand", "LJ-SR:Protein-Ligand", "Coul-14:Protein-Ligand", "LJ-14:Protein-Ligand"]

got = [w for w in want if w in avail]
sys.stderr.write("  energy terms selected: " + ", ".join(got) + "\n")
print(" ".join(avail[w] for w in got))
