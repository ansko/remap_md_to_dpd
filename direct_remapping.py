'''
***********************************************
*** MD structure -> DPD structure converter ***
***********************************************
---
Perform a direct fair remapping of the atomistic structure onto the
mesosocpic structure based on the beads geometrical centers positions.
---
MD structure desctiption is made of 48420 atoms:
    6480 atoms of mmt
    7560 atoms of modifier (70 atoms x 108 molecules)
    34380 atoms of polymer (20 monomers x 90 chains)
A big MD system is made by replication (3x3x1 times along x, y, and z,
correspondingly) of the small system (subsystem), which indices are:
    1-720 mmt
    721-1560 modifier
    1561-5380 polymer
---
DPD structure:
    ??? MMT beads
    3 beads x 108 molecules -> 324 beads
    20 beads x 90 molecules -> 1800 beads
    2124 soft beads

    3*12 + 20*10 = 36 + 200 = 236 beads in subsystem
'''


from src.datafile_content import DatafileContent
from src.get_atoms_info_from_bead_idx import get_atoms_info_from_bead_idx
from src.get_bead_geometry_center import get_bead_geometry_center
from src.convert_mmt import convert_mmt
from src.convert_soft import convert_soft

from options import *


if __name__ == '__main__':
    atomistic_input_fname = 'L_d.12500000.data'
    mesoscopic_outuput_fname = 'test.data'

    atomistic_dfc = DatafileContent(atomistic_input_fname)
    meso_dfc = DatafileContent()

    meso_dfc.comment = '# Mesostructure made from atomistic structure'

    # cell
    meso_dfc.xlo = atomistic_dfc.xlo / BEAD_R_C
    meso_dfc.xhi = atomistic_dfc.xhi / BEAD_R_C
    meso_dfc.ylo = atomistic_dfc.ylo / BEAD_R_C
    meso_dfc.yhi = atomistic_dfc.yhi / BEAD_R_C
    meso_dfc.zlo = atomistic_dfc.zlo / BEAD_R_C
    meso_dfc.zhi = atomistic_dfc.zhi / BEAD_R_C

    # beads

    mmt_data = convert_mmt(atomistic_dfc,
        BEAD_R_C=BEAD_R_C)

    meso_dfc.atoms = mmt_data['beads']
    meso_dfc.bonds = mmt_data['bonds']

    for atom in meso_dfc.atoms:
        atom['x'] /= BEAD_R_C
        atom['y'] /= BEAD_R_C
        atom['z'] /= BEAD_R_C

    convert_soft(atomistic_dfc, meso_dfc,
        BEAD_R_C=BEAD_R_C,
        MODIFIER_TAIL=MODIFIER_TAIL,
        MODIFIERS_IN_SUBSYSTEM=MODIFIERS_IN_SUBSYSTEM,
        POLYMERIZATION=POLYMERIZATION,
        POLYMERS_IN_SUBSYSTEM=POLYMERS_IN_SUBSYSTEM,
        SUBSYSTEMS_COUNT=SUBSYSTEMS_COUNT,
        MODIFIER_HEAD_SIZE=MODIFIER_HEAD_SIZE,
        POLYMER_TAIL_BEAD_ATOMS=POLYMER_TAIL_BEAD_ATOMS,
        POLYMER_MIDDLE_BEAD_ATOMS=POLYMER_MIDDLE_BEAD_ATOMS,
        bead_types=bead_types,
        bead_charges=bead_charges,
        bond_types=bond_types
    )

    meso_dfc.atoms_count = len(meso_dfc.atoms)
    meso_dfc.atom_types = 4
    meso_dfc.bonds_count = len(meso_dfc.bonds)
    meso_dfc.bond_types = 5
    meso_dfc.write(mesoscopic_outuput_fname)
