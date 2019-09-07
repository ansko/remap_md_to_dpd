MODIFIER_TAIL = 2    # number of beads in the modifier tail
MODIFIERS_IN_SUBSYSTEM = 12    # number of modifier molecules in subsystem
POLYMERS_IN_SUBSYSTEM = 10    # number of polymer molecules in subsystem
POLYMERIZATION = 20
SUBSYSTEMS_COUNT = 9
POLYMER_TAIL_BEAD_ATOMS = 20    # number of atoms in terminating monomers
POLYMER_MIDDLE_BEAD_ATOMS = 19    # number of atoms in the middle monomers
CLAY_IN_SUBSYSTEM = 720    # MMT atoms in subsystem

MODIFIER_HEAD_SIZE = 21    # number of atoms in the modifier head bead

BEAD_R_C = 8.3

bead_types = {
        'clay': 1,
        'modifier_head': 2,
        'modifier_tail': 3,
        'polymer': 4,
        'polymer_head': 4,
        'polymer_tail': 4
}
bead_charges = {
        'clay': -1,    # sometimes
        'modifier_head': 1,
        'modifier_tail': 0,
        'polymer': 0,
        'polymer_head': 0,
        'polymer_tail': 0
}
bond_types = {
        'mmt_edge': 1,
        'mmt_diagonal': 2,
        'modifier_head_tail': 3,
        'modifier_tail_tail': 4,
        'polymer': 5
}
