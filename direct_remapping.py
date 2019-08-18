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


from datafile_content import DatafileContent


MODIFIER_TAIL = 2    # number of beads in the modifier tail
MODIFIERS_IN_SUBSYSTEM = 12    # number of modifier molecules in subsystem
POLYMERS_IN_SUBSYSTEM = 10    # number of polymer molecules in subsystem
POLYMERIZATION = 20
SUBSYSTEMS_COUNT = 9
POLYMER_TAIL_BEAD_ATOMS = 20    # number of atoms in terminating monomers
POLYMER_MIDDLE_BEAD_ATOMS = 19    # number of atoms in the middle monomers
CLAY_IN_SUBSYSTEM = 720    # MMT atoms in subsystem

MODIFIER_HEAD_SIZE = 21    # number of atoms in the modifier head bead


def get_atoms_info_from_bead_idx(bead_idx):
    '''Get the index of a soft bead and return
       indices of correspoding atoms.
    '''
    atoms_idcs = []
    phase = None

    beads_in_subsystem = (1 + MODIFIER_TAIL) * MODIFIERS_IN_SUBSYSTEM \
        + POLYMERS_IN_SUBSYSTEM * POLYMERIZATION

    subsystem_idx = bead_idx // beads_in_subsystem
    subsystem_addition = subsystem_idx * 5380    # atoms, not beads!

    in_subsystem_bead_idx = bead_idx % beads_in_subsystem

    if in_subsystem_bead_idx < MODIFIERS_IN_SUBSYSTEM * (MODIFIER_TAIL + 1):
        # this is a bead of the modifier
        phase = 'modifier'
        modifier_molecule_idx = in_subsystem_bead_idx // (MODIFIER_TAIL + 1)
        modifier_addition = modifier_molecule_idx * 70    # atoms, not beads!
        in_subsystem_bead_idx -= modifier_molecule_idx * (MODIFIER_TAIL + 1)
        in_modifier_idx = in_subsystem_bead_idx % (MODIFIER_TAIL + 1)

         # Head is always the same
        if in_modifier_idx == 0:
            atoms_idcs_tmp = list(range(MODIFIER_HEAD_SIZE))
            phase += '_head'
        else:
            phase += '_tail'
            if MODIFIER_TAIL == 2:
                if in_modifier_idx == 1:
                    atoms_idcs_tmp = list(
                        range(MODIFIER_HEAD_SIZE, MODIFIER_HEAD_SIZE + 24))
                elif in_modifier_idx == 2:
                    atoms_idcs_tmp = list(range(MODIFIER_HEAD_SIZE + 24, 70))
                else:
                    raise IndexError('Modifier tail bead index out of range')

            elif MODIFIER_TAIL == 3:
                if in_modifier_idx == 1:
                    atoms_idcs_tmp = list(
                        range(MODIFIER_HEAD_SIZE, MODIFIER_HEAD_SIZE + 18))
                elif in_modifier_idx == 2:
                    atoms_idcs_tmp = list(
                        range(MODIFIER_HEAD_SIZE + 18, MODIFIER_HEAD_SIZE + 33))
                elif in_modifier_idx == 3:
                    atoms_idcs_tmp = list(range(MODIFIER_HEAD_SIZE + 33, 70))
                else:
                    raise IndexError('Modifier tail bead index out of range')

            else:
                raise NotImplementedError

        atoms_idcs = [idx + modifier_addition + subsystem_addition
            for idx in atoms_idcs_tmp]
    else:
        # this is a bead of polymer
        phase = 'polymer'
        in_subsystem_bead_idx -= MODIFIERS_IN_SUBSYSTEM * (MODIFIER_TAIL + 1)
        modifier_addition = MODIFIERS_IN_SUBSYSTEM * 70    # atoms, not beads!
        polymer_molecule_idx = in_subsystem_bead_idx // POLYMERIZATION
        polymer_addition = polymer_molecule_idx * 382    # atoms, not beads!
        in_polymer_idx = in_subsystem_bead_idx\
                         - polymer_molecule_idx * POLYMERIZATION
        if in_polymer_idx == 0:
            phase += '_head'
            atoms_idcs_tmp = list(range(POLYMER_TAIL_BEAD_ATOMS))
        elif in_polymer_idx == POLYMERIZATION - 1:
            start = POLYMER_TAIL_BEAD_ATOMS +\
                    POLYMER_MIDDLE_BEAD_ATOMS * (POLYMERIZATION - 2)
            stop = start + POLYMER_TAIL_BEAD_ATOMS
            atoms_idcs_tmp = list(range(start, stop))
        else:
            phase += '_tail'
            polymer_addition += POLYMER_TAIL_BEAD_ATOMS    # atoms, not beads
            polymer_addition += (in_polymer_idx - 1) * POLYMER_MIDDLE_BEAD_ATOMS
            atoms_idcs_tmp = list(range(POLYMER_MIDDLE_BEAD_ATOMS))

        atoms_idcs = [
            idx + subsystem_addition + modifier_addition + polymer_addition
                for idx in atoms_idcs_tmp]

    info = {
        'atoms_idcs': atoms_idcs,
        'phase': phase
    }
    return info


def bead_idx_from_atom_idx(atom_idx):
    beads_in_subsystem = ((1 + modifier_tail_length ) * modifiers_in_subsystem
        + monomers_count*polymers_in_subsystem)
    subsystem_idx = ( atomic_idx - 1 ) // subsystem_size


def get_bead_geometry_center(dfc, atoms_idcs_list):
    '''Get list of indices of atoms, return theirs geometrical center
    '''
    # Since indices in list start from 0, while atom id-s start from 1
    atom_ids = [idx + 1 for idx in atoms_idcs_list]

    chosen_atoms = []
    for atom in dfc.atoms:
        if atom['atom_id'] in atom_ids:
            chosen_atoms.append(atom)
    chosen_atoms.sort(key=lambda atom: atom['atom_id'])

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo

    geom_center = {'x': 0, 'y': 0, 'z': 0}
    for atom in chosen_atoms:
        geom_center['x'] += atom['x'] + lx * atom['nx']
        geom_center['y'] += atom['y'] + ly * atom['ny']
        geom_center['z'] += atom['z'] + lz * atom['nz']
    geom_center['x'] /= len(chosen_atoms)
    geom_center['y'] /= len(chosen_atoms)
    geom_center['z'] /= len(chosen_atoms)

    if geom_center['x'] <= dfc.xlo:
        geom_center['x'] += lx
    if geom_center['y'] <= dfc.ylo:
        geom_center['y'] += ly
    if geom_center['z'] <= dfc.zlo:
        geom_center['z'] += lz
    return geom_center


if __name__ == '__main__':
    bead_types = {
        'clay': 1
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

    atomistic_input_fname = 'L_d.12500000.data'
    mesoscopic_outuput_fname = 'test.data'
    bead_r_c = 8.3

    '''
    chosen_atoms = set()
    for bead_idx in range(2124):
        new_atoms = get_atoms_idcs_from_bead_idx(bead_idx)
        print(bead_idx, len(new_atoms), new_atoms)
        chosen_atoms.update(new_atoms)
    print(len(chosen_atoms))
    '''

    atomistic_dfc = DatafileContent(atomistic_input_fname)
    meso_dfc = DatafileContent()

    meso_dfc.comment = '# Mesostructure made from atomistic structure'

    # cell
    meso_dfc.xlo = atomistic_dfc.xlo / bead_r_c
    meso_dfc.xhi = atomistic_dfc.xhi / bead_r_c
    meso_dfc.ylo = atomistic_dfc.ylo / bead_r_c
    meso_dfc.yhi = atomistic_dfc.yhi / bead_r_c
    meso_dfc.zlo = atomistic_dfc.zlo / bead_r_c
    meso_dfc.zhi = atomistic_dfc.zhi / bead_r_c

    '''
    meso_dfc.masses = [
        {'mass_id': 1, 'mass': 1, 'comment': 'mh'},
        {'mass_id': 2, 'mass': 1, 'comment': 'mt'},
        {'mass_id': 3, 'mass': 1, 'comment': 'p'}
    ]
    '''

    # beads
    beads_count = ((MODIFIER_TAIL + 1) * MODIFIERS_IN_SUBSYSTEM
                   + POLYMERIZATION * POLYMERS_IN_SUBSYSTEM) * SUBSYSTEMS_COUNT

    meso_dfc.atoms = []
    bead_id = 1
    molecule_tag = 0
    for bead_idx in range(beads_count):
        info = get_atoms_info_from_bead_idx(bead_idx)
        c = get_bead_geometry_center(atomistic_dfc, info['atoms_idcs'])
        if info['phase'] in ['polymer_head', 'modifier_head']:
            molecule_tag += 1
        meso_dfc.atoms.append({
            'atom_id': bead_id,
            'molecule-tag': molecule_tag,
            'atom_type_id': bead_types[info['phase']],
            'q': bead_charges[info['phase']],
            'x': c['x'] / bead_r_c,
            'y': c['y'] / bead_r_c,
            'z': c['z'] / bead_r_c,
            'nx': 0,
            'ny': 0,
            'nz': 0,
            'comment': None
        })
        bead_id += 1
    meso_dfc.atoms_count = len(meso_dfc.atoms)
    meso_dfc.atom_types = 4

    meso_dfc.write(mesoscopic_outuput_fname)
