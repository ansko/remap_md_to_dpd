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

BEAD_R_C = 8.3


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
            phase += '_tail'
            start = POLYMER_TAIL_BEAD_ATOMS +\
                    POLYMER_MIDDLE_BEAD_ATOMS * (POLYMERIZATION - 2)
            stop = start + POLYMER_TAIL_BEAD_ATOMS
            atoms_idcs_tmp = list(range(start, stop))
        else:
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


def get_mmt(atomistic_dfc):
    '''Make mmt and return its atoms and bonds
       (without appending to structure etc.)
    '''
    density_coeff = 2  # distance between beads = 2*bead_r_c / density_coeff

    lx = atomistic_dfc.xhi - atomistic_dfc.xlo
    ly = atomistic_dfc.yhi - atomistic_dfc.ylo

    nx = int(lx / (2 * BEAD_R_C / density_coeff))
    ny = int(ly / (2 * BEAD_R_C / density_coeff))

    count = 0
    mmt_zlo = atomistic_dfc.zhi
    mmt_zhi = atomistic_dfc.zlo
    for atom in atomistic_dfc.atoms:
        if 0 < atom['atom_id'] % 5380 < 721:    # 0 for non 0th subsystem
            count += 1
            mmt_zlo = min(mmt_zlo, atom['z'])
            mmt_zhi = max(mmt_zhi, atom['z'])

    bead_id = 1
    bond_id = 1
    beads = []
    bonds = []

    z_center = (mmt_zlo + mmt_zhi) / 2
    mmt_zlo = z_center - BEAD_R_C / density_coeff
    mmt_zhi = z_center + BEAD_R_C / density_coeff

    beads_idcs_map = {
        idx_x: {idx_y: {'bottom': None, 'top': None} for idx_y in range(ny)}
            for idx_x in range(nx)
    }

    nz = 0
    for idx_x in range(nx):
        x = idx_x * 2 * BEAD_R_C / density_coeff
        nx = 0
        for idx_y in range(ny):
            y = idx_y * 2 * BEAD_R_C / density_coeff
            ny = 0
            beads.append({
                'atom_id': bead_id,
                'molecule-tag': 1,
                'atom_type_id': 1,
                'q': 0,
                'x': x,
                'y': y,
                'z': mmt_zlo,
                'nx': nx, 'ny': ny, 'nz': nz,
                'comment': None
            })
            beads_idcs_map[idx_x][idx_y]['bottom'] = bead_id
            bead_id += 1
            beads.append({
                'atom_id': bead_id,
                'molecule-tag': 1,
                'atom_type_id': 1,
                'q': 0,
                'x': x,
                'y': y,
                'z': mmt_zhi,
                'nx': nx, 'ny': ny, 'nz': nz,
                'comment': None
            })
            beads_idcs_map[idx_x][idx_y]['top'] = bead_id
            bead_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads[-1]['atom_id'],
                'atom_two_id': beads[-2]['atom_id'],
            })
            bond_id += 1

    for idx_x in range(nx):
        bot_idx_x = idx_x - 1 if idx_x > 0 else nx - 1
        for idx_y in range(ny):
            bot_idx_y = idx_y - 1 if idx_y > 0 else ny - 1
            # edge bonds
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[bot_idx_x][idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[bot_idx_x][idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1

            # diagonal bonds
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[bot_idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[bot_idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[bot_idx_x][idx_y]['top']
            })
            bond_id += 1
            bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[bot_idx_x][idx_y]['bottom']
            })
            bond_id += 1

    return {'beads': beads, 'bonds': bonds}


if __name__ == '__main__':
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

    mmt_data = get_mmt(atomistic_dfc)

    meso_dfc.atoms = mmt_data['beads']
    meso_dfc.bonds = mmt_data['bonds']

    for atom in meso_dfc.atoms:
        atom['x'] /= BEAD_R_C
        atom['y'] /= BEAD_R_C
        atom['z'] /= BEAD_R_C

    bead_id = len(meso_dfc.atoms) + 1
    molecule_tag = 1

    soft_beads_count = ((MODIFIER_TAIL + 1) * MODIFIERS_IN_SUBSYSTEM
                   + POLYMERIZATION * POLYMERS_IN_SUBSYSTEM) * SUBSYSTEMS_COUNT

    bond_lengths = {
        'mod_head_tail': {'count': 0, 'value': 0},
        'mod_tail_tail': {'count': 0, 'value': 0},
        'polymer': {'count': 0, 'value': 0},
    }

    for bead_idx in range(soft_beads_count):
        info = get_atoms_info_from_bead_idx(bead_idx)
        c = get_bead_geometry_center(atomistic_dfc, info['atoms_idcs'])
        if info['phase'] in ['polymer_head', 'modifier_head']:
            molecule_tag += 1
        meso_dfc.atoms.append({
            'atom_id': bead_id,
            'molecule-tag': molecule_tag,
            'atom_type_id': bead_types[info['phase']],
            'q': bead_charges[info['phase']],
            'x': c['x'] / BEAD_R_C,
            'y': c['y'] / BEAD_R_C,
            'z': c['z'] / BEAD_R_C,
            'nx': 0,
            'ny': 0,
            'nz': 0,
            'comment': info['phase']
        })
        try:
            dx = abs(meso_dfc.atoms[bead_id - 1]['x']
                     -meso_dfc.atoms[bead_id]['x'])
            dy = abs(meso_dfc.atoms[bead_id - 1]['y']
                     -meso_dfc.atoms[bead_id]['y'])
            dz = abs(meso_dfc.atoms[bead_id - 1]['z']
                     -meso_dfc.atoms[bead_id]['z'])
            dr = (dx**2 + dy**2 + dz**2)**0.5
            print(info['phase'])
            if (set((info['phase'], meso_dfc.atoms[bead_id - 1]['comment']))
                == set(('modifier_tail', 'modifier_head'))):
                    bond_lengths['mod_head_tail']['value'] += dr
                    bond_lengths['mod_head_tail']['count'] += 1
            elif (set(info['phase'], meso_dfc.atoms[bead_id - 1]['comment'])
                  == set(('modifier_tail'))):
                    bond_lengths['mod_tail_tail']['value'] += dr
                    bond_lengths['mod_tail_tail']['count'] += 1
            elif (set(info['phase'], meso_dfc.atoms[bead_id - 1]['comment'])
                  in [set(('polymer', 'polymer')),
                      set(('polymer_head', 'polymer')),
                      set(('polymer_tail', 'polymer'))]):
                    bond_lengths['polymer']['value'] += dr
                    bond_lengths['polymer']['count'] += 1
        except Exception:
            pass

    bead_id += 1

    for k in ['mod_head_tail', 'mod_tail_tail', 'polymer']:
        try:
            print(k, bond_lengths[k]['value'] / bond_lengths[k]['count'])
        except ZeroDivisionError:
            print(k)

    meso_dfc.atoms_count = len(meso_dfc.atoms)
    meso_dfc.atom_types = 4
    meso_dfc.bonds_count = len(meso_dfc.bonds)
    meso_dfc.bond_types = 5
    meso_dfc.write(mesoscopic_outuput_fname)
