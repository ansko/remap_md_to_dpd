from .get_atoms_info_from_bead_idx import get_atoms_info_from_bead_idx
from .get_bead_geometry_center import get_bead_geometry_center


def convert_soft(atomistic_dfc, meso_dfc, **kwargs):
    BEAD_R_C = kwargs['BEAD_R_C']
    MODIFIER_TAIL = kwargs['MODIFIER_TAIL']
    MODIFIERS_IN_SUBSYSTEM = kwargs['MODIFIERS_IN_SUBSYSTEM']
    POLYMERIZATION = kwargs['POLYMERIZATION']
    POLYMERS_IN_SUBSYSTEM = kwargs['POLYMERS_IN_SUBSYSTEM']
    SUBSYSTEMS_COUNT = kwargs['SUBSYSTEMS_COUNT']
    MODIFIER_HEAD_SIZE = kwargs['MODIFIER_HEAD_SIZE']
    POLYMER_TAIL_BEAD_ATOMS = kwargs['POLYMER_TAIL_BEAD_ATOMS']
    POLYMER_MIDDLE_BEAD_ATOMS = kwargs['POLYMER_MIDDLE_BEAD_ATOMS']
    bead_types = kwargs['bead_types']
    bead_charges = kwargs['bead_charges']
    bond_types = kwargs['bond_types']

    bead_id = len(meso_dfc.atoms) + 1
    molecule_tag = 1

    soft_beads_count = ((MODIFIER_TAIL + 1) * MODIFIERS_IN_SUBSYSTEM
                   + POLYMERIZATION * POLYMERS_IN_SUBSYSTEM) * SUBSYSTEMS_COUNT

    for bead_idx in range(soft_beads_count):
        info = get_atoms_info_from_bead_idx(bead_idx,
            MODIFIER_TAIL=MODIFIER_TAIL,
            MODIFIERS_IN_SUBSYSTEM=MODIFIERS_IN_SUBSYSTEM,
            POLYMERS_IN_SUBSYSTEM=POLYMERS_IN_SUBSYSTEM,
            POLYMERIZATION=POLYMERIZATION,
            MODIFIER_HEAD_SIZE=MODIFIER_HEAD_SIZE,
            POLYMER_TAIL_BEAD_ATOMS=POLYMER_TAIL_BEAD_ATOMS,
            POLYMER_MIDDLE_BEAD_ATOMS=POLYMER_MIDDLE_BEAD_ATOMS)
        #if bead_idx == 28:
        #    for idx in info['atoms_idcs']:
        #        print(idx, atomistic_dfc.atoms[idx]['x'],
        #              atomistic_dfc.atoms[idx]['y'],
        #              atomistic_dfc.atoms[idx]['z'])
        #    1/0
        c = get_bead_geometry_center(atomistic_dfc, info['atoms_idcs'])
        print(info['phase'])
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
        bead_id += 1

    # analyze bond lengths
    bond_lengths = {
        'mod_head_tail': {'count': 0, 'value': 0},
        'mod_tail_tail': {'count': 0, 'value': 0},
        'polymer': {'count': 0, 'value': 0},
    }
    meso_lx = meso_dfc.xhi - meso_dfc.xlo
    meso_ly = meso_dfc.yhi - meso_dfc.ylo
    meso_lz = meso_dfc.zhi - meso_dfc.zlo

    bond_id = meso_dfc.bonds[-1]['bond_id'] + 1
    for idx, atom in enumerate(meso_dfc.atoms):
        if atom['comment'] is None:
            continue
        if atom['molecule-tag'] != meso_dfc.atoms[idx - 1]['molecule-tag']:
           continue
        atom2 = meso_dfc.atoms[idx - 1]
        dx = abs(atom2['x'] - atom['x'])
        dx = min(dx, meso_lx - dx)
        dy = abs(atom2['y'] - atom['y'])
        dy = min(dy, meso_ly - dy)
        dz = abs(atom2['z'] - atom['z'])
        dz = min(dz, meso_lz - dz)
        dr = (dx**2 + dy**2 + dz**2)**0.5
        k = None
        if atom['comment'] == atom2['comment'] == 'modifier_tail':
            k = 'mod_tail_tail'
            meso_dfc.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': bond_types['modifier_tail_tail'],
                'atom_one_id': atom['atom_id'],
                'atom_two_id': atom2['atom_id'],
            })
            bond_id += 1
        elif atom['comment'] == 'modifier_tail':
            k = 'mod_head_tail'
            meso_dfc.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': bond_types['modifier_head_tail'],
                'atom_one_id': atom['atom_id'],
                'atom_two_id': atom2['atom_id'],
            })
            bond_id += 1
        elif atom['comment'].startswith('polymer') and atom2['comment'].startswith('polymer'):
            k = 'polymer'
            meso_dfc.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': bond_types['polymer'],
                'atom_one_id': atom['atom_id'],
                'atom_two_id': atom2['atom_id'],
            })
            bond_id += 1

        if k is not None and dr > 5:
            print('---', atom['comment'], atom2['comment'])
            print(idx, k, dr, dx, dy, dz)
        bond_lengths[k]['value'] += dr
        bond_lengths[k]['count'] += 1

    for k in ['mod_head_tail', 'mod_tail_tail', 'polymer']:
        try:
            print(k, bond_lengths[k]['value'] / bond_lengths[k]['count'])
        except ZeroDivisionError:
            print(k)
