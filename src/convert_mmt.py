def convert_mmt(atomistic_dfc, **kwargs):
    '''Make mmt and return its atoms and bonds
       (without appending to structure etc.)
    '''
    BEAD_R_C = kwargs['BEAD_R_C']

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
