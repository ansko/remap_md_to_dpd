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
    if geom_center['x'] >= dfc.xhi:
        geom_center['x'] -= lx
    if geom_center['y'] <= dfc.ylo:
        geom_center['y'] += ly
    if geom_center['y'] >= dfc.yhi:
        geom_center['y'] -= ly
    if geom_center['z'] <= dfc.zlo:
        geom_center['z'] += lz
    if geom_center['z'] >= dfc.zhi:
        geom_center['z'] -= lz

    return geom_center
