def get_atoms_info_from_bead_idx(bead_idx, **kwargs):
    '''Get the index of a soft bead and return
       indices of correspoding atoms.
    '''
    MODIFIER_TAIL = kwargs['MODIFIER_TAIL']
    MODIFIERS_IN_SUBSYSTEM = kwargs['MODIFIERS_IN_SUBSYSTEM']
    POLYMERS_IN_SUBSYSTEM = kwargs['POLYMERS_IN_SUBSYSTEM']
    POLYMERIZATION = kwargs['POLYMERIZATION']
    MODIFIER_HEAD_SIZE = kwargs['MODIFIER_HEAD_SIZE']
    POLYMER_TAIL_BEAD_ATOMS = kwargs['POLYMER_TAIL_BEAD_ATOMS']
    POLYMER_MIDDLE_BEAD_ATOMS = kwargs['POLYMER_MIDDLE_BEAD_ATOMS']

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
