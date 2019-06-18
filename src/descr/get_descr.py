from collections import defaultdict
import numpy as np
import pandas as pd

from descr import contacts, dihedrals, hbonds

def get_descr(ATOM, HETATM, hb2, dsr_snos, cid):
    """
    Original form have these lines: "ATOM", "ANISOU", "HETATM", "TER", but
    ANISOU and TER are not used, so only keeping ATOM and HETATM.
    """
    res, C, CA, N = _from_considered_elements(ATOM, dsr_snos, cid)
    pept_bonds = _get_pept_bonds(CA, dsr_snos)

    # For filling descr df
    res_CA = _get_res_CA(res, CA, dsr_snos)

    angles, CA = dihedrals.get_descr_dihedrals(C, CA, N, dsr_snos)

    hbond_descr = hbonds.get_descr_hb(hb2, ATOM, HETATM, dsr_snos)

    hetatom_contacts, hetatom_covalent, heavy_atom_contacts = \
        contacts.get_contacts(ATOM, HETATM, cid, dsr_snos)

    descr = _assemble_descr(hetatom_contacts, hetatom_covalent,
                                 heavy_atom_contacts, angles, hbond_descr,
                                 res_CA, pept_bonds)
    return descr

def _get_pept_bonds(CA, dsr_snos):
    """
    This consider that there exist a peptide bond, if two residues with same
    cid, aname is CA, and resi in resi_list, and are adjacent to each other in
    CA, are at a distance of less than 4, in x/y/z coord units.
    Return a set() of indices (i, i+1) indicating the pair of linked atoms.
    Relative to position along dsr_snos.
    """
    peptide_pairs = set()
    for i in range(len(CA) - 1):
        a = CA[i]
        b = CA[i + 1]
        c = a-b
        if np.sqrt(np.einsum('i,i', c, c)) < 4:
            peptide_pairs.add((i, i + 1))

    matrix = [[True if (i, j) in peptide_pairs else False
               for j in range(len(dsr_snos))] for i in range(len(dsr_snos))]

    peptide_bonds = dict()
    peptide_bonds['sno'] = dsr_snos
    peptide_bonds['pept_bonds'] = matrix

    return peptide_bonds

def _get_res_CA(res, CA, dsr_snos):
    res_CA = defaultdict(list)
    for i in range(len(res)):
        res_CA['sno'].append(dsr_snos[i])
        res_CA['res'].append(res[i])
        res_CA['CA'].append(CA[i])
    return res_CA

def _from_considered_elements(ATOM, dsr_snos, cid):
    _ATOM = ATOM.filter(['cid', 'sno', 'aname', 'coord', 'res'])
    _ATOM = _ATOM[(_ATOM.sno.isin(dsr_snos)) &
                  (_ATOM.aname.isin(("N", "C", "CA"))) &
                  (_ATOM.cid == cid)]

    _ATOM = _ATOM.set_index(['sno', 'aname'], drop=False)
    _ATOM = _ATOM.sort_index(axis=0, sort_remaining=True)
    np.testing.assert_array_equal(_ATOM.aname.values[:3],
                                  np.array(['C', 'CA', 'N']))
    np.testing.assert_array_equal(_ATOM.aname.values[3:6],
                                  np.array(['C', 'CA', 'N']))
    res = np.array([i for i in _ATOM.res.values[::3]])
    coords = _ATOM.coord.values
    C = np.array([i for i in coords[::3]])
    CA = np.array([i for i in coords[1::3]])
    N = np.array([i for i in coords[2::3]])

    assert len(res) == len(dsr_snos)
    assert len(C) == len(dsr_snos)

    return res, C, CA, N

def _assemble_descr(hetatom_contacts, hetatom_covalent,
                                 heavy_atom_contacts, angles, hbond_descr,
                                 res_CA, pept_bonds):
    descr = dict()
    descr.update(hetatom_contacts)
    descr.update(hetatom_covalent)
    descr.update(heavy_atom_contacts)
    descr.update(angles)
    descr.update(hbond_descr)
    descr.update(res_CA)
    descr.update(pept_bonds)

    if __debug__:
        ref_length = len(descr['sno'])
        for key, val in descr.items():
            assert len(val) == ref_length
    descr = pd.DataFrame.from_dict(descr)

    return descr

