import logging
import pandas as pd

import config
from utils import exceptions
from descr import dihedrals, contacts, hbonds

def calculate(pdb_info_data_map):
    descrs = pd.DataFrame()
    for pdb_info, pdb_data in pdb_info_data_map.items():
        filename, sno_marker, cid = pdb_info
        ATOM, HETATM, hb2 = pdb_data
        dsr_snos = _get_sno_range(ATOM, cid, sno_marker)
        if not dsr_snos:
            continue

        res, C, CA, N = _from_considered_elements(ATOM, dsr_snos, cid)
        pept_bonds = _get_pept_bonds(CA, dsr_snos)
        # For filling descr df
        res_CA = _get_res_CA(res, CA, dsr_snos)
        angles, CA = dihedrals.get_descr_dihedrals(C, CA, N, dsr_snos)
        if hb2 is None:
            hbond_descr = hbonds.build_empty_hb_descr(dsr_snos)
        else:
            assert not hb2.empty
            hbond_descr = hbonds.get_descr_hb(hb2, ATOM, HETATM, dsr_snos)

        heavy_atom_contacts, hetatom_contacts, hetatom_covalent =
        contacts.get_contacts(ATOM, HETATM, cid, dsr_snos)

        # heavy_atom_contacts = contacts.get_heavy_atom_contacts(
        #     ATOM, cid, dsr_snos)
        #
        # if HETATM is None:
        #     hetatom_contacts = contacts.build_empty_cov_or_con(dsr_snos)
        #     hetatom_covalent = contacts.build_empty_cov_or_con(dsr_snos)
        # else:
        #     hetatom_contacts, hetatom_covalent = \
        #         contacts.get_HETATOM_con_cov(ATOM, HETATM, cid, dsr_snos)

        descr = _assemble_descr(hetatom_contacts, hetatom_covalent,
                                heavy_atom_contacts, angles, hbond_descr,
                                res_CA, pept_bonds)

        # descr = get_descr.get_descr(ATOM, HETATM, hb2, dsr_snos, cid)
        full_descr = _add_columns(descr, filename, sno_marker, cid)
        descrs = descrs.append(full_descr, ignore_index=True)
    return descrs

def _get_param_to_consider(ATOM, marker, cids):
    param_to_consider = []
    for cid in cids:
        dsr_snos = _get_sno_range(ATOM, cid, marker)
        if not dsr_snos:
            continue
        param_to_consider.append((cid, dsr_snos, marker))
    return param_to_consider

def _get_sno_range(ATOM, cid, seq_marker):
    start_sno = seq_marker + config.offsets[0]
    end_sno = seq_marker + config.offsets[1]
    while ATOM[(ATOM.cid == cid) & (ATOM.sno == start_sno)].empty:
        start_sno += 1
        if start_sno == end_sno:
            msg = "No ATOM lines found in dsr_snos range {}-{}.".format(
                seq_marker + config.offsets[0], seq_marker + config.offsets[1])
            logging.warning(msg)
            return None
    while ATOM[(ATOM.cid == cid) & (ATOM.sno == end_sno)].empty:
        end_sno -= 1
    dsr_snos = range(start_sno, end_sno)
    return dsr_snos

def _add_columns(descr, filename, seq_marker, cid):
    length = len(descr['sno'])
    descr['filename'] = [filename for __ in range(length)]
    descr['seq_marker'] = [seq_marker for __ in range(length)]
    descr['cid'] = [cid for __ in range(length)]
    descr['relative_sno'] = descr.sno.values - descr.seq_marker.values
    descr = descr.reindex(sorted(descr.columns), axis=1)
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
        c = a - b
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

    ref_length = len(descr['sno'])
    for key, val in descr.items():
        assert len(val) == ref_length
    descr = pd.DataFrame.from_dict(descr)

    return descr