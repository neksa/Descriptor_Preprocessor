import logging
import numpy as np

import config
from utils import exceptions

# def get_contacts(ATOM, HETATM, cid, dsr_snos):
#     heavy_atoms_1 = ATOM[(ATOM.aname != "H")
#                        & (ATOM.cid == cid)
#                        & (ATOM.sno.isin(dsr_snos))]
#
#     heavy_atoms_2_ATOM = ATOM[(ATOM.aname != "H")
#                               & (ATOM.res != "HOH")
#                               & (ATOM.sno.isin(dsr_snos))]
#
#     a1_sno = heavy_atoms_1.sno.values
#
#     a2_ATOM_sno = heavy_atoms_2_ATOM.sno.values
#
#     a1_coord_flat = np.concatenate(heavy_atoms_1.coord.values)
#     a1_coord2 = a1_coord_flat.reshape(len(heavy_atoms_1), 3)
#
#     # np.concatenate to convert pd.df.values into normal np.array
#     a2_ATOM_coord_flat = np.concatenate(heavy_atoms_2_ATOM.coord.values)
#     a2_ATOM_coord2 = a2_ATOM_coord_flat.reshape(len(heavy_atoms_2_ATOM), 3)
#
#     heavy_atom_contacts = _get_heavy_atom_contacts(
#         a1_coord2, a1_sno, a2_ATOM_coord2, a2_ATOM_sno, dsr_snos)
#
#     if HETATM.empty:
#         hetatom_covalent = dict()
#         hetatom_covalent['sno'] = dsr_snos
#         hetatom_covalent['covalent'] = np.zeros(len(dsr_snos), dtype='int64')
#
#         hetatom_contacts = dict()
#         hetatom_contacts['sno'] = dsr_snos
#         hetatom_contacts['contact'] = np.zeros(len(dsr_snos), dtype='int64')
#         return hetatom_contacts, hetatom_covalent, heavy_atom_contacts
#     else:
#         heavy_atoms_2_HETATM = HETATM[(HETATM.aname != "H")
#                                       & (HETATM.res != "HOH")]
#     if heavy_atoms_2_HETATM.empty:
#         hetatom_covalent = dict()
#         hetatom_covalent['sno'] = dsr_snos
#         hetatom_covalent['covalent'] = np.zeros(len(dsr_snos), dtype='int64')
#
#         hetatom_contacts = dict()
#         hetatom_contacts['sno'] = dsr_snos
#         hetatom_contacts['contact'] = np.zeros(len(dsr_snos), dtype='int64')
#         return hetatom_contacts, hetatom_covalent, heavy_atom_contacts
#     else:
#         a2_HETATM_coord = heavy_atoms_2_HETATM.coord.values
#         a2_HETATM_coord2 = np.concatenate([i for i in a2_HETATM_coord]).reshape(len(a2_HETATM_coord), 3)
#         hetatom_contacts, hetatom_covalent = _get_HETATOM_contacts(
#             a1_coord2, a1_sno, a2_HETATM_coord2, dsr_snos)
#     return hetatom_contacts, hetatom_covalent, heavy_atom_contacts


def get_contacts(ATOM, HETATM, cid, dsr_snos):
    """
    Variable names taken from original code. I've no idea what's going on,
    so they are kept as-is.

    For con/cov calculation, code is equivalent to:
    hetatom_contacts = set()
    hetatom_covalent = set()
    for a1_sno, a1_coord in zip(a1_snos, a1_coords):
        for a2_coord in a2_HETATM_coords:
            if dist(a1_coord, a2_coord) < HC_THRESHOLD:
                hetatom_contacts.add(a1_sno)
            if dist(a1_coord, a2_coord) < CC_THRESHOLD:
                hetatom_covalent.add(a1_sno)
    return hetatom_contacts, hetatom_covalent

    Speeds up by a factor of 4 when tested on len(a1_coords) = 258
    and len(a2_HETATM_coords)=5
    """
    heavy_atoms_1 = ATOM[(ATOM.aname != "H")
                         & (ATOM.cid == cid)
                         & (ATOM.sno.isin(dsr_snos))]

    heavy_atoms_2_ATOM = ATOM[(ATOM.aname != "H")
                              & (ATOM.res != "HOH")
                              & (ATOM.sno.isin(dsr_snos))]

    a1_sno = heavy_atoms_1.sno.values

    a2_ATOM_sno = heavy_atoms_2_ATOM.sno.values

    a1_coord_flat = np.concatenate(heavy_atoms_1.coord.values)
    a1_coord2 = a1_coord_flat.reshape(len(heavy_atoms_1), 3)

    # np.concatenate to convert pd.df.values into normal np.array
    a2_ATOM_coord_flat = np.concatenate(heavy_atoms_2_ATOM.coord.values)
    a2_ATOM_coord2 = a2_ATOM_coord_flat.reshape(len(heavy_atoms_2_ATOM), 3)

    heavy_atom_contacts = _get_heavy_atom_contacts(
        a1_coord2, a1_sno, a2_ATOM_coord2, a2_ATOM_sno, dsr_snos)

    if HETATM is None:
        hetatom_contacts = build_empty_cov_or_con(dsr_snos)
        hetatom_covalent = build_empty_cov_or_con(dsr_snos)
        return heavy_atom_contacts, hetatom_contacts, hetatom_covalent

    heavy_atoms_2_HETATM = HETATM[(HETATM.aname != "H") & (HETATM.res != "HOH")]

    if heavy_atoms_2_HETATM.empty:
        hetatom_contacts = build_empty_cov_or_con(dsr_snos)
        hetatom_covalent = build_empty_cov_or_con(dsr_snos)
        return heavy_atom_contacts, hetatom_contacts, hetatom_covalent

    a2_HETATM_coord = heavy_atoms_2_HETATM.coord.values
    a2_HETATM_coord2 = np.concatenate([i for i in a2_HETATM_coord]).reshape(
        len(a2_HETATM_coord), 3)

    a1_coords = a1_coord2
    a1_snos = a1_sno
    a2_HETATM_coords = a2_HETATM_coord2

    # Get distance between points
    # a1_coords.shape = (258, 3)
    # a2_HETATM_coords.shape = (5, 3)
    # squared_dist.shape = (258, 5, 3)
    coordinates = np.array([a1 - a2_HETATM_coords for a1 in a1_coords])

    # Mean-squared distance
    # distance.shape = (258, 5)
    distance = np.sqrt(np.einsum('ijk,ijk->ij', coordinates, coordinates))

    hetatom_contacts = dict()
    hetatom_contacts['sno'] = dsr_snos
    hetatom_contacts['contact'] = _retrieve_num_bonds(
        distance, config.HC_THRESHOLD, a1_snos, dsr_snos)

    hetatom_covalent = dict()
    hetatom_covalent['sno'] = dsr_snos
    hetatom_covalent['covalent'] = _retrieve_num_bonds(
        distance, config.CC_THRESHOLD, a1_snos, dsr_snos)

    return heavy_atom_contacts, hetatom_contacts, hetatom_covalent

def get_heavy_atom_contacts(ATOM, cid, dsr_snos):
    heavy_atoms_1 = ATOM[(ATOM.aname != "H")
                         & (ATOM.cid == cid)
                         & (ATOM.sno.isin(dsr_snos))]

    heavy_atoms_2_ATOM = ATOM[(ATOM.aname != "H")
                              & (ATOM.res != "HOH")
                              & (ATOM.sno.isin(dsr_snos))]

    a1_sno = heavy_atoms_1.sno.values

    a2_ATOM_sno = heavy_atoms_2_ATOM.sno.values

    a1_coord_flat = np.concatenate(heavy_atoms_1.coord.values)
    a1_coord2 = a1_coord_flat.reshape(len(heavy_atoms_1), 3)

    # np.concatenate to convert pd.df.values into normal np.array
    a2_ATOM_coord_flat = np.concatenate(heavy_atoms_2_ATOM.coord.values)
    a2_ATOM_coord2 = a2_ATOM_coord_flat.reshape(len(heavy_atoms_2_ATOM), 3)

    heavy_atom_contacts = _get_heavy_atom_contacts(
        a1_coord2, a1_sno, a2_ATOM_coord2, a2_ATOM_sno, dsr_snos)
    return heavy_atom_contacts

def get_HETATOM_con_cov(ATOM, HETATM, cid, dsr_snos):
    """

    """
    heavy_atoms_1 = ATOM[(ATOM.aname != "H")
                         & (ATOM.cid == cid)
                         & (ATOM.sno.isin(dsr_snos))]

    a1_sno = heavy_atoms_1.sno.values

    a1_coord_flat = np.concatenate(heavy_atoms_1.coord.values)
    a1_coord2 = a1_coord_flat.reshape(len(heavy_atoms_1), 3)

    heavy_atoms_2_HETATM = HETATM[(HETATM.aname != "H")
                                  & (HETATM.res != "HOH")]
    if heavy_atoms_2_HETATM.empty:
        hetatom_contacts = build_empty_cov_or_con(dsr_snos)
        return hetatom_contacts
    else:
        a2_HETATM_coord = heavy_atoms_2_HETATM.coord.values
        a2_HETATM_coord2 = np.concatenate([i for i in a2_HETATM_coord]).reshape(
            len(a2_HETATM_coord), 3)

    a1_coords = a1_coord2
    a1_snos = a1_sno
    a2_HETATM_coords = a2_HETATM_coord2

    # Get distance between points
    # a1_coords.shape = (258, 3)
    # a2_HETATM_coords.shape = (5, 3)
    # squared_dist.shape = (258, 5, 3)
    coordinates = np.array([a1 - a2_HETATM_coords for a1 in a1_coords])

    # Mean-squared distance
    # distance.shape = (258, 5)
    distance = np.sqrt(np.einsum('ijk,ijk->ij', coordinates, coordinates))

    hetatom_contacts = dict()
    hetatom_contacts['sno'] = dsr_snos
    hetatom_contacts['contact'] = _retrieve_num_bonds(
        distance, config.HC_THRESHOLD, a1_snos, dsr_snos)

    hetatom_covalent = dict()
    hetatom_covalent['sno'] = dsr_snos
    hetatom_covalent['covalent'] = _retrieve_num_bonds(
        distance, config.CC_THRESHOLD, a1_snos, dsr_snos)

    return hetatom_contacts, hetatom_covalent

def _retrieve_num_bonds(distance, threshold, a1_snos, dsr_snos):
    bool_mask = distance < threshold
    indice_mask = np.argwhere(bool_mask)[:, 0]
    considered_snos = a1_snos[indice_mask]
    bond_pos = np.unique([dsr_snos.index(i) for i in considered_snos],
                      return_counts=True)

    if bond_pos[0].size == 0:
        num_bonds = np.zeros(len(dsr_snos), dtype='int64')
    else:
        num_bonds = np.zeros(len(dsr_snos), dtype='int')
        for sno, count in list(zip(*bond_pos)):
            num_bonds[sno] = count
        num_bonds = np.array(num_bonds, dtype='int64')
    return num_bonds

def build_empty_cov_or_con(dsr_snos):
    empty_cov_con = dict()
    empty_cov_con['sno'] = dsr_snos
    empty_cov_con['covalent'] = np.zeros(len(dsr_snos), dtype='int64')
    return empty_cov_con


def _get_heavy_atom_contacts(a1_coords, a1_snos, a2_ATOM_coord, a2_ATOM_sno,
                            dsr_snos):
    """
    Equivalent except set() to:
    heavy_atom_contacts = []
    for a1_sno, a1_coord in zip(a1_snos, a1_coords):
            if dist(a1_coord, a2_coord) < HC_THRESHOLD:
                heavy_atom_contacts.append(
                    (dsr_snos.index(a1_sno), dsr_snos.index(a2_sno)))
    return heavy_atom_contacts
    """
    assert len(a1_coords) == len(a1_snos)
    assert len(a2_ATOM_sno) == len(a2_ATOM_coord)

    # coords have an additional (,3) behind shape
    dsr_sno_map = dict([(sno, i) for i, sno in enumerate(dsr_snos)])
    value = np.array([a1 - a2_ATOM_coord for a1 in a1_coords])

    # Mean-squared distance
    # value.shape = (258, 258, 3)
    # distance.shape = (258, 258)
    distance = np.sqrt(np.einsum('ijk,ijk->ij', value, value))
    contact_mask = distance < config.HC_THRESHOLD

    # remove duplicate contacts
    # masked_array mask if True, hence need to reverse, to mask the False ones
    contacts = set()
    for a1_sno, contact_for_each_sno in zip(a1_snos, contact_mask):
        matched_snos = np.ma.compressed(
            np.ma.masked_array(a2_ATOM_sno, mask=~contact_for_each_sno))
        for a2_sno in matched_snos:
            contacts.add((dsr_sno_map[a1_sno], dsr_sno_map[a2_sno]))
    heavy_atom_contacts = dict()
    heavy_atom_contacts['sno'] = dsr_snos
    heavy_atom_contacts['h_contacts'] = [[True if (i, j) in contacts else False
                                        for j in range(len(dsr_snos))]
                                        for i in range(len(dsr_snos))]
    return heavy_atom_contacts
