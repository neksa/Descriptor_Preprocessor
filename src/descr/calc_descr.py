import logging
import pandas as pd

from .exceptions import BreakToNextFile
from .get_descr import get_descr
from .descr_config import offsets
# from .loaders import load_all, load_all_after, load_specific
from .write_descr import write_descr

def calc_descrs(data):
    descrs = pd.DataFrame()
    for (filename, sno_markers, cids), dfs in data.items():
        print(filename)
        logging.info(filename)
        try:
            ATOM, HETATM, hb2 = dfs
            if hb2.empty:
                logging.warning("{} has empty hbonds.".format(filename))
                raise BreakToNextFile
            # cids = ATOM[ATOM.aname == "CA"].cid.unique()
            param_to_consider = _get_param_to_consider(ATOM, sno_markers, cids)
            for param in param_to_consider:
                logging.info(param)
                cid, dsr_snos, seq_marker = param
                descr = get_descr(ATOM, HETATM, hb2, dsr_snos, cid)
                full_descr = _add_columns(descr, filename, seq_marker, cid)
                descrs = descrs.append(full_descr, ignore_index=True)
        except BreakToNextFile:
            logging.warning("BreakToNextFile called on {}.".format(filename))
            continue
        except:
            print(f"Flagged: {filename}")
            continue
    return descrs

# def calc_descrs_wrapped(filename=False, after=False, to_write=False):
#     """
#     Convenience legacy method.
#     """
#     data = load_pointer_file(ptr_file)
#     descrs = calc_descrs(data)
#     if to_write:
#         for __, descr in descrs.groupby(['filename', 'cid', 'seq_marker']):
#             write_descr(descr)
#     return descrs

def _get_param_to_consider(ATOM, sno_markers, cids):
    param_to_consider = []

    for cid in cids:
        for marker in sno_markers:
            dsr_snos = _get_sno_range(ATOM, cid, marker)
            if not dsr_snos:
                continue
            param_to_consider.append((cid, dsr_snos, marker))
    return param_to_consider

def _get_sno_range(ATOM, cid, seq_marker):
    start_sno = seq_marker + offsets[0]
    end_sno = seq_marker + offsets[1]
    while ATOM[(ATOM.cid == cid) & (ATOM.sno == start_sno)].empty:
        start_sno += 1
        if start_sno == end_sno:
            msg = "No ATOM lines found in dsr_snos range {}-{}.".format(
                seq_marker + offsets[0], seq_marker + offsets[1])
            logging.warning(msg)
            return False
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