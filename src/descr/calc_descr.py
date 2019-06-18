import logging
import pandas as pd

import config
from descr import get_descr
from utils import exceptions

def calc_descrs(data):
    descrs = pd.DataFrame()
    for (filename, sno_markers, cids), dfs in data.items():
        try:
            ATOM, HETATM, hb2 = dfs
            if hb2.empty:
                logging.warning("{} has empty hbonds.".format(filename))
                raise exceptions.BreakToNextFile
            # cids = ATOM[ATOM.aname == "CA"].cid.unique()
            param_to_consider = _get_param_to_consider(ATOM, sno_markers, cids)
            for param in param_to_consider:
                logging.info(param)
                cid, dsr_snos, seq_marker = param
                descr = get_descr.get_descr(ATOM, HETATM, hb2, dsr_snos, cid)
                full_descr = _add_columns(descr, filename, seq_marker, cid)
                descrs = descrs.append(full_descr, ignore_index=True)
        except exceptions.BreakToNextFile:
            logging.info("BreakToNextFile called on {}.".format(filename))
            continue
        except:
            print(f"Flagged: {filename}")
            continue
    return descrs

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
    start_sno = seq_marker + config.offsets[0]
    end_sno = seq_marker + config.offsets[1]
    while ATOM[(ATOM.cid == cid) & (ATOM.sno == start_sno)].empty:
        start_sno += 1
        if start_sno == end_sno:
            msg = "No ATOM lines found in dsr_snos range {}-{}.".format(
                seq_marker + config.offsets[0], seq_marker + config.offsets[1])
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