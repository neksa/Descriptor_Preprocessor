import os

from global_config import output_dir


def write_descr(d):
    d = d.drop(columns=['a_aname', 'a_cid', 'a_d_dd', 'a_res', 'a_sno',
                        'd_a_aa', 'd_a_dist', 'd_aname', 'd_cid', 'd_res',
                        'd_sno', 'pdb_id', 'planar1', 'planar2'])
    d = d.reindex(sorted(d.columns), axis=1)
    assert tuple(d.keys()) == ('CA', 'acc', 'category', 'cid', 'contact',
                               'covalent', 'donor', 'ext', 'filename',
                               'h_contacts', 'pept_bonds', 'phi', 'psi',
                               'region', 'relative_sno', 'res', 'role',
                               'seq_marker', 'sno', 'ss')
    values = tuple([value[1].values for value in d.items()])
    CA, acc, category, cid, contact, covalent, donor, ext, filename, \
    h_contacts, pept_bonds, phi, psi, region, relative_sno, res, role, \
    seq_marker, sno, ss = values

    assert len(set([len(i) for i in values])) == 1
    assert len(set(filename)) == 1
    assert len(set(cid)) == 1
    assert len(set(seq_marker)) == 1

    filename = filename[0]
    cid = cid[0]
    seq_marker = seq_marker[0]
    fname = os.path.join(output_dir,
                         "DES_{}_{}_{}".format(filename, cid, seq_marker))
    with open(fname + ".txt", "w") as f:
        f.write("#DESCRIPTOR V0.1\n")
        f.write("HEADER %s %s %d %d\n\n" % (filename, cid, min(sno),
                                            max(sno)))
        for i in range(len(sno)):
            f.write("NODE {}\n".format(sno[i]))
            f.write("NODE.RESN {}\n".format(res[i]))
            f.write("NODE.ROLE \n") # chemical_role
            f.write("NODE.CA {:.2f}, {:.2f}, {:.2f}\n".format(*CA[i]))
            f.write("NODE.DIHEDRAL {:.2f}, {:.2f}, {}, {}\n"
                    .format(phi[i], psi[i], region[i], ss[i]))

            if contact[i] == 0:
                f.write("NODE.HET_CONTACT \n")
            else:
                f.write("NODE.HET_CONTACT {}\n".format(contact[i]))
            if covalent[i] == 0:
                f.write("NODE.HET_COVALENT \n")
            else:
                f.write("NODE.HET_COVALENT {}\n".format(covalent[i]))
            for _ext, _role, _category, _donor, _acc \
                    in zip(ext[i], role[i], category[i], donor[i], acc[i]):
                f.write("NODE.{}.HBOND {}; {}; {}, {}, {}; {}, {}, {}\n"
                    .format(_ext, _role, _category, *_donor, *_acc))
            f.write("\n")
        f.write("\n")

        f.write("MATRIX.PEPT \t{}\n".format("\t".join([str(i) for i in sno])))

        for i, row in enumerate(pept_bonds):
            _row = _convert_to_str(row)
            f.write("MATRIX.DATA \t{}\t{}\n".format(sno[i], _row))

        f.write("\n\n")

        f.write("MATRIX.CONT \t{}\n".format("\t".join([str(k) for k in sno])))
        for i, row in enumerate(h_contacts):
            _row = _convert_to_str(row)
            f.write("MATRIX.DATA \t{}\t{}\n".format(sno[i], _row))

    with open(fname + ".contacts", "w") as f_cont:
        for row in h_contacts:
            _row = _convert_to_str(row)
            f_cont.write("{}\n".format(_row))

def _convert_to_str(row):
    string = "\t".join([str(i) for i in row]).replace("True", "1")\
                                             .replace("False", ".")
    return string