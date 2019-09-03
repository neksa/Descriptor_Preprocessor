import find_cid_from_pname
from utils import logs

from config import paths
import preprocess
import os
from collections import Counter, OrderedDict

# todo: 1. From this: https://prosite.expasy.org/PS00164#TP, Get the aligned
#  sequences (the short pattern) in .fasta file, run meme on it to construct
#  a meme.txt, use mast on the meme txt, cut out the left-and-right portion
#  to make a length-30, and use that as my motif_pos.

# todo: meme takes forever. An alternative is to convert the aligned seqs
#  into a matrix, then call matrix2meme.

# todo: need to build a composition file using the full seq first.

import numpy as np
import pickle

from utils import generic, build_composition, build_meme_from_aligned, \
    get_pname_seq
"""
1. Take aligned profile, build matrix
2. Download uniprot...?
3. Converge profile on uniprot...?
4. Collect sequences that match
5. Download .pdb for those that match => can we therefore run converge on 
.pdb list only?
6. Somehow cluster the sequences
7. Build for each cluster, the descr plot.

Pain points:
1. How long will running conv on uniprot take?
2. Suggestion to do a c++ rewrite of converge?

"""
# todo: preprocess, in run_prosite_aligned_cropped_flat(), settle todo.
# 2. send matrix in as seed alignment
# 3

# def derive_motif_pos_from_prosite():
#     seq_file = "something"
#     aligned_seq_file = "something"
#
#     composition_file = "something"
#     build_composition.build(seq_file, composition_file)
#
#     meme_output = "meme.txt"
#     build_meme_from_aligned.build(aligned_seq_file, composition_file,
#                                  meme_output)
#
#     pname_seq_map = get_pname_seq.parse(seq_file)
#     download_pdb_files.download(pname_seq_map)
#     download_pdb_files.trim_pname_cid(pname_seq_map, paths.PDB_FOLDER)
#     pname_cid_map = find_cid_from_pname.find(pname_seq_map, paths.PDB_FOLDER)
#     with open(pname_cid_map, "wb") as file:
#
#
#     motif_pos_path = "something"
#
#     find_motifs_mast(pname_cid_path, seq_file, meme_output, 30, motif_pos_path, meme_folder=paths.MEME_MAST_FOLDER)
#     with open(motif_pos_path, 'rb') as file:
#       print(pickle.load(file))
#
#     print("\n")
#     load_pdb_info(motif_pos_path, output)
#     if storage_path is None:
#         shutil.move(conv_meme_file, paths.TRASH)
#         shutil.move(meme_folder, paths.TRASH)
#         shutil.move(motif_pos_path, paths.TRASH)


def matrix_builder(aligned_path):
    alphabets = set(generic.AA3_to_AA1.values())
    AA_to_index = {AA: i for i, AA in enumerate(sorted(alphabets))}
    matrix_counter = None
    with open(aligned_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            if matrix_counter is None:
                matrix_counter = np.zeros((len(line), len(alphabets)),
                                          dtype=int)
            for i, char in enumerate(line):
                try:
                    AA_index = AA_to_index[char]
                except KeyError:

                    # Key not found for some reason
                    continue
                matrix_counter[i, AA_index] += 1

    # matrix_counter = None
    # with open(aligned_path, 'r') as file:
    #     for line in file:
    #         if line.startswith(">"):
    #             continue
    #         line = line.strip().upper()
    #         if matrix_counter is None:
    #             matrix_counter = [Counter() for __ in range(len(line))]
    #         for i, char in enumerate(line):
    #             matrix_counter[i].update([char])
    # matrix_ordered = []
    # for counter in matrix_counter:
    #     matrix_ordered.append(OrderedDict(sorted(counter.items())))
    return matrix_counter

def write_matrix_file(matrix_ordered, output):
    """
    For meme_suite matrix2meme
    """
    output_lines = []
    for AA_counts in matrix_ordered:
        output_lines.append(" ".join(str(i) for i in AA_counts))
    single_str_line = "\n".join(output_lines)
    generic.warn_if_exist(output)
    with open(output, 'w') as file:
        file.write(single_str_line)


def convert_nbdb_matrix_to_conv_encoder(nbdb_file, num_seqs):
    matrix = [[] for __ in range(50)]
    with open(nbdb_file, 'r') as file:
        for i, line in enumerate(file):
            probs = line.split(" ")[:20]
            for prob in probs:
                matrix[i].append(round(float(prob) * num_seqs))
    for i in matrix:
        if sum(i) > num_seqs:
            if i[0] > 0:
                i[0] -= sum(i) - num_seqs
            else:
                i[1] -= sum(i) - num_seqs
        elif sum(i) < num_seqs:
            i[0] += num_seqs - sum(i)
    for i in matrix:
        assert sum(i) == num_seqs
    return matrix


def main():
    """
    Two possible workflows at this point:
    1. Prosite
        1. Download aligned_seqs
        2. Run run_prosite_aligned_cropped, output as number of Kmatches,
        and input matrix for converge.
        3. Encode into output_matrix_binary using the output_matrix.
        4. Run converge using that output_matrix_binary as seed_alignment, on uniprot.
        Manually update Kmatches in converge code.
        5. With converge output and composition.txt, run conv_to_meme to
        convert to meme.
        6. Run ceqlogo on meme to get profile.

    2. NBDB:
        1. Download NBDB matrix file
        2. Manually get number of Kmatches, from length of matching seqs.
        3. Convert this NBDB matrix file to converge seed matrix file,
        using convert_nbdb_matrix_to_conv_encoder()
        4. See 3-6 above.

    3. todo: remember that igor also wants one that's wider, at length 50.
    One for flat pssm for the outer region, and one for included region.
    """
    # For GxGGxG: 2602
    # For GxGxxG: 4077
    # For GxxGxG: 4872
    # todo: check how dependent Kmatches here is. Also, since pdb with
    #  different cid but same pdb_id are probably the same thing, should the
    #  Kmatches then be based on number of unique pdb_id rather than all...?
    # logs.set_logging_level()
    # matrix = convert_nbdb_matrix_to_conv_encoder(os.path.join(
    #     paths.ROOT, "GxGxxG_pssm.txt"), 4077)
    # with open("output_matrix.txt", "w") as file:
    #     for line in matrix:
    #         for value in line:
    #             file.write(str(value) + " ")
    #         file.write("\n")


    file_folder = os.path.join(paths.ROOT, "Reports", "19082019", "PDOC00706")

    full_seq_path = os.path.join(file_folder, "full_seq.fasta")
    aligned_seq_path = os.path.join(file_folder, "aligner.txt")
    output_matrix_path = os.path.join(file_folder, "output_matrix.txt")
    # preprocess.run_prosite_aligned_cropped_flat(full_seq_path, aligned_seq_path,
    #                                        output_matrix_path)
    # preprocess.run_prosite_aligned_cropped(full_seq_path, aligned_seq_path,
    #                                             output_matrix_path)
    # import sys
    # sys.exit()

    # 1. move output_matrix.txt to converge_encoder. Run converge_encoder.
    # 2. Move input_matrix_binary into conv_rewrite. Take note of Kmatches
    # here, and edit in single_Kmatches in conv_rewrite. Run conv_rewrite.
    # 3. Move output.1.matrix and composition.txt back to Reports folder,
    # run the rest below.




    # After that
    from converge.conv_interface import convert_conv_to_meme_full_num
    converge_output_path = os.path.join(file_folder, "output.1.matrix")
    composition_path = os.path.join(file_folder, "composition.txt")
    output_meme_path = os.path.join(file_folder, "meme.txt")
    convert_conv_to_meme_full_num(converge_output_path, composition_path,
                                  output_meme_path)

    logo_path = os.path.join(file_folder, "logo.eps")
    import subprocess
    ceqlogo_exec = os.path.join(paths.ROOT, "src", "meme_suite", "meme",
                                "ceqlogo")
    command = f"{ceqlogo_exec} -i1 {output_meme_path} -o {logo_path} -f EPS"
    subprocess.run(command, shell=True)

    tmp_composition_path = os.path.join(file_folder, "compos_align.txt")
    aligned_meme_path = os.path.join(file_folder, "aligned_meme.txt")
    build_composition.build(aligned_seq_path, tmp_composition_path)
    build_meme_from_aligned.build(aligned_seq_path, aligned_meme_path,
                                  tmp_composition_path)
    aligned_logo_path = os.path.join(file_folder, "aligned_logo.eps")
    command = f"{ceqlogo_exec} -i1 {aligned_meme_path} -o {aligned_logo_path} -f EPS"
    subprocess.run(command, shell=True)


    # # with open(os.path.join(paths.ROOT, "motif_pos.pkl"), 'rb') as file:
    # #     print(len(pickle.load(file)))
    # preprocess.run_prosite_aligned_cropped(paths.PROSITE_ENOLASE_SEQS,
    #                                paths.PROSITE_ALIGNED_SEQS,
    #                                os.path.join(paths.ROOT,
    #                                             "output_motif_pos_new.pkl"))



    # preprocess.run_prosite_mast(paths.PROSITE_EXTRACT, 30,
    #                             paths.REF_MEME_TXT,
    #                  os.path.join(paths.ROOT, "prosite_mast_motif_pos.pkl"))
    # preprocess.run_prosite_aligned(paths.PROSITE_ENOLASE_SEQS,
    #                                paths.PROSITE_ALIGNED_SEQS,
    #                                os.path.join(paths.ROOT,
    #                                             "output_motif_pos.pkl"))

    # seq_path = os.path.join(paths.USER_INPUT, "aligner.txt")
    # matrix_ordered = matrix_builder(seq_path)
    # write_matrix_file(matrix_ordered,
    #                   os.path.join(paths.ROOT, "mine__.txt"))
    # print(matrix_ordered)
    # import sys
    # sys.exit()
    # logs.set_logging_level()
    #
    # from utils import get_pname_seq
    #
    # # pname_seq_map = get_pname_seq.parse(seq_path)
    # from meme_suite import meme_interface
    # meme_interface.create_meme_from_aligned(seq_path, 14, paths.MEME_MAST_FOLDER, num_p=7)
    # meme_txt = os.path.join(paths.MEME_MAST_FOLDER, "meme.txt")


    # logs.set_logging_level()
    # extract_path = paths.PROSITE_EXTRACT
    # motif_len = 13
    # output = paths.PID_PDB_MAP
    # num_p = 7
    # preprocess.run_prosite_meme(extract_path, motif_len, output, num_p)

    # extract_path = paths.PROSITE_EXTRACT
    # motif_len = 13
    # ref_meme_txt = paths.REF_MEME_TXT
    # print(ref_meme_txt)
    # output = paths.PID_PDB_MAP
    # preprocess.run_prosite_mast(extract_path, motif_len, ref_meme_txt, output)
    # with open(output, 'rb') as file:
    #     print(pickle.load(file))
    #
    # extract_path = paths.IONCOM_EXTRACT
    # motif_len = 13
    # ref_meme_txt = paths.REF_MEME_TXT
    # output = paths.PID_PDB_MAP
    # preprocess.run_ioncom_mast(extract_path, motif_len, ref_meme_txt, output)
    pass


if __name__ == "__main__":
    main()