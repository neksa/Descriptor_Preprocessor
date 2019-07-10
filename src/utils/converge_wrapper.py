from collections import OrderedDict
import logging
import os
import re
import shutil
import subprocess

def run_conv(seed_seqs_path, seqs_path, conv_folder, output_path, bash_exec,
             num_p=7):
    # Input as conv seed seqs, and seqs file paths
    # ioncom_binding_sites.txt to seed need to run separately
    # Output as conv, in meme format.
    # Also need dir for conv folder, where I can assert presence of
    # conv_exec, and blosum file, and assert absence of composition,
    # conv output file, and final output file dir.

    # Check, set up paths
    assert os.path.isfile(seed_seqs_path)
    assert os.path.isfile(seqs_path)
    assert os.path.isdir(conv_folder)
    if os.path.isfile(output_path):
        logging.warning(f"In run_conv(), output <{output_path}> is not "
                        f"empty. Deleting.\n")
        os.remove(output_path)
    conv_exec = os.path.join(conv_folder, 'converge')
    assert os.path.isfile(conv_exec)
    blosum_path = os.path.join(conv_folder, 'BLOSUM62')
    assert os.path.isfile(blosum_path)
    composition_path = os.path.join(conv_folder, 'composition.txt')
    output_1_path = os.path.join(conv_folder, 'output.1.matrix.0')
    if os.path.isfile(output_1_path):
        logging.warning(f"In run_conv(), conv_output_1 <{output_1_path}> "
                        f"is not empty. Deleting.\n")
        os.remove(output_1_path)
    output_4_path = os.path.join(conv_folder, 'output.4.matrix.0')
    if os.path.isfile(output_4_path):
        logging.warning(f"In run_conv(), conv_output_4 <{output_4_path}> "
                        f"is not empty. Deleting.\n")
        os.remove(output_4_path)
    if os.path.isfile(composition_path):
        logging.warning(f"In run_conv(), composition <{composition_path}>"
                        f" is not empty. Deleting.\n")
        os.remove(composition_path)
    # Converge requires files to be in converge folder (I think).
    seed_filename = seed_seqs_path.rsplit("/", maxsplit=1)[-1]
    new_seed_path = os.path.join(conv_folder, seed_filename)
    seqs_filename = seqs_path.rsplit("/", maxsplit=1)[-1]
    new_seqs_path = os.path.join(conv_folder, seqs_filename)
    # End

    command = f"cd {conv_folder} && " \
        f"mpirun -np {num_p} ./converge -B -E 1 -r 1 " \
        f"-f 1 -c ./composition.txt -i ./{seed_filename} -p ./{seqs_filename}"
    shutil.move(seed_seqs_path, new_seed_path)
    shutil.move(seqs_path, new_seqs_path)
    try:
        subprocess.run(command, shell=True, executable=bash_exec,
                       stdout=subprocess.DEVNULL)
    except Exception as e:
        raise e
    finally:
        shutil.move(new_seed_path, seed_seqs_path)
        shutil.move(new_seqs_path, seqs_path)
    if (not os.path.isfile(output_1_path) or
        not os.path.isfile(output_4_path) or
        not os.path.isfile(composition_path)):
        logging.error(f"In folder {conv_folder}, converge failed to produce "
                      f"output files. ")
        raise Exception

    _converge_to_minimal(output_4_path, composition_path, output_path)
    assert os.path.isfile(output_path)
    return

def _converge_to_minimal(input_conv_path, composition_path, output_path):
    # Converge output to minimal
    # Converts converge motif format to minimal meme format
    # see http://meme-suite.org/doc/examples/sample-protein-motif.meme

    # input_conv=''output.4.matrix.0''
    # composition='composition.txt'
    # output="meme_format.txt"
    alphabets, matrices = _parse_converge_output(input_conv_path)
    composition_map = _parse_converge_composition(composition_path)
    _format_minimal_from_conv(alphabets, composition_map, matrices, output_path)
    return

def _parse_converge_output(filename):
    alphabets = ""
    length = 30
    matrices = OrderedDict()
    matrix = []
    nsite = 0
    matrix_count = 0
    with open(filename, "r") as file:
        for line in file:
            if line.startswith("BEGIN") and matrix_count != 0:
                assert len(matrix) == length, len(matrix)
                motif_name = f"MEME-{matrix_count}"
                matrices[motif_name] = (nsite, matrix)
                assert nsite != 0
                matrix = []
                nsite = 0
                continue
            if line.startswith("MATRIX"):
                matrix_count += 1
                match = re.search(r"K=([0-9]+)", line)
                if match is None:
                    raise AssertionError
                nsite = int(match[1])
                continue
            if (line.startswith("50") or line.startswith("30")):
                if not alphabets:
                    matched_alphabets = re.findall("[A-Z]", line)
                    alphabets = "".join(matched_alphabets)
                continue
            if re.match(" [0-9]", line) or re.match("[0-9]+", line):
                probs = re.findall(r"[0-1]\.[0-9]+", line)
                assert len(probs) == len(alphabets)
                matrix.append(probs)
                continue
    return alphabets, matrices

def _parse_converge_composition(filename):
    composition_map = dict()
    with open(filename, "r") as file:
        for line in file:
            if re.match("[A-Z]", line):
                alphabet = line[0]
                composition = line[2:]
                composition_map[alphabet] = float(composition)
                continue
    summed_composition = sum(composition_map.values())
    for key, value in composition_map.items():
        composition_map[key] = value / summed_composition
    return composition_map

def _format_minimal_from_conv(alphabets, composition_map, matrices, output):
    m_to_write = list(range(len(matrices)))
    with open(output, 'w') as file:
        file.write("MEME version 4\n")
        file.write("\n")
        file.write("ALPHABET= " + alphabets + "\n")
        file.write("\n")
        file.write("Background letter frequencies\n")
        for i, alphabet in enumerate(alphabets):
            composition = round(composition_map[alphabet], 4)
            file.write(f"{alphabet} {composition} ")
            if (i != 0) and (i % 9 == 0):
                file.write("\n")
        file.write("\n")
        file.write("\n")
        m_count = 0
        while matrices:
            motif_name, (nsite, matrix) = matrices.popitem(last=False)
            if m_count not in m_to_write:
                m_count += 1
                continue
            m_count += 1
            file.write(f"MOTIF {motif_name}")
            file.write("\n")
            file.write(f"letter-probability matrix: alength= 20 w= 30 "
                       f"nsites= {nsite} E= 0.000")  # alength = len(alphabets)
            # E is just some random number, filled in by subsequent eval calc.
            # w = width of motif
            file.write("\n")
            for line in matrix:
                to_write = ""
                for prob in line:
                    to_write += prob + " "
                file.write(to_write)
                file.write("\n")
            file.write("\n")
    return

from config import paths

seed_seqs_path = paths.ROOT + '/data/seed_seqs.fasta'
seqs_path = paths.ROOT + '/data/input/mg_full/mg_50.fasta'
conv_folder = paths.ROOT + '/external/converge'

output_path = paths.ROOT + '/output_mine.txt'
bash_exec = "/bin/bash"
run_conv(seed_seqs_path, seqs_path, conv_folder, output_path, bash_exec,
         num_p=7)