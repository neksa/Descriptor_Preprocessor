from collections import OrderedDict
import logging
import os
import re
import shutil
import subprocess

from preprocessing.converge import conv_to_meme

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

    conv_to_meme.convert(output_4_path, composition_path, output_path)
    assert os.path.isfile(output_path)
    os.remove(output_1_path)
    os.remove(output_4_path)
    os.remove(composition_path)