import logging
import os
import shutil
import subprocess

from converge import conv_to_meme

def encode_proteome(proteome_fname, output, conv_folder, bash_exec):
    assert os.path.isfile(proteome_fname)
    assert os.path.isdir(conv_folder)
    if os.path.isfile(output):
        logging.warning(f"In encode_proteome(), output <{output}> is not "
                        f"empty. Deleting.\n")
        os.remove(output)
    conv_exec = os.path.join(conv_folder, 'converge_encoder')
    assert os.path.isfile(conv_exec)
    command = f"{conv_exec} -pi {proteome_fname} -po {output} -silent"
    return_code = subprocess.run(command, shell=True, executable=bash_exec)
    if return_code != 0:
        raise Exception


def encode_blosum(output, conv_folder, bash_exec):
    assert os.path.isdir(conv_folder)
    if os.path.isfile(output):
        logging.warning(f"In encode_blosum(), output <{output}> is not "
                        f"empty. Deleting.\n")
        os.remove(output)
    conv_exec = os.path.join(conv_folder, 'converge_encoder')
    assert os.path.isfile(conv_exec)
    command = f"{conv_exec} -bo {output} -silent"
    return_code = subprocess.run(command, shell=True, executable=bash_exec)
    if return_code != 0:
        raise Exception


def encode_matrix(input_matrix, output, conv_folder, bash_exec):
    assert os.path.isfile(input_matrix)
    assert os.path.isdir(conv_folder)
    if os.path.isfile(output):
        logging.warning(f"In encode_matrix(), output <{output}> is not "
                        f"empty. Deleting.\n")
        os.remove(output)
    conv_exec = os.path.join(conv_folder, 'converge_encoder')
    assert os.path.isfile(conv_exec)
    command = f"{conv_exec} -mi {input_matrix} -mo {output}"
    return_code = subprocess.run(command, shell=True, executable=bash_exec).returncode
    if return_code != 0:
        print(return_code)
        raise Exception


def converge_calculate(profile_length, kmatches, input_matrix_b, proteome_b,
                       output_matrix, output_composition, conv_folder,
                       bash_exec, iteration=5):
    assert iteration >= 1
    assert profile_length >= 1
    assert kmatches >= 1
    assert os.path.isfile(input_matrix_b)
    assert os.path.isfile(proteome_b)
    assert os.path.isdir(conv_folder)
    if os.path.isfile(output_matrix):
        logging.warning(f"In converge_calculate(), output <{output_matrix}> is not "
                        f"empty. Deleting.\n")
        os.remove(output_matrix)
    if os.path.isfile(output_composition):
        logging.warning(
            f"In converge_calculate(), output <{output_composition}> is not "
            f"empty. Deleting.\n")
        os.remove(output_composition)
    conv_exec = os.path.join(conv_folder, 'converge_calculator')
    assert os.path.isfile(conv_exec)
    command = f"{conv_exec} {profile_length} {kmatches} {input_matrix_b} " \
              f"{proteome_b} {output_matrix} {output_composition}"
    while iteration != 0:
        return_code = subprocess.run(command, shell=True, executable=bash_exec).returncode
        if return_code == 2:
            iteration -= 1
            continue
        if return_code == 0:
            break
        else:
            raise Exception