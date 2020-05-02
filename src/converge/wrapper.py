import logging
import os
import shutil
import subprocess

from converge import conv_to_meme
from config import paths

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


def converge_calculate(profile_length, kmatches, iteration=5,
                       maxS_start_factor=1, min_match_percent=0.7,
                       maxS_factor_decrement=0.7):
    assert iteration >= 1
    assert profile_length >= 1
    assert kmatches >= 1
    assert os.path.isfile(paths.CONV_INPUT_MATRIX)
    if os.path.isfile(paths.CONV_OUTPUT):
        logging.warning(f"In converge_calculate(), output <{paths.CONV_OUTPUT}> is not "
                        f"empty. Deleting.\n")
        os.remove(paths.CONV_OUTPUT)
    if os.path.isfile(paths.CONV_COMPOSITION):
        logging.warning(
            f"In converge_calculate(), output <{paths.CONV_COMPOSITION}> is not "
            f"empty. Deleting.\n")
        os.remove(paths.CONV_COMPOSITION)
    assert os.path.isfile(paths.CONV_EXEC)
    min_matches = min_match_percent * kmatches
    # todo: need to call converge_encoder first, then converge

    command = f"{paths.CONV_EXEC} {profile_length} {kmatches} {min_matches} " \
              f"{maxS_start_factor}"
    while iteration != 0:
        return_code = subprocess.run(command, shell=True,
                                     executable=paths.BASH_EXEC).returncode
        if return_code == 2:
            iteration -= 1
            maxS_start_factor *= maxS_factor_decrement
            command = f"{paths.CONV_EXEC} {profile_length} {kmatches} " \
                      f"{min_matches} {maxS_start_factor}"
            continue
        if return_code == 0:
            break
        else:
            raise Exception