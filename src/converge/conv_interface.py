import os

from config import paths
from converge import wrapper, conv_to_meme


def encode_proteome(proteome_fname, output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_proteome(proteome_fname, output, conv_folder,
                                   bash_exec)


def encode_blosum(output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_blosum(output, conv_folder,
                                   bash_exec)


def encode_matrix(input_matrix, output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_matrix(input_matrix, output, conv_folder,
                                 bash_exec)


def run(input_matrix, profile_length, kmatches, output_meme,
        proteome_binary=paths.UNIPROT_BINARY, iteration=5):
    assert os.path.isfile(proteome_binary)
    matrix_binary = paths.TMP_FILE_TEMPLATE.format(0)
    output_matrix = paths.TMP_FILE_TEMPLATE.format(1)
    output_composition = paths.TMP_FILE_TEMPLATE.format(2)
    encode_matrix(input_matrix, matrix_binary)
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    wrapper.converge_calculate(profile_length, kmatches, matrix_binary,
                               proteome_binary, output_matrix,
                               output_composition, conv_folder,
                               paths.BASH_EXEC, iteration=iteration)
    conv_to_meme.convert(output_matrix, output_composition, output_meme)