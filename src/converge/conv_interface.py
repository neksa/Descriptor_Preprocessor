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

from converge import meme_to_conv

def convert_meme_to_conv(meme, composition, matrix, minimal=False):
    if minimal:
        ret_code = meme_to_conv.convert_minimal(meme, composition, matrix)
    else:
        ret_code = meme_to_conv.convert_full(meme, composition, matrix)
    return ret_code

def encode_input_matrix(matrix, filename):
    with open(filename, 'w') as file:
        for line in matrix:
            for term in line[:-1]:
                file.write(str(term)+",")
            file.write(str(line[-1]))
        file.write("\n")


def run(input_matrix, profile_length, kmatches, output_meme,
        proteome_binary=paths.UNIPROT_BINARY, iteration=5):
    assert os.path.isfile(proteome_binary)
    input_matrix_file = paths.CONV_INPUT_MATRIX
    output_matrix = paths.CONV_OUTPUT
    output_composition = paths.CONV_COMPOSITION

    # converge_encoder can be removed from the whole loop
    encode_input_matrix(input_matrix, input_matrix_file)
    # encode_matrix(input_matrix, matrix_binary)
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    wrapper.converge_calculate(profile_length, kmatches, iteration=iteration)
    conv_to_meme.convert(output_matrix, output_composition, output_meme)