import logging
import os

from config import paths
from preprocessing.converge import make_seed_seqs, wrapper, meme_to_conv, \
    conv_to_meme
from utils import generic


def run(seqs_path, output_path,
        binding_site_file=None, bash_exec=paths.BASH_EXEC, num_p=7,
        seed_seq_path=paths.CONV_SEED_SEQS,
        conv_folder="./"):
    if binding_site_file:
        os.remove(seed_seq_path)
        make_seed_seqs.make(binding_site_file, seed_seq_path)
        generic.quit_if_missing(seed_seq_path)
    else:
        if not os.path.isfile(seed_seq_path):
            logging.error(f"Converge needs either a binding_site file, "
                          f"or a seed_seq file. No seed_seq file in "
                          f"<{seed_seq_path}>.")
            raise Exception
    wrapper.run_conv(seed_seq_path, seqs_path, conv_folder, output_path,
                     bash_exec, num_p=num_p)
    return

def convert_meme_to_conv(meme, composition, matrix, minimal=False):
    generic.quit_if_missing(meme)
    generic.warn_if_exist(composition)
    generic.warn_if_exist(matrix)
    if minimal:
        meme_to_conv.convert_minimal(meme, composition,
                                     matrix)
    else:
        meme_to_conv.convert_full(meme, composition,
                                  matrix)
    generic.quit_if_missing(composition)
    generic.quit_if_missing(matrix)


def convert_conv_to_meme(matrix, composition, meme):
    generic.quit_if_missing(composition)
    generic.quit_if_missing(matrix)
    generic.warn_if_exist(meme)
    conv_to_meme.convert(matrix, composition, meme)
    generic.quit_if_missing(meme)
