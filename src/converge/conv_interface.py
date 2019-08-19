import logging
import os

from config import paths
from converge import make_seed_seqs, wrapper, meme_to_conv, conv_to_meme
from utils import generic


def run(seqs_path, output_path, bash_exec=paths.BASH_EXEC, num_p=7,
        seed_seq_path=paths.CONV_SEED_SEQS,
        conv_folder=paths.CONV_FOLDER):
    generic.quit_if_missing(seed_seq_path)
    wrapper.run_conv(seed_seq_path, seqs_path, conv_folder, output_path,
                     bash_exec, num_p=num_p)
    return

def make_seed_seq(binding_site_path, seed_seq_path):
    return make_seed_seqs.make(binding_site_path, seed_seq_path)

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

def convert_conv_to_meme_full_num(matrix, composition, meme):
    """
    When using output.1.matrix instead, so we get full numbers instead of probs.
    """
    generic.quit_if_missing(composition)
    generic.quit_if_missing(matrix)
    generic.warn_if_exist(meme)
    conv_to_meme.convert_fullnum(matrix, composition, meme)
    generic.quit_if_missing(meme)

# if __name__ == "__main__":
