import subprocess

import numpy as np

from utils import generic
from config import paths

def build(aligned, output, composition=None):
    generic.quit_if_missing(aligned)
    generic.warn_if_exist(output)
    matrix_file = paths.TMP_FILE_TEMPLATE.format("matrix")
    matrix_counter = _matrix_builder(aligned)
    _write_matrix_file(matrix_counter, matrix_file)
    generic.quit_if_missing(matrix_file)
    if composition is not None:
        generic.quit_if_missing(composition)
        command = f"{paths.MATRIX_2_MEME_EXEC} -protein -bg" \
                  f" {composition} < {matrix_file} > {output}"
    else:
        command = f"{paths.MATRIX_2_MEME_EXEC} -protein < {matrix_file} > " \
                  f"{output}"
    subprocess.run(command, shell=True)
    generic.quit_if_missing(output)


def _matrix_builder(aligned_path):
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
    return matrix_counter


def _write_matrix_file(matrix_ordered, output):
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