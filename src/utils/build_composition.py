from collections import Counter
from utils import generic
from config import paths
import os
import shutil
import subprocess

import decimal


def build(seq_path, output):
    generic.quit_if_missing(seq_path)
    generic.warn_if_exist(output)
    if os.path.isfile(output):
        try:
            shutil.move(output, paths.TRASH)
        except shutil.Error:
            os.remove(output)
    command = f"{paths.FASTA_2_MARKOV_EXEC} -protein < {seq_path} > {output}"
    subprocess.run(command, shell=True)
    generic.quit_if_missing(output)


    # with open(seq_path, 'r') as file:
    #     lines = file.readlines()
    # counter_obj = Counter()
    # for line in lines:
    #     if line.startswith(">"):
    #         pass
    #     counter_obj.update(line.strip())
    # total_count = sum(counter_obj.values())
    # alphabets = sorted(generic.AA3_to_AA1.values())
    # print(counter_obj)
    # with open(output, 'w') as file:
    #     for letter in alphabets:
    #         percentage = "%.5f" % (counter_obj[letter] / total_count)
    #         file.write(percentage + " ")

# from tests.src import paths_test
# build(paths_test.UNIPROT_SEQ, "./composition.txt")