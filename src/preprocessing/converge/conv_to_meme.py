from collections import OrderedDict
import re

def convert(input_conv_path, composition_path, output_path):
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
                nsite_match = re.search(r"K = ([0-9]+)", line)
                assert nsite_match is not None
                nsite = int(nsite_match[1])
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
        if matrix:
            motif_name = f"MEME-{matrix_count}"
            matrices[motif_name] = (nsite, matrix)
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

# if __name__ == "__main__":
#     import os
#     from config import paths
#
#     input_conv_path = os.path.join(paths.ROOT, "output.4.matrix.0")
#     composition_path = os.path.join(paths.ROOT, "composition.txt")
#     output_path = os.path.join(paths.ROOT, "output_meme.txt")
#     convert(input_conv_path, composition_path, output_path)