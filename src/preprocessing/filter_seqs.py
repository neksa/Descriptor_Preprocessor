def delete_short_seqs(fasta_fname, threshold=30):
    """
    Delete sequences with len < threshold.
    """
    output_str = []
    tmp = ""
    with open(fasta_fname, 'r') as file:
        for line in file:
            if line.startswith(">"):
                tmp = line
                continue
            if len(line) < threshold:
                continue
            else:
                output_str.append(tmp + line)
                continue
    output_str = ''.join(output_str)
    with open(fasta_fname, 'w') as file:
        file.write(output_str)