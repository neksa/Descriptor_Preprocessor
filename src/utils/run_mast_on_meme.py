import logging
import os
import shutil
import subprocess

import config

def create_mast_output(meme_path=None, output_dir=None):
    # mast_binary_path = "/home/yincp/meme/bin/mast"
    if meme_path is None:
        meme_file_name = "meme_mg_50_30_all.txt"
        meme_path = os.path.join(config.ROOT, meme_file_name)
    if output_dir is None:
        output_dir = os.path.join(config.ROOT, 'mast')
    if os.path.isdir(output_dir):
        logging.warning(f"Output dir <{output_dir}> is not empty. Deleting.")
        shutil.rmtree(output_dir)
    # mast creates its own output_dir so there's no need to mkdir().
    seq_path = os.path.join(config.ROOT, 'data', 'input',
                            'fasta_template.fasta')

    command = f'mast -o {output_dir} {meme_path} {seq_path}'
    subprocess.run(command, shell=True)
    assert os.path.isdir(output_dir)

if __name__ == "__main__":
    create_mast_output()
