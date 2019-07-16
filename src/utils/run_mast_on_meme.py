import logging
import os
import shutil
import subprocess

from config import paths

def create(meme_path=None, output_dir=None,
                       seq_file=paths.TEMPLATE_SEQFILE):
    # mast_binary_path = "/home/yincp/meme/bin/mast"
    if meme_path is None:
        meme_path = paths.REF_MEME_TXT
    if output_dir is None:
        output_dir = os.path.join(paths.ROOT, 'mast')
    if os.path.isdir(output_dir):
        logging.warning(f"Output dir <{output_dir}> is not empty. Deleting.")
        shutil.rmtree(output_dir)
    # mast creates its own output_dir so there's no need to mkdir().
    command = f'mast -o {output_dir} {meme_path} {seq_file}'
    subprocess.run(command, shell=True)
    assert os.path.isdir(output_dir)

if __name__ == "__main__":
    # import os
    # create_mast_output(os.path.join(paths.ROOT, "output_meme.txt"),
    #                    os.path.join(paths.ROOT, "mast_folder"))
    create(os.path.join(paths.ROOT, "conv_meme.txt"),
           os.path.join(paths.ROOT, "output_dir"))