import os
import contextlib
from urllib import request

def download_pdb_files(seq_cid_map,
                        output_folder,
                        replace_existing=True,
                        file_suffix='.pdb',
                        url_template='https://files.rcsb.org/view/{}.pdb'):
    # This downloads the .pdb files listed in pdb_list, from rcsb server.
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    stored_pdb_files = set(os.listdir(output_folder))
    for pname in seq_cid_map.keys():
        pname = pname.lower()
        url = url_template.format(pname.strip())
        output_path = os.path.join(output_folder, pname+file_suffix)
        if replace_existing or pname+file_suffix not in stored_pdb_files:
            with contextlib.closing(request.urlopen(url)) as contents:
                with open(output_path, 'w') as output_file:
                    output_file.write(contents.read().decode("utf-8"))
                    # for line in contents:
                    #     output_file.write(line.decode("utf-8"))