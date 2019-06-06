import os

src = os.path.dirname(__file__)
root = src.rsplit(os.sep, 1)[0]

##############################################################################
# Common File Directories
# input_dir = os.path.join(root, 'data', 'input')
pdb_files_dir = os.path.join(root, 'pdb_files')
input_dir = os.path.join(root, 'data', 'input')
store_dir = os.path.join(root, 'data', 'store')
output_dir = os.path.join(root, 'data', 'output')

##############################################################################

_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

AA3_to_AA1 = dict(_aa_index)

