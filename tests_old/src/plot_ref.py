import os
import pickle

import matplotlib.pyplot as plt

from config import paths
import plots

def plot_prosite_meme():
    ref_ = os.path.join(paths.ROOT, 'tests', 'data', 'ref')
    prosite_meme_descr_path = os.path.join(ref_, 'descr_prosite_meme.pkl')
    with open(prosite_meme_descr_path, 'rb') as file:
        prosite_meme_descr = pickle.load(file)
    plots.plot_signature_logo(prosite_meme_descr, title='prosite_meme_descr')
    plots.plot_dihedral(prosite_meme_descr)

def plot_prosite_mast():
    ref_ = os.path.join(paths.ROOT, 'tests', 'data', 'ref')
    prosite_mast_descr_path = os.path.join(ref_, 'descr_prosite_mast.pkl')
    with open(prosite_mast_descr_path, 'rb') as file:
        prosite_mast_descr = pickle.load(file)
    plots.plot_signature_logo(prosite_mast_descr, title='prosite_mast_descr')

def plot_ioncom_mast():
    ref_ = os.path.join(paths.ROOT, 'tests', 'data', 'ref')
    ioncom_mast_descr_path = os.path.join(ref_, 'descr_ioncom_mast.pkl')
    with open(ioncom_mast_descr_path, 'rb') as file:
        ioncom_mast_descr = pickle.load(file)
    plots.plot_signature_logo(ioncom_mast_descr, title='ioncom_mast_descr')

if __name__ == "__main__":
    plot_prosite_meme()
    # plot_prosite_mast()
    # plot_ioncom_mast()
    plt.show()