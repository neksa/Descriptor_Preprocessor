"""
To use: Logo(to_logo, min_sno); plt.show()
to_logo of form [[(res, proportion),... for 0th position] for nth position]
min_sno = first relative position of sno, for xaxis labels.
Call _to_align to align each individual text character.

Adapted from https://stackoverflow.com/questions/42615527/sequence-logos-in-
matplotlib-aligning-xticks, with color properties taken from weblogo.
"""

from enum import Enum, auto
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from matplotlib.text import TextPath
import matplotlib.pyplot as plt

# Class definition kept in global so pickle works fine.
class Res(Enum):
    polar=auto()
    neutral=auto()
    basic=auto()
    acidic=auto()
    hydrophobic=auto()

class Logo:
    def __init__(self, data, init_pos, figsize=(20, 3), convert_AA3=True,
                 title=None):
        self._title = title
        self._init_pos = init_pos
        self._convert_AA3 = convert_AA3
        self._AA3_to_AA1 = self.get_AA3_to_AA1()
        self._colorcode = self.get_colorcode()
        self._res_grp_map = self.get_res_grp_map()
        self._letter_to_glyph = self.get_letter_to_glyph()
        self._figsize = figsize
        self.plot = self.get_plot(data)

    def get_AA3_to_AA1(self):
        AA3_to_AA1 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
                          HIS='H', HSE='H', HSD='H', ILE='I', LYS='K', LEU='L',
                          MET='M', MSE='M', ASN='N', PRO='P', GLN='Q', ARG='R',
                          SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
        return AA3_to_AA1

    def get_colorcode(self):
        colorcode = dict()
        colorcode[Res.acidic] = '#ff0000'
        colorcode[Res.basic] = '#0000ff'
        colorcode[Res.hydrophobic] = '#000000'
        colorcode[Res.neutral] = '#800080'
        colorcode[Res.polar] = '#008000'
        return colorcode

    def get_res_grp_map(self):
        res_grp_map = dict(
            C=Res.polar, G=Res.polar, S=Res.polar, T=Res.polar, Y=Res.polar,
            N=Res.neutral, Q=Res.neutral,
            H=Res.basic, K=Res.basic, R=Res.basic,
            D=Res.acidic, E=Res.acidic,
            A=Res.hydrophobic, F=Res.hydrophobic, I=Res.hydrophobic,
            L=Res.hydrophobic, M=Res.hydrophobic, P=Res.hydrophobic,
            V=Res.hydrophobic, W=Res.hydrophobic)
        return res_grp_map

    def get_letter_to_glyph(self):
        # preferably Arial, but that requires separate installation
        fp = FontProperties(family="sans", weight="bold")
        letters = dict(
            A=TextPath((-0.35, 0), "A", size=1, prop=fp),
            C=TextPath((-0.35, 0), "C", size=1, prop=fp),
            D=TextPath((-0.36, 0), "D", size=1, prop=fp),
            E=TextPath((-0.34, 0), "E", size=1, prop=fp),
            F=TextPath((-0.32, 0), "F", size=1, prop=fp),
            G=TextPath((-0.378, 0), "G", size=1, prop=fp),
            H=TextPath((-0.36, 0), "H", size=1, prop=fp),
            I=TextPath((-0.14, 0), "I", size=1, prop=fp),
            K=TextPath((-0.39, 0), "K", size=1, prop=fp),
            L=TextPath((-0.33, 0), "L", size=1, prop=fp),
            M=TextPath((-0.41, 0), "M", size=1, prop=fp),
            N=TextPath((-0.36, 0), "N", size=1, prop=fp),
            P=TextPath((-0.345, 0), "P", size=1, prop=fp),
            Q=TextPath((-0.392, 0), "Q", size=1, prop=fp),
            R=TextPath((-0.385, 0), "R", size=1, prop=fp),
            S=TextPath((-0.32, 0), "S", size=1, prop=fp),
            T=TextPath((-0.31, 0), "T", size=1, prop=fp),
            V=TextPath((-0.33, 0), "V", size=1, prop=fp),
            W=TextPath((-0.47, 0), "W", size=1, prop=fp),
            Y=TextPath((-0.33, 0), "Y", size=1, prop=fp))
        return letters

    def get_plot(self, data):
        fig, ax = plt.subplots(figsize=self._figsize)
        x = self._init_pos
        for per_position in data:
            y = 0
            for AA3, proportion in per_position:
            # for AA3, proportion in per_position[0], per_position[1]:
                ax = self.add_patch(AA3, x, y, proportion, ax)
                y += proportion
            x += 1
        plt.xticks(range(self._init_pos, x))
        plt.xlim((self._init_pos-0.75, x-0.25))
        plt.ylim((0, 1))
        for name, value in self._colorcode.items():
            plt.scatter(0, 2, label=name.name, s=100, color=value)
        plt.legend()
        if self._title:
            plt.title(self._title)
        return ax

    def add_patch(self, AA3, x, y, yscale, ax):
        if self._convert_AA3:
            AA1 = self._AA3_to_AA1[AA3]
        else:
            AA1 = AA3
        glyph = self._letter_to_glyph[AA1]
        res_group = self._res_grp_map[AA1]
        globscale = 1.35
        t = Affine2D().scale(globscale, yscale*globscale) + \
            Affine2D().translate(x, y) + ax.transData
        p = PathPatch(glyph, lw=0, color=self._colorcode[res_group],
                      transform=t)
        ax.add_patch(p)
        return ax


def _to_align():
    AA3_to_AA1 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
                      HIS='H', HSE='H', HSD='H', ILE='I', LYS='K', LEU='L',
                      MET='M', MSE='M', ASN='N', PRO='P', GLN='Q', ARG='R',
                      SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
    data = [[]]
    length = len(AA3_to_AA1)
    per_count = 1/length
    for i in AA3_to_AA1.keys():
        data[0].append((i, per_count))
    plot = Logo(data, -1, (5, 20)).plot
    plot.legend().remove()
    plot.axvline(-1)
    plot.axvline(-1.478)
    plot.axvline(-0.499)
    plt.show()