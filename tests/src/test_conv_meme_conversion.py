import os
import unittest

from preprocessing.converge import conv_interface
from tests.src import paths_test
from utils import generic


class TestMemeToConv(unittest.TestCase):
    """
    Because they use the same REF_MEME_FROM_CONV, it's possible in
    theory for conv_to_meme and meme_to_conv to both be wrong, and yet
    pass.
    """
    def setUp(self):
        self.meme_full = paths_test.ORIG_MEME_FOR_CONVERT
        self.meme_minimal = paths_test.REF_MEME_FROM_CONV
        self.ref_matrix = paths_test.REF_CONV_MATRIX
        self.ref_composition = paths_test.REF_CONV_COMPOSITION

        self.tmp_1 = paths_test.TMP_FILE_TEMPLATE.format(1)
        self.tmp_2 = paths_test.TMP_FILE_TEMPLATE.format(2)

        generic.warn_if_exist(self.tmp_1)
        generic.warn_if_exist(self.tmp_2)

    def test_convert_meme_to_conv_full(self):
        conv_interface.convert_meme_to_conv(self.meme_full,
                                            self.tmp_1,
                                            self.tmp_2,
                                            minimal=False)
        self.assertTrue(os.path.isfile(self.meme_full))
        self.assertTrue(os.path.isfile(self.tmp_1))
        self.assertTrue(os.path.isfile(self.tmp_2))

        with open(self.tmp_1, "r") as file:
            act_lines = file.readlines()
        with open(self.ref_composition, "r") as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)

        with open(self.tmp_2, "r") as file:
            act_lines = file.readlines()
        with open(self.ref_matrix, "r") as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)

        os.remove(self.tmp_1)
        os.remove(self.tmp_2)

    def test_convert_meme_to_conv_minimal(self):
        conv_interface.convert_meme_to_conv(self.meme_minimal,
                                            self.tmp_1,
                                            self.tmp_2,
                                            minimal=True)
        self.assertTrue(os.path.isfile(self.meme_minimal))
        self.assertTrue(os.path.isfile(self.tmp_1))
        self.assertTrue(os.path.isfile(self.tmp_2))

        with open(self.tmp_1, "r") as file:
            act_lines = file.readlines()
        with open(self.ref_composition, "r") as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)
        with open(self.tmp_2, "r") as file:
            act_lines = file.readlines()
        with open(self.ref_matrix, "r") as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)
        os.remove(self.tmp_1)
        os.remove(self.tmp_2)


class TestConvToMeme(unittest.TestCase):
    def setUp(self):
        self.matrix = paths_test.REF_CONV_MATRIX
        self.composition = paths_test.REF_CONV_COMPOSITION
        self.ref_meme = paths_test.REF_MEME_FROM_CONV
        self.tmp_1 = paths_test.TMP_FILE_TEMPLATE.format(1)
        generic.warn_if_exist(self.tmp_1)

    def test_convert_meme_to_conv_full(self):
        conv_interface.convert_conv_to_meme(self.matrix,
                                            self.composition,
                                            self.tmp_1)
        self.assertTrue(os.path.isfile(self.ref_meme))
        self.assertTrue(os.path.isfile(self.tmp_1))
        with open(self.tmp_1, "r") as file:
            act_lines = file.readlines()
        with open(self.ref_meme, "r") as file:
            ref_lines = file.readlines()
        self.assertListEqual(act_lines, ref_lines)
        os.remove(self.tmp_1)