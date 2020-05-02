from config import paths

import os
import numpy as np
import random

random.seed(10)

# with open(os.path.join(paths.ROOT, 'test_input.txt'), 'w') as file:
#     g_y = 5
#     g_x = (7, 9, 12)
#     total = 200
#     per_pos = total // 20
#     var_percent = 0.2
#     for i in range(30):
#         output_str = ""
#         output_array = []
#         num_to_dist = int(total * var_percent)
#         if i in g_x:
#             for j in range(20):
#                 if j == g_y:
#                     output_array.append(int(total - num_to_dist))
#                 else:
#                     output_array.append(0)
#
#             for __ in range(num_to_dist):
#                 output_array[random.randint(0, 19)] += 1
#         else:
#             for j in range(20):
#                 output_array.append(int(total - num_to_dist) // 20)
#
#             for __ in range(num_to_dist):
#                 output_array[random.randint(0, 19)] += 1
#
#         output_str = ""
#         for j in range(20):
#             output_str += str(int(output_array[j])) + ","
#         output_str = output_str[:-1] + "\n"
#         file.write(output_str)

# def get_aligned_seq(filename):
#     AA3_to_AA1 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
#                       HIS='H', HSE='H', HSD='H', ILE='I', LYS='K', LEU='L',
#                       MET='M', MSE='M', ASN='N', PRO='P', GLN='Q', ARG='R',
#                       SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
#     counts =

matrix = []
with open(os.path.join(paths.ROOT, "converged_matrix2.txt"), 'r') as file:
    row = []
    for line in file:
        split_num = [int(i) for i in line.strip().split(" ")]
        summed = sum(split_num)
        for i in range(len(split_num)):
            split_num[i] = split_num[i] / summed
        matrix.append(split_num)
print(matrix)
with open(os.path.join(paths.ROOT, "output_tmp.txt"),
          'w') as file:
    for row in matrix:
        line = ""
        for prob in row:
            prob_in_str = "{:.6f}".format(prob)
            line += f" {prob_in_str}"
        line += "\n"
        file.write(line)