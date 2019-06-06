import numpy as np
import os
import pickle as pkl

from global_config import store_dir
from matchers.matcher import Matcher
import sys

with open(os.path.join(store_dir, "current.pkl"),
          "rb") as pklfile:
    df = pkl.load(pklfile)
df = df.sort_values(['filename', 'seq_marker', 'cid'])
matcher = Matcher()
matcher.load(df)
df_ = df.groupby(['filename', 'seq_marker', 'cid'])
count = 0
for key, df_per_file in df_:
    p_all_sno = matcher.query(df_per_file)
    for sno, p in p_all_sno.items():
        print(sno)
        print(np.array(p))









