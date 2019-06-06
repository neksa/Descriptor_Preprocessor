self.weights_sig = self._weight_based_on_signature(df)

def _weight_based_on_signature(self, df):
    param = list()
    df = df.groupby('relative_sno')
    total = 0
    for relative_sno, df_values in df:
        raw_res = df_values.res
        res = [AA3_to_AA1[x] for x in raw_res]
        names, counts = np.unique(res, return_counts=True)
        top_3 = sum(sorted(counts)[:3])
        param.append([relative_sno, top_3])
        total += top_3
    weights = defaultdict(list)
    for relative_sno, top_3 in param:
        weights['relative_sno'].append(relative_sno)
        weights['weight'].append(top_3 / total)
    to_return = pd.DataFrame.from_dict(weights)
    to_return = to_return.set_index('relative_sno')
    assert np.isclose(sum(to_return.values), 1)
    return to_return


class _PhipsiMatcher:
    def __init__(self):
        # only accept from a single relative_sno, only values.
        self.to_skip = False

    def load(self, phipsis):
        if np.allclose(phipsis, np.full(phipsis.shape, 360)):
            self.to_skip = True
            self.probs = np.zeros(len(phipsis))
            self.weight = 0.
            return
        phipsis = (phipsis + 180) / 2
        assert np.amin(phipsis) > 0
        phipsis = np.sin(phipsis)
        assert isinstance(phipsis, np.ndarray)
        distances = np.linalg.norm(phipsis, axis=1)
        self.stdev = np.std(distances)
        max_stdev = 147 # temp magic normalisation factor
        self.weight = (max_stdev - self.stdev) / max_stdev
        self.distance = distances
        self.phipsi = phipsis
        return

    def query(self, q_phipsi):
        if self.to_skip:
            return self.weight, self.probs
        if np.allclose(q_phipsi, np.full(q_phipsi.shape, 360)):
            return 0., np.zeros(len(self.phipsi))
        # 1D gaussian cdf
        q_phipsi = (q_phipsi + 180) / 2
        assert np.amin(q_phipsi) > 0
        q_phipsi = np.sin(q_phipsi)
        assert isinstance(q_phipsi, np.ndarray)
        dists = np.linalg.norm(self.phipsi - q_phipsi, axis=1)
        erf_arg = dists / (self.stdev * math.sqrt(2))
        probs = np.array([1 - math.erf(i) for i in erf_arg])
        return self.weight, probs



# 1D gaussian cdf
# q_phipsi = (q_phipsi + 180) / 2
# assert np.amin(q_phipsi) > 0
# q_phipsi = np.sin(q_phipsi)
# assert isinstance(q_phipsi, np.ndarray)
# dists = np.linalg.norm(self.phipsi - q_phipsi, axis=1)
# erf_arg = dists / (self.stdev * math.sqrt(2))
# probs = np.array([1 - math.erf(i) for i in erf_arg])