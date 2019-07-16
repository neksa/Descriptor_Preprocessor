# import urllib.parse
# import urllib.request
#
# url = 'https://www.uniprot.org/uploadlists/'
#
# params = {
# 'from': 'ACC','to': 'PDB_ID',
# 'format': 'tab',
# 'query': 'RPOZ_THET8 P40926 O43175 Q9UM73 P97793'}
#
# data = urllib.parse.urlencode(params)
# data = data.encode('utf-8')
# req = urllib.request.Request(url, data)
# with urllib.request.urlopen(req) as f:
#     response = f.read()
# print(response.decode('utf-8'))

# from collections import defaultdict
# #
# # d1 = defaultdict(lambda: defaultdict(list))
# # d1['rgg']['dsf'].append(3)
# # print(d1)
# #
# # d2 = {key:dict(value) for key, value in d1.items()}
# # print(d2)
# from preprocessing.converge import conv_to_meme
# from config import paths
# import os
#
# conv_to_meme.convert(os.path.join(paths.ROOT, "output.4.matrix.0"),
#                      os.path.join(paths.ROOT, "composition.txt"),
#                      os.path.join(paths.ROOT, "conv_meme.txt"))