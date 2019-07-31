import urllib.parse
import urllib.request


def convert(from_id, to_id, input_ids):
    """
    Assumption that new_id is unique, and that for multiple matches to
    old_id, the first is selected. Works well when new_id is PDB_ID.
    :param from_id: "ACC"
    :param to_id: "PDB_ID"
    :param input_ids:
    """
    assert isinstance(input_ids, list)
    query_line = " ".join(input_ids)
    url = 'https://www.uniprot.org/uploadlists/'

    params = {'from': from_id, 'to': to_id, 'format': 'tab',
              'query': query_line}
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    parsed_response = response.decode('utf-8')
    old_new_id_str = parsed_response.strip().split("\n")
    if len(old_new_id_str) == 1:  # Only header line
        return dict()
    old_new_id_str = old_new_id_str[1:]
    old_new_id_map = dict()
    for line in old_new_id_str:
        # ['MSHA_SALTO\tmshA', 'NADE_STRCO\tnadE']
        old_id, new_id = line.split("\t")
        if old_id not in old_new_id_map:
            old_new_id_map[old_id] = new_id.lower()
    return old_new_id_map