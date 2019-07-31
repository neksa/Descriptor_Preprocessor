from utils import logs

from config import paths
import preprocess


def main():
    logs.set_logging_level()
    extract_path = paths.PROSITE_EXTRACT
    motif_len = 13
    output = paths.PID_PDB_MAP
    num_p = 7
    preprocess.run_prosite_meme(extract_path, motif_len, output, num_p)

    extract_path = paths.PROSITE_EXTRACT
    motif_len = 13
    ref_meme_txt = paths.REF_MEME_TXT
    output = paths.PID_PDB_MAP
    preprocess.run_prosite_mast(extract_path, motif_len, ref_meme_txt, output)

    extract_path = paths.IONCOM_EXTRACT
    motif_len = 13
    ref_meme_txt = paths.REF_MEME_TXT
    output = paths.PID_PDB_MAP
    preprocess.run_ioncom_mast(extract_path, motif_len, ref_meme_txt, output)


if __name__ == "__main__":
    main()
