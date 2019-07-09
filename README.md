# Descriptor of Elementary Function





Tests
- Generate Reference Output
    - `/tests/src/setup_ref.py`

- Visualise Reference Output
    - `/tests/src/plot_ref.py`

- Checks against reference output
    - `/tests/src/test_motif_finder.py`
    - `/tests/src/test_descr_main.py`


Data files
- `/data`
    - `/tmp`: created during runtime, should be deleted at end of run, 
            except for debugging. Does not get deleted for tests that fail.
    - `/input`
        - `/ioncom`
           - `allsulfate.txt`: Raw sequence-binding_site match, for mg, in 
             dataIonCom.zip, downloaded from 
             https://zhanglab.ccmb.med.umich.edu/IonCom/ >> 
             download dataset used to...
           - `ioncom.txt`: allid_reso3.0_len50_nr40.txt in dataIonCom, shows 
             list of sequences. (deprecated eventually)
        - `/mg_full`
            - `mg_50.fasta`: From uniprot, uniref50 for seqs with MG as 
            co-factor/ligand.
            - `mg_100.fasta`: uniref100 for MG cofactor seqs
        - `/pdb_files`: Stored pdb_files. Both tests and main should call this,
         since downloading takes a while. Automatically downloaded from rscb 
         server, via link https://files.rcsb.org/view/{1ABC}.pdb
         - `/prosite`
            - `prosite_extract.txt`: Copy-pasted from html (inspect source 
            code) from prosite website (https://prosite.expasy.org/cgi-bin/pdb/pdb_structure_list.cgi?src=PS00018).
        - `/internal`
            - `fasta_template.fasta`: Used for running mast, when we only want 
            the seqlogo and doesn't actually care about matching for motifs.
            - `meme.txt`: Motif file for Calcium EF-hand. 

Linter
- `pylint`, mostly following google style guide with some additional disabled
 clauses. 
