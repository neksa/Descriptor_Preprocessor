# Descriptor of Elementary Function


Workflow

Preprocessing

Unless otherwise stated, functions come from `src/preprocessing/preprocess.py`.

1. Extract seq-name and chain-ID from source extracts<br>
    - Input: <br>
        1. Website extract files (e.g. `data/user/input/prosite_extract.txt`)<br>
    - Output: <br>
        1. `data/internal/pname_cid_map.pkl`<br>
    - Description: <br>
        1. Place prosite or ioncom extract in `/data/user/input/`. <br>
        2. Run `parse_extracts(source, filename)` in `preprocessing/preprocess.py`, 
        specifying the source (`ioncom` or `prosite`) and name of extract file. This 
        extracts the sequence names and chain-ids to be processed. 

2. (Optional) Download relevant `.pdb` files from rscb server<br>
    - Input: <br>
        1. `data/internal/pname_cid_map.pkl`<br>
        2. Internet connection<br>
    - Output: <br>
        1. Populated `data/internal/pdb_files/`<br>
    - Description: <br>
        Downloads corresponding `.pdb` files from rscb server. Delete entries
         in `pname_cid_map` if `.pdb` files are not in folder.
        1. Run `download_pdb()`
        2. Run `trim_pnames_based_on_pdb()`
        
3. Create sequence `.fasta` file
    - Input: <br>
        1. `data/internal/pname_cid_map.pkl`<br>
    - Output: <br>
        1. `data/internal/seqs.fasta`<br>
    - Description: <br>
         The motif-finding binaries require the sequences to be in a `.fasta`
          file. 
        1. Run `create_seq()`

4. Filter short sequences
    - Input: <br>
        1. `data/internal/seqs.fasta`<br>
        2. Populated `data/internal/pdb_files/`<br>
    - Output: <br>
        1. Updated `data/internal/seqs.fasta`<br>
    - Description: <br>
        Sequences shorter than the desired motif length (30 residues) can lead
        to errors when performing the motif search, and need to be dropped.
        1. Run `filter_seq_file()`
         
5. (Optional) Create seed sequence file for `converge`
    - Input: <br>
        1. `data/user/input/ioncom_binding_sites.txt`<br>
    - Output: <br>
        1. `data/internal/seed_seqs.fasta`<br>
    - Description: <br>
        The motif-finding binary `converge` requires seed sequences from 
        which it generates its initial set of motifs.
        1. Place ioncom binding-site file in `/data/user/input/`. <br>
        2. Run `make()` in `src/preprocessing/make_conv_seed_seqs.py`.<br>
        
5. Run motif-search binary to find motif positions<br>
    - Input: <br>
        1. `data/internal/seqs.fasta`<br>
        2. (Optional) Populated `data/internal/pdb_files/`<br>
        3. (Optional) Provided motif file (e.g. `data/user/input/meme.txt`)
        4. (Optional) `data/internal/seed_seqs.fasta`<br>
    - Output: <br>
        1. `data/internal/motif_pos.pkl`<br>
    - Description: <br>
        This finds the positions of the desired motif for each sequence-chain. 
        There are three implemented ways of running this locally:
        1. Motifs can be derived from scratch, using `meme`. This 
        generates both the motif file and the motif positions. Run `find
        (process='meme', num_p=<num_processors>)` in 
        `src/preprocessing/motif_finder.py`.
        2. Motifs can be found using a given motif file. First, put the motif 
        file (in MEME format) in `data/user/input/<filename>`. Then, run `find
        (process='mast', motif_fname=<filename>, num_p=<num_processors>)` in 
        `src/preprocessing/motif_finder`. 
        3. Motifs can be derived from scratch using `converge`, which also 
        provides the motif file and positions. Run `make
          (input_fname=<filename>, num_p=<num_processors>)` in 
          `src/preprocessing/make_conv_seed_seqs.py`.<br>
          
        Because of long run-time for the motif-finding process, it is recommended
        to run this step in a server. Instructions for doing so are in [1] below. 

Descriptor Generation

6. Calculate descriptor properties<br>
    - Input: <br>
        1. `data/internal/motif_pos.pkl`<br>
        2. Populated `data/internal/pdb_files/`<br>
    - Output: <br>
        1. `data/internal/descrs.pkl`<br>
    - Description: <br>
        This calculates the descriptor properties, for each motif. Run 
        `calculate()` in `src/descr/descr_main.py`.
    
7. Visualise properties<br>
    - Input: <br>
        1. `data/internal/descrs.pkl`<br>
    - Output: <br>
        1. (Optional) `data/user/output/`<br>
    - Description: <br>
        Plots for different descriptor properties can be generated via 
        `src/utils/plots.py`. Run each `plot_<something>(save=False)` as needed, 
        and set `save=True` to keep the generated plots in the output folder.
    

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


TODO:
1. pdb_list from prosite need to be extracted too...?