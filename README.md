# Descriptor of Elementary Function, Preprocessor

Parse various raw sources to produce input data for descriptor calculator. 

Takes in as input either from prosite or ioncom. To start, we need the pdb code
 and chain id for identified sequences. 

For prosite, this can be obtained using the html extract and a manually 
copy-pasted pdb list. The latter provides the pdb code, while the former 
allows for matching between pdb code and chain id. 

For ioncom, the pdb code-chain id matching is derived directly from the 
downloaded test dataset. 

Requires compiling of converge, and meme-suite, both in `src`. Call `make 
converge` in converge folder, and `./configure --prefix=$PWD 
--with-url=http://meme-suite.org/ --enable-build-libxml2 
--enable-build-libxslt` `make` `make install` in meme-suite. Binaries for 
both cannot be used because absolute path-ing done in compile stage (not sure
 for converge). meme-suite requires `mpicc` for parallel runs. 

Workflow

Unless otherwise stated, functions come from `src/preprocess.py`.

1. Extract seq-name and chain-ID from source extracts<br>
    - Input: <br>
        1. Website extract files (e.g. `data/input/prosite_extract.txt`)<br>
    - Output: <br>
        1. `data/internal/pname_cid_map.pkl`<br>
    - Description: <br>
        1. Place prosite or ioncom extract in `/data/input/`. <br>
        2. Run `parse_extract_prosite()` or `parse_extract_ioncom()`.

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
        1. `data/input/ioncom_binding_sites.txt`<br>
    - Output: <br>
        1. `data/internal/seed_seqs.fasta`<br>
    - Description: <br>
        The motif-finding binary `converge` requires seed sequences from 
        which it generates its initial set of motifs.
        1. Place ioncom binding-site file in `/data/input/`. <br>
        2. Run `make()` in `src/make_conv_seed_seqs.py`.<br>
        
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
        generates both the motif file and the motif positions. Run 
        `find_motifs_meme()`
        2. Motifs can be found using a given motif file. First, put the motif 
        file (in MEME format) in `data/input/<filename>`. Then, run 
        `find_motifs_mast()`
        3. Motifs can be derived from scratch using `converge`, which also 
        provides the motif file and positions. Run `make
          (input_fname=<filename>, num_p=<num_processors>)` in 
          `src/make_conv_seed_seqs.py`.<br>
          
        The motif-finding process takes a while. Instructions to run it 
        separately are in [1] below. 

Tests
- Generate Reference Output
    - `/tests/src/setup_ref.py`

- Visualise Reference Output
    - `/tests/src/plot_ref.py`

- Checks against reference output
    - `/tests/src/test_preprocessing.py`

Data files
- `/data`
    - `/debug`: created during runtime, should be deleted at end of run, 
            except for debugging.
    - `/input`
        - `ioncom_extract.txt`: Raw sequence-binding_site match, for mg, in 
        dataIonCom.zip, downloaded from 
        https://zhanglab.ccmb.med.umich.edu/IonCom/ >> 
        download dataset used to...
        - `ioncom_binding_sites.txt`: allid_reso3.0_len50_nr40.txt in 
        dataIonCom, shows list of sequences.
        - `mg_50.fasta`: From uniprot, uniref50 for seqs with MG as 
        co-factor/ligand.
        - `mg_100.fasta`: uniref100 for MG cofactor seqs
        - `pdb_files`: Stored pdb_files. Both tests and main should call this,
        since downloading takes a while. Automatically downloaded from rscb 
        server, via link https://files.rcsb.org/view/{1ABC}.pdb
        - `prosite_extract.txt`: Copy-pasted from html (inspect source 
        code) from prosite website (https://prosite.expasy
        .org/cgi-bin/pdb/pdb_structure_list.cgi?src=PS00018).
        - `ref_meme.txt`: Motif file for Calcium EF-hand. 

Linter
- `pylint`, mostly following google style guide with some additional disabled
 clauses. 
 
To build sequence logo for final UI, take the descr.pkl from
 Descriptor_Calculator, and see build() in seq_logo.py. After that is done, run 
 `./src/meme_suite/meme/src/ceqlogo -i ./gxggxg_descr.txt -m motif_name 
 -o gxggxg_logo.eps` to generate the logo figure in .eps format. Finally, go to 
 https://www.epsconverter.com/ to convert it into .png form with white
  background. Move to Descriptor_Calculator/src/ui/static with the
   appropriate filename (see IndividualFigure for that) to load it in. 
 
 Todo:
 1. Add compilation instructions for `converge` and `meme-suite`. 