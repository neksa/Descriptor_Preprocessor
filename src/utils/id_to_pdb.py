"""
Given a ACC+ID, query uniprot to retrieve the corresponding pdb_id

From https://www.uniprot.org/help/fasta-headers:
UniRef
>UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier
RepID=RepresentativeMember
Where:

UniqueIdentifier is the primary accession number of the UniRef cluster.
ClusterName is the name of the UniRef cluster.
Members is the number of UniRef cluster members.
TaxonName is the scientific name of the lowest common taxon shared by all
UniRef cluster members.
TaxonIdentifier is the NCBI taxonomy identifier of the lowest common taxon
shared by all UniRef cluster members.
RepresentativeMember is the entry name of the representative member of the
UniRef cluster.

##########################

UniProtKB
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName
OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
Where:

db is ‘sp’ for UniProtKB/Swiss-Prot and ‘tr’ for UniProtKB/TrEMBL.
UniqueIdentifier is the primary accession number of the UniProtKB entry.
EntryName is the entry name of the UniProtKB entry.
ProteinName is the recommended name of the UniProtKB entry as annotated in
the RecName field. For UniProtKB/TrEMBL entries without a RecName field,
the SubName field is used. In case of multiple SubNames, the first one is
used. The ‘precursor’ attribute is excluded, ‘Fragment’ is included with the
name if applicable.
OrganismName is the scientific name of the organism of the UniProtKB entry.
OrganismIdentifier is the unique identifier of the source organism, assigned
by the NCBI.
GeneName is the first gene name of the UniProtKB entry. If there is no gene
name, OrderedLocusName or ORFname, the GN field is not listed.
ProteinExistence is the numerical value describing the evidence for the
existence of the protein.
SequenceVersion is the version number of the sequence.
"""

import urllib.parse
import urllib.request

from utils import uniprot_id_converter

def query_genename(genenames):
    id_pdb_map = _query(genenames, "GENENAME")
    return id_pdb_map


def query_acc(acc_ids):
    id_pdb_map = _query(acc_ids, "ACC")
    return id_pdb_map


def _query(input_names, tag):
    id_pdb_map = uniprot_id_converter.convert(tag, "PDB_ID", input_names)
    return id_pdb_map


if __name__ == "__main__":
    genenames = ["MSHA_SALTO", "NADE_STRCO"]
    x = query_genename(genenames)

    print(x.strip().split("\n")[1:])

