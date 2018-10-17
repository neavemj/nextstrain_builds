from Bio import Entrez
from Bio import SeqIO

# download BTV genome segments from NCBI with required metadata
# Matthew J. Neave 17.10.2018

# get NCBI sequence IDs from spreadsheet from John White
# need to enter an email so users don't overload their servers

Entrez.email = "matthewjneave1@gmail.com"


def retrieve_ncbi_record(ncbi_id):
    new_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    # extracting several bits of information from different parts of the record
    print(seq_record.features[0])
    collection_date = seq_record.features[0].qualifiers["collection_date"][0]
    host =seq_record.features[0].qualifiers["host"][0]
    authors = seq_record.annotations['references'][0].authors.split(",")[0] + " et al"
    title = seq_record.annotations['references'][0].title
    journal = seq_record.annotations['references'][0].journal
    metadata = [seq_record.id, collection_date, host, authors, title, journal]
    return(metadata, seq_record)


# use this function to retrieve all the records and write to fasta file
# with corresponding metadata file correctly formatted for Nextstrain

#mahar_fasta = open("mahar_RHDV.fasta", "w")

print(retrieve_ncbi_record("JN881985"))

#for record_id in mahar_ids:
#    record = retrieve_ncbi_record(record_id)
    #mahar_fasta.write(record[0] + "\n" + str(record[1].seq) + "\n")



