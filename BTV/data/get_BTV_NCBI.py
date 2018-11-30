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
    collection_date = seq_record.features[0].qualifiers["collection_date"][0]
    host = seq_record.features[0].qualifiers["host"][0]
    try:
        lat_long = seq_record.features[0].qualifiers["lat_lon"][0]
    except:
        lat_long = "unknown"
    authors = ",".join(seq_record.annotations['references'][0].authors.split(",")[0:6]) + " et al"
    title = seq_record.annotations['references'][0].title
    journal = seq_record.annotations['references'][0].journal
    metadata = [seq_record.id, collection_date, host, lat_long, authors, title, journal]
    return(metadata, seq_record)


# use this function to retrieve all the records and write to fasta file
# with corresponding metadata file correctly formatted for Nextstrain
# think I'm going to need separate files for each of the segments

btv_segment1_fasta = open("btv_segment1.fasta", "w")
btv_segment1_meta = open("btv_segment1.meta.tsv", "w")
btv_segment2_fasta = open("btv_segment2.fasta", "w")
btv_segment2_meta = open("btv_segment2.meta.tsv", "w")
btv_segment3_fasta = open("btv_segment3.fasta", "w")
btv_segment3_meta = open("btv_segment3.meta.tsv", "w")
btv_segment4_fasta = open("btv_segment4.fasta", "w")
btv_segment4_meta = open("btv_segment4.meta.tsv", "w")
btv_segment5_fasta = open("btv_segment5.fasta", "w")
btv_segment5_meta = open("btv_segment5.meta.tsv", "w")
btv_segment6_fasta = open("btv_segment6.fasta", "w")
btv_segment6_meta = open("btv_segment6.meta.tsv", "w")
btv_segment7_fasta = open("btv_segment7.fasta", "w")
btv_segment7_meta = open("btv_segment7.meta.tsv", "w")
btv_segment8_fasta = open("btv_segment8.fasta", "w")
btv_segment8_meta = open("btv_segment8.meta.tsv", "w")
btv_segment9_fasta = open("btv_segment9.fasta", "w")
btv_segment9_meta = open("btv_segment9.meta.tsv", "w")
btv_segment10_fasta = open("btv_segment10.fasta", "w")
btv_segment10_meta = open("btv_segment10.meta.tsv", "w")

fasta_file_list = [btv_segment1_fasta, btv_segment2_fasta, btv_segment3_fasta, btv_segment4_fasta,
                   btv_segment5_fasta, btv_segment6_fasta, btv_segment7_fasta, btv_segment8_fasta,
                   btv_segment9_fasta, btv_segment10_fasta]

meta_file_list = [btv_segment1_meta, btv_segment2_meta, btv_segment3_meta, btv_segment4_meta, btv_segment5_meta,
                  btv_segment6_meta, btv_segment7_meta, btv_segment8_meta, btv_segment9_meta, btv_segment10_meta]

# want to write a header line into each meta file

for meta_file in meta_file_list:
    meta_file.write("\t".join(["strain", "virus", "accession", "segment", "type", "state", "date", "host", "lat_long", "authors", "title", "journal"]) + "\n")


# need to go through each isolate, then each segment
with open("AUS_BTV_ISOLATES.txt") as f:
    header = next(f) # skip header
    count = 0
    for lines in f:
        lines = lines.strip()
        cols = lines.split("\t")
        strain = cols[0].strip()
        count += 1
        print("analysing virus: {}, number: {}".format(strain, str(count)))
        virus = cols[1].strip()
        btv_type = "type_" + cols[2].strip()
        state = cols[3].strip()
        date = cols[4].strip() + "-XX-XX"
        # get record for each segment
        # should be 10 segments
        for index, segment_id in enumerate(cols[5:]):
            record = retrieve_ncbi_record(segment_id)
            segment = index + 1
            print("segment", str(segment))
            fasta_file_list[index].write(">" + strain + "\n" + str(record[1].seq) + "\n")
            meta_file_list[index].write("\t".join([strain, virus, segment_id, str(segment), btv_type, state, date,
                                                   record[0][2], record[0][3], record[0][4],
                                                   record[0][5], record[0][6]
                                                   ]) + "\n")



