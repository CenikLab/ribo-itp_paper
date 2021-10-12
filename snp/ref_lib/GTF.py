
import gzip
from collections import OrderedDict, defaultdict

# Pleasse note that there are differences between GTF an d BED file conventions.
# We had to adjust for these correspondences by adding or subtracting one
# to coordinates when necessary. 
# In particualr:
# NOTE: GTF is 1-based
# start and end coordinates are INCLUSIVE
# BED is 0 based start is inclusive and end is EXCLUSIVE

# IMPORTANT NOTE:
# IN GENCODE GTF FILE
# For - strand transcripts
# The exons are given in DECREASING position
# This is because the "first" exon is the rightmost exon (as opposed to left-pmost in the + case)
# Exons are numbered accordingly and this needs to be taken into account
# For finding CDS, UTRs etc.

class GTFEntry:

    def __init__(self , gtf_line_contents ):
        
        self.contents = gtf_line_contents
        
        #name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
        self.seqname = gtf_line_contents[0]
        
        # name of the program that generated this feature, or the data source (database or project name)
        self.source = gtf_line_contents[1]
        
        # feature type name, e.g. Gene, Variation, Similarity
        self.feature = gtf_line_contents[2]
        
        # Start position of the feature, with sequence numbering starting at 1.
        self.start = int(gtf_line_contents[3])
        
        # End position of the feature, with sequence numbering starting at 1.
        self.end = int(gtf_line_contents[4])
        
        # A floating point value.
        self.score = gtf_line_contents[5] 
        
        # defined as + (forward) or - (reverse).
        self.strand = gtf_line_contents[6]
        
        # frame
        # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
        self.frame = gtf_line_contents[7]
        
        # - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
        self.attribute_str = gtf_line_contents[8]
        
        # Holds parsed attributes
        self.attribute_dict = OrderedDict()
        
        self.__extract_attributes()
        

    def __extract_attributes(self):
        att_contents = self.attribute_str.split(";")
        for p in att_contents:
            
            pairs = p.strip().split()
            if len(pairs) < 2:
                continue
            self.attribute_dict[pairs[0]] = pairs[1].strip("\"")

    def __str__(self ):
        return self.feature + " && "+ self.attribute_dict["gene_id"] + \
                 " && " + self.attribute_dict["transcript_id"]


class GTFfile:
    '''
    TO BE COMPLETED
    '''

    def __init__(self , file):
        myopen = open
        if file.endswith(".gz"):
            myopen = gzip.open

        if(file):
            self.f = myopen(file , "rt")
        else:
            self.f = stdin


    #####################################################

    def __enter__(self):
        return self

    #####################################################

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


    ######################################################

    def __getitem__(self, index):
        line = self.f.readline().strip()
        while line.startswith("#"):
            line = self.f.readline().strip()
        
        if line == "":
            raise IndexError
        line_contents = line.split("\t")
        if len(line_contents) < 9:
            raise IndexError
        return GTFEntry(line_contents)
                
    #########################################################

    def __del__(self):
        self.f.close()



def find_transcript_lengths(gtf_dict):
    """
    Find the total number of nucleotides in each transcript
    """

    for gene_name, transcripts in gtf_dict.items():
        for t_name, t_contents in transcripts.items():
            length = 0
            for start , end in t_contents["exons"]:
                length += (( end - start ) + 1)
            t_contents["length"] = length
    return gtf_dict


def get_gtf_contents(gtf_file):
    """
    Reads a gtf file into a dictionary. 
    
    NOTE 1:
    We tested this on encode gtf file.
    Some unexpected behavior might occur on GTF files from other sources.
    
    NOTE 2:
    For compatibility reasons, in gene_id and transcript_id, we are removing the part after
    ".". So, for example, ENSG00000178199.13 becomes ENSG00000178199
    
    It returns a dictionary where its keys are gene ids and values are again (sub)dictionaries
    where each (sub dictionary is of the form)
        transcript_id -> { "exons": list(), "CDS": list(), "strand": k.strand,
                           "start": k.start, "end": k.end, length: int}
                           
    Later when we want to determine relative UTR, CDS coordinates  etc. 
    We can place them in these (sub)dictionaries
    """
    gtf_contents = dict()
    with GTFfile(gtf_file) as gtf:
        for k in gtf:
            g_contents = k.attribute_dict["gene_id"].strip().split(".")
            
            # Discard Y chromosome paralogs
            if len(g_contents) >=2 and "PAR_Y" in g_contents[1]:
                continue

            gene_id = g_contents[0]

            
            if k.feature == "gene":
                gtf_contents[gene_id] = dict()
                
            if k.feature == "transcript":
                thranscript_id = k.attribute_dict["transcript_id"].strip().split(".")[0]
                gtf_contents[gene_id][thranscript_id] =\
                  { "exons": list(), "CDS": list(), "strand": k.strand,
                    "start": k.start, "end": k.end,\
                    "gene_type": k.attribute_dict["gene_type"],
                    "chr": k.seqname}
                
            if k.feature == "exon":
                thranscript_id = k.attribute_dict["transcript_id"].strip().split(".")[0]
                gtf_contents[gene_id][thranscript_id]["exons"].append( (k.start, k.end) )
                
            if k.feature == "CDS":
                # The third element gives us the corresponding exon of the CDS
                # we subtract 1 to adjust for array index
                thranscript_id = k.attribute_dict["transcript_id"].strip().split(".")[0]
                gtf_contents[gene_id][thranscript_id]["CDS"].append( \
                (k.start, k.end, int( k.attribute_dict["exon_number"] ) - 1 ) )
                
    result = find_transcript_lengths(gtf_contents)
    return result
