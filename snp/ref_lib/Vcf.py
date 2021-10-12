
import gzip

VCF_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO",    "FORMAT",  "CAST_EiJ"]

class VcfEntry:

    def __init__(self , vcf_line_contents ):
        assert len(vcf_line_contents) >= len(VCF_FIELDS)
        
        self.fields = { VCF_FIELDS[i] : vcf_line_contents[i] for i in range( len(VCF_FIELDS) ) }
        
        
    def __str__(self ):
        """
        This needs to be rewritten 
        """
        return "\t".join( [self.fields[f] for f in VCF_FIELDS] )


############################################################################
    
class VcfFile:
    '''
    This is a reader for 
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
        #line_contents = line.split("\t")
        line_contents = line.split()
        if len(line_contents) < 9:
            raise IndexError
        return VcfEntry(line_contents)
                
    #########################################################

    def __del__(self):
        self.f.close()

