usage: catPiss_v.1_0.py [-h] [-o OUT] [-u UPSTREAM] [-d DOWNSTREAM]
                        [-s SORTBY] [-n DENOVOVAR] [-f REFTABLEFILE]
                        [-t {Cas9,Cas12A}] [-c CASFILEDATA] [-p CLONEDEST]
                        geneNames speciesName {activation,interference,custom}

The 'catPiss' program identifies a list of guide RNAs for CRISPR-based gene
attenuation, promotion, or custom mode for a given list of genes and
associated species. The three required inputs are a gene (or a list of genes),
the species of interest (see pattern requirements below), and CRISPR mode. In
addition, the user must also indicate either deNovo scoring or use of a
reference scoring table. The output is a .bed file that contains a list of
gRNAs for each gene and the associated primers for those gRNAs. If you want to
make deNovo gRNA (-n) from a fasta, -t and -c must be set. If you want to use
your own list of gRNAs downloaded from UCSC, -f must be set to the file path.

positional arguments:
  geneNames             The specific gene(s) of interest entered as a string
                        or the name of the file with the list of genes.
  speciesName           The name of the target species. Uses underscore for
                        spaces (e.g., homo_sapiens)
  {activation,interference,custom}
                        User determines if goal is gene activation, gene
                        interference, or a custom mode. This argument sets the
                        upstream and downstream boundaries relative to TSS. If
                        "custom" is selected, then "-u" and "-d" also need to
                        be defined

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --output OUT  Optional: supply output file directory name, otherwise
                        write to program default
  -u UPSTREAM, --upstream UPSTREAM
                        Optional: supply distance upstream of the TSS site to
                        check for gRNAs
  -d DOWNSTREAM, --downstream DOWNSTREAM
                        Optional: supply distance downstream of the TSS site
                        to check for gRNAs
  -s SORTBY, --sortby SORTBY
                        Optional: if 'de novo' is 'True' sort by is an option
                        for sorting by score using either 'Doench2016_perc',
                        'Doench2016_score', 'Moreno_Matos_perc',
                        'Moreno_Matos_score', 'MIT_specificity'
  -n DENOVOVAR, --deNovo DENOVOVAR
                        Choice of using a de novo table from this program
                        ('True')
  -f REFTABLEFILE, --refTable REFTABLEFILE
                        User provided reference table of scored gRNAs.
  -t {Cas9,Cas12A}, --casType {Cas9,Cas12A}
                        Specify type of cas: Cas9 or Cas12A. This must be
                        entered if deNovo = 'True'
  -c CASFILEDATA, --casFile CASFILEDATA
                        Identify file name with cas protein information
                        detailed. This should be in the working directory.
                        This must be entered if deNovo = 'True'
  -p CLONEDEST, --cloning CLONEDEST
                        Identify the cloning strategy for generation of
                        primers. Default is pX330. Other options for cloning
                        strategies are in 'cloningstrategy_primers.txt, which
                        MUST be in your working directory.


##### Examples:

Command line for a deNovo run on Homo sapiens:
./catPiss_v.1_0.py genelist.txt homo_sapiens activation -n True -t Cas9 -c cas_prot_list.tsv

Command line for a run using an external gRNA 
./catPiss_v.1_0.py genelist.txt homo_sapiens activation -f UCSC_DownloadGrnas.txt


#### Dependent files
Required Python files:
catPiss_v.1_0.py
genelist_to_TSS_bed.py
gff3_to_TSSbed.py
grna_gen.py
overlap.py
targetsite_to_primers.py

Required TSV:
cat_prot_list.tsv --- place in working directory

Required TXT:
UCSC_DownloadGrnas.txt --- for '-f' option --- place in working directory
genelist.txt --- place in working directory
cloningstrategy_primers.txt --- place in working directory

Optional TXT:
cat.txt
