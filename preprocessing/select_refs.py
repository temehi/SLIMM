#!/usr/bin/env python
from helper_methods import *

parser = argparse.ArgumentParser(description =
''' Download reference genomes of microorganisms
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-wd', '--work_dir', type=str, required=True,
                    help = 'The path of working directory where (intermediate) results will be saved')
parser.add_argument('-g', '--groups',  type=str, default = "AB",
                    help = '''Which group of microbes to consider any combination of the letters [A], [B] and [V]
                                where B =  Bacteria, A = Archaea and V = Viruses and Viroids (default: AB)''')
parser.add_argument('-s', '--species_only', dest='species_lv', action='store_true',
                    help = 'download one reference per species.')
parser.add_argument('-c', '--complete', dest='complete', action='store_true',
                    help = 'download only complete genomes.')
parser.add_argument('-t', '--taxa_ids',  type=str, default = "",
                    help = '''comma separated list of taxonomic ids to be included (in addition to --groups) into
                                the reference database. This way you might even add the genome of Eukaryotes.
                                 e.g. the host genome''')
parser.add_argument('-d', '--database', type=str, choices = ['refseq', 'genbank'], default = 'refseq',
                    help = 'From which database references should be downloaded  (default: refseq)')


args = parser.parse_args()

parallel = args.threads
working_dir = args.workdir
groups = args.groups
only_species = args.species_lv
only_complete = args.complete
db_choice = args.database
taxid_list = args.taxa_ids


##############################################################################
# For KNIME workflow only
##############################################################################
# parallel        = flow_variables['threads']
# working_dir     = flow_variables['workdir']
# groups          = flow_variables['groups']
# only_species    = flow_variables['species_lv']
# db_choice       = flow_variables['database']
# testing         = flow_variables['testing']

subset_taxids   = []
if len(taxid_list) >= 1:
    subset_taxids   = map(int, taxid_list.split(','))

groups_name     = groups
if len(groups_name) < 1:
    groups_name = "CUSTOM"
elif len(subset_taxids) >= 1:
    groups_name += "_CUSTOM"

if not os.path.isdir(working_dir):
    os.makedirs(working_dir)
# else:
#     empty = os.listdir(working_dir) == [];
#     if not empty:
#         print ("[ERROR!] Working directory [" + working_dir + "] should be empty!")
#         sys.exit(0)

today_string = (datetime.datetime.now()).strftime("%d%m%Y")
genomes_dir = working_dir + "/genomes_" + today_string
if not os.path.isdir(genomes_dir):
    os.makedirs(genomes_dir)

slimmDB_dir = working_dir + "/slimmDB_" + today_string
if not os.path.isdir(slimmDB_dir):
    os.makedirs(slimmDB_dir)

##############################################################################
# Get the path of Names, nodes, catagories file from the extracted folder
##############################################################################

taxdmp_extract_dir = taxonomy_download("taxdump", working_dir, today_string)
taxcat_extract_dir = taxonomy_download("taxcat", working_dir, today_string)


names_path = taxdmp_extract_dir + "/names.dmp"
nodes_path = taxdmp_extract_dir + "/nodes.dmp"
catagories_path = taxcat_extract_dir + "/categories.dmp"
reduced_names_path = slimmDB_dir + "/names.dmp"
reduced_nodes_path = slimmDB_dir + "/nodes.dmp"
taxaid2sp_path=nodes_path.replace("nodes.dmp", "taxaid2sp_" + groups_name + ".dmp")


assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
if db_choice == "genbank":
    assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

assembly_summary_file = working_dir + "/assembly_summary_" + db_choice + "_" + today_string + ".txt"

genomes_to_download_path = working_dir + "/" + groups_name + "_genomes_to_download.txt"
ncbi_get_script_path =  working_dir + "/" + groups_name + "_genomes_ncbi_download.sh"

print "Downloading assembly_summary file ..."
print assembly_summary_url
urllib.urlretrieve(assembly_summary_url, assembly_summary_file)


################################################################################################
# From the assembly summary file, identify, download and annotate the genomes to download.
# Unique genomes per taxon will be downloaded. Genomes are selected in the following order:
# 1. Complete genomes (longest first)
# 2. Chromosomes (longest first)
################################################################################################
taxid_col = 5
if only_species :
    taxid_col = 6

intial_taxids = {}
inpf = open(catagories_path, 'r')
for line in inpf:
    l = line.split('\t')
    if l[0] in groups or int(l[1]) in subset_taxids or  int(l[2]) in subset_taxids:
        intial_taxids[int(l[1])] = 1
        intial_taxids[int(l[2])] = 1
inpf.close()
taxid_genomes =  {}
inpf = open(assembly_summary_file, 'r')
firstLine = True

valid_assembly_levels=["Complete Genome"]
if only_complete:
    print "Downloading only complete genomes ..."
else:
    valid_assembly_levels = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]
# 0   # assembly_accession
# 1   bioproject
# 2   biosample
# 3   wgs_master
# 4   refseq_category
# 5   taxid
# 6   species_taxid
# 7   organism_name
# 8   infraspecific_name
# 9   isolate
# 10  version_status
# 11  assembly_level
# 12  release_type
# 13  genome_rep
# 14  seq_rel_date
# 15  asm_name
# 16  submitter
# 17  gbrs_paired_asm
# 18  paired_asm_comp
# 19  ftp_path
# 20  excluded_from_refseq
# 21  relation_to_type_material

for line in inpf:
    if line[0] == "#":
        continue
    # if not firstLine:
    if has_valid_genome_type(line, valid_assembly_levels):
        l = line.split('\t')
        taxid = int(l[taxid_col])
        path = l[19].replace("\n", "")
        if (taxid in intial_taxids and (path != 'na') and l[10] == "latest") :
            refseq_category = l[4] #[0]
            species_taxid = l[6] #[1]
            assembly_level = l[11] #[2]
            genome_rep = l[13] #[3]
            seq_rel_date = l[14] #[4]
            ftp_path = path  #[5]
            fna_file_name = path[path.rfind('/')+1:] + "_genomic.fna.gz" #[6]
            if taxid in taxid_genomes :
                taxid_genomes[taxid].append([refseq_category, species_taxid, assembly_level, genome_rep, seq_rel_date, ftp_path, fna_file_name])
            else:
                taxid_genomes[taxid] = []
                taxid_genomes[taxid].append([refseq_category, species_taxid, assembly_level, genome_rep, seq_rel_date, ftp_path, fna_file_name])
    firstLine = False
inpf.close()
outf = open(genomes_to_download_path, 'w')
count = 0
smallest_date = "0001/01/01"
for taxid in taxid_genomes:
    found = False
    candidate_op = taxid_genomes[taxid]
    selected_op = []
    best_date=smallest_date
    for options in candidate_op: # try to look for the latest reference genome
        if (options[0] == "reference genome" and best_date < options[4]) :
            selected_op = options
            best_date=options[4]
            found = True
    if not found :
        best_date=smallest_date
        for options in candidate_op :# try to look for the latest representative genome
            if (options[0] == "representative genome" and best_date < options[4]) :
                selected_op = options
                best_date=options[4]
                found = True
                break;
    if not found : # try with different assembly levels
        for genome_type in valid_assembly_levels :
            best_date=smallest_date
            for options in candidate_op:
                if (options[2] == genome_type and best_date < options[4]) :
                    selected_op = options
                    best_date=options[4]
                    found = True
            if found:
                break
    if found:
        outf.write(str(taxid)+ '\t' + '\t'.join(selected_op) + "\n")
        count += 1
outf.close()

print str(count) + " Genomes matching the criteria were found!"
print "Their FTP path for downloading them are written at:"
print genomes_to_download_path