#!/usr/bin/env python
from helper_methods import *

parser = argparse.ArgumentParser(description =
''' Download reference genomes of microorganisms
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help = 'The path of working directory where (intermediate) results will be saved')
parser.add_argument('-g', '--groups',  type=str, default = "AB",
                    help = '''Which group of microbes to consider any combination of the letters [A], [B], [V] and [F]
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

output_dir = args.output_dir
groups = args.groups
only_species = args.species_lv
only_complete = args.complete
db_choice = args.database
taxid_list = args.taxa_ids

organism_group = {'A':'archaea', 'B':'bacteria', 'V':'viral', 'F':'fungi'}
subset_taxids   = []
if len(taxid_list) >= 1:
    subset_taxids   = map(int, taxid_list.split(','))

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

sp_string = "ST"
taxid_col = 5
if only_species :
    taxid_col = 6
    sp_string = "SP"

complete_string = "ALL"
valid_assembly_levels=["Complete Genome"]
if only_complete:
    complete_string = "CMP"
    print "Downloading only complete genomes ..."
else:
    valid_assembly_levels = ["Complete Genome", "Chromosome", "Scaffold", "Contig"]

taxid_genomes =  {}
for g in groups:
    ################################################################################################
    # Get the assembly file for the current specified group.
    # If a file exists from the same day skip downloading.
    ################################################################################################
    assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/" + db_choice + "/" + organism_group[g] + "/assembly_summary.txt"
    assembly_summary_file = output_dir + "/assembly_summary_" + db_choice + "_"  + organism_group[g] + ".txt"
    if (os.path.isfile(assembly_summary_file) ) :
        print "Using an existing assembly_summary file ..."
        print assembly_summary_file
        print "[WARNING] delete this file if you want to update to a new set of assemblies!"
    else:
        print "Downloading assembly_summary file ..."
        print assembly_summary_url
        urllib.urlretrieve(assembly_summary_url, assembly_summary_file)
    ############################################################################################
    # From the assembly summary file, identify the genomes to download based on criterion given.
    # Unique genomes per taxon will be downloaded. Genomes are selected in the following order:
    # 1. Complete genomes (longest first)
    # 2. Chromosomes (longest first)
    ############################################################################################
    # 0   # assembly_accession          # 11  assembly_level
    # 1   bioproject                    # 12  release_type
    # 2   biosample                     # 13  genome_rep
    # 3   wgs_master                    # 14  seq_rel_date
    # 4   refseq_category               # 15  asm_name
    # 5   taxid                         # 16  submitter
    # 6   species_taxid                 # 17  gbrs_paired_asm
    # 7   organism_name                 # 18  paired_asm_comp
    # 8   infraspecific_name            # 19  ftp_path
    # 9   isolate                       # 20  excluded_from_refseq
    # 10  version_status                # 21  relation_to_type_material
    inpf = open(assembly_summary_file, 'r')
    for line in inpf:
        if line[0] == "#":
            continue
        # if not firstLine:
        if has_valid_genome_type(line, valid_assembly_levels):
            l = line.split('\t')
            taxid = int(l[taxid_col])
            path = l[19].replace("\n", "")
            if ((path != 'na') and l[10] == "latest") :
                refseq_category = l[4] #[0]
                species_taxid = l[6] #[1]
                assembly_level = l[11] #[2]
                genome_rep = l[13] #[3]
                seq_rel_date = l[14] #[4]
                ftp_path = path  #[5]
                fna_file_name = path[path.rfind('/')+1:] + "_genomic.fna.gz" #[6]
                if taxid in taxid_genomes :
                    taxid_genomes[taxid].append([refseq_category, species_taxid, assembly_level, genome_rep, seq_rel_date, ftp_path, fna_file_name, g])
                else:
                    taxid_genomes[taxid] = []
                    taxid_genomes[taxid].append([refseq_category, species_taxid, assembly_level, genome_rep, seq_rel_date, ftp_path, fna_file_name, g])
    inpf.close()

################################################################################################
# Write the selected genomes to file  select_refs_refseq_SP_ALL_B.log
################################################################################################
selected_genomes_tsv_path = output_dir + "/selected_refs_" + db_choice + "_" + sp_string + "_" + complete_string + "_" + groups  + ".tsv"
outf = open(selected_genomes_tsv_path, 'w')

organism_group_count = {'A':0, 'B':0, 'V':0, 'F':0}
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
        organism_group_count[selected_op[7]] += 1
outf.close()

detail_count = ""
total_count = 0
for g in groups:
    total_count += organism_group_count[g]
    detail_count += organism_group[g] + " = " + str (organism_group_count[g]) + ", "
print str(total_count) + " Genomes matching the criteria were found!"
print "(" + detail_count[:len(detail_count)-2] + ")"
print "Their FTP path for downloading them are written at:"
print selected_genomes_tsv_path