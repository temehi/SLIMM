#!/usr/bin/env python
from helper_methods import *

parser = argparse.ArgumentParser(
    description='Download reference genomes using a list compiled by select_refs.py ')

parser.add_argument('-t', '--tsv_file', type=str, required=True,
                    help = 'Path to a TSV file containing genomes to download and their corresponding FTP path. '
                            '(a result from running select_refs.py)')
parser.add_argument('-tr', '--threads',  type=int, choices=xrange(1, 11), default = 1,
                    help = 'number of threads for downloading in parallel in the range 1..10 (default: 1)')
parser.add_argument('-bd', '--buffer_dir', type=str, default = "",
                    help = 'directory containing already downloaded files to avoid repetitive downloads. '
                            'Symbolic link will be created to already existing files in this directory')
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help = 'Path of a directory where (intermediate) results will be saved')


args = parser.parse_args()

output_dir = args.output_dir
tsv_file = args.tsv_file
buffer_dir = args.buffer_dir
threads = args.threads
testing = False

use_buffer = True;

if not os.path.isdir(buffer_dir):
    use_buffer = False
    print "the provided buffer directory: '"+ buffer_dir + "' doesn't exist!"
    print "all files will be downloaded!"

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

to_download = {}
to_symlink = {}

inpf = open(tsv_file, 'r')
for line in inpf:
    if line[0] == "#":
        continue
    # if not firstLine:
    l = line.replace("\n", "").split('\t')
    taxid = l[0]
    species_taxid = l[2]
    ftp_path = l[6]
    fna_file_name = l[7]
    if (use_buffer and os.path.isfile(buffer_dir + "/" + fna_file_name) ) :
        to_symlink[fna_file_name] = os.path.abspath(buffer_dir + "/" + fna_file_name)
    else :
        to_download[fna_file_name] = ftp_path + '/' + fna_file_name
inpf.close()

print str(len(to_symlink)) + " Genomes will be symbolically linked to locally existing files!"
print str(len(to_download)) + " Genomes will be downloaded!"

for fna_file_name in to_symlink:
    dst = output_dir + '/' + fna_file_name
    os.symlink(to_symlink[fna_file_name], dst)
print str(len(to_symlink)) + " symbolic links created"


download_queue = Queue.Queue()
# DELETE FOR TESTING ONLY
# for line in inpf:
count = 0
for fna_file_name in to_download:
    if count == 10 and testing:
        break
    count += 1
    download_item = [fna_file_name, to_download[fna_file_name]]
    download_queue.put(download_item) # produce

total_to_download = download_queue.qsize()
while not download_queue.empty():
    if threading.activeCount() <= threads and not download_queue.empty():
        threading.Thread(target=download_one, args=(download_queue, output_dir)).start()
        sys.stdout.write('\r')
        sys.stdout.write("%d of %d remaining ..." % (download_queue.qsize(), total_to_download))
        sys.stdout.flush()
sys.stdout.write('\n')
