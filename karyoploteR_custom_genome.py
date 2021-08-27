import subprocess
import warnings
import sys, getopt
import os
import re
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def print_help():
    print ('karyoploteR_custom_genome.py -i REFERENCE <optional>\n')
    print('Dependencies are: minimap2, bwa, seqkit, bioawk\n')
    print('This script is used to create bedfiles, genome size file, and karyoploteR genome file from reference.\nReference is indexed using minimap2, unless set to false. Can also index with bwa.\n')
    print('Requires gzip, pandas, and re libraries for python\n')
    print('-i      Input reference file in fastX format\n')
    print('-o      Change the output file name Default is REFERENCE_karyoploteR_genome_file.txt\n')
    print('-d      desitination folder\n')


BAI_INDEX_SIZE = 536870912
class Sequence:
    # def __init__(self):
    #     self.name = ""
    #     self.seq = []
    #
    def __init__(self,n,l):
        self.name = n
        self.len = l


def karyomerge(names, startpos, lengths):
    merged_list = [(names[i], startpos[i], lengths[i]) for i in range(0, len(names))]
    return merged_list

def create_karyoploteR(output_name, seq_list, seq_count):
    print('Creating karyoploteR genome file:')

    readnames = []
    readlengths = []
    for seq in seq_list:
        readnames.append(seq.name)
        readlengths.append(seq.len+1)

    start_position = []
    for z in range (0, seq_count):
        start_position.append(1)

    karyotuples = karyomerge(readnames,start_position,readlengths)

    tab = pd.DataFrame(karyotuples)
    tab.columns = ['chr', 'start', 'end']
    print(tab)
    tab.to_csv(output_name, sep='\t', index=False, header=False)


def output_name_arg(output_name,input_name, destination_folder):
    if output_name == '':
        path_name = os.path.basename(input_name)
        if path_name.endswith('.gz'):
            path_name = os.path.splitext(path_name)[0]
        path_name = os.path.splitext(path_name)[0]
        output_name = '%s_karyoploteR_genome_file.txt' % path_name

    if destination_folder != '':
        output_name = os.path.join(destination_folder, output_name)

    return output_name

def get_seq_records(input_name):
    if input_name.endswith('.gz'):
        fin=gzip.open(input_name, "rt")
    else:
        fin=open(input_name, "rt")

    seq_list = []
    seq_count = 0
    print('Loading reference')
    for record in SeqIO.parse(fin, 'fasta'):
        seq_count += 1
        seq_list.append(Sequence(record.id, len(record.seq)))

    return seq_count, seq_list

def max_len(seq_list):
    seqlist_len = seq_list.len
    return seqlist_len

def sort_seq_list(seq_list):
    seqlist_name = seq_list.name
    lower_case = seqlist_name.lower()
    return lower_case

def main(argv):
    input_name = ''
    output_name= ''
    destination_folder = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:k:d:",["help=", "inputfilename=", "outputfilename=","destinationfolder="])
    except getopt.GetoptError:
        print ('Reference.bed.karyoploteR.py -i <inputfilename>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print_help()
            sys.exit()
        elif opt in ("-i", "--inputfilename") :
            input_name = arg
        elif opt in ("-o", "--outputfilename"):
            output_name = arg
        elif opt in ("-d", "--desitination folder"):
            destination_folder = arg

    output_name = output_name_arg(output_name,input_name,destination_folder)
    seq_count, seq_list = get_seq_records(input_name)
    print(type(seq_list))

### checking number of chromosomes present in file or list

    print('Number of chromosomes/scaffolds in reference:', seq_count)

    yes = {'yes','y', 'ye', '',}
    no = {'no','n'}

    while True:
        if seq_count > 50:
            print('Scaffolds in reference exceed 50. It is recommended to merge the scaffolds. \nDo you still wish to continue? [Y/N]')
            choice = input().lower()
            if choice in yes:
               print('Will proceed using all scaffolds')
               break
            elif choice in no:
               print('Program aborted')
               quit()
            else:
               print("Please respond with 'Y','N'")
               continue
        else:
            break

    sorted_seq_list = sorted(seq_list, key = sort_seq_list)

    create_karyoploteR(output_name, sorted_seq_list, seq_count)

if __name__ == "__main__":
   main(sys.argv[1:])
