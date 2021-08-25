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
    print ('prep_reference.py -i <inputfilename> <optional>\n')
    print('Dependencies are: minimap2, bwa, seqkit, bioawk\n')
    print('This script is used to create bedfiles, genome size file, and karyoploteR genome file from reference.\nReference is indexed using minimap2, unless set to false. Can also index with bwa.\n')
    print('Requires gzip, pandas, and re libraries for python\n')
    print('-i      Input reference file in fastX format\n')
    print('-c      list of chromosomes you wish to extract and use for file creation\n')
    print('-r      list of replacement names which will rename the chromosomes provided with -c\n')
    print('-o      Change the output file name for for filtered or replaced reference file. Using -c the output is <inputfilename>.filtered.fastq. Using -c and -r the output is <inputfilename>.replaced.fastq\n')
    print('-m      Minimap2 ONT indexing of reference sequence. Default true. Use false for no indexing\n')
    print('-b      bwa index of reference genome. Default false. Use true for indexing\n')
    print('-a      Create bedfiles. Default is true. Use false to change setting\n')
    print('-g      Create genome size file. Default is true. Use false to change setting\n')
    print('-k      Create genome file for karyoploteR. Default is true. Use false to change setting.\n')
    print('-s      Create file containing size of longest reference chromosome. Default is true. Use false to change setting\n')
    print('-d      desitination folder\n')


def karyomerge(names, startpos, lengths):
    merged_list = [(names[i], startpos[i], lengths[i]) for i in range(0, len(names))]
    return merged_list

def create_karyoploteR(longest_chromo,output_name, seq_list, seq_count):
    if longest_chromo > 536870912:
        print('One or more reference chromosomes exceed the length limit of 536,870,912 bp for SAMtools bai index\nbai index is necessary to create karyoploteR plots\nkaryoploteR genome file will not be output.')
    print('Creating karyoploteR genome file:')

    readnames = []
    readlengths = []
    for readname, readlength in seq_list:
        readnames.append(readname)
        readlengths.append(readlength+1)

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
        seq_list.append((record.id, len(record.seq)))

    longest_chromo = max(seq_list, key=lambda item:item[1])[1]
    print('Longest chromo:', longest_chromo)

    return seq_count, seq_list, longest_chromo


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
    seq_count, seq_list, longest_chromo = get_seq_records(input_name)
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

    create_karyoploteR(longest_chromo, output_name, seq_list, seq_count)

if __name__ == "__main__":
   main(sys.argv[1:])
