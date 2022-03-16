"""
Motivation: CRISPRCas discovery and subtyping
This script will run CRISPRDetect on a FASTA or GFF for a single genome.
A FASTA file of unique spacers will be created with CD-HIT-EST.


Author: Ian Rambo
Contact: ian.rambo@utexas.edu

Thirteen... that's a mighty unlucky number... for somebody!
"""
version = 'v1.0.0.1'

from Bio import SeqIO
import argparse
import logging
import os
import numpy as np
import pandas as pd
import subprocess
import magic
import gzip
import re
from datetime import datetime
import shell_tools
from gff3 import gff3_to_pddf
#==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', type=str, dest='fasta_file', action='store',
help='input nucleotide FASTA file. Required.')
parser.add_argument('-o', '--out_root', type=str, dest='out_root', action='store',
help='root output directory. Required.')
parser.add_argument('-s', '--tmp_dir', type=str, dest='tmp_dir', action='store',
default='/tmp', help='temporary file directory. Default = /tmp')
parser.add_argument('-t', '--threads', type=int, dest='threads', action='store', default=1,
help='number of threads for each command. Default = 1')
parser.add_argument('-j', '--joblog', type=str, dest='joblog', action='store',
nargs = '?', help='Path to logging joblog.')

parser.add_argument('-p', '--prefix', type=str, dest='prefix', action='store', nargs='?',
help='optional prefix for output files. Uses the input nucleotide fasta basename if not supplied.')

#CRISPRDetect 2.4 Options
parser.add_argument('-c', '--crispr_detect_dir', type=str, dest='CRISPRDetectDir', action='store',
help='directory containing CRISPRDetect.pl', default='/build/CRISPRDetect_2.4/CRISPRDetect_2.4')
parser.add_argument('-v', '--crispr_qual_cutoff', type=int, dest='crispr_qual_cutoff', action='store',
default=3, help='exclude CRISPR arrays with CRISPRDetect quality score less than this value. Default = 3')
parser.add_argument('-a', '--minimum_repeat_length', type=int, dest='minimum_repeat_length', action='store',
default=23, help='minimum length of CRISPR repeats. Default = 23')
parser.add_argument('-n', '--minimum_no_of_repeats', type=int, dest='minimum_no_of_repeats', action='store',
default=3, help='minimum number of CRISPR repeats. Default = 3')
parser.add_argument('-i', '--left_flank_length', type=int, dest='left_flank_length', action='store',
default=500, help="length of the 5' (upstream) region of the CRISPRs. Default = 500")
parser.add_argument('-d', '--right_flank_length', type=int, dest='right_flank_length', action='store',
default=500, help="length of the 3' (downstream) region of the CRISPRs. Default = 500")
parser.add_argument('-q', '--quiet', type=int, dest='quiet', action='store',
default=0, help="run CRISPRDetect in quiet mode (1). Default = 0")

#CD-HIT-EST Options
parser.add_argument('-k', '--seq_id_thresh', type=float, dest='seq_id_thresh', action='store',
default=1.0, help="CD-HIT-EST sequence identity threshold for spacers. Default = 1.0")
parser.add_argument('-b', '--band_width', type=int, dest='band_width', action='store',
default=20, help="CD-HIT-EST aligment bandwidth. Default = 20")
parser.add_argument('-l', '--desc_length', type=int, dest='desc_length', action='store',
default=50, help="CD-HIT-EST description length for .clstr file. Default = 50")
parser.add_argument('-w', '--word_length', type=int, dest='word_length', action='store',
default=11, help="CD-HIT-EST word length. Default = 11")

opts = parser.parse_args()
#==============================================================================
#Set up logger
logging_format = '%(levelname)s %(asctime)s - $(message)s'
if opts.joblog:
    #Create the joblog directory if not specified
    if not os.path.exists(os.path.dirname(opts.joblog)):
        os.makedirs(os.path.dirname(opts.joblog))
    else:
        pass
    logging.basicConfig(filename = opts.joblog, level = logging.DEBUG, format = logging_format)
else:
    logging.basicConfig(filename = os.path.join(opts.out_root, 'crisprdetect_cluster.%s.joblog' % str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))),
    level = logging.DEBUG, format = logging_format)

logger = logging.getLogger()
#==============================================================================
#Create the output directories
try:
    output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect']}

    for key, value in output_paths.items():
        if not os.path.isdir(value):
            os.makedirs(value)
        else:
            pass
except:
    print('no --out_root specified')
    logging.error('Cannot make output directories, no --out_root specified')


#==============================================================================
nt_fasta = opts.fasta_file
#Get the file basename to name output files
nt_fasta_basename = shell_tools.get_basename(nt_fasta)
#If --prefix is not set, use the name of the nucleotide fasta
prefix = opts.prefix
if not prefix:
    prefix = nt_fasta_basename
#==============================================================================
#If the nucleotide fasta input is gzipped, gunzip it for use with CRISPRDetect
gzip = False
if shell_tools.is_gzipped(nt_fasta):
    gzip = True
    print('%s is gzip compressed, gunzip file...' % nt_fasta)
    subprocess.run(['gunzip', nt_fasta], shell=False)
    nt_fasta = os.path.splitext(opts.fasta_file)[0]

else:
    pass
#==============================================================================
###---BEGIN CRISPRDetect---###
###---Run CRISPRDetect---###
logger.info('Run CRISPRDetect')
print('Begin CRISPRDetect...')

###---CRISPRDetect---###
#Pattern for CRISPRDetect output
crispr_detect_outpatt = os.path.join(output_paths['CRISPRDetect'], prefix + '_CRISPRDetect')
crispr_detect_log = os.path.join(output_paths['CRISPRDetect'], prefix + '_CRISPRDetect.log')

#CRISPRDetect options
crispr_detect_optdict = {'-f':nt_fasta,
                        '-o':crispr_detect_outpatt,
                        '-T':opts.threads,
                        '-minimum_repeat_length':opts.minimum_repeat_length,
                        '-array_quality_score_cutoff':opts.crispr_qual_cutoff,
                        '-tmp_dir':opts.tmp_dir,
                        '-logfile':crispr_detect_log,
                        '-left_flank_length':opts.left_flank_length,
                        '-right_flank_length':opts.right_flank_length,
                        '-minimum_no_of_repeats':opts.minimum_no_of_repeats,
                        '-quiet':opts.quiet}

#Path to CRISPRDetect executable
crispr_detect_exec = os.path.join(opts.CRISPRDetectDir, 'CRISPRDetect.pl')

crispr_detect_cmd = shell_tools.exec_cmd_generate(crispr_detect_exec, crispr_detect_optdict)
#Run CRISPRDetect
subprocess.run(crispr_detect_cmd, shell=False)

###---Read the GFF file produced by CRISPRDetect---###
crispr_detect_gff = crispr_detect_outpatt + '.gff'

#if os.path.exists(crispr_detect_gff) and os.stat(crispr_detect_gff).st_size != 0:
crispr_array_df = pd.DataFrame()
crispr_spacer_df = pd.DataFrame()

if os.stat(crispr_detect_gff).st_size > 0:
    try:
        #Convert the GFF to a pandas data frame, selecting full CRISPR arrays coords
        crispr_array_df = gff3_to_pddf(gff = crispr_detect_gff, ftype = 'repeat_region', index_col=False)
    except FileNotFoundError:
        logger.error('CRISPRDetect GFF file %s not found' % crispr_detect_gff)

    if not crispr_array_df.empty:
        #Split up attributes for CRISPR arrays into new columns
        crispr_array_df[['ID', 'Repeat', 'Dbxref', 'OntologyTerm', 'ArrayQualScore']] = crispr_array_df['attributes'].str.replace('[A-Za-z]+\=', '', regex=True).str.split(pat = ";", expand=True)
        #Select entries for spacers
        crispr_spacer_df = gff3_to_pddf(gff = crispr_detect_gff, ftype = 'binding_site', index_col=False)
        #Split up attributes for CRISPR spacers into new columns
        crispr_spacer_df[['ID', 'Name', 'Parent', 'Spacer', 'Dbxref', 'OntologyTerm', 'ArrayQualScore']] = crispr_spacer_df['attributes'].str.replace('[A-Za-z]+\=', '', regex=True).str.split(pat = ";", expand=True)

        #####=====SPACERS=====#####
        #Write the CRISPR spacers to an output nucleotide FASTA
        crispr_spacer_fna = os.path.join(output_paths['CRISPRDetect'], '%s_crispr_spacers.fna' % prefix)
        with open(crispr_spacer_fna, 'w') as spacer_fa:
            logger.info('writing spacers to %s' % crispr_spacer_fna)
            for index, row in crispr_spacer_df.iterrows():
                spacer_fasta_record = '>%s_____%s' % (row['seqid'], row['ID']) + '\n' + row['Spacer'] + '\n'
                spacer_fa.write(spacer_fasta_record)

        ###---Cluster unique spacers @ 100% identity with CD-HIT-EST---###
        crispr_spacers_cluster = os.path.join(output_paths['CRISPRDetect'], '%s_crispr_spacers_cd-hit-est_cluster100.fna' % prefix)
        ###---CD-HIT-EST---###
        cdhit_est_opts = {'-i':crispr_spacer_fna,
                          '-o':crispr_spacers_cluster,
                          '-c':opts.seq_id_thresh,
                          '-b':opts.band_width,
                          '-d':opts.desc_length,
                          '-T':opts.threads,
                          '-n':opts.word_length}
        cdhit_est_spc_cmd = shell_tools.exec_cmd_generate('cd-hit-est', cdhit_est_opts)
        subprocess.run(cdhit_est_spc_cmd)

    else:
        print('CRISPRDetect GFF file has entries', \
              'CRISPR pandas df is empty')
else:
    print('No putative CRISPRs found with CRISPRDetect for', opts.fasta_file)

#Compress the input if it was gzipped originally
if gzip:
    print('re-gzip compressing file %s' % opts.fasta_file)
    subprocess.run(['gzip', nt_fasta], shell=False)
