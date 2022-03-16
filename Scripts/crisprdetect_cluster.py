"""
Motivation: CRISPRCas discovery and subtyping
This script will run CRISPRDetect on a FASTA or GFF for a single genome.
A FASTA file of unique spacers will be created with CD-HIT-EST.


Author: Ian Rambo
Contact: ian.rambo@utexas.edu, imrambo@lbl.gov

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
import gene_clusters
#==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument('--fasta', type=str, dest='fasta_file', action='store',
help='input nucleotide FASTA file. Required.')
parser.add_argument('--out_root', type=str, dest='out_root', action='store',
help='root output directory. Required.')
parser.add_argument('--tmp_dir', type=str, dest='tmp_dir', action='store',
default='/tmp',
help='temporary file directory. Default = /tmp')
parser.add_argument('--threads', type=int, dest='threads', action='store', default=1,
help='number of threads for each command. Default = 1')
parser.add_argument('--joblog', type=str, dest='joblog', action='store',
nargs = '?', help='Path to logging joblog.')

parser.add_argument('--prefix', type=str, dest='prefix', action='store', nargs='?',
help='optional prefix for output files. Uses the input nucleotide fasta basename if not supplied.')
parser.add_argument('--prodigal_mode', type=str, dest='prodigal_mode', action='store',
default='single', help='Prodigal 2.6.3 search mode - choose single or meta. Default = single')
parser.add_argument('--clst_window_extent', type=int, dest='clst_window_extent', action='store', default=10000,
help='Integer bp extension to search for ORFs near a CRISPR array. Default = 10000.')

#CRISPRDetect 2.4 Options
parser.add_argument('--crispr_detect_dir', type=str, dest='CRISPRDetectDir', action='store',
help='Directory for CRISPRDetect.pl', default='/build/CRISPRDetect_2.4/CRISPRDetect_2.4')
parser.add_argument('--crispr_qual_cutoff', type=int, dest='crispr_qual_cutoff', action='store',
default=3, help='Exclude CRISPR arrays with CRISPRDetect quality score less than this value. Default = 3')
parser.add_argument('--minimum_repeat_length', type=int, dest='minimum_repeat_length', action='store',
default=23, help='Minimum length of CRISPR repeats. Default = 23')
parser.add_argument('--minimum_no_of_repeats', type=int, dest='minimum_no_of_repeats', action='store',
default=3, help='Minimum number of CRISPR repeats. Default = 3')
parser.add_argument('--left_flank_length', type=int, dest='left_flank_length', action='store',
default=500, help="length of the 5' (upstream) region of the CRISPRs. Default = 500")
parser.add_argument('--right_flank_length', type=int, dest='right_flank_length', action='store',
default=500, help="Length of the 3' (downstream) region of the CRISPRs. Default = 500")

# parser.add_argument('--database', type=str, dest='database', action='store',
# nargs='?', help='Comma-separated list of HMM database paths to use for CRISPR-proximity gene hmmsearch. You can also specify paths to certain profiles, e.g. /path/to/K00001.hmm')
# parser.add_argument('--profile_suffix', type=str, dest='profile_suffix', action='store', default='.hmm',
# help='suffix for HMM gene profiles. Default is ".hmm"')
# parser.add_argument('--ccs_typing', type=str, dest='ccs_typing', action='store', default='sub-typing',
# help='Level of CRISPR-Cas system typing specificity. Choose from: [general, typing, sub-typing]')
# parser.add_argument('--macsy_dbtype', type=str, dest='macsy_dbtype', action='store',
# default='ordered_replicon', help='Specify dataset type for MacSyFinder. [unordered_replicon, unordered, ordered_replicon, gembase]')
# parser.add_argument('--macsy_eval', type=float, dest='macsy_eval', action='store',
# default=1e-4, help='E-value cutoff for MacSyFinder. Default = 1e-4')
# parser.add_argument('--macsy_coverage', type=float, dest='macsy_coverage', action='store',
# default=0.4, help='Minimum profile coverage for MacSyFinder. Default = 0.4')
# parser.add_argument('--macsy_systems', type=str, dest='macsy_systems', action='store',
# default='all', help='Systems to search for with MacSyFinder, e.g. CasIF. Default = all')
# parser.add_argument('--profile_dir', type=str, dest='profile_dir', action='store',
# default='/build/data/profiles', help='Path to directory containing HMMs for MacSyFinder. Default = /build/data/profiles')


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
    logging.basicConfig(filename = os.path.join(opts.tmp_dir, 'CRISPRFlow.%s.joblog' % str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))),
    level = logging.DEBUG, format = logging_format)

logger = logging.getLogger()
#==============================================================================
#Create the output directories
try:
    output_paths = {p: os.path.join(opts.out_root, p) for p in ['CRISPRDetect','Prodigal','Cluster']}

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
                        '-minimum_no_of_repeats':opts.minimum_no_of_repeats}

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
                          '-c':'1.0',
                          '-b':'20',
                          '-d':'50',
                          '-T':opts.threads,
                          '-n':'11'}
        cdhit_est_spc_cmd = shell_tools.exec_cmd_generate('cd-hit-est', cdhit_est_opts)
        subprocess.run(cdhit_est_spc_cmd)
        #==============================================================================
        ###---BEGIN Prodigal---###
        """
        Predict genes in contigs containing a putative CRISPR array.
        """
        prodigal_outfmt = 'gff'
        prodigal_out = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.%s' % prodigal_outfmt)
        prodigal_aa = os.path.join(output_paths['Prodigal'], prefix + '_prodigal.faa')

        prodigal_opts = {'-o':prodigal_out, '-a':prodigal_aa, '-p':opts.prodigal_mode, '-i':nt_fasta}

        #Generate and run the Prodigal command
        prodigal_cmd = shell_tools.exec_cmd_generate('prodigal', prodigal_opts)
        subprocess.run(prodigal_cmd, shell = False)


        prodigal_aa_dict = dict()

        if os.path.exists(prodigal_aa) and os.stat(prodigal_aa).st_size != 0:
            prodigal_aa_dict = gene_clusters.make_seqdict(prodigal_aa, format='fasta')

        ###---Fetch the proximal ORFs for each putative CRISPR array---###
        #Loop through the CRISPR array DataFrame and pull out Prodigal ORF entries
        #that fall within the window coordinate extent
        #cluster_seq_paths = []

        for index, row in crispr_array_df.iterrows():
            cluster_orfs = []
            seqid = re.escape(row['seqid'])
            pattern = r'%s_\d+$' % seqid
            #Are ORFs within range of the CRISPR array? If so, gather them.
            for orfid in prodigal_aa_dict.keys():
                if re.match(pattern, orfid):
                    orf_coord = [int(a.strip()) for a in prodigal_aa_dict[orfid].description.split('#')[1:3]]
                    #if orf_coord[0] >= int(row['start']) - opts.clst_window_extent or orf_coord[1] <= int(row['end']) + opts.clst_window_extent:
                    if orf_coord[0] >= int(row['start']) - opts.clst_window_extent and orf_coord[1] < int(row['start']):
                        cluster_orfs.append(prodigal_aa_dict[orfid])
                    elif orf_coord[0] > int(row['end']) and orf_coord[1] <= int(row['end']) + opts.clst_window_extent:
                        cluster_orfs.append(prodigal_aa_dict[orfid])
                    else:
                        pass
            print('CRISPR', row['ID'], row['start'], ':', \
                 row['end'], '----', len(cluster_orfs), 'proximal ORFs', \
                 'within', str(opts.clst_window_extent), 'bp')
            #write FASTA amino acid file of translated ORFs within window extent of CRISPR
            cluster_seqs = os.path.join(output_paths['Cluster'], '%s_%s_orfclust_%d.faa' % (row['seqid'], row['ID'], opts.clst_window_extent))
            #cluster_seq_paths.append(cluster_seqs)
            with open(cluster_seqs, 'w') as clustseq:
                print('writing CRISPR-proximal translated ORFs to %s' % cluster_seqs)
                SeqIO.write(cluster_orfs, cluster_seqs, 'fasta')
    else:
        print('CRISPRDetect GFF file has entries', \
              'CRISPR pandas df is empty')
else:
    print('No putative CRISPRs found with CRISPRDetect for', opts.fasta_file)
###---END Prodigal---###
#Compress the input if it was gzipped originally
if gzip:
    print('re-gzip compressing file %s' % opts.fasta_file)
    subprocess.run(['gzip', nt_fasta], shell=False)
#==============================================================================
###---BEGIN MacSyFinder---###
#(Sub)type the groups of CRISPR-neighboring genes
# macsyfinder_opts = {'--db-type':opts.macsy_dbtype,
# '--e-value-search':opts.macsy_eval, '--i-evalue-select':opts.macsy_eval, '--coverage-profile':opts.macsy_coverage,
# '--def':'/build/data/definitions/%s' % opts.ccs_typing,
# '--res-search-suffix':'hmmout', '--res-extract-suffix':'res_hmm_extract',
# '--profile-suffix':'.hmm', '--profile-dir':opts.profile_dir,
# '--worker':opts.threads}
#
# for csp in cluster_seq_paths:
#     macsyfinder_opts['--sequence-db'] = csp
#     macsyfinder_opts['--out-dir'] = os.path.join(output_paths['MacSyFinder'], '%s_%s' % (prefix, shell_tools.get_basename(csp)))
#     macsyfinder_cmd = shell_tools.exec_cmd_generate('macsyfinder', macsyfinder_opts)
#     #Search the CRISPR-Cas (sub)types systems of choice
#     macsyfinder_cmd.append(opts.macsy_systems)
#     #Run MacSyFinder
#     subprocess.run(macsyfinder_cmd, shell=False)
#     logger.info('%s with MacSyFinder performed for %s' % (opts.ccs_typing, csp))
###---END MacSyFinder---###
# #==============================================================================
