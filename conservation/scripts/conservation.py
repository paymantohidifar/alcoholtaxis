import matplotlib.pyplot as plt
from re import split
import os


def read_fasta_file(filename):
    """ Reads in file in fasta format and returns sequences of tuples containing
    header name and the sequence"""
    with open(filename, 'r') as file:
        return tuple((name, seq.replace('\n', ''))
                     for name, ignore, seq in [entry.partition('\n')
                                               for entry in file.read().split('>')[1:]])


def blastp_parser(filename):
    """ Reads in blastp output file (-outfmt 6 & 7 ) and returns content of table as sequences of
    keys and values, where keys holds the attributes in the table"""
    with open(filename, 'r') as file:
        table = file.read().split('#')
        attribs = table[4].replace(' Fields: ', '').rstrip('\n').split(', ')
        return tuple({attribs[n]: value for n, value in enumerate(hit.split('\t'))}
                     for hit in table[5].split('\n')[1:-1])


def skip_intro(fil):
    """Skip description lines at the beginning of files"""
    line = fil.readline()
    while line.startswith('#'):
        line = fil.readline()


def hmmscan_parser(filename):
    """ Reads in hmmscan output file (domtblout format) and returns a tuple of tuples
    each containing the query name and end position of the HAMP domain in that query"""
    with open(filename, 'r') as file:
        domains = tuple()
        skip_intro(file)
        line = file.readline()
        while line:
            if line.startswith('HAMP'):
                parts = split('\s+', line, 22)
                domains += ((parts[3], int(parts[18])+1),)
            line = file.readline()
    return domains


def extract_seqs(hits, start, end=None):
    """ Takes in the hits and the desired range (start, end) and update sequence field
    for all hits"""
    for n in range(0, len(hits)):
        hits[n]['subject seq'] = hits[n]['subject seq'][start:end]


def collapse_duplicates(hits, key):
    """ Collapses hits with 100% similar sequences into one hit based on hits order"""
    duplicates = set()
    return tuple(hit for hit in hits
                 if not (hit.get(key) in duplicates or duplicates.add(hit.get(key))))


def collect_loci(hits, nhelix_range, chelix_range):
    """ Takes in hits and updates N helix ocus and C helix locus sequence
    fields using loci ranges provided as arguments """
    for n in range(0, len(hits)):
        hits[n]['N helix'] = hits[n]['subject seq'][nhelix_range[0]:nhelix_range[1]]
        hits[n]['C helix'] = hits[n]['subject seq'][chelix_range[0]:chelix_range[1]][::-1]


def process_domains(pident):
    """ Runs hmmscan to obtain the HAMP domain and returns start position of
    signaling domain for each hit"""

    # run hmmcan to identify the best domain matches for each hit
    print('Running hmmscan to search each hit against Pfam 33.1 hmm profile database to determine HAMP domain')
    os.system('hmmscan --cut_ga --cpu 4 --domtblout {} {} {} > {}'
              .format(resultdir+'domains_'+pident+'.tab', pfamdb, resultdir+'_processed_hits_'+pident+'.fa',
                      logsdir+'hmmscan'+pident+'.log'))
    return hmmscan_parser(resultdir + 'domains_'+pident+'.tab')


def alignment_file_parser(filename):
    """ Reads alignment file (in fasta format) and returns sequence of header name
      and the corresponding sequence."""
    aln_hits = tuple()
    for entry in read_fasta_file(filename):
        hit = dict()
        hit['subject acc.'] = entry[0]
        hit['subject seq'] = entry[1]
        aln_hits += (hit,)
    return aln_hits


def process_hits():
    """ Processes blastp output """

    print('Parsing hits obtained from blastp search ...')

    # parse hit table obtained form blastp search
    hits = blastp_parser(datadir + 'hits.tab')
    print('Total {} of hits were parsed from blastp search output.'.format(len(hits)))
    hits = tuple(hit for hit in hits if '1386' not in hit.get('subject tax ids'))

    # remove gaps from sequences
    for n in range(0, len(hits)):
        hits[n]['subject seq'] = hits[n]['subject seq'].replace('-', '')

    # collect the hits that have at least 98% coverage of the query sequence (650 aa).
    hits = tuple(hit for hit in hits if len(hit['subject seq']) > 0.98*662)
    print('Total {} of hits remained after removing the hits that had less than 98% coverage [excluding gaps].'
          .format(len(hits)))

    # generate % identity histogram plot
    identity = [float(hit.get('% identity')) for hit in hits]
    fig, ax = plt.subplots()
    ax.hist(identity, bins=100)
    ax.set_xlabel("Percent identity")
    ax.set_ylabel("Count")
    plt.savefig(resultdir+'pident.eps', dpi=None, facecolor='w', edgecolor='w', format='eps', pad_inches=0.1)
    # plt.show()
    return hits


def get_domain_seq(hits, pident):
    """Write amino acid sequences signaling domains of the hits with certain percent identity
    into a FASTA file"""
    hits = tuple(hit for hit in hits if float(hit.get('% identity')) > pident)
    print('Total {} of hits with % identity greater than {}.'.format(len(hits), pident))

    # remove duplicate sequences
    hits = collapse_duplicates(hits, 'subject seq')
    print('Total {} of hits remained after removing duplicate sequences.'.format(len(hits)))

    # write processed hits in a fasta file used for domain identification
    with open(resultdir + '_processed_hits_'+str(pident)+'.fa', 'w') as file:
        file.writelines('>' + hit['subject acc.'] + '\n' + hit['subject seq'] + '\n' for hit in hits)
    print(100*'-')

    # identify domains for each hit
    domains = process_domains(str(pident))
    assert len(domains) == len(hits)

    # write sequences after HAMP domain (350-end) hits in a fasta file used for domain identification
    with open(resultdir + 'domain_seqs_'+str(pident)+'.fa', 'w') as file:
        file.writelines('>'+hit['subject acc.']+'\n'+hit['subject seq'][domains[n][1]:]+'\n'
                        for n, hit in enumerate(hits))
    print(100 * '-')


def align_domain_seqs(pident):
    """ Runs multiple-sequence alignments on processed hits"""
    print('Running Muscle for accurate MSA ...')
    os.system('muscle3.8.31_i86linux64 -in {} -out {} -log {} > /dev/null'
              .format(resultdir + 'domain_seqs_'+pident+'.fa', resultdir + 'domain_seqs_'+pident+'.afa',
                      logsdir + 'muscle_'+pident+'.logs'))
    print(100 * '-')


def process_aln_seqs(pident):
    """ Processes alignment sequences and write sequences of ethanol sensing loci into a file"""
    print('Collecting putative ethanol sensing loci')
    aln_seqs = alignment_file_parser(resultdir + 'domain_seqs_'+pident+'.afa')
    # collect amino acids residues within ethanol sensing region on N helix and C helix
    nhelix = max([aln_seq['subject seq'].find('QTSSDITKASIQST') for aln_seq in aln_seqs])
    chelix = max([aln_seq['subject seq'].find('MTNEIAGKLQTMNS') for aln_seq in aln_seqs])
    collect_loci(aln_seqs, (nhelix, nhelix + 14), (chelix, chelix + 14))

    with open(resultdir + 'nhelix_logo_'+pident+'.fa', 'w') as file1, \
            open(resultdir + 'chelix_logo_'+pident+'.fa', 'w') as file2:
        file1.writelines('>' + hit['subject acc.'] + '\n' + hit['N helix'] + '\n' for hit in aln_seqs)
        file2.writelines('>' + hit['subject acc.'] + '\n' + hit['C helix'] + '\n' for hit in aln_seqs)


def clear_intermediate_files():
    """ Remove all intermediate files generated during analysis"""
    os.system('rm -f {}_*'.format(resultdir))


def main():
    # process hits generated from blastp output
    hits = process_hits()

    # generate putative ethanol-sensing loci with all McpB homologs in Bacillus subtilis strains
    get_domain_seq(hits, 50)
    align_domain_seqs('50')
    process_aln_seqs('50')

    # generate putative ethanol-sensing loci with McpB-like receptors in Bacillus subtilis strains
    get_domain_seq(hits, 90)
    align_domain_seqs('90')
    process_aln_seqs('90')

    # remove intermediate files
    clear_intermediate_files()


if __name__ == '__main__':
    datadir = '../data/'
    resultdir = '../results/'
    pfamdb = '$HOME/hmmer-3.3/db/Pfam-A.hmm'
    logsdir = '../logs/'
    main()
