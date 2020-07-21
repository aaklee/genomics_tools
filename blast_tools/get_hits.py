#!/usr/bin/python

import sys
import os
import argparse
import pandas
from Bio import SeqIO

def generate_contig_hits(blastresults, evalf, lenf, total_num, pull):
    
    hits = pandas.read_csv(blastresults, sep=',')
    hits.columns = ['qseqid', 'sseqid', 'pident', 'qlen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    hits = hits.sort_values(by=['qseqid', 'evalue']) # sort by evalue first, then qseqid to group
    
    dup_hits = []
    to_retrieve = {}

    num_ct = 0
    curr_qseqid = ''
    for hit in hits.itertuples():
        evalue = float(hit.evalue)
        sseqid = str(hit.sseqid)
        qseqid = str(hit.qseqid)
        sstart = int(hit.sstart)
        send = int(hit.send)
        qlen = int(hit.qlen)
        
        # only include a certain number of blast hits per query sequence
        if num_ct == 0:
            curr_qseqid = qseqid
            
        if total_num and num_ct == total_num:
            if curr_qseqid == qseqid:
                continue
            elif curr_qseqid != qseqid:
                num_ct = 0
                pass
        
        # evalue filter
        if evalf:
            if evalue > evalf:
                pass
            else:
                continue
            
        # add hit to to_retrieve dictionary
        if sseqid not in to_retrieve.keys():
            to_retrieve[sseqid] = []
        else:
            dup_hits.append(sseqid)
            
        # get positions of hit sequence
        if send > sstart:
            send = send + pull
            sstart = sstart - pull
            if sstart < 0:
                sstart = 0
            
            # do not include if does not meet length requirement
            if lenf:
                if (send - sstart) > lenf:
                    pass
                else:
                    continue

            to_retrieve[sseqid].append((sstart, send, qseqid))
            num_ct += 1
            
        elif send < sstart:
            send = send - pull
            sstart = sstart + pull
            if send < 0:
                send = 0

            # do not include if does not meet length requirement
            if lenf:
                if (sstart - send) > lenf:
                    pass
                else:
                    continue

            to_retrieve[sseqid].append((send, sstart, qseqid))
            num_ct += 1
        

    with open('contig_hits.txt', 'w+') as outf:
        for hit in to_retrieve.keys():
            for pos in to_retrieve[hit]:
                direction = ''
                if pos[0] > pos[1]:
                    direction = 'minus'
                    outf.write('%s %i-%i %s\n' % (hit, pos[0], pos[1], direction))
                if pos[1] > pos[0]:
                    direction = 'plus'
                    outf.write('%s %i-%i %s\n' % (hit, pos[0], pos[1], direction))
                    
    return to_retrieve
    
    
    
    
def append_query_sequences(to_retrieve, csv):
    fasta = []
    try:
        with open('%s_hits.out' % csv, 'r') as inf:
            for record in SeqIO.parse(inf, 'fasta'):
                new_desc = [i[2] for i in to_retrieve[record.id]]
            
                record.id = record.description.replace(' ', '_')
                record.description = ' '.join(new_desc)
                
                fasta.append(record)

    except FileNotFoundError as e:
        print("Is blast+ in the PATH?")
   
   
    ### HITS.FA STICKS ON THE BLAST QUERY SEQUENCE THAT RECOVERED THE BLASTDB SEQUENCE ###
    outfn = '%s_hits.fa' % csv
    with open(outfn, 'w+') as outf:
        SeqIO.write(fasta, outf, 'fasta')

    print('blast hits in %s_hits.fa' % csv)
    print('done!')
    print()






def main():
    parser = argparse.ArgumentParser(description='Use blastdbcmd to batch recover BLAST results from a provided BLASTable database.')
    parser.add_argument('-d', '--database', help='path to BLASTable database')
    parser.add_argument('-r', '--results', help='BLAST results in csv format')
    parser.add_argument('-p', '--pull', default=500, help='bps to pull up/downstream of the BLAST result (default=500 bp)')
    parser.add_argument('-l', '--len', help='len filter')
    parser.add_argument('-e', '--eval', help='eval filter (float)')
    parser.add_argument('-n', '--total_num', help='number of blast hits to get')

    # parse arguments
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    blastdb = args.database
    blastresults = os.path.abspath(args.results)
    pull = int(args.pull)
    csv = blastresults.strip().split('/')[-1]

    # parse filtration parameters
    if args.len:
        lenf = int(args.len)
    else:
        lenf = False

    if args.eval:
        evalf = float(args.eval)
    else:
        evalf = False

    if args.total_num:
        total_num = int(args.total_num)
    else:
        total_num = False

    # generate contig_hits.txt intermediate file with information needed to retrieve hits from blastdb
    to_retrieve = generate_contig_hits(blastresults, evalf, lenf, total_num, pull)

    os.system('blastdbcmd -db %s -entry_batch contig_hits.txt -outfmt %s -out %s_hits.out' % (blastdb, '%f', csv))

    # append query sequences for each hit to produce final hits file
    append_query_sequences(to_retrieve, csv)





if __name__ == '__main__':
    main()
