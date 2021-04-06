#!/usr/bin/python

import sys
import os
import argparse
import pandas
from Bio import SeqIO


def recursive_combine(hits, index=0):
    if hits and len(hits) > 1:
        for i in range(index, len(hits)-1):  
            if hits[i][1] > hits[i+1][0]:
                hits[i] = (hits[i][0], hits[i+1][1], hits[i][2], hits[i][3])
                del hits[i+1]
                #print('RECURSE')
                return recursive_combine(hits, i)
    return hits


def generate_contig_hits(separate_hits, evalf, lenf, total_num, pull):
    # find hits corresponding to one query sequence
    dup_hits = []
    to_retrieve = {}
    retrieve_dict = {}
    retrieve_list = []

    for hit in separate_hits.itertuples():
        evalue = float(hit.evalue)
        sseqid = str(hit.sseqid)
        qseqid = str(hit.qseqid)
        sstart = int(hit.sstart)
        send = int(hit.send)
        qlen = int(hit.qlen)

        # e-value filter
        if evalf and evalue > evalf:
            continue

        # check whether fwd or rev hit
        fwd = True
        if send < sstart:
            fwd = False

        # add positions of hit sequence
        if fwd:
            send = send + pull
            sstart = sstart - pull
            if sstart < 0:
                sstart = 0
            
            if lenf and (send - sstart) < lenf:
                continue

            if sseqid not in to_retrieve.keys():
                to_retrieve[sseqid] = []
            else:
                dup_hits.append(sseqid)

            to_retrieve[sseqid].append((sstart, send, qseqid, 'plus'))

        elif not fwd:
            send = send - pull
            sstart = sstart + pull
            if send < 0:
                send = 0

            if lenf and (sstart - send) < lenf:
                continue

            if sseqid not in to_retrieve.keys():
                to_retrieve[sseqid] = []
            else:
                dup_hits.append(sseqid)

            to_retrieve[sseqid].append((send, sstart, qseqid, 'minus'))



    # combine hits by coordinates if they overlap on the same scaffold/chromosome
    new_hits = []
    for hit in to_retrieve.keys():
        for direction in ['plus', 'minus']:
            hits = [i for i in to_retrieve[hit] if direction in i[3]]
            hits = sorted(hits, key = lambda x: x[0])
            new_hits.extend(recursive_combine(hits))
        to_retrieve[hit] = new_hits

    for hit in to_retrieve.keys():
        for pos in to_retrieve[hit]:
            sub_sseqid = '{}_{}-{}'.format(hit, pos[0], pos[1])
            if 'minus' in pos[3]:
                sub_sseqid = '{}_{}-{}'.format(hit, pos[1], pos[0])

            if sub_sseqid not in retrieve_dict.keys():
                retrieve_dict[sub_sseqid] = []

            retrieve_dict[sub_sseqid].append(pos[2])
            retrieve_list.append(sub_sseqid)


    with open('contig_hits.txt', 'w+') as outf:
        ct = 0
        for hit in to_retrieve.keys():
            for pos in to_retrieve[hit]:
                #if 'minus' in pos[3]:
                #    outf.write('{} {}-{} {}\n'.format(hit, pos[1], pos[0], pos[3]))
                #else:
                #    outf.write('{} {}-{} {}\n'.format(hit, pos[0], pos[1], pos[3]))
                outf.write('{} {}-{} {}\n'.format(hit, pos[0], pos[1], pos[3]))
                
                ct += 1
                if total_num and ct > total_num:
                    return to_retrieve, retrieve_dict, retrieve_list
                    
    return to_retrieve, retrieve_dict, retrieve_list



def append_query_sequences(to_retrieve, retrieve_list, retrieve_dict, outputf):
    directions = []
    for hit in to_retrieve.keys():
        for pos in to_retrieve[hit]:
            directions.append(pos[3])

    fasta = []
    ct = 0
    with open(outputf, 'r') as inf:
        for record in SeqIO.parse(inf, 'fasta'):
            new_desc = retrieve_dict[retrieve_list[ct]]
        
            record.id = '{}-'.format(ct) + record.description.replace(' ', '_')
            record.description = ' '.join(new_desc)
            
            fasta.append(record)
            ct += 1

    with open(outputf, 'w+') as outf:
        SeqIO.write(fasta, outf, 'fasta')

    print('{} blast hits in {}'.format(ct, outputf))



def separate_hits(blastresults, evalf, lenf, total_num, pull, csv, blastdb):
    hits = pandas.read_csv(blastresults, sep=',')
    hits.columns = ['qseqid', 'sseqid', 'pident', 'qlen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    hits = hits.sort_values(by=['qseqid', 'evalue']) # sort by evalue first, then qseqid to group

    for i in hits.qseqid.unique():
        # create input file for blastdbcmd retrieval
        to_retrieve, retrieve_dict, retrieve_list = generate_contig_hits(hits.loc[hits['qseqid'] == i], evalf, lenf, total_num, pull)
        
        # blastdbcmd retrieval of blast hits
        outputf = '{}-{}-hits.fa'.format(i, csv)
        os.system('blastdbcmd -db {} -entry_batch contig_hits.txt -outfmt %f -out {}'.format(blastdb, outputf))

        # append query sequence name to fasta headers
        append_query_sequences(to_retrieve, retrieve_list, retrieve_dict, outputf)




def main():
    parser = argparse.ArgumentParser(description='Use blastdbcmd to batch recover BLAST results from a provided BLASTable database.')
    parser.add_argument('-d', '--database', help='path to BLASTable database')
    parser.add_argument('-r', '--results', help='BLAST results in csv format')
    parser.add_argument('-p', '--pull', default=500, help='bps to pull up/downstream of the BLAST result (default=500 bp)')
    parser.add_argument('-l', '--len', help='len filter')
    parser.add_argument('-e', '--eval', help='e-val filter (float)')
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

    
    separate_hits(blastresults, evalf, lenf, total_num, pull, csv, blastdb)
    


if __name__ == '__main__':
    main()
