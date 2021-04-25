#!/usr/bin/python

"""


DESCRIPTION:
BLAST search a list or multiple lists of databases with a list of query sequences.



INPUTS (REQUIRED):
    -dbl        db1,db2,...     List of the paths to BLASTable databases (use makeblastdb)
    -query      query.fasta     FASTA formatted list of query sequences
    -csv        t/f             default: false
                                Format output as a CSV file (probably should have used -parse_seqids
                                when running makeblastdb)
    -search     blastn,tblastx,blastx


"""

import sys
import os
import argparse


def blast(search, database, query, outfmt):
    q = query.split('/')[-1]
    db = database.split('/')[-1]
    if outfmt:
        print('{} -db {} -query {} -out {}.{}.{}.csv -outfmt \"10 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore\"'.format(search, database, query, db, search, q))
        os.system('{} -db {} -query {} -out {}.{}.{}.csv -outfmt \"10 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore\"'.format(search, database, query, db, search, q))
    else:
        print('{} -db {} -query {} -out {}.{}.{}.out'.format(search, database, query, db, search, q))
        os.system('{} -db {} -query {} -out {}.{}.{}.out'.format(search, database, query, db, search, q))



def blast_db(database, query, search, outfmt):
    print(search, database.split('/')[-1], query.split('/')[-1])
    blast(search, database, query, outfmt)


def blast_dbl(database_list, query, search, outfmt):
    for i in database_list:
        blast_db(i, query, search, outfmt)


def blast_dbf(database_file, query, search, outfmt):
    with open(database_file, 'r') as inf:
        for line in inf:
            blast_db(line.strip(), query, search, outfmt)



def main():
    parser = argparse.ArgumentParser(description='BLAST search a list or multiple lists of databases with a list of query sequences.')
    parser.add_argument('-d', '--database', help='path to BLASTable database, do not use with -dbl/-dbf')
    parser.add_argument('-l', '--database_list', help='comma-separated list of paths to BLASTable databases, do not use with -db/-dbf')
    parser.add_argument('-f', '--database_file', help='file containing list of paths to BLASTable databases, do not use with -db/-dbl')
    parser.add_argument('-q', '--query', help='FASTA formatted list of query sequences', required=True)
    parser.add_argument('-s', '--search', help='blastn,blastx,tblastx', required=True)
    parser.add_argument('-c', '--csv', action='store_true', help='format output as a CSV file; -parse_seqids in makeblastdb!!')


    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit()

    # parse arguments
    args = parser.parse_args()

    query = os.path.abspath(args.query)
    search = args.search

    if args.csv:
        outfmt = True
    else:
        outfmt = False

    if args.database:
        database = args.database
        blast_db(database, query, search, outfmt)
    elif args.database_list:
        database_list = args.database_list.split(',')
        blast_dbl(database_list, query, search, outfmt)
    elif args.database_file:
        database_file = os.path.abspath(args.database_file)
        blast_dbf(database_file, query, search, outfmt)



if __name__ == '__main__':
    main()
