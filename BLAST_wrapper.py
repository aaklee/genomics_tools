#!/usr/bin/python

"""

DESCRIPTION
    use this to complete all the things you need to do for a BLAST search.
    runs your BLAST query (blastn), retrieves your hits (blastdbcmd), and 
    filters them for uniqueness (cd-hit-est).
    
OPTIONS
    REQUIRED:
    -db     path_to_blastdb     path to BLASTable database
    -q      query.fa            name of query file (or path to file)
    -o      BLASTout.csv        name of BLAST output file
    -kmers  comma-separated-ks  e.g. 21,27,33,51,55,63

    OPTIONAL:
    -pull   #bps                number of base pairs to pull up/downstream your
                                BLAST hit, default is 500
    -c      0.8-1.0             sequence identity threshold for cd-hit-est,
                                default is 0.8 (least sensitive)
    -n      3-10                word size for cd-hit-est, default is 4
    
    
"""

import sys, os

def main():

    if len(sys.argv) <= 1:
        print(__doc__)
        exit()

    db, q, o, pull, c, n, kmers = '', '', '', 500, 0.8, 4, []

    print(sys.argv)

    for i in range(1, len(sys.argv)):
        if '-db' in sys.argv[i].lower():
            db = os.path.abspath(sys.argv[i+1])
        elif '-q' in sys.argv[i].lower():
            q= os.path.abspath(sys.argv[i+1])
        elif '-o' in sys.argv[i].lower():
            o = sys.argv[i+1]
        elif '-pull' in sys.argv[i].lower():
            pull = int(sys.argv[i+1])
        elif '-c' in sys.argv[i].lower():
            c = float(sys.argv[i+1])
        elif '-n' in sys.argv[i].lower():
            n = int(sys.argv[i+1])
        elif '-kmers' in sys.argv[i].lower():
            kmers = sys.argv[i+1].split(',')

    print(kmers)
    try:
        for k in range(0,len(kmers)+1):
            os.system('blastn -outfmt "10 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -db {} -out {} -query {}'.format(db, o, q))
            os.system('python /projects/clement-lab/resources/blast_tools/get_hits.py -db {} -res {} -pull {}'.format(db, o, pull))
            of = o.strip('.csv')
            os.system('/home/hpc/leea33/github/cdhit/cd-hit-est -c 0.8 -n 4 -i {}_hits.fa -o {}_cd-hit-est.fa'.format(of, of))
            db = db.replace(kmers[k-1],'{}'.format(kmers[k]))
            o = o.replace(kmers[k-1],'{}'.format(kmers[k]))
            print(db,o,k)
    except FileNotFoundError as e:
        print('Did you "module add blast+"?')
    except IndexError as e:
        pass

    print('\n\n~ ~ ~ ~ ~ DONE! ~ ~ ~ ~ ~\n\n')
    
if __name__ == '__main__':
    main()
                                
