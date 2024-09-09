import argparse
import os

def get_arguments():
    parser=argparse.ArgumentParser()
    parser.add_argument("--lib_key",required=True,type=str,help="Targets key")
    parser.add_argument("--tumor_list",required=True,type=str,help="File containing tumor list")
    args=parser.parse_args()
    return args

def main(argv=None):
    argv=get_arguments()
    with open(argv.tumor_list,'r') as file:
        TUMORS=file.read().splitlines()

    #purecn
    with open('purecn_purity.tsv','w') as outfile0,open('purecn_ploidy.tsv','w') as outfile1:
        for tumor in TUMORS:
            file_check=f'data/work/{argv.lib_key}/{tumor}/purecn/{tumor}.csv'
            if os.path.isfile(file_check):
                with open(file_check,'r') as file:
                    line=file.readline()#header
                    line=file.readline()
                    purity,ploidy=line.replace('"','').rstrip().split(',')[1:3]
                    purity=f'{float(purity):.2f}'
                    ploidy=f'{float(ploidy):.2f}'
                outfile0.write(f'{tumor}\t{purity}\n')
                outfile1.write(f'{tumor}\t{ploidy}\n')
            else:
                outfile0.write(f'{tumor}\t1.0\n')
                outfile1.write(f'{tumor}\tNA\n')

    #facets
    with open('facets_purity.tsv','w') as outfile0,open('facets_ploidy.tsv','w') as outfile1:
        for tumor in TUMORS:
            file_check=f'data/work/{argv.lib_key}/{tumor}/facets/purity_ploidy.csv'
            if os.path.isfile(file_check):
                with open(file_check,'r') as file:
                    line=file.readline()#header
                    line=file.readline()
                    purity,ploidy=line.replace('"','').rstrip().split(',')
                    if purity!="NA":
                        purity=f'{float(purity):.2f}'
                    else:
                        purity="1.0"
                    ploidy=f'{float(ploidy):.2f}'
                outfile0.write(f'{tumor}\t{purity}\n')
                outfile1.write(f'{tumor}\t{ploidy}\n')
            else:
                outfile0.write(f'{tumor}\t1.0\n')
                outfile1.write(f'{tumor}\tNA\n')

if __name__=="__main__":
    main()
