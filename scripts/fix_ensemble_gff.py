#!/usr/bin/env python

import argparse
import re
import sys

'''
Usage:
./fix_ensemble_gff.py [broken.gff3] > [fixed_name.gff3]
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('out',  nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    
    ens_source_pattern = re.compile(r'\[Source\:.*?\]')
    ens_name_pattern = re.compile(r'Name=ENSFHE.*?\b(.*?) \;Target=')
    ens_target_pattern = re.compile(r'Target=ENSFHE.*?\b(.*?) [\d]*? [\d]*? [\+\-]\;database=')

    fun_name_pattern = re.compile(r'Name=Funhe2EK.*?\b(.*?\;)\;Target=')
    fun_target_pattern = re.compile(r'Target=Funhe2EK.*?\b(.*?) [\d]*? [\d]*? [\+\-]\;database=')

    for line in args.gff3:
        if ens_source_pattern.search(line) is not None:
            tokens = re.split(ens_source_pattern, line)
            line = ''.join(tokens)
        
        match = ens_name_pattern.search(line)
        if match is not None:
            line = line[:match.start(1)] + line[match.end(1):]
        match = ens_target_pattern.search(line)
        if match is not None:
            line = line[:match.start(1)] + line[match.end(1):]

        match = fun_name_pattern.search(line)
        if match is not None:
            line = line[:match.start(1)] + line[match.end(1):]
        match = fun_target_pattern.search(line)
        if match is not None:
            line = line[:match.start(1)] + line[match.end(1):]


        print(line, file=args.out)


if __name__ == '__main__':
    main()