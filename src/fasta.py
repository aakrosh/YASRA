#!/usr/bin/env python

"""
Read and write fasta files
"""

__author__ = "Aakrosh Ratan"
__email__  = "ratan@bx.psu.edu"

from sys import argv, exit

class fastarecord:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    @staticmethod
    def prettyprintdna(seq, size):
        str = ""
        for i in range(0, len(seq), size):
            str += seq[i:i+size]
            str += "\n"
        return str[:-1]
       
    def __str__(self):
        str =  ">%s\n" % self.header
        str += self.prettyprintdna(self.sequence, 60)
        return str

    def __len__(self):
        return len(self.sequence)

class fastafile:
    def __init__(self, filename):
        self.filename = filename
        self.file     = open(self.filename, "r")
        self.cache    = None
                
    def __iter__(self):
        return self

    def next(self):
        if self.cache != None: 
            line = self.cache
        else:                  
            line = self.file.readline()

        # are we at the end of the file?
        if not line:
            self.file.close()
            raise StopIteration

        assert line[0] == ">", "the first letter on the header should be >"
        header = line.strip()[1:]
        sequence = ""
       
        line = self.file.readline()
        while line[0] != ">":
            sequence += line.strip()    

            if line[0].isdigit():
                # if this is a fasta file of quality values, we need an extra
                # space here
                sequence += " "

            line = self.file.readline()
            if not line:
                break

        self.cache = line
        return fastarecord(header, sequence)

    def close(self):
        self.file.close()
    
if __name__ == "__main__":
    if len(argv) != 2:
        print >> stderr, "please input a test fasta file"
        exit(2)

    for entry in fastafile(argv[1]):
        print entry

    fastafile.close()
