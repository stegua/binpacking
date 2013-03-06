#!/data/gualandi/LIBS/python/bin/python

import sys

def printAMPL(F):
    for i in range(36):
        F.readline()
        F.readline()
        line = F.readline().replace(" ","")
        out = open(str(line[1:].strip()), "w")
        out.write("data;\nparam m := 6;\nparam n := 18;\nparam k := 6;\n")
        out.write("param A : 1 2 3 4 5 6 :=\n")
        for l in range(7):
            F.readline()
        for l in range(18):
            line = F.readline().split()
            out.write("\t")
            for n in line:
                out.write(n)
                out.write(" ")
            out.write("\n")
        out.write(";\n")
        for l in range(4):
            F.readline()
        line = F.readline().split()
        out.write("param b := \n")
        for i,n in enumerate(line[1:]):
            out.write("\t"+str(i+1)+" "+n+"\n")
        for l in range(7):
            F.readline()
        out.write(";\n")
        out.close()

# MAIN PART
F = open(sys.argv[1], 'r')
printAMPL(F)
