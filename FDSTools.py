# -*- coding: utf-8 -*-
"""
Created on Wed Oct 08 10:00:58 2014

@author: tstopi
"""

import numpy as np
import sys
import os
import struct


# Read fortran record, return payload as bytestring
def readFRec(infile,fmt):
    len1   = infile.read(4)
    if not len1:
        return None
    len1   = struct.unpack('@I', len1)[0]

    if len1==0:
        len2=infile.read(4)
        return None
    num    = int(len1/struct.calcsize(fmt))
    fmt2   = str(num)+fmt
    if num>1:
        result = struct.unpack(fmt2,infile.read(len1))
    else:
        result = struct.unpack(fmt2,infile.read(len1))[0]
    len2   = struct.unpack('@I', infile.read(4))[0]
    if fmt == 's':
        result=result[0].rstrip()
    return result


def read_bndf(fname,max_time=np.Inf):
    fin = open(fname,'rb')
    quantity = readFRec(fin,'s')
    short_name = readFRec(fin,'s')
    units = readFRec(fin,'s')
    n_patch = readFRec(fin,'I')
    patch_extents = []
    for n in range(0,n_patch):
        temp = np.array(readFRec(fin,'I'))
        patch_extents.append(temp)

    Q=[]
    T  = []
    while True:
        Time  = readFRec(fin,'f')
        if Time == None or Time>max_time:
            break
        T.append(Time)
        Q.append([])
        for n in range(0,n_patch):
            temp = np.array(readFRec(fin,'f'))
            Q[-1].append(temp)
    fin.close()
    return (quantity, np.array(T),Q, np.hstack(patch_extents))



#Assumes single precision
def read_prt(fname,max_time=np.Inf):
    fin = open(fname,'rb')
    one_integer=readFRec(fin,'I')
    version=readFRec(fin,'I')
    n_part=readFRec(fin,'I')
    q_labels = []
    q_units  = []
    for npc in range(0,n_part):
        n_quant,zero_int=readFRec(fin,'I')
        for nq in range(0,n_quant):
            smv_label=readFRec(fin,'s')
            units    =readFRec(fin,'s')
            q_units.append(units)
            q_labels.append(smv_label)
    Q=[]
    T  = []
    diam =[]
    while True:

        Time  = readFRec(fin,'f')
        if Time == None or Time>max_time:
            break
        nplim = readFRec(fin,'I')
        if nplim>0:
            xyz  = np.array(readFRec(fin,'f'))
            tag  = np.array(readFRec(fin,'I'))
            q    = np.array(readFRec(fin,'f'))
            xyz.shape = (3,nplim)
            q.shape   = (n_quant,nplim)
            # process timestep data
            T.append(Time)
            Q.append(q)
        else:
            tmp = fin.read(24)
    fin.close()
    return (np.array(T),np.hstack(Q),q_labels,q_units)

def parse_smv(fname):
    lines = open(fname).readlines()

    def find_keyword(lines,keyword):
        return [n for n,line in enumerate(lines) if line.startswith(keyword)]

    # findgrid definitions
    igrid = find_keyword(lines,"GRID")
    ipdim = find_keyword(lines,"PDIM")
    itrnx = find_keyword(lines,"TRNX")
    itrny = find_keyword(lines,"TRNY")
    itrnz = find_keyword(lines,"TRNZ")

    grids = []
    for n,i in enumerate(igrid):
        ibar,jbar,kbar,_ =map(int,lines[igrid[n]+1].split())
        xmin, xmax, ymin, ymax, zmin, zmax, r, g, b = map(float,lines[ipdim[n]+1].split())

        trnx = [float(num.split()[1]) for num in lines[itrnx[n]+2:itrnx[n]+ibar+2]]
        trny = [float(num.split()[1]) for num in lines[itrny[n]+2:itrny[n]+jbar+2]]
        trnz = [float(num.split()[1]) for num in lines[itrnz[n]+2:itrnz[n]+kbar+2]]
        grid = {'ibar': ibar,
                'jbar': jbar,
                'kbar': kbar,
                'xmin': xmin,
                'xmax': xmax,
                'ymax': ymax,
                'ymin': ymin,
                'zmin': zmin,
                'zmax': zmax,
                'trnx': trnx,
                'trny': trny,
                'trnz': trnz
                }
        grids.append(grid)
        
    return grids

if __name__ == "__main__":
    if len(sys.argv)<2:
        sys.exit("Give a file as argument")
    if sys.argv[1].endswith(".prt5"):
        T,Q,labels,units =  read_prt(sys.argv[1])
        print(labels,units)
        print(len(T),Q.shape)
    elif sys.argv[1].endswith(".bf"):    
        quantity, T, Q, patches = read_bndf(sys.argv[1])
        print(quantity)
        print(patches)
    elif sys.argv[1].endswith(".smv"):    
        grids = parse_smv(sys.argv[1])
        print(grids)