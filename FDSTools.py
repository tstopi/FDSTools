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


def read_bndf(fname,max_time=np.Inf, patches_only=False,statistics=None):
    fin = open(fname,'rb')
    quantity = readFRec(fin,'s')
    short_name = readFRec(fin,'s')
    units = readFRec(fin,'s')
    n_patch = readFRec(fin,'i')
    patch_extents = []
    for n in range(0,n_patch):
        temp = readFRec(fin,'i')
        patch_extents.append(temp)

    if patches_only:
        return patch_extents
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
    if statistics is not None:
        if statistics='max':
            Q = np.amax(Q,axis=0)
        elif statistics='min':
            Q = np.amin(Q,axis=0)
        elif statistics='mean':
            Q = np.mean(Q,axis=0)
    return quantity, T,Q, patch_extents



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

    # find grid definitions
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


    # find output files
    ibndf = find_keyword(lines,"BNDF")
    islcf = find_keyword(lines,"SLCF")

    bndfs = []
    slcfs = []

    for n,ind in enumerate(islcf):
        # SLCF     1 # STRUCTURED &     0    46     0   109     7     7
        slcf = {}
        nm, imin, imax, jmin, jmax, kmin, kmax = [int(num) for num in lines[ind].split() if num.isnumeric()] 
        slcf['fname'] = lines[ind+1].strip()
        slcf['quantity'] = lines[ind+2].strip()
        slcf['short_name'] = lines[ind+3].strip()
        slcf['units'] = lines[ind+4].strip()
        slcf['nm'] = nm
        slcf['imin'] = imin
        slcf['imax'] = imax 
        slcf['jmin'] = jmin
        slcf['jmax'] = jmax
        slcf['kmin'] = kmin
        slcf['kmax'] = kmax
        slcfs.append(slcf)

    for n,ind in enumerate(ibndf):
        # BNDF     1     1
        bndf = {}
        nm, _ = [int(num) for num in lines[ind].split() if num.isnumeric()] 
        bndf['fname'] = lines[ind+1].strip()
        bndf['quantity'] = lines[ind+2].strip()
        bndf['short_name'] = lines[ind+3].strip()
        bndf['units'] = lines[ind+4].strip()
        bndf['nm'] = nm
        bndf['patches'] = read_bndf(bndf['fname'],patches_only=True)
        bndfs.append(bndf)

    return grids,slcfs,bndfs


    

def get_data_from_slcfs():
    yield

def proc_commandline():
    import argparse
    parser = argparse.ArgumentParser(description='FDS Tools')
    parser.add_argument('fname', help='Smokeview filename')
    parser.add_argument('-o', help='Output filename')
    parser.add_argument('-t','--type',help='Type of query [slcf,bndf,prt5]')
    parser.add_argument('--pbx',help='X-value',type=float)
    parser.add_argument('--pby',help='Y-value',type=float)
    parser.add_argument('--pbz',help='Z-value',type=float)
    parser.add_argument('-s','--statistics',help='Type of statistics [max,min,mean] default: None')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args =proc_commandline()
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
        grids,slcfs,bndfs = parse_smv(sys.argv[1])
        print(grids,slcfs,bndfs)