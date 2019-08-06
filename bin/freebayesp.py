#!/usr/bin/env python
# #-*- coding:utf-8 -*-
import sys
import tempfile
import os
import shutil
from subprocess import Popen, PIPE
import atexit

sys.setrecursionlimit(20000)

def quicksort(array, index):#快速排序
    if len(array) < 2:
        return array
    else:
        less = []
        greater = []
        pivot = array[0]
        for i in array[1:]:
            if i[index] <= pivot[index]:
                less.append(i)
        for i in array[1:]:
            if i[index] > pivot[index]:
                greater.append(i)
        return quicksort(less, index) + [pivot] + quicksort(greater, index)
def overlap(array):#检查是否有重叠
    for i in array:
        if int(i[1]) >= int(i[2]):
            return i[0],i[1] + ' >= ' + i[2]
    if len(array) < 2:
        return 'n'
    a = array[0]
    less = array[1:]
    for i in less:
        #print a[0],a[1],a[2],i[0],i[1],i[2]
        if a[0] == i[0]:
            if int(a[2]) > int(i[2]):
                if int(a[1]) < int(i[2]):
                    return 'y'
            if int(a[2]) < int(i[2]):
                if int(a[2]) > int(i[1]):
                    return 'y'
            if int(a[2]) == int(i[2]):
                return 'y'
    return overlap(less)
                    
def unit(array): #重叠区域合并
    if overlap(array) == 'n':
        return array
    a = array[0]
    less = array[1:]
    c = 0
    for index, b in enumerate(less):
        #print a[0],a[1],a[2],b[0],b[1],b[2]
        if a[0] == b[0]:
            if int(a[2]) > int(b[2]): #a(5 10) b(4 6)
                if int(a[1]) < int(b[2]):
                    less.append([a[0], b[1], a[2], int(a[2]) - int(b[1])])
                    del less[index]
                    c += 1
            if int(a[2]) < int(b[2]):#a(4 6) b(5 10)
                if int(a[2]) > int(b[1]):
                    less.append([a[0], a[1], b[2], int(b[2]) - int(a[1])])
                    del less[index]
                    c += 1
            if int(a[2]) == int(b[2]):
                if int(a[1]) < int(b[1]): #a(4 10) b(5 10)
                    less.append([a[0], a[1], b[2], int(b[2]) - int(a[1])])
                    del less[index]
                    c += 1
                if int(a[1]) > int(b[1]): #a(6 10) b(5 10)
                    less.append([a[0], b[1], b[2], int(b[2]) - int(b[1])])
                    del less[index]
                    c += 1
                if int(a[1]) == int(b[1]): #a(5 10) b(5 10)
                    less.append([a[0], b[1], b[2], int(b[2]) - int(b[1])])
                    del less[index]
                    c += 1
    if c == 0:
        less.append(a)
    return unit(less)
def cutSplit(array):
    length = 0
    for i in array:
        length += i[3]
    avg = int(length/len(array))
    #avg = 500
    out = []
    for index, i in enumerate(array):
        if i[3] <= avg:
            out.append(i)
        if i[3] > avg:
            a, b = i[3]//avg, i[3]%avg
            for m in range(0, a + 1):
                if m == a and b != 0:
                    c = [i[0], str(int(i[1]) + avg*m), i[2], int(i[2]) - (int(i[1]) + avg*m)]
                    out.append(c)
                elif m != a:
                    c = [i[0], str(int(i[1]) + avg*m), str(int(i[1]) + avg*(m+1)), (int(i[1]) + avg*(m+1)) - (int(i[1]) + avg*m)]
                    out.append(c)
    return out
    

if len(sys.argv) < 3:
    sys.exit('Usage: ' + sys.argv[0] +' threadNumber freebayes parameters')
p = int(sys.argv[1])
bed = ""
bedindex = 0
for i in range(len(sys.argv)):
    if sys.argv[i] == "-t":
        bed = sys.argv[i+1]
        bedindex = i+1
if bed == "":
    sys.exit('This script only work for target based parallel!')
bedfile = open(bed).readlines()
bedline = len(bedfile)
cmd = []
r = []
f = []
bedpath = []
bedsize = []
dirpath = tempfile.mkdtemp()

########ytx
raw = []
for i in bedfile:
    l = i.strip().split("\t")
    l.append(int(l[2])-int(l[1]))
    raw.append(l)

rawUnit = unit(raw)  #将 重叠区域合并
rawSplit = cutSplit(rawUnit) #将 较长区域切分
rawSort = quicksort(rawSplit, 3) #快速排序

pList = range(p)
for i in range(len(rawSort)/p + 1):
    pList = pList + range(p)
pList = pList[:len(rawSort)]

listP = []
for i in range(p):
    bedpath.append(os.path.join(dirpath,str(i) + ".bed"))
    f.append(open(os.path.join(dirpath,str(i) + ".bed"),"w"))
    f[i].write("chrM\t1\t2\n") # hack for first line loss bug.
    listP.append([])
    
for l, i in zip(rawSort, pList):
    x = l[0].split('chr')[1]
    try:
        x = int(x)
    except:
        x = 0
    l.append(x)
    listP[i].append(l)
    
for i in range(p):
    for l in quicksort(quicksort(listP[i], 1), 4):
        f[i].write(l[0] + '\t' + l[1] + '\t' + l[2] + "\n")
########
'''
for i in range(p):
    bedsize.append(0)
    bedpath.append(os.path.join(dirpath,str(i) + ".bed"))
    f.append(open(os.path.join(dirpath,str(i) + ".bed"),"w"))
    f[i].write("chrM\t1\t2\n") # hack for first line loss bug.
    if i == 0:
        r.append(range(0, bedline / p))
    elif i == p - 1:
        r.append(range(bedline * i/p,bedline))
    else:
        r.append(range(bedline * i/p,bedline * (i+1)/p))
for i,v in enumerate(bedfile):
    for index,j in enumerate(r):
        if i in j:
            bedsize[index] += int(v.strip().split("\t")[2]) - int(v.strip().split("\t")[1])
            f[index].write(v)
'''
            
            
            
for i in f:
    i.close()
    
    
    
    
for i in range(p):
    sys.argv[bedindex] = bedpath[i]
    cmd.append(sys.argv[2:])
out = []
outfile = []
pids=set()
vcftmp = []


def clean():
    shutil.rmtree(dirpath)
atexit.register(clean)


for i in range(p):
    f = open(os.path.join(dirpath,str(i) + ".vcf"),'w+')
    vcftmp.append(f)
    process = Popen(' '.join(cmd[i]),stdout=f, shell=True)
    out.append(process)


bedsub = False
for index,i in enumerate(out):
    i.wait()
    vcftmp[index].close
    for j in open(os.path.join(dirpath,str(index) + ".vcf"),'r'):
        if j[0] == '#' and index != 0:
            continue
        if not bedsub:
            if j.find(bedpath[0]) != -1:
                bedsub = True
                j = j.replace(bedpath[0],bed)
        outfile.append(j.strip())
for i in outfile:
    print i
