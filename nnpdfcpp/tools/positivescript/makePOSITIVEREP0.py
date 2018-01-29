#!/usr/bin/python
import os
import sys
import shutil

if len(sys.argv) < 2:
    print "\nusage: ./buildgrid [prior set name] \n"
    exit()

dirb  = sys.argv[1]

pdfdir = "./"

rep = 100
base= dirb
out = base + "_mc"
print "output:", out
print "replicas:", rep

try:
    os.mkdir(out)
except:
    pass

index = 1
for i in range(1,rep+1):    

    id = int(i)
    if id < 10:
        a = "000" + str(id)
    elif id < 100:
        a = "00" + str(id)
    elif id < 1000:
        a = "0" + str(id)
    else:
        a = str(id)

    if index < 10:
        b = "000" + str(index)
    elif index < 100:
        b = "00" + str(index)
    elif index < 1000:
        b = "0" + str(index)
    else:
        b = str(index)

    print "copy ", a, "\t->\t", b
    src = pdfdir + base + "/" + base + "_" + a + ".dat"
    dest= out + "/" + out + "_" + b + ".dat"
    shutil.copy(src,dest)
    index+=1

src = pdfdir + base + "/" + base + ".info"
dest= out + "/" + out + ".info"
f = open(src,'rb')
o = open(dest,'wb')

for i in f.readlines():
    o.write(i)

o.close()
f.close()

    
# replica 0
xpdf = []
xgrid = []
qgrid = []
fgrid = []
for i in range(1,rep+1):

    print "- reading replica", i
    if i < 10:
        a = "000" + str(i)
    elif i < 100:
        a = "00" + str(i)
    elif i < 1000:
        a = "0" + str(i)
    else:
        a = str(i)
    
    w = open(out + "/" + out + "_" + a + ".dat",'rb')    
    xpdf.append([])
    for j in range(0,2): w.readline()        
    
    s = 0
    while True: 
        w.readline()        

        xs = w.readline()
        nx = len(xs.split())

        qs = w.readline()
        nq = len(qs.split())

        fs = w.readline()
        nfl = len(fs.split())            
        
        xpdf[i-1].append([])
        
        if nx == 0:
            break
        
        if i == 1:
            xgrid.append(xs)
            qgrid.append(qs)
            fgrid.append(fs)        
        
        for ix in range(0,nx):
            xpdf[i-1][s].append([])
            for iq in range(0,nq):
                xpdf[i-1][s][ix].append([])
                line = w.readline().split()
                for ifl in range(0,nfl):                    
                    number = float(line[ifl]) 
                    if number < 0:
                        xpdf[i-1][s][ix][iq].append(1e-10)
                    else:
                        xpdf[i-1][s][ix][iq].append(float(line[ifl]))
                    
        s+=1

    w.close()
    

print "computing replica 0"
w = open(out + "/" + out + "_0000.dat",'wb')    
w.write("PdfType: central\n")
w.write("Format: lhagrid1\n---\n")

for s in range(len(qgrid)):
    w.write(xgrid[s])
    w.write(qgrid[s])
    w.write(fgrid[s])
    for ix in range(len(xgrid[s].split())):
        for iq in range(len(qgrid[s].split())):
            w.write(" ")
            for ifl in range(len(fgrid[s].split())):
                sum = 0
                for irep in range(0,rep):
                    sum += xpdf[irep][s][ix][iq][ifl]
                sum /= rep
                print >> w, "%14.7E" % sum,
            w.write("\n")
    w.write("---\n")
w.close()
    

