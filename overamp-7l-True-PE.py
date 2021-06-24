#! /usr/bin/env python2
import sys, math, os, time
from optparse import OptionParser
from collections import defaultdict

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2011
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This program looks at a sam file and picks out the unique reads.
#It can be used on a single or paired ended alignment.
#
#INPUT:
#THis program takes a samfile as input mapped SE
#
#OUTPUT:
#This program outputs a unique read sma file, and prints alignment statistics to the screen.
#The output columns are [Filename, #unique, #aligned, #unique, total reads].
#
#NOTE:
#This program is designed to be used in conjuction wit bwa-samtools-do-all.py, 
#with its printed output being redirected to a combined all lib files, 
#that is why there are no headers for the printed output.
#Also, overamp generates the unique files while doing the counts, not independently.
#
#PARAMETERS:
#1. REQUIRED:
#   i. -f or --samfile, The input.sam file
#   ii.  -o or--outfile, The output unique .sam file name
#2. OPTIONAL:
#   i. -s or --sort, sort sam by chrom and start position
#   ii. -p or --paired, Switches to pair-end sam file. Default is single ended

usage = "USAGE: overAmp.py -f samfile.sam -o outputfile [-p] [-s]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--samfile", dest="f", help="Input sorted .sam file.")
parser.add_option("-o", "--outfile", dest="o", help="Output file.")
parser.add_option("-b", "--blocksize", dest="b", type = 'int', default = 10000, help="Blocksize before squish unpaired reads")
parser.add_option("-m", "--maxdist", dest="m", type = 'int', default = 800, help="largest maximum normal insert")

(opt, args) = parser.parse_args()

squishtrig = opt.b
maxdist = opt.m

actsize = 150000
iact = actsize+0
goodflags = ['99-147','147-99', '83-163', '163-83']

allreads = defaultdict(int)
outreads = defaultdict(int)


print "CURRENT VERSION ONLY FOR PE MEM MAPPED SE, WITH SORTED REFERENCES!"

def outblock(newblock, ofile):
   for item in newblock:
      b1 = item[0]
      b2 = item[1]
      b1[3] = str(b1[3])
      b2[3] = str(b2[3])
      o1 = '\t'.join(b1)+'\n'
      o2 = '\t'.join(b2)+'\n'
      ofile.write(o1+o2)   

def collapse(block):
   newblock = []
   bcpy = block[:]
   bcpy.sort(key = lambda x: len(x[0][9]), reverse = True)
   bcpy.sort(key = lambda x: int(x[0][4]), reverse = True)
   newblock.append(bcpy.pop(0))
   while len(bcpy) > 0:
      top = bcpy.pop(0)
      flag = 1
      for item in newblock:
         if top[0][9] in item[0][9]:
            flag = 0
         if top[1][9] in item[1][9]:
            flag = 0
   
      if flag == 1:
         newblock.append(top)
   return newblock

pretime = time.time()
f = open(opt.f)
o = open(opt.o, 'w')
chrom = ''
pos = ''
hold = ''
flen = 0


unpaired = defaultdict(list)
blocks = defaultdict(list)
oldunpaired = defaultdict(list)

activereads = defaultdict(lambda: defaultdict(list))
sereads = defaultdict(list)
twochrom = defaultdict(lambda: defaultdict(list))


starttime = time.time()
seen = []

f.seek(0)
i = 0
headerlist = []
sflag = 0
while 1:
   i+=1
   l = f.readline()
   
   if l == '':
      break
   if l[0] == '@':
      if "@SQ" in l and 'LN:0' not in l:
         o.write(l)
         headerlist.append(l.split()[1][3:])
      continue

   x = l.split()
   name, flag, ref, pos = x[:4]
   x[3] = int(x[3])
   pos = int(pos)
   #if auxilary score / read
   if int(flag) > 200:
      continue

   #if seen before
   if name in unpaired:
      #get other pai of read, remove from unpaired
      otherhalf = unpaired.pop(name)
      flags = [flag, otherhalf[1]]
      flagtype = '-'.join(flags)
      allreads[flagtype] +=1
      #if one half unmapped
      if "*" in flags:
         pass
      else:   
         tup = '^'.join([ref, str(pos), otherhalf[2], str(otherhalf[3])])
      if flagtype not in goodflags:
         continue
      activereads[flagtype][tup].append([x,otherhalf])
      if ref == otherhalf[2]:
         activereads[flagtype][tup].append([x,otherhalf])
      else:
         twochrom[flagtype][tup].append([x,otherhalf])  
   #if new read name
   else:
      unpaired[name] = x
   #if too many unpaired, compare to old lsit
   if len(unpaired) > squishtrig:
      #print "Squish", len(unpaired), len(oldunpaired)
      removelist = []
      rem = 0
      for item in unpaired:
         ti = unpaired[item]
         if pos - ti[3] > maxdist or ref != ti[2] or ti[2] == '*':
            if ti[0] in oldunpaired:
               #something!!!
               firsthalf = unpaired[ti[0]]
               removelist.append(ti[0])
               otherhalf = oldunpaired.pop(ti[0])
               flags = [firsthalf[1], otherhalf[1]]
               tup = '^'.join([firsthalf[2], str(firsthalf[3]), otherhalf[2], str(otherhalf[3])])
               flagtype = '-'.join(flags)
               allreads[flagtype]+=1
               if flagtype in goodflags:
                  if firsthalf[2] == otherhalf[2]:
                     activereads[flagtype][tup].append([firsthalf,otherhalf])
                  else:
                     twochrom[flagtype][tup].append([firsthalf,otherhalf])   
               rem +=1
            else:
               oldunpaired[ti[0]] = ti
               removelist.append(ti[0])
      for sub in removelist:
         temp = unpaired.pop(sub)
      #print "Squished", len(unpaired), len(oldunpaired), rem, round(time.time() - starttime,1)
      if len(unpaired) > squishtrig*.9:
         squishtrig = int(squishtrig * 1.5)
      if len(unpaired) < 100:
         squishtrig = opt.b
      #print i, allreads, sum(allreads.values())
   #if many positions found and ready to output  
   if len(activereads["147-99"])+len(activereads["83-163"]) > actsize:
      #FLUSH
      if sflag == i-1:
         actsize = actsize*2.0
      else:
         actsize = iact+0
      sflag = i
      print "Flushing", i, len(activereads["147-99"]), len(activereads["83-163"]), len(twochrom["163-83"]), len(twochrom["99-147"]), actsize
      for flaggroup in activereads:
         ck = activereads[flaggroup].keys()
         if len(ck) == 0:
            continue
         ck.sort(key = lambda x: int(x.split('^')[3]))
         ck.sort(key = lambda x: headerlist.index(x.split('^')[2]))
         ck.sort(key = lambda x: int(x.split('^')[1]))
         ck.sort(key = lambda x: headerlist.index(x.split('^')[0]))

         flushchrom = ck[-1].split('^')[0]
         flushpos = min(int(ck[-1].split('^')[1]),int(ck[-1].split('^')[3]))-maxdist
         #print flushchrom, flushpos, ck[0], ck[-1], headerlist.index(flushchrom)
         if flushchrom not in seen:
            seen.append(flushchrom)
            sflag = 1
            print flushchrom, round(time.time() - starttime,1)
         for item in ck:
            ch1, pos1, ch2, pos2 = item.split('^')
            pos1 = int(pos1)
            pos2 = int(pos2)
            if (ch1 == flushchrom and max(pos1,pos2) < flushpos) or ch1 != flushchrom:
               #outit
               block = activereads[flaggroup].pop(item)
               if len(block) == 1:
                  outblock(block, o)
                  outreads[flaggroup]+=1             
               else:
                  #collapse!!!!
                  newblock = collapse(block)
                  #for item in newblock:
                  outreads[flaggroup]+=len(newblock)
                  outblock(newblock, o)
#      print allreads
#      print outreads                    
      #print len(activereads["16-0"]), len(activereads["0-16"]), '\n'

       
      ###### USE ACTIVE RECENTLY READ IN CHROM? IF NOT THAT?
      ########  ADD A BROKEN ACTIVE FOR BEWTEEN TWO CHROMS
      ######## ADD ONE FOR no chrom * ref unmapped
#      x[32332]
      if sflag != 1:
         continue
      print opt.f, "2Chr Flushing", i, len(twochrom["147-99"]), len(twochrom["83-163"]), round(time.time() - starttime,1)
      for flaggroup in twochrom:
         ck = twochrom[flaggroup].keys()
         ck.sort(key = lambda x: int(x.split('^')[3]))
         ck.sort(key = lambda x: headerlist.index(x.split('^')[2]))
         ck.sort(key = lambda x: int(x.split('^')[1]))
         ck.sort(key = lambda x: headerlist.index(x.split('^')[0]))
         for item in ck:
#            x[323223]
            ch1, pos1, ch2, pos2 = item.split('^')
            pos1 = int(pos1)
            pos2 = int(pos2)
            if headerlist.index(ch1) < headerlist.index(flushchrom) and headerlist.index(ch2) < headerlist.index(flushchrom):
               #outit
               block = twochrom[flaggroup].pop(item)
               if len(block) == 1:
                  outblock(block, o)
                  outreads[flaggroup]+=1            
               else:
                  #collapse!!!!
                  newblock = collapse(block)
                  outreads[flaggroup]+=len(newblock)
                  outblock(newblock, o)                   
      #print len(twochrom["16-0"]), len(twochrom["0-16"]),'\n'
#      print allreads, round(time.time() - starttime,1)
#      print outreads, round(time.time() - starttime,1)

         
   
#REDO THIS SECTION
removelist = []
rem = 0
for item in unpaired:
   ti = unpaired[item]
   if ti[0] in oldunpaired:
      #something!!!
      firsthalf = unpaired[ti[0]]
      removelist.append(ti[0])
      otherhalf = oldunpaired.pop(ti[0])
      flags = [firsthalf[1], otherhalf[1]]
      tup = '^'.join([firsthalf[2], str(firsthalf[3]), otherhalf[2], str(otherhalf[3])])
      flagtype = '-'.join(flags)
      allreads[flagtype] +=1
      if flagtype in goodflags:
         if firsthalf[2] == otherhalf[2]:
            activereads[flagtype][tup].append([firsthalf,otherhalf])
         else:
            twochrom[flagtype][tup].append([firsthalf,otherhalf])   
      rem +=1
   else:
      oldunpaired[ti[0]] = ti
      removelist.append(ti[0])
for sub in removelist:
   temp = unpaired.pop(sub)






#output remainers
for flaggroup in activereads:
   ck = activereads[flaggroup].keys()
   for item in ck:
      block = activereads[flaggroup].pop(item)
      if len(block) == 1:
         outreads[flaggroup]+=1 
         outblock(block, o)             
      else:
         #collapse!!!!
         newblock = collapse(block)
         outreads[flaggroup]+=len(newblock)
         outblock(newblock, o)                     

for flaggroup in twochrom:
   ck = twochrom[flaggroup].keys()
   for item in ck:
      block = twochrom[flaggroup].pop(item)
      if len(block) == 1:
         outreads[flaggroup]+=1 
         outblock(block, o)            
      else:
         #collapse!!!!
         newblock = collapse(block)
         outreads[flaggroup]+=len(newblock)
         outblock(newblock, o)          


f.close()
o.close()
o2 = open("stats-"+opt.f, 'w')

good = 0
all = 0
#print opt.f, ctout, ctin, weird, float(ctout+weird)/ctin
for key in goodflags:
   o2line = [opt.f, key, outreads[key], round(float(outreads[key])/float(allreads[key]),2), allreads[key]]
   o2line = map(lambda z: str(z), o2line)
   o2line = '\t'.join(o2line)+'\n'
   good+=outreads[key]
   all +=allreads[key]
   o2.write(o2line)



o2line = [opt.f, good, all, good/float(all), sum(allreads.values())]
o2line = map(lambda z: str(z), o2line)
o2line = '\t'.join(o2line)+'\n'
o2.write(o2line)
















