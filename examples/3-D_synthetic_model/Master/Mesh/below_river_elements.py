import string,sys,math,os,numpy, time

# ========= READ RIVERBED ELEMENTS from Background Zone file ==============

Elements_list = []

currentfile = open('./Zones','r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    elem,zone = line.split()
    if float(zone) > 15:
        Elements_list.append(elem)

currentfile.close()

# ========= WRITE RIVERBED ELEMENTS ==============

outputfile = './Belowriver_elements.dat'
currentfile = open(outputfile,'w')

for elem in Elements_list:
    currentfile.write(str(elem) +'\n')

currentfile.close()
