import string,sys,math,os,numpy, time

# ========= READ 2D injection nodes file ==============

Nodes2D_list = []

currentfile = open('./Injection_nodes.txt','r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    node_id = float(line)
    Nodes2D_list.append(node_id)


sheet = 21
nodes_2d = 7015

Nodes3D_list = np.array(Nodes2D_list) * (nodes_2d * sheet-1)


    if float(zone) == 72:
        Elements_list.append(elem)
    if float(zone) == 73:
        Elements_list.append(elem)
    if float(zone) == 74:
        Elements_list.append(elem)
    if float(zone) == 75:
        Elements_list.append(elem)
    if float(zone) == 76:
        Elements_list.append(elem)
    if float(zone) == 77:
        Elements_list.append(elem)
    if float(zone) == 78:
        Elements_list.append(elem)
    if float(zone) == 79:
        Elements_list.append(elem)
    if float(zone) == 80:
        Elements_list.append(elem)

currentfile.close()

# ========= WRITE RIVERBED ELEMENTS ==============

outputfile = './Riverbed_elements.dat'
currentfile = open(outputfile,'w')

for elem in Elements_list:
    currentfile.write(str(elem) +'\n')

currentfile.close()
