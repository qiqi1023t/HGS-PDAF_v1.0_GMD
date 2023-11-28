import string,sys,math,struct,array,csv
import numpy
import time

Subsurface_vsl_path = '../'
Gridname = 'SuperGrid.dat'

tecplot_prism_filename = Subsurface_vsl_path + Gridname

#READ X and Y COORDINATES for SURFACE
olffile = '../Flowo.olf.dat'

currentfile = open(olffile,'r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '# x':
        break

# Read X coordinates

X = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# y':
        break
    temp = line.split()
    for i in range(len(temp)):
        X.append(float(temp[i]))

Xsurf = X
X = []

# Read Y coordinates

Y = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# z':
        break
    temp = line.split()
    for i in range(len(temp)):
        Y.append(float(temp[i]))

Ysurf = Y
Y = []

# Read Z coordinates

Z = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# zone (cell centred':
        break
    temp = line.split()
    for i in range(len(temp)):
        Z.append(float(temp[i]))

Zsurf = Z
Z = []


xmin = 1e38
xmax = -1e38
ymin = 1e38
ymax = -1e38
zmin = 1e38
zmax = -1e38

#baricentric_coor = ((1./3.,1./3.,1./3.),(.5,.25,.25),(.25,.5,.25),(.25,.25,.5),(2./3.,1./6.,1./6.),(1./6.,2./3.,1./6.),(1./6.,1./6.,2./3.))
baricentric_coor = ((1./3.,1./3.,1./3.),(2./3.,1./6.,1./6.),(1./6.,2./3.,1./6.),(1./6.,1./6.,2./3.))

nsubdiv = 5  
if nsubdiv < 3:   # on subdivise le prism en nsubdiv sections mais on elimine l,inferieure et la superieur
  nsubdiv = 3     # d<ou la necessite d<avoir au moins 3 subdivions pour en obtenir au mins 1

# ------------------ DEFINITION OF FUNTIONS ------------------

def getnodesinelem(nodes_number_list):  # recoit [node1_numbre,node2_numbre,node3_numbre,node4_numbre,node5_numbre,node6_numbre]
                                        #retourne [[x1,y1,z1],[x2,y2,z2].....[xn,yn,zn]] 
                                       # les coordonnees des triangles inferieurs et superieurs sont les memes
  pos_node1 = nodes_number_list[0]-1   # triangle inferieur
  pos_node2 = nodes_number_list[1]-1
  pos_node3 = nodes_number_list[2]-1
  pos_node4 = nodes_number_list[3]-1   # triangle superieur
  pos_node5 = nodes_number_list[4]-1
  pos_node6 = nodes_number_list[5]-1
  coorlist_inf = []
  for each_baricentric_coor in baricentric_coor:
     # calcule uniquement pour l'inferieur car idem pour le superieur
     xnew = each_baricentric_coor[0] * xcoor[pos_node1] +each_baricentric_coor[1] * xcoor[pos_node2] +each_baricentric_coor[2] * xcoor[pos_node3]
     ynew=   each_baricentric_coor[0] * ycoor[pos_node1] +each_baricentric_coor[1] * ycoor[pos_node2] +each_baricentric_coor[2] * ycoor[pos_node3]
     znew=   each_baricentric_coor[0] * zcoor[pos_node1] +each_baricentric_coor[1] * zcoor[pos_node2] +each_baricentric_coor[2] * zcoor[pos_node3]
     coorlist_inf.append([xnew,ynew,znew])
  coorlist_sup = []
  for each_baricentric_coor in baricentric_coor:
     xnew = each_baricentric_coor[0] * xcoor[pos_node4] +each_baricentric_coor[1] * xcoor[pos_node5] +each_baricentric_coor[2] * xcoor[pos_node6]
     ynew=   each_baricentric_coor[0] * ycoor[pos_node4] +each_baricentric_coor[1] * ycoor[pos_node5] +each_baricentric_coor[2] * ycoor[pos_node6]
     znew=   each_baricentric_coor[0] * zcoor[pos_node4] +each_baricentric_coor[1] * zcoor[pos_node5] +each_baricentric_coor[2] * zcoor[pos_node6]
     coorlist_sup.append([xnew,ynew,znew])
  #print coorlist_inf
  #print coorlist_sup
  coorlist = []
  for each_node_twin in zip(coorlist_inf,coorlist_sup):
    # each_node_twin  = ([xinf,yinf,zinf],[xsup,ysup,zsup]) 
    #   |x|x|x|x|x|     on ne conserver que les points internes, donc nsubdiv - 2
    #                   valeurs
    #print each_node_twin
    #print nsubdiv
    varx = (each_node_twin[1][0] - each_node_twin[0][0] ) / nsubdiv
    vary = (each_node_twin[1][1] - each_node_twin[0][1] ) / nsubdiv
    varz = (each_node_twin[1][2] - each_node_twin[0][2] ) / nsubdiv
    #print varx,vary,varz
    for isub in range(1,nsubdiv):
      xnew = each_node_twin[0][0] + isub * varx
      ynew = each_node_twin[0][1] + isub * vary
      znew = each_node_twin[0][2] + isub * varz
      coorlist.append([xnew,ynew,znew])
  #print coorlist
  #sys.exit(0)
  return coorlist 


#################################### Lecture des prism ########################################

print('Reading prism...')


""" file format

TITLE = "Ile_1st_mesh                                                "
TITLE = "Tobi, Donnerstag, Oct 05, 2011 at 10:22                     "
TITLE = "Beispielgrid als Grundlage zur Darstellung in TecPlot       "
VARIABLES = "X", "Y", "Z", "Zone"                                                                                                                                                                                                                               
ZONE  T="Grid only", DATAPACKING=BLOCK, N= 364375, E= 656280, ZONETYPE=FEBRICK, VARLOCATION=([      4]=CELLCENTERED)                                                                                                                                                                                                                                                                                                                                                                                                            
# x
   282219.00000000000        282288.00000000000        282275.00000000000        2
# y
   5209341.0000000000        5209159.0000000000        5208964.0000000000        5208776.0000000000        5208610.0000000000     
   5208445.0000000000        5208304.0000000000        5208172.0000000000        5207985.0000000000        5207787.0000000000    

# z
  -42.479701995849609       -44.415435791015625       -46.597087860107422       -46.930252075195312       -46.850318908691406     

# zone (cell-centered)
   1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000     
   1.

# element node lists 
       205     28818     29903     29903     33330     61943
     63028     6302

"""

tecplot_prism = open(tecplot_prism_filename)

nnodes = 0
nelem = 0
xcoor = []
ycoor= []
zcoor = []

connect =[]


# READ GENERAL FILE PROPERTIES

while nnodes == 0 or nelem == 0:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne specifiant nombre noeuds et nombre elem introuvable")
    sys.exit(0)
  if "N=" and "E=" and "ZONETYPE" and "DATAPACKING"   in lineread:
    allvals = lineread.split(',')
    for eachval in allvals:
      if "N= " in eachval:     
        nnodes = int(eachval.replace("N= ",""))
      if "E= " in eachval:
        nelem = int(eachval.replace("E= ",""))
      if "ZONETYPE= " in eachval:
        if not "FEBRICK" in eachval:
          print("Erreur type fichier: il doit contenir des elements avec une connectivite a 8 noeuds")
          sys.exit(0)
      if "DATAPACKING= " in eachval:
        if not "BLOCK" in eachval:
          print("Erreur type fichier: le datapacking doit etre en mode BLOCK")
          sys.exit(0)

# GET X GRID COORDINATES

while 1:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne contenant '# x' introuvable")
    sys.exit(0)
  if "# x" in lineread:
    break

while 1:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne contenant '# y' introuvable")
    sys.exit(0)
  if "# y" in lineread:
     break
  xvals = lineread.split()
  for eachval in xvals:
    xcoor.append(float(eachval))
    if xmin > float(eachval):
      xmin = float(eachval)
    if xmax < float(eachval):
      xmax = float(eachval)
      

if len(xcoor) != nnodes:
  print("fichier des prism incorrect: nombre de coordonnees x trouvees: " + str(len(xcoor)))
  print("nombre de coordonnees x attendues: " + ncoors)
  sys.exit(0)

# GET Y GRID COORDINATES

while 1:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne contenant '# z' introuvable")
    sys.exit(0)
  if "# z" in lineread:
     break
  yvals = lineread.split()
  for eachval in yvals:
    ycoor.append(float(eachval))
    if ymin > float(eachval):
      ymin = float(eachval)
    if ymax < float(eachval):
      ymax = float(eachval)


if len(ycoor) != nnodes:
  print("fichier des prism incorrect: nombre de coordonnees y trouvees: " + str(len(ycoor)))
  print("nombre de coordonnees y attendues: " + ncoors)
  sys.exit(0)


# GET Z GRID COORDINATES

while 1:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne contenant '# zone' introuvable")
    sys.exit(0)
  if "# propnum" in lineread:
     break
  zvals = lineread.split()
  for eachval in zvals:
    zcoor.append(float(eachval))
    if zmin > float(eachval):
      zmin = float(eachval)
    if zmax < float(eachval):
      zmax = float(eachval)


if len(zcoor) != nnodes:
  print("fichier des prism incorrect: nombre de coordonnees z trouvees: " + str(len(zcoor)))
  print("nombre de coordonnees z attendues: " + ncoors)
  sys.exit(0)

# GET MEAN BARYCENTRIC COORDINATES FOR ALL PRISM

baryelem = open('./baryelem.txt','w')

while 1:
  lineread = tecplot_prism.readline()
  if lineread == '':
    print("fichier des prism incorrect: ligne contenant '# element node lists' introuvable")
    sys.exit(0)
  if "# element node lists" in lineread:
     break

ielem = 0

#connectvals = lineread.split()
#connect_temp=[]

#for eachval in connectvals:
#    intval = int(eachval)
#    if intval not in connect_temp:
#        connect_temp.append(intval)

while 1:
    lineread = tecplot_prism.readline()
    if lineread == '':
        break
    connectvals = lineread.split()
    connect_temp=[]
    for eachval in connectvals:
        intval = int(eachval)
        if intval not in connect_temp:
            connect_temp.append(intval)
    ielem = ielem + 1
    coorlist2interpol = numpy.array(getnodesinelem(connect_temp))
    xr = numpy.mean(coorlist2interpol[:,0])
    yr = numpy.mean(coorlist2interpol[:,1])
    zr = numpy.mean(coorlist2interpol[:,2])
    # Get the three top nodes associated to the barycentric position
    dist2 = numpy.square(xr-Xsurf) + numpy.square(yr-Ysurf)
    posminA= int(numpy.argsort(dist2)[0])
    posminB= int(numpy.argsort(dist2)[1])
    posminC= int(numpy.argsort(dist2)[2])
    # Compute the top mean elevation associated to the three top nodes
    z_list = []
    z_list.append(Zsurf[posminA])
    z_list.append(Zsurf[posminB])
    z_list.append(Zsurf[posminC])
    ztop = numpy.mean(z_list)
    baryelem.write(str(ielem) + ';' + str(xr) + ';' + str(yr) + ';' + str(zr) + ';' + str(ztop) + '\n')

baryelem.close()
