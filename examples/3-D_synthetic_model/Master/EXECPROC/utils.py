import string,sys,math,os,time
from math import *
import numpy as np
from os import listdir
import glob

# FUNCTIONS

def chunks(l, n):
    for i in range(0, len(l),n):
        yield l[i:i + n]

def get_ages(decay,Ceq,Ct):
    t = 1/decay * log(Ceq/max(1,(Ceq-val)))
    return t


def P_Rn(porosity_aq,Co_Rn, decay_cst_rn):
    '''
    This is the production function for 222Rn as deascribed in Guillon et al. 2016.
    UNIT : kg cm-3 s-1
    '''
    p = Co_Rn * decay_cst_rn * porosity_aq # Bq.m-3.d-1
    return p

def frange(start, stop, numelements):
    """range function for floats"""
    incr = (stop - start) / numelements
    return (start + x * incr for x in range(numelements))

def exprange(start, stop, numelements):
    """exponential range - each element is a fixed factor bigger than the previous"""
    return (exp(x) for x in frange(log(start), log(stop), numelements))

def P_func(minP,maxP,n_layer):
    P = np.array([x for x in exprange(minP,int(maxP),int(n_layer))])
    return P

def vg_func(P,alpha,beta,Swr,lp):
    n = (1-(1/beta))
    if P >= 0:
        Sw = 1
    else:
        Sw = Swr + (1 - Swr)*(1+abs(alpha*P)**beta)**(-n)
    Se = (Sw - Swr) / (1-Swr)
    Kr = Se**lp * ( (1 - (1-Se**(1/n))**n)**2 )
    return Sw, Kr, Se

def C0_37Ar_background(depth,l,p_surf,p_situ,porosity_aq,decay_cst):
    Co = p_surf * np.exp(-depth*(1/l))+p_situ # Bq m-3
    p = Co*decay_cst*porosity_aq # Bq m-3 d-1
    return Co, p

def get_P_Ar_background(Depths,l,p_surf,p_situ,porosity_aq,decay_cst_Ar):
    # Argon 37
    P = []
    Co = []
    for depth in Depths:
        co,p = C0_37Ar_background(depth,l,p_surf,p_situ,porosity_aq,decay_cst_Ar)
        P.append(p)
        Co.append(co)
    return P, Co

def get_P_Ar_background_river(Depths,l,p_surf,p_situ,porosity_aq,decay_cst_Ar):
    # Argon 37
    P = []
    Co = []
    # Production below river
    #h = 1.0 # m
    for depth in Depths:
        co,p = C0_37Ar_background(depth,l,p_surf,p_situ,porosity_aq,decay_cst_Ar)
        #att_fact = np.exp( -(float(h)/(float(l))) * (0.997) )
        #p_river = p * att_fact
        P.append(p)
        Co.append(co)
    return P, Co


def transport_BC_inc_firsttype(n_layer,Co_Rn,Co_Ar,Co_He):
    outputfile = './HGS/Grokfiles/transport_BC.inc'
    currentfile = open(outputfile,'w')
    for i in range(n_layer-1):
        co_Ar = Co_Ar[::-1][i]
        currentfile.write('Choose nodes y plain' + '\n' )
        currentfile.write('500.0'  + '\n' )
        currentfile.write('0.5'  + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('Specified concentration' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(Co_Rn) + ' ' +  str(co_Ar) + ' ' +  str(Co_He) + ' ' +  str(200) + ' ' + str(0.0) + ' ' + str(100.0) + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    # Riverbed production
    for i in range(n_layer-1):
        co_Ar = Co_Ar[::-1][i]
        currentfile.write('Choose faces vertical from am nodes' + '\n' )
        currentfile.write('../Mesh/Selections/River_inflow'  + '\n' )
        currentfile.write(str(i+1) + ',' + str(i+2) + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('Specified third-type concentration' + '\n' )
        currentfile.write('True' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(Co_Rn) + ' ' +  str(co_Ar) + ' ' +  str(Co_He) + ' ' +  str(200) + ' ' + str(0.0) + ' ' + str(100.0) +  '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    currentfile.close()

def transport_BC_inc(n_layer,Co_Rn,Co_Ar,Co_He):
    outputfile = './HGS/Grokfiles/transport_BC.inc'
    currentfile = open(outputfile,'w')
    for i in range(n_layer-1):
        co_Ar = Co_Ar[::-1][i]
        currentfile.write('Choose faces vertical from am nodes' + '\n' )
        currentfile.write('../Mesh/Selections/Upstream_BC'  + '\n' )
        currentfile.write(str(i+1) + ',' + str(i+2) + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('Specified third-type concentration' + '\n' )
        currentfile.write('True' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(Co_Rn) + ' ' +  str(co_Ar) + ' ' +  str(Co_He) + ' ' +  str(200) + ' ' + str(0.0) + ' ' + str(100.0) + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    # Riverbed production
    for i in range(n_layer-1):
        co_Ar = Co_Ar[::-1][i]
        currentfile.write('Choose faces vertical from am nodes' + '\n' )
        currentfile.write('../Mesh/Selections/River_inflow'  + '\n' )
        currentfile.write(str(i+1) + ',' + str(i+2) + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('Specified third-type concentration' + '\n' )
        currentfile.write('True' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(Co_Rn) + ' ' +  str(co_Ar) + ' ' +  str(Co_He) + ' ' +  str(200) + ' ' + str(0.0) + ' ' + str(100.0) +  '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    currentfile.close()


def rch_inc(rch):
    outputfile = './HGS/Grokfiles/recharge.inc'
    currentfile = open(outputfile,'w')
    currentfile.write('Choose nodes top am' + '\n' )
    currentfile.write('../Mesh/Selections/Floodplain_elements' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('create face set' + '\n' )
    currentfile.write('rch_faces' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('boundary condition' + '\n' )
    currentfile.write('type' + '\n' )
    currentfile.write('Flux' + '\n' )
    currentfile.write('name' + '\n' )
    currentfile.write('Recharge' + '\n' )
    currentfile.write('face set' + '\n' )
    currentfile.write('rch_faces' + '\n' )
    currentfile.write('time value table' + '\n' )
    currentfile.write('0.0 ' + str(rch) + '\n' )
    currentfile.write('end' + '\n' )
    currentfile.write('end' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('clear chosen zones' + '\n' )
    currentfile.write('clear chosen nodes' + '\n' )
    currentfile.write('clear chosen elements' + '\n' )
    currentfile.write('clear chosen faces' + '\n' )
    currentfile.write('clear chosen segments' + '\n' )


def write_mprop(props_soil,P,Sw,Kr):
    outputpath = ('./HGS/Parameters/props_pm.mprops')
    outputfile = open(outputpath,'w')
    outputfile.write('Aquifer' + '\n')
    outputfile.write('Unsaturated tables' + '\n')
    outputfile.write('pressure-saturation' + '\n')
    for i in range(len(P)):
        p = P[i]
        sat = Sw[i]
        outputfile.write(str(p) + '       ' + str(sat) + '\n')
    outputfile.write('end ! pressure-saturation table' + '\n')
    outputfile.write('' + '\n')
    outputfile.write('saturation-relative k' + '\n')
    for i in range(len(P)):
        kr = Kr[i]
        sat = Sw[i]
        outputfile.write(str(sat) + '       ' + str(kr) + '\n')
    outputfile.write('end ! saturation-relative k table' + '\n')
    outputfile.write('end ! unsaturated tables' + '\n')
    outputfile.write('specific storage' + '\n')
    outputfile.write(str(props_soil['Ss']) + '\n')
    outputfile.write('end material' + '\n')
    outputfile.close()

def P_Ar(z,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w):
    '''
    This is the production function for 37Ar as deascribed in Guillon et al. 2016.
    UNIT : kg cm-3 s-1
    '''
    p = (Ca + 0.38 * K) * Xsc * p_surf * np.exp( -(float(z)/float(l)) * ( ( (1-porosity_aq) * rho_g ) + ( porosity_aq * rho_w * 1) ) )
    return p

def P_Ar_river(z,hp,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w,Sw):
    '''
    This is the production function for 37Ar as deascribed in Guillon et al. 2016.
    UNIT : kg cm-3 s-1
    '''
    p = (Ca + 0.38 * K) * Xsc * p_surf * np.exp( -(float(z)/float(l)) * ( ( (1-porosity_aq) * rho_g ) + ( porosity_aq * rho_w * Sw) ) )
    att_fact = np.exp( -(float(hp)/(float(l)/1000.0)) * (rho_w/1000.0) )
    p_river = p * att_fact
    #print(att_fact)
    return p_river

def get_P_Ar(Depths,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w,emmanation,Sw):
    P = []
    for depth in Depths:
        p = P_Ar(depth,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w,Sw)
        ac = (p*(emmanation/porosity_aq))/(365.0*24.0*3600.0) # Bq cm-3 d-1 
        #ac = (p*(emmanation/porosity_aq))/(365.0) # at cm-3 d-1 
        #ac_w = (ac / porosity_aq) # Bq/cm3eau
        ac_w_hgs = ac * 10**6 # Bq m-3 d-1
        #print(ac_w_hgs)
        ac_w_hgs_scaled = ac_w_hgs / 1000.0 *10**6 # scaled production in Bq.L-1.d-1
        P.append(ac_w_hgs_scaled)
    return P

def get_P_Ar_river(Depths,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w,emmanation,Sw):
    P = []
    # Production below river
    h = 100 # cm
    hp = h * 101.97 # g/cm2
    hp = hp/1000.0 # kg/cm2
    for depth in Depths:
        #depth_corrigee = (depth_m + 3.0) * 101.97 # g/cm2
        #print(depth_corrigee)
        p = P_Ar_river(depth,hp,Ca,K,Xsc,p_surf,l,porosity_aq,rho_g,rho_w,Sw)
        #print(p)
        ac = (p*(emmanation/porosity_aq))/(365.0*24.0*3600.0) # Bq cm-3 d-1 
        #ac = (p*(emmanation/porosity_aq))/(365.0) # at cm-3 d-1 
        ac_w_hgs = ac * 10**6 # Bq m-3 d-1
        ac_w_hgs_scaled = ac_w_hgs / 1000.0 *10**6 # scaled production in Bq.L-1.d-1
        P.append(ac_w_hgs)
    return P

def transport_Source_Partition_Elements_inc(prod_Ar_list,prod_Rn_list,Storage_list):
    outputfile = './HGS/Grokfiles/transport_Source.inc'
    currentfile = open(outputfile,'w')
    for i in range(len(prod_Ar_list)):
        prod_Rn = prod_Rn_list[i]
        prod_Ar = prod_Ar_list[i]
        storage = Storage_list[i]
        currentfile.write('choose zone number' + '\n' )
        currentfile.write(str(i+1)  + '\n' )
        currentfile.write('Zero order source with partitioning' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(prod_Rn) + ' ' +  str(prod_Ar) + ' ' +  str(0.0) + ' ' + '0.0' +  ' ' + str(storage) + ' ' + '0.0' + ' ' + '0.35' + ' ' + '0.04182' + ' ' + '0.00' + ' ' + '0.00' + ' ' + '1.00' + ' ' + '1.00' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    currentfile.close()


def transport_Source_Partition_inc(n_layer,P,P_Ar_rb,prod_rate_Rn_Aq,prod_rate_Rn_Ch):
    outputfile = './HGS/Grokfiles/transport_Source.inc'
    currentfile = open(outputfile,'w')
    for i in range(n_layer-1):
        prod_rate_Ar = P[::-1][i]
        currentfile.write('choose zone number' + '\n' )
        currentfile.write(str(i+1)  + '\n' )
        currentfile.write('Zero order source with partitioning' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(prod_rate_Rn) + ' ' +  str(prod_rate_Ar) + ' ' +  str(0.0) + ' ' + '0.0' +  ' ' + '0.15' + ' ' + '0.0' + ' ' + '0.35' + ' ' + '0.04182' + ' ' + '0.00' + ' ' + '0.00' + ' ' + '1.00' + ' ' + '1.00' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    i = i+1
    # Riverbed production
    for n in range(n_layer-1):
        zone_rb = i+n+1
        prod_rate_Ar = P_Ar_rb[::-1][n]
        currentfile.write('choose zone number' + '\n' )
        currentfile.write(str(zone_rb)  + '\n' )
        currentfile.write('Zero order source with partitioning' + '\n' )
        currentfile.write('1' + '\n' )
        currentfile.write('0 1.0d20 ' + str(prod_rate_Rn) + ' ' +  str(prod_rate_Ar)  + ' ' +  str(0.0) + ' ' + '0.0' +  ' ' + '0.15' + ' ' + '0.0' + ' ' + '0.35' + ' ' + '0.04182' + ' ' + '0.00' + ' ' + '0.00' + ' ' + '1.00' + ' ' + '1.00' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('clear chosen zones' + '\n' )
        currentfile.write('clear chosen nodes' + '\n' )
        currentfile.write('clear chosen elements' + '\n' )
        currentfile.write('clear chosen faces' + '\n' )
        currentfile.write('clear chosen segments' + '\n' )
        currentfile.write('' + '\n' )
        currentfile.write('' + '\n' )
    currentfile.close()

def write_GHB(u_bc,d_bc):
    outputpath = ('./HGS/Grokfiles/bc.inc')
    currentfile = open(outputpath,'w')
    currentfile.write('clear chosen zones' + '\n' )
    currentfile.write('clear chosen nodes' + '\n' )
    currentfile.write('clear chosen elements' + '\n' )
    currentfile.write('clear chosen faces' + '\n' )
    currentfile.write('clear chosen segments' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('Choose nodes am' + '\n' )
    currentfile.write('../Mesh/Selections/Upstream_BC' + '\n' )
    currentfile.write('1,34' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('create face set' + '\n' )
    currentfile.write('Upstream' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('boundary condition' + '\n' )
    currentfile.write('  type' + '\n' )
    currentfile.write('    fluid transfer' + '\n' )
    currentfile.write('  name' + '\n' )
    currentfile.write('    Upstream_Specified_head' + '\n' )
    currentfile.write('  face set' + '\n' )
    currentfile.write('    Upstream' + '\n' )
    currentfile.write('  time value table' + '\n' )
    currentfile.write('    0.0  ' + str(u_bc) + '\n' )
    currentfile.write('  end' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('fluid transfer coefficients' + '\n' )
    currentfile.write(' 500' + '\n' )
    currentfile.write(' 1e3' + '\n' )
    currentfile.write('end' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('clear chosen zones' + '\n' )
    currentfile.write('clear chosen nodes' + '\n' )
    currentfile.write('clear chosen elements' + '\n' )
    currentfile.write('clear chosen faces' + '\n' )
    currentfile.write('clear chosen segments' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('Choose nodes am' + '\n' )
    currentfile.write('../Mesh/Selections/Downstream_BC' + '\n' )
    currentfile.write('1,34' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('create face set' + '\n' )
    currentfile.write('Downstream' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('boundary condition' + '\n' )
    currentfile.write('  type' + '\n' )
    currentfile.write('    fluid transfer' + '\n' )
    currentfile.write('  name' + '\n' )
    currentfile.write('    Downstream_Specified_head' + '\n' )
    currentfile.write('  face set' + '\n' )
    currentfile.write('    Downstream' + '\n' )
    currentfile.write('  time value table' + '\n' )
    currentfile.write('    0.0  ' + str(d_bc) + '\n' )
    currentfile.write('  end' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('fluid transfer coefficients' + '\n' )
    currentfile.write(' 500' + '\n' )
    currentfile.write(' 1e3' + '\n' )
    currentfile.write('end' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.close()
