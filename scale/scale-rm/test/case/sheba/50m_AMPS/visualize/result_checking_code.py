import struct
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import netCDF4 as nc4
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transform
from mpl_toolkits.axes_grid1 import make_axes_locatable


poly_exist = True

iceaxis_version = 3


#->->-> BIN INFORMATION

# maximum number of liquid bins
max_binr = 30

# define liquid bins
massbin_ = []
massbin_bnd = []
c_min = 0.000000000000004188790205
c_max = 0.0000000654498
c_maxmass = 0.52359870
temp1 = c_min
dbin = (c_max/c_min)**(0.05) # 0.1 = 1/10, 0.05 = 1/20
temp2 = temp1*dbin
massbin_bnd.append(temp1)
i = 1
while (i <= 20):
    massbin_bnd.append(temp2)
    massbin_.append(0.5*(temp1+temp2)) # liquid haze bins
    temp1 = temp2
    temp2 = temp2*dbin
    i = i + 1
dbin = (c_maxmass/c_max)**(0.05)
temp2 = temp1*dbin
i = 1
while (i <= 20):
    massbin_bnd.append(temp2)
    massbin_.append(0.5*(temp1+temp2)) # liquid drop bins
    temp1 = temp2
    temp2 = temp2*dbin
    i = i + 1
    
massbin = []
diameterbin = []
diameterbin_bnd = []
diameterbin_bnd.append((massbin_bnd[0]/(3.1416/6.0))**(1.0/3.0)*10000.0)
print("liquid diameter: ",0, "---", (massbin_bnd[0]/(3.1416/6.0))**(1.0/3.0)*10000.0)
i = 1
while (i <= max_binr):
    massbin.append(massbin_[i-1])
    diameterbin.append((massbin_[i-1]/(3.1416/6.0))**(1.0/3.0)*10000.0) # unit is \mu m
    diameterbin_bnd.append((massbin_bnd[i]/(3.1416/6.0))**(1.0/3.0)*10000.0)
    print("liquid diameter: ",i, (massbin_[i-1]/(3.1416/6.0))**(1.0/3.0)*10000.0, (massbin_bnd[i]/(3.1416/6.0))**(1.0/3.0)*10000.0)
    i = i + 1
    
# define ice bins
massbini = []
c_min = 4.18879020478639e-12
c_max = 1.0e-6
temp1 = c_min
dbin = (c_max/c_min)**(0.25)
temp2 = temp1*dbin
i = 1
while (i <= 4): 
    massbini.append(0.5*(temp1+temp2)) # ice deposition bins
    temp1 = temp2
    temp2 = temp2*dbin
    i = i + 1
c_min = c_max
c_max = 1.0e-2
dbin = (c_max/c_min)**(1.0/10.0)
temp2 = temp1*dbin
i = 1
while (i <= 10): 
    massbini.append(0.5*(temp1+temp2)) # ice collection bins
    temp1 = temp2
    temp2 = temp2*dbin
    i = i + 1
c_min = c_max
c_max = 1.0e+1
dbin = (c_max/c_min)**(1.0/6.0)
temp2 = temp1*dbin
i = 1
while (i <= 6): 
    massbini.append(0.5*(temp1+temp2)) # ice riming bins
    temp1 = temp2
    temp2 = temp2*dbin
    i = i + 1

diameteribin_bnd = []
diameteribin = [] 
diameteribin_bnd.append(0.0)
diameteribin_bnd.append(50.0)
diameteribin.append(25.0)
i = 1
while(True):
    diameteribin_bnd.append(diameteribin_bnd[i]*1.25)
    diameteribin.append(0.5*(diameteribin_bnd[i] + diameteribin_bnd[i]*1.25))
    if (diameteribin_bnd[i]*1.25 > 10000.0):
        break
    i = i + 1

for i, temp in enumerate(diameteribin):
    print("Ice diameter: ", temp, diameteribin_bnd[i])
    i = i + 1
print("Ice diameter end: ", diameteribin_bnd[i])



input_file_dir = '../history_one.pe000000.nc'
plot_dir = './'

### OPEN NETCDF DATA
# create a file object to handle the data inside the file
f = nc4.Dataset(input_file_dir,'r')

print('\n \n \n')
print("LOG -- START READING NETCDF DIMENSIONS")

if (f.file_format == 'NETCDF4'):
    #print(f)
    
    ### READ DIMENSION VARIABLES
    # time
    time = f.variables['time'][:]
    # y
    y = f.variables['y'][:]
    # z
    z = f.variables['z'][:]
    # fz
    fz = f.variables['FZ'][:]

    tlim = len(time)
    ylim = len(y)
    zlim = len(z)

    #-- dimensions are time, z, y, x
    print('tlim, zlim, ylim: ', tlim, zlim, ylim)

else:
    print("ERROR -- NETCDF VERSION NOT SUPPORTED")
    sys.exit(0)

#->->-> special method
def quickSorti(inputArray,start,end):

    if (start >= end):
        return
        
    loop_hdl = True
    #outer_loop = 0
    temp = end - start + 1
    #while (loop_hdl):
    #    temp = temp/2
    #    temp = int(temp)
    #    outer_loop = outer_loop + 1
    #    if (temp < 1):
    #        loop_hdl = False
    outer_loop = int(temp/2) + 1
    
    pivots = []
    pivots.append(start)
    pivots.append(end)
    k = 1
    segmentNumber = 1
    while (k <= outer_loop):
        kk = 1
        pivot_array = []
        temp1 = 0
        kk_pivot = 1
        while (kk <= segmentNumber):
            pivot_index = quickPartitioni(inputArray,pivots[kk_pivot-1],pivots[kk_pivot])
            if (pivot_index-1 > pivots[kk_pivot-1]):
                pivot_array.append(pivots[kk_pivot-1])
                pivot_array.append(pivot_index-1)
                temp1 = temp1 + 1
            if (pivot_index+1 < pivots[kk_pivot]):
                pivot_array.append(pivot_index+1)
                pivot_array.append(pivots[kk_pivot])
                temp1 = temp1 + 1
            kk = kk + 1
            kk_pivot = kk_pivot + 2
        
        segmentNumber = temp1
        
        if (temp1 == 0):
            break
        
        pivots = []
        for temp in pivot_array:
            pivots.append(temp)
        k = k + 1
    
def quickPartitioni(inputArray,start,end):
    
    if (start >= end):
        return start
        
    k = start
    pivot = inputArray[end]
    i = start
    while (k < end):
        if (inputArray[k] < pivot):
            temp = inputArray[k]
            inputArray[k] = inputArray[i]
            inputArray[i] = temp
            i = i + 1
        k = k + 1
    inputArray[end] = inputArray[i]
    inputArray[i] = pivot
    return i

def Iice_habit(icemass,crymass,aggmass,rimmass):
    if (rimmass > 0.1 * icemass):
        if (aggmass < 0.1 * icemass):
            if (rimmass < crymass):
                return 5 # rimed crystals
            else:
                return 4 # graupel
        else:
            if (rimmass > aggmass):
                return 4 # graupel
            else:
                return 3 # rimed aggregates
    else:
        if (aggmass < 0.1 * icemass):
            return 2 # pristine crystals
        else:
            return 1 # aggregates

def Iice_shape(a_axis,c_axis,d_axis):
    if (c_axis > a_axis):
        return 3 # columnar crystals
    else:
        if (d_axis > 2.0/3.0 * a_axis):
            return 2 # dendrites
        else:
            return 1 # hexagonal plates

def Iice_shapePoly(a_axis,c_axis,d_axis,ag_axis,cg_axis,ex_cry):
    if (ex_cry < 0.5):
        if (c_axis > a_axis):
            return 3 # columnar crystals
        else:
            if (d_axis > 2.0/3.0 * a_axis):
                return 2 # dendrites
            else:
                return 1 # hexagonal plates
    else:
        if (a_axis == 0.0 and c_axis == 0.0):
            return -1
        elif (a_axis == 0.0):
            temp = 1.0
        elif (c_axis == 0.0):
            temp = -1.0
        else:
            temp = ag_axis/a_axis - cg_axis/c_axis
        if (temp >= 0.5):
            return 4 # planar polycrystals
        elif (temp <= -0.5):
            return 5 # columnar polycrystals
        else:
            return 6 # irregular polycrystals


print("LOG -- FINISH READING NETCDF DIMENSIONS")


### INITIALIZE GLOBAL VARIABLES
# liquid water path
LWP = []
# ice water path
IWP = []
# liquid number concentration
Nliq_timeseries = []
# ice number concentration
Nice_timeseries = []
# liquid effective diameter
Dliq_timeseries = []
# ice effective diameter (maximum dimension)
Dice_timeseries = []
# surface downwelling longwave flux
LW = []
# surface downwelling shortwave flux
SW = []
# precipitation
prec =[]
# temperature profile
T1d = np.zeros(zlim)
ini_T1d = np.zeros(zlim)
# qliq profile
qliq1d = np.zeros(zlim)
ini_qliq1d = np.zeros(zlim)
# qice profile
qice1d = np.zeros(zlim)
ini_qice1d = np.zeros(zlim)
# Nliq profile
Nliq1d = np.zeros(zlim)
ini_Nliq1d = np.zeros(zlim)
# Nice profile
Nice1d = np.zeros(zlim)
ini_Nice1d = np.zeros(zlim)
# Dliq profile
Dliq1d = np.zeros(zlim)
# Dice profile
Dice1d = np.zeros(zlim)
# ice type profile
type_agg1d = np.zeros(zlim)
type_rag1d = np.zeros(zlim)
type_gra1d = np.zeros(zlim)
type_cry1d = np.zeros(zlim)
type_rcr1d = np.zeros(zlim)
# ice habit profile
habit_pla1d = np.zeros(zlim)
habit_col1d = np.zeros(zlim)
habit_den1d = np.zeros(zlim)
habit_ppl1d = np.zeros(zlim)
habit_pcl1d = np.zeros(zlim)
habit_pir1d = np.zeros(zlim)
# qv profile
qv1d = np.zeros(zlim)
ini_qv1d = np.zeros(zlim)
# wind speed  profile
v1d = np.zeros(zlim)
ini_v1d = np.zeros(zlim)
# RH profile
rh1d = np.zeros(zlim)
ini_rh1d = np.zeros(zlim)
# TKE profile
TKE1d = np.zeros(zlim)

print("LOG -- START READING NETCDF VARIABLES AND ANALYSIS")

for t in (range(tlim)):
    print("LOG -- PROGRESS",t)

    ### READ THERMODYNAMICS VARIABLES
    # density
    dens = f.variables['DENS'][t,:,:,0]

    # temperature
    temperature = f.variables['T'][t,:,:,0]

    # horizontal wind speed
    v = f.variables['V'][t,:,:,0]

    # verticle wind speed
    w = f.variables['W'][t,:,:,0]
    
    # vapor mixing ratio
    qv = f.variables['QV'][t,:,:,0]

    # RH
    rh = f.variables['RH'][t,:,:,0]

    # liquid mixing ratio
    qliq = f.variables['QLIQ'][t,:,:,0]

    # all mixing ratio
    qhyd = f.variables['QHYD'][t,:,:,0]

    # ice mixing ratio
    qice = qhyd - qliq

    # liquid number concentration
    bin_liq_con = np.zeros((len(massbin),zlim,ylim))
    for j in range(len(massbin)):
        variable_name = 'con_lmass'+str(j+1)
        bin_liq_con[j] = f.variables[variable_name][t,:,:,0]

    # ice number concentration
    bin_ice_con = np.zeros((len(massbini),zlim,ylim))
    for j in range(len(massbini)):
        variable_name = 'con_imass'+str(j+1)
        bin_ice_con[j] = f.variables[variable_name][t,:,:,0]

    # ice maximum dimension
    bin_ice_Dmax = np.zeros((len(massbini),zlim,ylim))
    for j in range(len(massbini)):
        variable_name = 'vol_imass'+str(j+1)
        bin_ice_Dmax[j] = f.variables[variable_name][t,:,:,0]

    Nliq = np.zeros((zlim,ylim))
    Nice = np.zeros((zlim,ylim))
    Dliq = np.zeros((zlim,ylim))
    Dice = np.zeros((zlim,ylim))

    if (t == tlim - 1):
        # a axis
        bin_ice_alen = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'a_axis'+str(j+1)
            bin_ice_alen[j] = f.variables[variable_name][t,:,:,0]
            
        # c axis
        bin_ice_clen = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'c_axis'+str(j+1)
            bin_ice_clen[j] = f.variables[variable_name][t,:,:,0]

        # d axis
        bin_ice_dlen = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'd_axis'+str(j+1)
            bin_ice_dlen[j] = f.variables[variable_name][t,:,:,0]

        if (poly_exist):
            # ag axis
            bin_ice_aglen = np.zeros((len(massbini),zlim,ylim))
            for j in range(len(massbini)):
                variable_name = 'ag_axis'+str(j+1)
                bin_ice_aglen[j] = f.variables[variable_name][t,:,:,0]
            
            # cg axis
            bin_ice_cglen = np.zeros((len(massbini),zlim,ylim))
            for j in range(len(massbini)):
                variable_name = 'cg_axis'+str(j+1)
                bin_ice_cglen[j] = f.variables[variable_name][t,:,:,0]

            # nex
            bin_ice_nex = np.zeros((len(massbini),zlim,ylim))
            for j in range(len(massbini)):
                variable_name = 'ex_cry'+str(j+1)
                bin_ice_nex[j] = f.variables[variable_name][t,:,:,0]

        # ice mass
        bin_ice_mass = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'ice_imass'+str(j+1)
            bin_ice_mass[j] = f.variables[variable_name][t,:,:,0]
            
        # crystal mass
        bin_ice_mcry = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'cry_imass'+str(j+1)
            bin_ice_mcry[j] = f.variables[variable_name][t,:,:,0]

        # agregate mass
        bin_ice_magg = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'agg_imass'+str(j+1)
            bin_ice_magg[j] = f.variables[variable_name][t,:,:,0]

        # riming mass
        bin_ice_mrim = np.zeros((len(massbini),zlim,ylim))
        for j in range(len(massbini)):
            variable_name = 'rim_imass'+str(j+1)
            bin_ice_mrim[j] = f.variables[variable_name][t,:,:,0]

    # precipitation
    sfc_prec = f.variables['PREC'][t,:,0]

    # surface downwelling longwave flux
    sfc_LW_dn = f.variables['SFLX_LW_dn'][t,:,0]

    # surface downwelling shortwave flux
    sfc_SW_dn = f.variables['SFLX_SW_dn'][t,:,0]


    ### COMPUTE MEAN VARIABLES
    # liquid and ice water paths
    temp1, temp2 = 0.0, 0.0
    for i in range(ylim):
        temp3, temp4 = 0.0, 0.0
        for k in range(zlim):
            temp3 = temp3 + qliq[k][i]*dens[k][i]*(fz[k+3]-fz[k+2])*1000.0 # kg -> g
            temp4 = temp4 + qice[k][i]*dens[k][i]*(fz[k+3]-fz[k+2])*1000.0 # kg -> g
        temp1 = temp1 + temp3
        temp2 = temp2 + temp4
    LWP.append(temp1/float(ylim))
    IWP.append(temp2/float(ylim))
    
    # liquid and ice number concentrations and effective diameter
    temp1, temp2, temp3 = 0.0, 0.0, 0.0
    temp4, temp5, temp6 = 0.0, 0.0, 0.0
    for i in range(ylim):
        for k in range(zlim):
            D2, D3 = 0.0, 0.0
            for j in range(len(massbin)):
                bin_liq_con[j][k][i] = bin_liq_con[j][k][i]*dens[k][i]
                if (qliq[k][i] > 0.000001): # if in clouds
                    temp1 = temp1 + bin_liq_con[j][k][i]
                Nliq[k][i] = Nliq[k][i] + bin_liq_con[j][k][i]
                D3 = D3 + bin_liq_con[j][k][i]*diameterbin[j]**3
                D2 = D2 + bin_liq_con[j][k][i]*diameterbin[j]**2
            if (D2 != 0.0):
                Dliq[k][i] = D3/D2
                if (qliq[k][i] > 0.000001): # if in clouds
                    temp2 = temp2 + D3/D2
            if (qliq[k][i] > 0.000001): # if in clouds
                temp3 = temp3 + 1.0

            D3, D2 = 0.0, 0.0
            for j in range(len(massbini)):
                bin_ice_con[j][k][i] = bin_ice_con[j][k][i]*dens[k][i]
                if (qice[k][i] > 1.e-11):
                    temp4 = temp4 + bin_ice_con[j][k][i]*1000.0
                Nice[k][i] = Nice[k][i] + bin_ice_con[j][k][i]*1000.0
                if (bin_ice_con[j][k][i] > 1.e-11): # 10e-8 L-1 limit
                    local_Dmax = (bin_ice_Dmax[j][k][i]/bin_ice_con[j][k][i]*6.0/math.pi)**(1.0/3.0)*10000.0 # cm -> um
                else:
                    local_Dmax = 0.0
                D3 = D3 + bin_ice_con[j][k][i]*local_Dmax**3
                D2 = D2 + bin_ice_con[j][k][i]*local_Dmax**2
            if (D2 != 0.0):
                Dice[k][i] = D3/D2
                if (qice[k][i] > 1.e-11):
                    temp5 = temp5 + D3/D2
            if (qice[k][i] > 1.e-11):
                temp6 = temp6 + 1.0
    if (temp3 != 0.0):
        Nliq_timeseries.append(temp1/temp3)
        Dliq_timeseries.append(temp2/temp3)
    else:
        Dliq_timeseries.append(0.0)
        Nliq_timeseries.append(0.0)
    if (temp6 != 0.0):
        Nice_timeseries.append(temp4/temp6)
        Dice_timeseries.append(temp5/temp6)
    else:
        Nice_timeseries.append(0.0)
        Dice_timeseries.append(0.0)

    temp1, temp2, temp3 = 0.0, 0.0, 0.0
    for i in range(ylim):
        temp1 = temp1 + sfc_prec[i]*3600.0*24.0 # 1kg = 1mm, mm / Day
        temp2 = temp2 + sfc_LW_dn[i]
        temp3 = temp3 + sfc_SW_dn[i]
    prec.append(temp1/float(ylim))
    LW.append(temp2/float(ylim))
    SW.append(temp3/float(ylim))

    # vertical profiles
    if (t == 0):
        for k in range(zlim):
            ini_T1d[k] = ini_T1d[k] + sum(temperature[k][:])/float(ylim) - 273.15
            ini_qliq1d[k] = ini_qliq1d[k] + sum(qliq[k][:])/float(ylim)*1000.0
            ini_qice1d[k] = ini_qice1d[k] + sum(qice[k][:])/float(ylim)*1000.0
            ini_Nliq1d[k] = ini_Nliq1d[k] + sum(Nliq[k][:])/float(ylim)
            ini_Nice1d[k] = ini_Nice1d[k] + sum(Nice[k][:])/float(ylim)*1000.0
            ini_qv1d[k] = ini_qv1d[k] + sum(qv[k][:])/float(ylim)*1000.0
            ini_v1d[k] = ini_v1d[k] + sum(v[k][:])/float(ylim)
            ini_rh1d[k] = ini_rh1d[k] + sum(rh[k][:])/float(ylim)

    if (t == tlim - 1):
        for k in range(zlim):
            T1d[k] = sum(temperature[k][:])/float(ylim) - 273.15
            qliq1d[k] = sum(qliq[k][:])/float(ylim)*1000.0
            qice1d[k] = sum(qice[k][:])/float(ylim)*1000.0
            qv1d[k] = sum(qv[k][:])/float(ylim)*1000.0
            v1d[k] = sum(v[k][:])/float(ylim)
            rh1d[k] = sum(rh[k][:])/float(ylim)

            Hagg, Hcry, Hrcr, Hgra, Hrag = 0.0, 0.0, 0.0, 0.0, 0.0
            Hpla, Hcol, Hden, Hppl, Hpcl, Hpir = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            HDliq, HDice = 0.0, 0.0
            HNliq, HNice = 0.0, 0.0
            for i in range(ylim):
                TKE1d[k] = TKE1d[k] + dens[k][i]*(w[k][i]**2)/float(ylim)

                # nliq and Dliq vertical profiles
                if (qliq[k][i] > 0.000001): # in clouds
                    HNliq = HNliq + Nliq[k][i]
                    HDliq = HDliq + Dliq[k][i]*Nliq[k][i]
                    
                for j in range(len(massbini)):
                    if (bin_ice_con[j][k][i] > 1.e-11):
                        # rescale axis lengths
                        bin_ice_alen[j][k][i] = (bin_ice_alen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_clen[j][k][i] = (bin_ice_clen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_dlen[j][k][i] = (bin_ice_dlen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        if (poly_exist):
                            if (iceaxis_version == 1):
                                bin_ice_aglen[j][k][i] = (bin_ice_aglen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                                bin_ice_cglen[j][k][i] = (bin_ice_cglen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                            elif (iceaxis_version == 3):
                                bin_ice_aglen[j][k][i] = (bin_ice_aglen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)*bin_ice_alen[j][k][i]
                                bin_ice_cglen[j][k][i] = (bin_ice_cglen[j][k][i]*dens[k][i]/bin_ice_con[j][k][i])**(1.0/3.0)*bin_ice_clen[j][k][i]
                            bin_ice_nex[j][k][i] = bin_ice_nex[j][k][i]*dens[k][i]/bin_ice_con[j][k][i]
                        habit = Iice_habit(bin_ice_mass[j][k][i],bin_ice_mcry[j][k][i],bin_ice_magg[j][k][i],bin_ice_mrim[j][k][i])
                        if (habit == 1):
                            Hagg = Hagg + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 2):
                            Hcry = Hcry + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 3):
                            Hrag = Hrag + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 4):
                            Hgra = Hgra + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 5):
                            Hrcr = Hrcr + bin_ice_con[j][k][i]*1000.0

                        if (poly_exist):
                            habit = Iice_shapePoly(bin_ice_alen[j][k][i],bin_ice_clen[j][k][i],bin_ice_dlen[j][k][i],bin_ice_aglen[j][k][i],bin_ice_cglen[j][k][i],bin_ice_nex[j][k][i])
                            if (habit == 1):
                                Hpla = Hpla + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 2):
                                Hden = Hden + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 3):
                                Hcol = Hcol + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 4):
                                Hppl = Hppl + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 5):
                                Hpcl = Hpcl + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 6):
                                Hpir = Hpir + bin_ice_con[j][k][i]*1000.0
                        else:
                            habit = Iice_shape(bin_ice_alen[j][k][i],bin_ice_clen[j][k][i],bin_ice_dlen[j][k][i])
                            if (habit == 1):
                                Hpla = Hpla + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 2):
                                Hden = Hden + bin_ice_con[j][k][i]*1000.0
                            elif (habit == 3):
                                Hcol = Hcol + bin_ice_con[j][k][i]*1000.0
                HNice = HNice + Nice[k][i]
                HDice = HDice + Dice[k][i]*Nice[k][i]

            Nliq1d[k] = HNliq/float(ylim)
            Nice1d[k] = HNice/float(ylim)
            if (HNliq != 0.0):
                Dliq1d[k] = HDliq/HNliq
            if (HNice != 0.0):
                Dice1d[k] = HDice/HNice
            type_agg1d[k] = Hagg/float(ylim)
            type_gra1d[k] = Hgra/float(ylim)
            type_cry1d[k] = Hcry/float(ylim)
            type_rag1d[k] = Hrag/float(ylim)
            type_rcr1d[k] = Hrcr/float(ylim)
            habit_pla1d[k] = Hpla/float(ylim)
            habit_col1d[k] = Hcol/float(ylim)
            habit_den1d[k] = Hden/float(ylim)
            if (poly_exist):
                habit_ppl1d[k] = Hppl/float(ylim)
                habit_pcl1d[k] = Hpcl/float(ylim)
                habit_pir1d[k] = Hpir/float(ylim)
            
                    

    # COMPUTE LIQUID AND ICE SIZE DISTRIBUTIONS AT 3 HOURS
    if (t == tlim - 1):
        ### OPEN FRIDLIND DATA
        fridlind_fdir1 = 'may7_2dcp_ice.mean_values.txt'
        fridlind_fdir2 = 'may7_2dcp_drops.mean_values.txt'

        i = 0
        fridlind_nidata = 0
        fridlind_idiam, fridlind_idata = [], []
        with open(fridlind_fdir1,'r') as f:
            for line in f:
                data = line.split()
                if (i == 0):
                    fridlind_nidata = int(data[0])
                elif(i <= fridlind_nidata):
                    fridlind_idiam.append(float(data[0]))
                    fridlind_idata.append(float(data[1]))
                i = i + 1

        i = 0
        fridlind_nldata = 0
        fridlind_ldiam, fridlind_ldata = [], []
        with open(fridlind_fdir2,'r') as f:
            for line in f:
                data = line.split()
                if (i == 0):
                    fridlind_nldata = int(data[0])
                elif(i <= fridlind_nldata):
                    fridlind_ldiam.append(float(data[0]))
                    fridlind_ldata.append(float(data[1]))
                i = i + 1

        for i in range(len(fridlind_ldata)):
            if (fridlind_ldata[i] < 1.e-8):
                fridlind_ldata[i] = float('nan')
        for i in range(len(fridlind_idata)):
            if (fridlind_idata[i] < 1.e-8):
                fridlind_idata[i] = float('nan')

        # liquid size distribution
        dNdlogD, NP = np.zeros((max_binr,zlim,ylim)), np.zeros((max_binr,zlim,ylim))
        mean_dNdlogD = np.zeros(max_binr)
        median_dNdlogD = np.zeros(max_binr)

        for k in range(zlim):
            for i in range(ylim):
                if (qliq[k][i] > 0.000001):
                    NP[0][k][i] = bin_liq_con[0][k][i]
                    for j in range(max_binr-1):
                        NP[j+1][k][i] = NP[j][k][i] + bin_liq_con[j+1][k][i]
        
        for k in range(zlim):
            for i in range(ylim):
                if (qliq[k][i] > 0.000001):
                    j = 0
                    dNdlogD[j][k][i] = (NP[j+1][k][i] - NP[j][k][i])/(math.log10(diameterbin[j+1]) - math.log10(diameterbin[j]))
                    j = j + 1
                    while (j < max_binr-1):
                        dNdlogD[j][k][i] = (NP[j+1][k][i] - NP[j-1][k][i])/(math.log10(diameterbin[j+1]) - math.log10(diameterbin[j-1]))
                        j = j + 1
                    dNdlogD[j][k][i] = (NP[j][k][i] - NP[j-1][k][i])/(math.log10(diameterbin[j]) - math.log10(diameterbin[j-1]))

        for j in range(max_binr):
            temp = 0.0
            median_array = []
            for k in range(zlim):
                for i in range(ylim):
                    if (qliq[k][i] > 0.000001):
                        mean_dNdlogD[j] = mean_dNdlogD[j] + dNdlogD[j][k][i]
                        median_array.append(dNdlogD[j][k][i])
                        temp = temp + 1.0
            if (temp != 0.0):
                mean_dNdlogD[j] = mean_dNdlogD[j]/temp
                quickSorti(median_array,0,len(median_array)-1)
                if (len(median_array)%2 == 0):
                    median_dNdlogD[j] = (0.5*(median_array[int(len(median_array)/2)-1] + median_array[int(len(median_array)/2)]))
                else:
                    median_dNdlogD[j] = (median_array[int((len(median_array)-1)/2)])
    
        # ice size distribution
        Dibin_ice_con = np.zeros((len(diameteribin),zlim,ylim))
        dNidlogD, NPi = np.zeros((len(diameteribin),zlim,ylim)), np.zeros((len(diameteribin),zlim,ylim))
        mean_dNidlogD, median_dNidlogD = np.zeros(len(diameteribin)), np.zeros(len(diameteribin))

        ice_con = np.zeros((zlim,ylim))
        for k in range(zlim):
            for i in range(ylim):
                if (qice[k][i] > 1.e-11):
                    for j in range(len(massbini)):
                        if (bin_ice_con[j][k][i] > 1.e-11):
                            ivol_temp = (bin_ice_Dmax[j][k][i]*dens[k][i]/bin_ice_con[j][k][i]*6.0/math.pi)**(1.0/3.0)*10000.0
                            for jj in range(len(diameteribin)-1):
                                if (ivol_temp >= diameteribin_bnd[jj] and ivol_temp < diameteribin_bnd[jj+1]):
                                    Dibin_ice_con[jj][k][i] = Dibin_ice_con[jj][k][i] + bin_ice_con[j][k][i]
                                    break
                        ice_con[k][i] = ice_con[k][i] + bin_ice_con[j][k][i]

        for k in range(zlim):
            for i in range(ylim):
                if (ice_con[k][i] > 1.e-8):
                    NPi[0][k][i] = Dibin_ice_con[0][k][i]
                    for j in range(len(diameteribin)-1):
                        NPi[j+1][k][i] = NPi[j][k][i] + Dibin_ice_con[j+1][k][i]

        for k in range(zlim):
            for i in range(ylim):
                if (ice_con[k][i] > 1.e-8 and z[k] <= 280.0):
                    j = 0
                    dNdlogD[j][k][i] = (NPi[j+1][k][i] - NPi[j][k][i])/(math.log10(diameteribin[j+1]) - math.log10(diameteribin[j]))
                    j = j + 1
                    while (j < len(diameteribin)-1):
                        dNdlogD[j][k][i] = (NPi[j+1][k][i] - NPi[j-1][k][i])/(math.log10(diameteribin[j+1]) - math.log10(diameteribin[j-1]))
                        j = j + 1
                    dNdlogD[j][k][i] = (NPi[j][k][i] - NPi[j-1][k][i])/(math.log10(diameteribin[j]) - math.log10(diameteribin[j-1]))

        for j in range(len(diameteribin)):
            temp = 0.0
            median_array = []
            for k in range(zlim):
                for i in range(ylim):
                    if (ice_con[k][i] > 1.e-8 and z[k] <= 280.0):
                        mean_dNidlogD[j] = mean_dNidlogD[j] + dNdlogD[j][k][i]
                        median_array.append(dNdlogD[j][k][i])
                        temp = temp + 1.0
            if (temp != 0.0):
                mean_dNidlogD[j] = mean_dNidlogD[j]/temp
                quickSorti(median_array,0,len(median_array)-1)
                if (len(median_array)%2 == 0):
                    median_dNidlogD[j] = (0.5*(median_array[int(len(median_array)/2)-1] + median_array[int(len(median_array)/2)]))
                else:
                    median_dNidlogD[j] = (median_array[int((len(median_array)-1)/2)])
        
        plt.close()
        plt.clf()
        fig = plt.figure(figsize=(10.5,6))
        plt.plot(diameteribin,mean_dNidlogD,color="red")
        #plt.plot(diameteribin,mean_dNidlogD,color="red",linestyle='dashed')
        plt.xlabel("Diameter ($\mu$m)",fontsize=18)
        plt.ylabel("dN/dlogD (cm$^{-3}$)",fontsize=18)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(0.1,10000.0)
        plt.ylim(1.e-8,10000.0)
        plt.title('Liquid and ice particle size distribution at 1.25 hr',fontsize=18)
        plt.tick_params(axis='both',which='major',labelsize=18)

        plt.plot(diameterbin,mean_dNdlogD,color="blue")
        #plt.plot(diameterbin,median_dNdlogD,color="blue",linestyle='dashed')
        
        plt.plot(fridlind_idiam,fridlind_idata,color="black",linestyle="solid")
        plt.plot(fridlind_ldiam,fridlind_ldata,color="black",linestyle="solid")

        plt.savefig(plot_dir+'figure_distr.png',format='png',dpi=300)

        
        ### plot 2d qliq and qice (g kg-1)
        bluemap = plt.get_cmap('Blues')
        redmap = plt.get_cmap('Reds')

        # contour plot function
        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        y2d = np.zeros((zlim,ylim))
        z2d = np.zeros((zlim,ylim))
        for k in range(zlim):
            for i in range(ylim):
                qliq[k][i] = qliq[k][i]*1000.0
                qice[k][i] = qice[k][i]*1000.0
                y2d[k][i] = y[i]/1000.0
                z2d[k][i] = z[k]/1000.0

        plt.close()
        fig, axs = plt.subplots(ncols=2,figsize=(18,9))
        fig.subplots_adjust(wspace=0.2)

        qliq_boundaries = np.linspace(0,0.06,100)
        qliq_tboundaries = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
        qliq_norm = colors.BoundaryNorm(qliq_boundaries, bluemap.N, clip=True)

        qice_boundaries = np.linspace(0,0.4e-2,100)
        qice_tboundaries = [0, 0.1e-2, 0.2e-2, 0.3e-2, 0.4e-2]
        qice_norm = colors.BoundaryNorm(qice_boundaries, redmap.N, clip=True)

        #cboundaries = np.linspace(-1,1,11)
        cboundaries = [-0.75,-0.25,0.25,0.75]

        cp1 = axs[0].contourf(y2d,z2d,qliq,cmap=bluemap,levels=qliq_boundaries,norm=qliq_norm)
        cp2 = axs[1].contourf(y2d,z2d,qice,cmap=redmap,levels=qice_boundaries,norm=qice_norm)
        axs[0].contour(y2d,z2d,w,colors='black',levels=cboundaries)
        #axs[1].contour(y,z,ww,colors='black',levels=cboundaries)

        fontsize = 18

        divider1 = make_axes_locatable(axs[0])
        cax1 = divider1.append_axes("right", size="5%", pad=0.1)
        cb1 = plt.colorbar(cp1,ticks=qliq_tboundaries,cax=cax1)
        divider2 = make_axes_locatable(axs[1])
        cax2 = divider2.append_axes("right", size="5%", pad=0.1)
        cb2 = plt.colorbar(cp2,ticks=qice_tboundaries,cax=cax2,format=ticker.FuncFormatter(fmt))
        #cb1 = plt.colorbar(cp1,format=ticker.FuncFormatter(fmt),ticks=nice_tboundaries)
        cb1.ax.tick_params(labelsize=fontsize)
        cb2.ax.tick_params(labelsize=fontsize)

        axs[0].set_xlim(0,2)
        axs[0].set_ylim(0,1)
        axs[0].set_xticks([0,0.5,1,1.5,2])
        axs[0].set_yticks([0,0.25,0.5,0.75,1])
        axs[0].set_xlabel("y (km)",fontsize=fontsize)
        axs[0].set_ylabel("z (km)",fontsize=fontsize)
        axs[0].tick_params(axis='both',which='major',labelsize=fontsize)
        axs[0].set_title("Liquid water content (g kg$^{-3}$)",fontsize=fontsize)
        
        axs[1].set_xlim(0,2)
        axs[1].set_ylim(0,1)
        axs[1].set_xticks([0,0.5,1,1.5,2])
        axs[1].set_yticks([0,0.25,0.5,0.75,1])
        axs[1].set_xlabel("y (km)",fontsize=fontsize)
        #axs[1].set_ylabel("z (km)",fontsize=fontsize)
        axs[1].axes.yaxis.set_ticklabels([])
        axs[1].tick_params(axis='both',which='major',labelsize=fontsize)
        axs[1].set_title("Ice water content (g kg$^{-3}$)",fontsize=fontsize)

        #axs[0].set_aspect(1.0/axs[0].get_data_ratio(),adjustable='box')
        #axs[1].set_aspect(1.0/axs[1].get_data_ratio(),adjustable='box')
        plt.savefig(plot_dir+'figure_qliq_qice2d.png',format='png',dpi=300)

        
# convert m to km and second to hour
for i, temp in enumerate(z):
    z[i] = z[i]/1000.0
for i, temp in enumerate(time):
    time[i] = time[i]/3600.0


f.close()
print("LOG -- FINISH READING NETCDF")

print("LOG -- START PLOTTING")

# read sonde data
sonde1115_temperature, sonde1115_rh ,sonde1115_wtot = [], [], []
sonde_z, sonde1730_temperature, sonde1730_rh ,sonde1730_wtot = [], [], [], []
sonde2335_temperature, sonde2335_rh ,sonde2335_wtot = [], [], []
with open('sonde_data.dat','r') as f:
    i = 0
    for line in f:
        data = line.split()
        sonde_z.append(float(data[0])/1000.0)
        if (float(data[1]) < -30.0):
            sonde1115_temperature.append(float('nan'))
        else:
            sonde1115_temperature.append(float(data[1]))
        sonde1115_rh.append(float(data[2])*100.0)
        if (float(data[0])/1000.0 < 0.2 or float(data[3]) < 0.0):
            sonde1115_wtot.append(float('nan'))
        else:
            sonde1115_wtot.append(float(data[3]))
            
        sonde1730_temperature.append(float(data[4]))
        sonde1730_rh.append(float(data[5])*100.0)
        if (float(data[0])/1000.0 < 0.2 or float(data[6]) < 0.0):
            sonde1730_wtot.append(float('nan'))
        else:
            sonde1730_wtot.append(float(data[6]))
        
        sonde2335_temperature.append(float(data[7]))
        sonde2335_rh.append(float(data[8])*100.0)
        if (float(data[0])/1000.0 < 0.2 or float(data[9]) < 0.0):
            sonde2335_wtot.append(float('nan'))
        else:
            sonde2335_wtot.append(float(data[9]))

plt.close()
plt.clf()
fig = plt.figure(figsize=(20,8.5))
gs1 = gridspec.GridSpec(1, 5)
gs1.update(wspace=0.4, hspace=0.1)

ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])
ax4 = plt.subplot(gs1[3])
ax5 = plt.subplot(gs1[4])

init_color = 'black'
init_style = 'dashed'
styles = 'solid'
colors = 'red'
labels = '1.25 hr'
linewidths = 2.5

ax1.plot(T1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
ax1.plot(ini_T1d,z,color=init_color,linestyle=init_style,label="Init.")
ax1.fill_betweenx(sonde_z,sonde1730_temperature,sonde2335_temperature,color="gray", alpha=0.4)
ax1.set_xlim(-22,-14)
ax1.set_ylim(0,1)
ax1.xaxis.set_ticks([-22,-18,-14])
ax1.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1.0])
ax1.grid(True)
ax1.set_xlabel("T (C)",fontsize=22)
ax1.set_ylabel("Height (km)",fontsize=22)
ax1.tick_params(axis='both',which='major',labelsize=22)
ax1.legend(bbox_to_anchor=(0.027,1),loc='upper left',prop={'size':20})

ax2.plot(qliq1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
ax2.plot(ini_qliq1d,z,color=init_color,linestyle=init_style,label="Init.")
ax2.set_xlim(-0.005,0.08)
ax2.set_ylim(0,1)
ax2.xaxis.set_ticks([0,0.05])
ax2.grid(True)
ax2.set_xlabel("QLIQ (g kg$^{-1}$)",fontsize=22)
ax2.axes.yaxis.set_ticklabels([])
ax2.tick_params(axis='both',which='major',labelsize=22)

ax3.plot(qice1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
ax3.plot(ini_qice1d,z,color=init_color,linestyle=init_style,label="Init.")
ax3.set_xlim(-0.000125,0.00205)
ax3.set_ylim(0,1)
ax3.xaxis.set_ticks([0,0.001,0.002])
ax3.grid(True)
ax3.set_xlabel("QICE (g kg$^{-1}$)",fontsize=22)
ax3.axes.yaxis.set_ticklabels([])
ax3.tick_params(axis='both',which='major',labelsize=22)
ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax3.xaxis.get_offset_text().set_fontsize(18)
ax3.get_xaxis().get_offset_text().set_position((1.2,0.05))

ax4.plot(rh1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
ax4.plot(ini_rh1d,z,color=init_color,linestyle=init_style,label="Init.")
ax4.fill_betweenx(sonde_z,sonde1730_rh,sonde2335_rh,color="gray", alpha=0.4)
ax4.set_xlim(40.0,120.0)
ax4.set_ylim(0,1)
ax4.grid(True)
ax4.set_xlabel("RH (liquid)",fontsize=22)
ax4.set_yticklabels([])
ax4.tick_params(axis='both',which='major',labelsize=22)

ax5.plot(TKE1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
ax5.set_ylim(0,1)
ax5.xaxis.set_ticks([0,0.5,1])
ax5.grid(True)
ax5.set_xlabel("TKE (kg m$^{-1}$ s$^{-2}$)",fontsize=22)
ax5.axes.yaxis.set_ticklabels([])
ax5.tick_params(axis='both',which='major',labelsize=22)


ax1.text(0.5,1.025, "(a)", size=25, ha="center",transform=ax1.transAxes)
ax2.text(0.5,1.025, "(b)", size=25, ha="center",transform=ax2.transAxes)
ax3.text(0.5,1.025, "(c)", size=25, ha="center",transform=ax3.transAxes)
ax4.text(0.5,1.025, "(d)", size=25, ha="center",transform=ax4.transAxes)
ax5.text(0.5,1.025, "(e)", size=25, ha="center",transform=ax5.transAxes)

plt.subplots_adjust(bottom=0.12,top=0.93,left=0.05,right=0.95)

plt.savefig(plot_dir+'Figure_vertical_profile.pdf',format='pdf',dpi=300)


fig.set_size_inches(6.5, 9)
plt.clf()
plt.plot(qv1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.plot(ini_qv1d,z,color=init_color,linestyle=init_style,label="Init.")
plt.xlabel("QV (g/kg)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(0.5,1.5)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Vapor mixing ratio profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_qv1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(v1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.plot(ini_v1d,z,color=init_color,linestyle=init_style,label="Init.")
plt.xlabel("Wind speed (m s$^{-1}$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(3.5,5.5)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal wind speed profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_u1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(Nliq1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.plot(ini_Nliq1d,z,color=init_color,linestyle=init_style,label="Init.")
plt.xlabel("Nliq (cm$^{-3}$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-10,250)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal liquid conc. profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_Nliq1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(Dliq1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.xlabel("Dliq ($\mu$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-1,10)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal liquid Dmax profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_Dliq1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(Nice1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.plot(ini_Nice1d,z,color=init_color,linestyle=init_style,label="Init.")
plt.xlabel("Nice (L$^{-1}$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-0.02,0.5)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal ice conc. profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_Nice1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(Dice1d,z,linestyle=styles,color=colors,label=labels,linewidth=linewidths)
plt.xlabel("Dice ($\mu$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-20,450)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal ice Dmax profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.savefig(plot_dir+'Figure_Dice1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(type_cry1d,z,linestyle=styles,color='blue',label='Crystal',linewidth=linewidths)
plt.plot(type_agg1d,z,linestyle=styles,color='red',label='Aggregate',linewidth=linewidths)
plt.plot(type_gra1d,z,linestyle=styles,color='black',label='Graupel',linewidth=linewidths)
plt.plot(type_rag1d,z,linestyle=styles,color='purple',label='Rimed gra.',linewidth=linewidths)
plt.plot(type_rcr1d,z,linestyle=styles,color='green',label='Rimed cry.',linewidth=linewidths)
plt.xlabel("Conc. (L$^{-1}$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-0.02,0.5)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal ice type profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.legend(prop={'size':18})
plt.tight_layout()
plt.savefig(plot_dir+'Figure_icetype1d.pdf',format='pdf',dpi=300)

plt.clf()
plt.plot(habit_pla1d,z,linestyle=styles,color='blue',label='Plate',linewidth=linewidths)
plt.plot(habit_col1d,z,linestyle=styles,color='red',label='Column',linewidth=linewidths)
plt.plot(habit_den1d,z,linestyle=styles,color='green',label='dendrite',linewidth=linewidths)
plt.plot(habit_ppl1d,z,linestyle='dashed',color='blue',label='Poly plate',linewidth=linewidths)
plt.plot(habit_pcl1d,z,linestyle='dashed',color='red',label='Poly column',linewidth=linewidths)
plt.plot(habit_pir1d,z,linestyle='dashed',color='green',label='Poly irregular',linewidth=linewidths)
plt.xlabel("Conc. (L$^{-1}$)",fontsize=18)
plt.ylabel("Height (km)",fontsize=18)
plt.xlim(-0.02,0.5)
plt.ylim(0,1)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.title('Horizontal ice habit profile',fontsize=18)
plt.tick_params(axis='both',which='major',labelsize=18)
plt.grid(True)
plt.tight_layout()
plt.legend(prop={'size':18})
plt.savefig(plot_dir+'Figure_icehabit1d.pdf',format='pdf',dpi=300)

obs_lwp1, obs_lwp2 = [], []
obs_nliq1, obs_nice2 = [], []
obs_nice1, obs_nliq2 = [], []
obs_sw1, obs_sw2 = [], []
obs_lw1, obs_lw2 = [], []
for temp in time:
    obs_lwp1.append(4.0)
    obs_lwp2.append(10.0)
    obs_nliq1.append(200.0)
    obs_nliq2.append(230.0)
    obs_nice1.append(0.17)
    obs_nice2.append(1.7)
    obs_sw1.append(480.4)
    obs_sw2.append(503.3)
    obs_lw1.append(210.2)
    obs_lw2.append(218.8)

plt.clf()
fig = plt.figure(figsize=(18.75,24))
gs1 = gridspec.GridSpec(5, 1)
gs1.update(wspace=0.0, hspace=0.5)

ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])
ax4 = plt.subplot(gs1[3])
ax5 = plt.subplot(gs1[4])

fontsize = 26
labelsize = 26
textsize = 26

ax1.plot(time,Nice_timeseries,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax1.fill_between(time,obs_nice1,obs_nice2,color="gray", alpha=0.4)
ax1.set_yscale("log")
ax1.set_xlim(0.0,1.25)
ax1.set_ylim(1.e-3,10)
ax1.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax1.set_ylabel("Nice (L$^{-1}$)",fontsize=fontsize)
ax1.axes.xaxis.set_ticklabels([])
ax1.tick_params(axis='both',which='major',labelsize=labelsize)

ax2.plot(time,Nliq_timeseries,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax2.set_xlim(0.0,1.25)
ax2.set_ylim(0,250)
ax2.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax2.axes.xaxis.set_ticklabels([])
ax2.set_ylabel("Nliq (cm$^{-3}$)",fontsize=fontsize)
ax2.tick_params(axis='both',which='major',labelsize=labelsize)

ax3.plot(time,LWP,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax3.fill_between(time,obs_lwp1,obs_lwp2,color="gray", alpha=0.4)
ax3.set_xlim(0.0,1.25)
ax3.set_ylim(0,13)
ax3.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax3.axes.xaxis.set_ticklabels([])
ax3.set_ylabel("LWP (g m$^{-2}$)",fontsize=fontsize)
ax3.tick_params(axis='both',which='major',labelsize=labelsize)

ax4.plot(time,IWP,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax4.set_xlim(0.0,1.25)
ax4.set_ylim(1.e-5,10)
ax4.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax4.axes.xaxis.set_ticklabels([])
ax4.set_yscale("log")
ax4.set_ylabel("IWP (g m$^{-2}$)",fontsize=fontsize)
ax4.tick_params(axis='both',which='major',labelsize=labelsize)

ax5.plot(time,prec,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax5.set_xlim(0.0,1.25)
ax5.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax5.set_ylabel("Prep. (g m$^{-2}$ hr$^{-1}$)",fontsize=20)
ax5.set_ylabel("Prep. (mm hr$^{-1}$)",fontsize=fontsize)
ax5.set_xlabel("Time (hr)",fontsize=fontsize)
ax5.tick_params(axis='both',which='major',labelsize=labelsize)
ax5.xaxis.set_label_coords(0.5,-0.25)

legend = ax5.legend(ncol=4,bbox_to_anchor=(0.16,-0.45),loc='upper left',prop={'size':textsize})
frame = legend.get_frame()
frame.set_edgecolor('black')

ax1.text(0.5,1.1, "(a)", size=textsize, ha="center",transform=ax1.transAxes)
ax2.text(0.5,1.1, "(b)", size=textsize, ha="center",transform=ax2.transAxes)
ax3.text(0.5,1.1, "(c)", size=textsize, ha="center",transform=ax3.transAxes)
ax4.text(0.5,1.1, "(d)", size=textsize, ha="center",transform=ax4.transAxes)
ax5.text(0.5,1.1, "(e)", size=textsize, ha="center",transform=ax5.transAxes)

plt.subplots_adjust(bottom=0.13,top=0.93,left=0.1,right=0.9)
plt.savefig(plot_dir+'Figure_time_series1.pdf',format='pdf',dpi=300)


plt.clf()
fig = plt.figure(figsize=(18.75,17))
gs1 = gridspec.GridSpec(3, 1)
gs1.update(wspace=0.0, hspace=0.5)

ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])

ax1.plot(time,SW,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax1.fill_between(time,obs_sw1,obs_sw2,color="gray", alpha=0.4)
ax1.set_xlim(0.0,1.25)
ax1.set_ylim(350.0,550.0)
ax1.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax1.axes.xaxis.set_ticklabels([])
ax1.set_ylabel("SW flux (W m$^{-2}$)",fontsize=fontsize)
ax1.tick_params(axis='both',which='major',labelsize=labelsize)

ax2.plot(time,LW,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax2.fill_between(time,obs_lw1,obs_lw2,color="gray", alpha=0.4)
ax2.set_xlim(0.0,1.25)
ax2.set_ylim(200.0,230.0)
ax2.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax2.axes.xaxis.set_ticklabels([])
ax2.set_ylabel("LW flux (W m$^{-2}$)",fontsize=fontsize)
ax2.tick_params(axis='both',which='major',labelsize=labelsize)

ax3.plot(time,Dliq_timeseries,color=colors,linestyle=styles,label=labels,linewidth=linewidths)
ax3.set_xlim(0.0,1.25)
ax3.set_ylim(0.0,15)
ax3.set_xticks([0.0,0.25,0.5,0.75,1.0,1.25])
ax3.set_ylabel("D$_{eff}$ ($\mu$m)",fontsize=fontsize)
ax3.set_xlabel("Time (hr)",fontsize=fontsize)
ax3.tick_params(axis='both',which='major',labelsize=labelsize)
ax3.xaxis.set_label_coords(0.5,-0.25)

legend = ax3.legend(ncol=4,bbox_to_anchor=(0.16,-0.45),loc='upper left',prop={'size':textsize})
frame = legend.get_frame()
frame.set_edgecolor('black')

ax1.text(0.5,1.1, "(a)", size=textsize, ha="center",transform=ax1.transAxes)
ax2.text(0.5,1.1, "(b)", size=textsize, ha="center",transform=ax2.transAxes)
ax3.text(0.5,1.1, "(c)", size=textsize, ha="center",transform=ax3.transAxes)

plt.subplots_adjust(bottom=0.19,top=0.93,left=0.1,right=0.9)

plt.savefig(plot_dir+'Figure_time_series2.pdf',format='pdf',dpi=300)


print("LOG -- END")



