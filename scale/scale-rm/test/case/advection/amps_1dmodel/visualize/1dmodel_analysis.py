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


def Iice_types(icemass,crymass,aggmass,rimmass):
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
        if (isinstance(temp,complex)):
            print("Something wrong: ", ag_axis, a_axis, cg_axis, c_axis)
            sys.exit(0)
        if (temp >= 0.5):
            return 4 # planar polycrystals
        elif (temp <= -0.5):
            return 5 # columnar polycrystals
        else:
            return 6 # irregular polycrystals



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


print("LOG -- FINISH READING NETCDF DIMENSIONS")


### INITIALIZE GLOBAL VARIABLES
INITppv_con    = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_mass   = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_aaxis  = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_caxis  = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_agaxis = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_cgaxis = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INITppv_nex    = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_monopl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_monocl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_monodd = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_polypl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_polycl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
INIThabit_polyir = [[0.0 for _ in range(ylim)] for _ in range(zlim)]

ppv_con    = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_mass   = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_aaxis  = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_caxis  = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_agaxis = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_cgaxis = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
ppv_nex    = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_monopl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_monocl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_monodd = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_polypl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_polycl = [[0.0 for _ in range(ylim)] for _ in range(zlim)]
habit_polyir = [[0.0 for _ in range(ylim)] for _ in range(zlim)]

print("LOG -- START READING NETCDF VARIABLES AND ANALYSIS")

axis_version = 1

for t in (range(tlim)):
    if (t != tlim - 1 and t != 0):
        continue
    print("LOG -- PROGRESS",t)

    ### READ THERMODYNAMICS VARIABLES
    # density
    dens = f.variables['DENS'][t,:,:,0]

    # wind speed
    u = f.variables['U'][t,:,:,0]

    # ice crystal mass mixing ratio
    bin_ice_cry = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'cry_imass'+str(j+1)
        bin_ice_cry[j] = f.variables[variable_name][t,:,:,0]
        
    # ice number concentration
    bin_ice_con = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'con_imass'+str(j+1)
        bin_ice_con[j] = f.variables[variable_name][t,:,:,0]

    # ice maximum dimension
    bin_ice_Dmax = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'vol_imass'+str(j+1)
        bin_ice_Dmax[j] = f.variables[variable_name][t,:,:,0]

    # ice a-axis
    bin_ice_aaxis = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'a_axis'+str(j+1)
        bin_ice_aaxis[j] = f.variables[variable_name][t,:,:,0]

    # ice c-axis
    bin_ice_caxis = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'c_axis'+str(j+1)
        bin_ice_caxis[j] = f.variables[variable_name][t,:,:,0]

    # ice c-axis
    bin_ice_daxis = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'd_axis'+str(j+1)
        bin_ice_daxis[j] = f.variables[variable_name][t,:,:,0]

    # ice ag-axis
    bin_ice_agaxis = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'ag_axis'+str(j+1)
        bin_ice_agaxis[j] = f.variables[variable_name][t,:,:,0]

    # ice cg-axis
    bin_ice_cgaxis = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'cg_axis'+str(j+1)
        bin_ice_cgaxis[j] = f.variables[variable_name][t,:,:,0]

    # ice nex
    bin_ice_nex = [[[]] for _ in range(len(massbini))]
    for j in range(len(massbini)):
        variable_name = 'ex_cry'+str(j+1)
        bin_ice_nex[j] = f.variables[variable_name][t,:,:,0]


    ### COMPUTE MEAN VARIABLES
    # liquid and ice water paths
    if (axis_version == 3):
        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    if (bin_ice_con[j][k][i] > 1.e-12):
                        bin_ice_aaxis[j][k][i] = (bin_ice_aaxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_caxis[j][k][i] = (bin_ice_caxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_daxis[j][k][i] = (bin_ice_daxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_agaxis[j][k][i] = (bin_ice_agaxis[j][k][i]*bin_ice_aaxis[j][k][i]**3.0/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_cgaxis[j][k][i] = (bin_ice_cgaxis[j][k][i]*bin_ice_caxis[j][k][i]**3.0/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_nex[j][k][i] = bin_ice_nex[j][k][i]/bin_ice_con[j][k][i]
                    else:
                        bin_ice_aaxis[j][k][i] = 0.0
                        bin_ice_caxis[j][k][i] = 0.0
                        bin_ice_daxis[j][k][i] = 0.0
                        bin_ice_agaxis[j][k][i] = 0.0
                        bin_ice_cgaxis[j][k][i] = 0.0
                        bin_ice_nex[j][k][i] = 0.0
    elif (axis_version == 1):
        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    if (bin_ice_con[j][k][i] > 1.e-12):
                        bin_ice_aaxis[j][k][i] = (bin_ice_aaxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_caxis[j][k][i] = (bin_ice_caxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_daxis[j][k][i] = (bin_ice_daxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_agaxis[j][k][i] = (bin_ice_agaxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_cgaxis[j][k][i] = (bin_ice_cgaxis[j][k][i]/bin_ice_con[j][k][i])**(1.0/3.0)
                        bin_ice_nex[j][k][i] = bin_ice_nex[j][k][i]/bin_ice_con[j][k][i]
                    else:
                        bin_ice_aaxis[j][k][i] = 0.0
                        bin_ice_caxis[j][k][i] = 0.0
                        bin_ice_daxis[j][k][i] = 0.0
                        bin_ice_agaxis[j][k][i] = 0.0
                        bin_ice_cgaxis[j][k][i] = 0.0
                        bin_ice_nex[j][k][i] = 0.0
    else:
        print("currently axis version should not be 2", axis_version)
        sys.exit(0)

    if (t == 0):
        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    if (bin_ice_con[j][k][i] > 1.e-12):
                        habit = Iice_shapePoly(bin_ice_aaxis[j][k][i],bin_ice_caxis[j][k][i],bin_ice_daxis[j][k][i],bin_ice_agaxis[j][k][i],bin_ice_cgaxis[j][k][i],bin_ice_nex[j][k][i])
                        if (habit < 0):
                            print("FATAL ERROR IN HABIT: ", j, k, i)
                        elif (habit == 1):
                            INIThabit_monopl[k][i] = INIThabit_monopl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 2):
                            INIThabit_monodd[k][i] = INIThabit_monodd[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 3):
                            INIThabit_monocl[k][i] = INIThabit_monocl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 4):
                            INIThabit_polypl[k][i] = INIThabit_polypl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 5):
                            INIThabit_polycl[k][i] = INIThabit_polycl[k][i] + bin_ice_con[j][k][i]*1000.0
                        else:
                            INIThabit_polyir[k][i] = INIThabit_polyir[k][i] + bin_ice_con[j][k][i]*1000.0

        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    INITppv_con[k][i] = INITppv_con[k][i] + bin_ice_con[j][k][i]
                    INITppv_mass[k][i] = INITppv_mass[k][i] + bin_ice_cry[j][k][i]
                    INITppv_aaxis[k][i] = INITppv_aaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_aaxis[j][k][i]*10000.0
                    INITppv_caxis[k][i] = INITppv_caxis[k][i] + bin_ice_con[j][k][i]*bin_ice_caxis[j][k][i]*10000.0
                    INITppv_agaxis[k][i] = INITppv_agaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_agaxis[j][k][i]*10000.0
                    INITppv_cgaxis[k][i] = INITppv_cgaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_cgaxis[j][k][i]*10000.0
                    INITppv_nex[k][i] = INITppv_nex[k][i] + bin_ice_con[j][k][i]*bin_ice_nex[j][k][i]

        print("MAX",max(max(x) for x in INITppv_con))
        print("MAX",max(max(x) for x in INITppv_aaxis))

        for k in range(zlim):
            for i in range(ylim):
                if (INITppv_con[k][i] > 1.e-8):
                    INITppv_aaxis[k][i] = INITppv_aaxis[k][i]/INITppv_con[k][i]
                    INITppv_caxis[k][i] = INITppv_caxis[k][i]/INITppv_con[k][i]
                    INITppv_agaxis[k][i] = INITppv_agaxis[k][i]/INITppv_con[k][i]
                    INITppv_cgaxis[k][i] = INITppv_cgaxis[k][i]/INITppv_con[k][i]
                    INITppv_nex[k][i] = INITppv_nex[k][i]/INITppv_con[k][i]
                else:
                    INITppv_aaxis[k][i] = 0.0
                    INITppv_caxis[k][i] = 0.0
                    INITppv_agaxis[k][i] = 0.0
                    INITppv_cgaxis[k][i] = 0.0

    if (t == tlim - 1):
        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    if (bin_ice_con[j][k][i] > 1.e-12):
                        habit = Iice_shapePoly(bin_ice_aaxis[j][k][i],bin_ice_caxis[j][k][i],bin_ice_daxis[j][k][i],bin_ice_agaxis[j][k][i],bin_ice_cgaxis[j][k][i],bin_ice_nex[j][k][i])
                        if (habit < 0):
                            print("FATAL ERROR IN HABIT: ", j, k, i)
                        elif (habit == 1):
                            habit_monopl[k][i] = habit_monopl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 2):
                            habit_monodd[k][i] = habit_monodd[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 3):
                            habit_monocl[k][i] = habit_monocl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 4):
                            habit_polypl[k][i] = habit_polypl[k][i] + bin_ice_con[j][k][i]*1000.0
                        elif (habit == 5):
                            habit_polycl[k][i] = habit_polycl[k][i] + bin_ice_con[j][k][i]*1000.0
                        else:
                            habit_polyir[k][i] = habit_polyir[k][i] + bin_ice_con[j][k][i]*1000.0

        for j in range(len(massbini)):
            for k in range(zlim):
                for i in range(ylim):
                    ppv_con[k][i] = ppv_con[k][i] + bin_ice_con[j][k][i]
                    ppv_mass[k][i] = ppv_mass[k][i] + bin_ice_cry[j][k][i]
                    ppv_aaxis[k][i] = ppv_aaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_aaxis[j][k][i]*10000.0
                    ppv_caxis[k][i] = ppv_caxis[k][i] + bin_ice_con[j][k][i]*bin_ice_caxis[j][k][i]*10000.0
                    ppv_agaxis[k][i] = ppv_agaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_agaxis[j][k][i]*10000.0
                    ppv_cgaxis[k][i] = ppv_cgaxis[k][i] + bin_ice_con[j][k][i]*bin_ice_cgaxis[j][k][i]*10000.0
                    ppv_nex[k][i] = ppv_nex[k][i] + bin_ice_con[j][k][i]*bin_ice_nex[j][k][i]

        for k in range(zlim):
            for i in range(ylim):
                if (ppv_con[k][i] > 1.e-8):
                    ppv_aaxis[k][i] = ppv_aaxis[k][i]/ppv_con[k][i]
                    ppv_caxis[k][i] = ppv_caxis[k][i]/ppv_con[k][i]
                    ppv_agaxis[k][i] = ppv_agaxis[k][i]/ppv_con[k][i]
                    ppv_cgaxis[k][i] = ppv_cgaxis[k][i]/ppv_con[k][i]
                    ppv_nex[k][i] = ppv_nex[k][i]/ppv_con[k][i]
                else:
                    ppv_aaxis[k][i] = 0.0
                    ppv_caxis[k][i] = 0.0
                    ppv_agaxis[k][i] = 0.0
                    ppv_cgaxis[k][i] = 0.0


print("LOG -- END READING NETCDF VARIABLES AND ANALYSIS")

f.close()


print("LOG -- START PLOTTING")
    
INITppv_con1d    = [0.0 for _ in range(ylim)]
INITppv_mass1d   = [0.0 for _ in range(ylim)]
INITppv_aaxis1d  = [0.0 for _ in range(ylim)]
INITppv_caxis1d  = [0.0 for _ in range(ylim)]
INITppv_agaxis1d = [0.0 for _ in range(ylim)]
INITppv_cgaxis1d = [0.0 for _ in range(ylim)]
INITppv_nex1d    = [0.0 for _ in range(ylim)]
INITratio_caaxis1d = [0.0 for _ in range(ylim)]
INITratio_aagaxis1d = [0.0 for _ in range(ylim)]
INITratio_ccgaxis1d = [0.0 for _ in range(ylim)]
INIThabit_monopl1d = [0.0 for _ in range(ylim)]
INIThabit_monocl1d = [0.0 for _ in range(ylim)]
INIThabit_monodd1d = [0.0 for _ in range(ylim)]
INIThabit_polypl1d = [0.0 for _ in range(ylim)]
INIThabit_polycl1d = [0.0 for _ in range(ylim)]
INIThabit_polyir1d = [0.0 for _ in range(ylim)]

ppv_con1d    = [0.0 for _ in range(ylim)]
ppv_mass1d   = [0.0 for _ in range(ylim)]
ppv_aaxis1d  = [0.0 for _ in range(ylim)]
ppv_caxis1d  = [0.0 for _ in range(ylim)]
ppv_agaxis1d = [0.0 for _ in range(ylim)]
ppv_cgaxis1d = [0.0 for _ in range(ylim)]
ppv_nex1d    = [0.0 for _ in range(ylim)]
ratio_caaxis1d = [0.0 for _ in range(ylim)]
ratio_aagaxis1d = [0.0 for _ in range(ylim)]
ratio_ccgaxis1d = [0.0 for _ in range(ylim)]
habit_monopl1d = [0.0 for _ in range(ylim)]
habit_monocl1d = [0.0 for _ in range(ylim)]
habit_monodd1d = [0.0 for _ in range(ylim)]
habit_polypl1d = [0.0 for _ in range(ylim)]
habit_polycl1d = [0.0 for _ in range(ylim)]
habit_polyir1d = [0.0 for _ in range(ylim)]

for i in range(ylim):
    for k in range(zlim):
        INITppv_con1d[i] = INITppv_con1d[i] + (INITppv_con[k][i])/zlim
        INITppv_mass1d[i] = INITppv_mass1d[i] + (INITppv_mass[k][i])/zlim
        INITppv_aaxis1d[i] = INITppv_aaxis1d[i] + (INITppv_aaxis[k][i])/zlim
        INITppv_caxis1d[i] = INITppv_caxis1d[i] + (INITppv_caxis[k][i])/zlim
        INITppv_agaxis1d[i] = INITppv_agaxis1d[i] + (INITppv_agaxis[k][i])/zlim
        INITppv_cgaxis1d[i] = INITppv_cgaxis1d[i] + (INITppv_cgaxis[k][i])/zlim
        INITppv_nex1d[i] = INITppv_nex1d[i] + (INITppv_nex[k][i])/zlim
        INIThabit_monopl1d[i] = INIThabit_monopl1d[i] + (INIThabit_monopl[k][i])/zlim
        INIThabit_monocl1d[i] = INIThabit_monocl1d[i] + (INIThabit_monocl[k][i])/zlim
        INIThabit_monodd1d[i] = INIThabit_monodd1d[i] + (INIThabit_monodd[k][i])/zlim
        INIThabit_polypl1d[i] = INIThabit_polypl1d[i] + (INIThabit_polypl[k][i])/zlim
        INIThabit_polycl1d[i] = INIThabit_polycl1d[i] + (INIThabit_polycl[k][i])/zlim
        INIThabit_polyir1d[i] = INIThabit_polyir1d[i] + (INIThabit_polyir[k][i])/zlim
    
        ppv_con1d[i] = ppv_con1d[i] + (ppv_con[k][i])/zlim
        ppv_mass1d[i] = ppv_mass1d[i] + (ppv_mass[k][i])/zlim
        ppv_aaxis1d[i] = ppv_aaxis1d[i] + (ppv_aaxis[k][i])/zlim
        ppv_caxis1d[i] = ppv_caxis1d[i] + (ppv_caxis[k][i])/zlim
        ppv_agaxis1d[i] = ppv_agaxis1d[i] + (ppv_agaxis[k][i])/zlim
        ppv_cgaxis1d[i] = ppv_cgaxis1d[i] + (ppv_cgaxis[k][i])/zlim
        ppv_nex1d[i] = ppv_nex1d[i] + (ppv_nex[k][i])/zlim
        habit_monopl1d[i] = habit_monopl1d[i] + (habit_monopl[k][i])/zlim
        habit_monocl1d[i] = habit_monocl1d[i] + (habit_monocl[k][i])/zlim
        habit_monodd1d[i] = habit_monodd1d[i] + (habit_monodd[k][i])/zlim
        habit_polypl1d[i] = habit_polypl1d[i] + (habit_polypl[k][i])/zlim
        habit_polycl1d[i] = habit_polycl1d[i] + (habit_polycl[k][i])/zlim
        habit_polyir1d[i] = habit_polyir1d[i] + (habit_polyir[k][i])/zlim
    if (INITppv_aaxis1d[i] != 0.0):
        INITratio_caaxis1d[i] = INITppv_caxis1d[i]/INITppv_aaxis1d[i]
        INITratio_aagaxis1d[i] = INITppv_agaxis1d[i]/INITppv_aaxis1d[i]
    if (INITppv_caxis1d[i] != 0.0):
        INITratio_ccgaxis1d[i] = INITppv_cgaxis1d[i]/INITppv_caxis1d[i]
    if (ppv_aaxis1d[i] != 0.0):
        ratio_caaxis1d[i] = ppv_caxis1d[i]/ppv_aaxis1d[i]
        ratio_aagaxis1d[i] = ppv_agaxis1d[i]/ppv_aaxis1d[i]
    if (ppv_caxis1d[i] != 0.0):
        ratio_ccgaxis1d[i] = ppv_cgaxis1d[i]/ppv_caxis1d[i]

for i, temp in enumerate(y):
    y[i] = y[i]/1000.0


plt.close()
#fig = plt.figure(figsize=(32,24.5))
fig = plt.figure(figsize=(24.5,24.5))
gs1 = gridspec.GridSpec(3, 3)
gs1.update(wspace=0.42, hspace=0.52)  # set the spacing between axes.
        
ax1 = plt.subplot(gs1[0,:-1]) # habits
ax2 = plt.subplot(gs1[3]) # ca ratio
ax3 = plt.subplot(gs1[4]) # aag ratio
ax4 = plt.subplot(gs1[5]) # ccg ratio
ax5 = plt.subplot(gs1[6]) # a axis
ax6 = plt.subplot(gs1[7]) # c axis
ax7 = plt.subplot(gs1[8]) # nex
        
p1, = ax1.plot(y,INIThabit_monopl1d,linestyle="solid",color="blue",label="Init. mono plate")
p2, = ax1.plot(y,habit_monopl1d,linestyle="dashed",color="blue",label="Final mono plate")
#ax1.plot(y,ini1_habit_monocl,linestyle="solid",color="magenta",label="Init mono column")
#ax1.plot(y,ini1_habit_polycl,linestyle="solid",color="red",label="Init. poly column")
#ax1.plot(y,ini1_habit_polypl,linestyle="solid",color="purple",label="Init. poly plate")
p3, = ax1.plot(y,INIThabit_polyir1d,linestyle="solid",color="green",label="Init. poly irregular")
p4, = ax1.plot(y,habit_polyir1d,linestyle="dashed",color="green",label="Final poly irregular")
p5, = ax1.plot(y,habit_monocl1d,linestyle="dashed",color="magenta",label="Final mono column")
p6, = ax1.plot(y,habit_polycl1d,linestyle="dashed",color="red",label="Final poly column")
#p5, = ax1.plot(y,fnl1_habit_polypl,linestyle="dashed",color="purple",label="Final poly plate")
ax1.set_ylim(1.e-5,1)
ax1.set_xlim(0,10)
ax1.set_yscale("log")
ax1.set_xticks([0,2,4,6,8,10])
ax1.grid(True)
#ax1.set_aspect()
ax1.set_title("(a) Habits",fontsize=30)
ax1.set_xlabel("x (km)",fontsize=30)
ax1.set_ylabel("Conc. (L$^{-1}$)",fontsize=30)
ax1.tick_params(axis='both',which='major',labelsize=30)
#ax1.legend(ncol=1,bbox_to_anchor=(1.125,1.07),loc='upper left',prop={'size':30})
        
p11, = ax2.plot(y,INITratio_caaxis1d,linestyle="solid",color="black",label="Init. distribution")
p12, = ax2.plot(y,ratio_caaxis1d,linestyle="dashed",color="black",label="Final distribution")
ax2.set_ylim(0,1.1)
ax2.set_xlim(0,10)
ax2.set_xticks([0,2,4,6,8,10])
ax2.grid(True)
#ax2.legend(prop={'size':24.5})
ax2.set_title("(b) $c/a$ axis ratio",fontsize=30)
ax2.set_xlabel("x (km)",fontsize=30)
ax2.tick_params(axis='both',which='major',labelsize=30)

lns = [p1, p2, p3, p4, p5, p6, p11, p12]
ax1.legend(handles=lns,ncol=1,bbox_to_anchor=(1.125,1.07),loc='upper left',prop={'size':30})

ax3.plot(y,INITratio_aagaxis1d,linestyle="solid",color="black")
ax3.plot(y,ratio_aagaxis1d,linestyle="dashed",color="black")
ax3.set_ylim(0,2)
ax3.set_xlim(0,10)
ax3.set_xticks([0,2,4,6,8,10])
ax3.grid(True)
#ax3.legend(prop={'size':24.5})
ax3.set_title("(c) $a_g/a$ axis ratio",fontsize=30)
ax3.set_xlabel("x (km)",fontsize=30)
ax3.tick_params(axis='both',which='major',labelsize=30)

ax4.plot(y,INITratio_ccgaxis1d,linestyle="solid",color="black")
ax4.plot(y,ratio_ccgaxis1d,linestyle="dashed",color="black")
ax4.set_ylim(0,2)
ax4.set_xlim(0,10)
ax4.set_xticks([0,2,4,6,8,10])
ax4.grid(True)
#ax4.legend(prop={'size':24.5})
ax4.set_title("(d) $c_g/c$ axis ratio",fontsize=30)
ax4.set_xlabel("x (km)",fontsize=30)
ax4.tick_params(axis='both',which='major',labelsize=30)

ax5.plot(y,INITppv_aaxis1d,linestyle="solid",color="black")
ax5.plot(y,ppv_aaxis1d,linestyle="dashed",color="black")
#ax5.set_ylim(0,250) # for 25 26
ax5.set_ylim(0,450) # for 29 30
ax5.set_xlim(0,10)
ax5.set_xticks([0,2,4,6,8,10])
ax5.grid(True)
#ax5.legend(prop={'size':24.5})
ax5.set_title("(e) $a$ axis length",fontsize=30)
ax5.set_xlabel("x (km)",fontsize=30)
ax5.tick_params(axis='both',which='major',labelsize=30)

ax6.plot(y,INITppv_caxis1d,linestyle="solid",color="black")
ax6.plot(y,ppv_caxis1d,linestyle="dashed",color="black")
ax6.set_ylim(0,200)
ax6.set_xlim(0,10)
ax6.set_xticks([0,2,4,6,8,10])
ax6.grid(True)
#ax6.legend(prop={'size':24.5})
ax6.set_title("(f) $c$ axis length",fontsize=30)
ax6.set_xlabel("x (km)",fontsize=30)
ax6.tick_params(axis='both',which='major',labelsize=30)

ax7.plot(y,INITppv_nex1d,linestyle="solid",color="black")
ax7.plot(y,ppv_nex1d,linestyle="dashed",color="black")
ax7.set_ylim(0,2)
ax7.set_xlim(0,10)
ax7.set_xticks([0,2,4,6,8,10])
ax7.grid(True)
#ax7.legend(prop={'size':24.5})
ax7.set_title("(g) $n_{ex}$",fontsize=30)
ax7.set_xlabel("x (km)",fontsize=30)
ax7.tick_params(axis='both',which='major',labelsize=30)


ax_dx = 0./72.; ax_dy=-8./72.
offset = transform.ScaledTranslation(ax_dx,ax_dy,fig.dpi_scale_trans)

for label in ax1.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)
            
for label in ax2.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax3.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax4.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax5.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax6.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax7.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

ax_dx = -8./72.; ax_dy=0/72.
offset = transform.ScaledTranslation(ax_dx,ax_dy,fig.dpi_scale_trans)

for label in ax1.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)
            
for label in ax2.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax3.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax4.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax5.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax6.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

for label in ax7.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

plt.savefig(plot_dir+'figure_advection.png',format='png',dpi=300)


print("LOG -- END")

