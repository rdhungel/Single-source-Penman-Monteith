#Model_3_Satellite_passing_date_soil_surface_resistance_correction_05172008

#This model computes soil moisture at surface, soil surface evaporation and soil surface resistance. It also refines the value of transpiration computed from previous model (Model_2). This model completes the inversion of METRIC ET at satellite overpass time. The final outputs soil surface resistance and soil surface moisture. The surface energy fluxes for this inversion process completes in this model.
#Data needed:
#NARR reanalysis: Data for satellite overpass time-Incoming longwave radiation, incoming shortwave radiation, windspeed, air temperature, specific humidity
#METRIC: Satellite overpass time-Evapotranspiration, Roughness length of momentum, board band emissivity, surface albedo, leaf area index, normalized difference vegetation index
#Model_2 data: Canopy resistance

# Import arcgisscripting library/module
import arcgisscripting
import os
import time
import sys

DEBUG = True

# Adding a test comment
# One more comment


def CleanUpTmp(tDir):
    #Cleanup (remove) temporary files and rasters in the scratch workspace.
    for root, dirs, files in os.walk(tDir, topdown=False):
        for name in files:
            os.remove(os.path.join(root,name))

        for name in dirs:
            os.rmdir(os.path.join(root,name))

    print 'Cleaned temporary raster'

def PrintCellValue(gp,raster,remark):
   
    rvalue=(gp.GetCellValue(raster,'2600290.913  1328482.369  1', '1').getoutput(0))
    print remark,rvalue

# Create a geoprocessor object, for example gp or rdh or something
gp = arcgisscripting.create(9.3)
# Check out ArcGIS Spatial Analyst extension license
gp.CheckOutExtension("Spatial")
# Overwrite existing outputs
gp.OverWriteOutput = 1
#____________________________________________________________________________


# Section 2- Defining cell size, extent
# Overwrite existing outputs
gp.CellSize = 'MINOF'
gp.Extent = 'MINOF'

#____________________________________________________________________________

# Section 3- Defining working directory, common files, folders

start = time.time()

print 'Start Index 1102'
print 'Date 05/17/2008'
print 'Time 11 am'

gp.ScratchWorkspace = r'J:\Ramesh\arcgis'

# Commons
ComDir = r'J:\Ramesh'

# Set a working project directory
workDir = r'I:\Ramesh\Combined_Ess_Ec_interpolation_05172008'
METRICDir= r'J:\Ramesh\Old_folders\METRICData_05172008'
ClipDir = r'J:\Ramesh\Ext_0517_NARR'

Model1Dir = r'J:\Ramesh\Combined_Ess_Ec_05172008'

UGRDDir=r'J:\Ramesh\UGRD_OUT\UGRD_'
VGRDDir=r'J:\Ramesh\VGRD_OUT\VGRD_'
TMPDir=r'J:\Ramesh\TMP_OUT\TMP_'
DLWRFDir=r'J:\Ramesh\DLWRF_OUT\DLWRF_'
DSWRFDir = r'J:\Ramesh\DSWRF_OUT\DSWRF_'
SPFHDir = r'J:\Ramesh\SPFH_OUT\SPFH_'

PrecipDir = r'J:\Ramesh\Final_output_precip'

# Defining Area of Interest
aoi=os.path.join(ComDir,'aoi_small_Idaho.shp')
aoi_extent = "2597181.4617845 1309037.436659 2615838.396443 1340785.540191"

# Check and apply reasonable limits on ETinst_hourly in mm/hr
low = '0.005' # For hourly
high = '1.2' # For hourly

# Limits on Surface Temperature in K ( be careful if we give big window in temperature, sometimes it won't converge)
tslo = '265.0'
tshi = '335.0'

# Limits on Sensible Heat in W/m2 ( be careful if we give big window in sensible heat flux, sometimes it won't converge)
shlo = '-100.0'
shhi = '500.0'

# Limits on ground Heat in W/m2
ghlo = '-100.0'
ghhi = '650.0'

# Limiting value aerodynamic resistance that should be iterated

# In this iteration 98% of rah is within a limit of -1 s/m to +1 s/m
iter_lo = '-1.0'
iter_hi = '1.0'

# Date and time
# Mountain time 05/17/2008 - 11am
# UTM time = 05/17/2008 - 6pm

# Limiting friciton velocity in m/s
u_fri_lo = '0.01'
u_fri_hi = '500.0'

# Limiting value of aerodynamic resistance in s/m
rahlo = '1.0'
rahhi = '200.0'

# UTM time day of the year
#Start_Day_of_the_year = 170

Start_Day_of_the_year = 138

# Hour start from 0, 3am, 6am, 9am, 15pm, 18pm, 21pm, 24pm
Start_hour = 18

Start_Index = ((Start_Day_of_the_year - 1 ) * 8 ) + (Start_hour / 3 )
print 'Done\n\tStart_Index: ',Start_Index ,

# Minimum and maximum fraction of cover
fc_min = '0.05'
fc_max = '1'

fc_full_veg = '1'

# Roughness of surface
Zos = '0.01'

# Mean boundary layer resistance
rb = '25'

# Partition of Albedo in soil and vegetation
Albedo_soil_min = '0.15'
Albedo_soil_max = '0.28'
Albedo_veg_max = '0.24'
Albedo_veg_min = '0.15'

# Soil moisture value for albedo computation
soilm_min = '0.0'

#Vegetation Part
Ref_ET = '0.0002389'
ET_min = '0.000003'

# Maximum and Minimum NDVI
NDVImax = '0.8'
NDVImin = '0.15'

#____________________________________________________________________________

# Section 4
Cp= '1013'

Stefan_Boltzamn= '0.0000000568'

# Emissivity of soil and vegetation
Emiss_soil = '0.95'
Emiss_veg = '0.98'

i= Start_Index

# METRIC data set
start = time.time()

DEM =os.path.join(METRICDir, 'DEM_clipped.img')
ETins =os.path.join(METRICDir, 'ET_ins_hourly_clipped.img')
Landuse =os.path.join(ComDir, 'Landuse_clipped.img')
NDVI=os.path.join(METRICDir, 'NDVI_clipped.img')
Zom_M=os.path.join(METRICDir, 'Zom_clipped.img')
Albedo=os.path.join(METRICDir, 'Albedo_clipped.img')
BB_Emissi=os.path.join(METRICDir, 'BBemmissivity_clipped.img')
LAI_M=os.path.join(METRICDir, 'LAI_clipped.img')
ETrF_M = os.path.join(METRICDir, 'etrf24_stretched_05172008_P39R30_L5_IDWR.img')

# Ground heat flux is for starting guess value, it will be updated in the iteration process
G_Flux=os.path.join(METRICDir, 'gflux_clipped.img')
H_Flux = os.path.join(METRICDir, 'hflux_clipped.img')

# Water bodies temperature is subtituted from METRIC tempeature (landuse 11)
Ts_METRIC = os.path.join(METRICDir, 'Ts_METRIC_clipped.img')

# Soil Map is developed
Soil_final = os.path.join(ComDir, 'Soil_map_clipped.img')

#****************************************************************
print' Generating values of given coordinates'

if DEBUG:
    PrintCellValue(gp,DEM,'\tResult -- DEM: ')
    PrintCellValue(gp,ETins,'\tResult -- ETins: ')
    PrintCellValue(gp,Landuse,'\tResult -- Landuse: ')
    PrintCellValue(gp,NDVI,'\tResult -- NDVI: ')
    PrintCellValue(gp,Zom_M,'\tResult -- Zom_M: ')
    PrintCellValue(gp,Albedo,'\tResult -- Albedo: ')
    PrintCellValue(gp,BB_Emissi,'\tResult -- BB_Emissi: ')
    PrintCellValue(gp,LAI_M,'\tResult -- LAI_M: ')
    PrintCellValue(gp,ETrF_M,'\tResult -- ETrF_M: ')
    PrintCellValue(gp,G_Flux,'\tResult -- G_Flux: ')
    PrintCellValue(gp,H_Flux,'\tResult -- H_Flux: ')
    PrintCellValue(gp,Ts_METRIC,'\tResult -- Ts_METRIC: ')
    PrintCellValue(gp,Soil_final,'\tResult -- Soil_final: ')

gp.SnapRaster = Albedo

Water_flag = os.path.join(workDir, 'Water_flag.img')
Water_flg_eqn= 'con(('+Landuse+' == 11) or ( '+NDVI+' < 0), 1, 0)'
gp.SingleOutputMapAlgebra_sa(Water_flg_eqn, Water_flag)

#************************************************************************
# NARR data
Tair_F = TMPDir +str(i)+'.img'
In_long =DLWRFDir +str(i)+'.img'
In_short_act =DSWRFDir +str(i)+'.img'
wind_u=UGRDDir +str(i)+'.img'
wind_v=VGRDDir +str(i)+'.img'
S_hum = SPFHDir +str(i)+'.img'

In_short =os.path.join(workDir, 'DSWRF'+ str(i)+'.img')
In_short_con = 'Con(IsNull('+In_short_act+'),0, '+In_short_act+')'
gp.SingleOutputMapAlgebra_sa(In_short_con,In_short)

# Precipitaton for the current time step
precip_att = os.path.join(PrecipDir, 'NARR_apcp_clip'+ str(i)+'.img')

# Precipitation for the current step without negatives
precip = os.path.join(workDir, 'NARR_apcp_new_'+ str(i)+'.img')
precip_con_eqn= 'CON( '+precip_att+' <= 0 , 0, '+precip_att+' )'
gp.SingleOutputMapAlgebra_sa(precip_con_eqn, precip)

irrigation = os.path.join(workDir, 'irrigation_'+ str(i)+'.img')
irrigation_eqn= precip_att+' * 0'
gp.SingleOutputMapAlgebra_sa(irrigation_eqn, irrigation)

if DEBUG:
    PrintCellValue(gp,Tair_F,'\tResult -- Tair_F: ')
    PrintCellValue(gp,In_long,'\tResult -- In_long: ')
    PrintCellValue(gp,In_short,'\tResult -- In_short: ')
    PrintCellValue(gp,wind_u,'\tResult -- wind_u: ')
    PrintCellValue(gp,wind_v,'\tResult -- wind_v: ')
    PrintCellValue(gp,S_hum,'\tResult -- S_hum: ')

################################################################################################################################
# Results from previous model

eoveg_M1=os.path.join(Model1Dir, 'eoveg_'+ str(i)+'.img')
Lambda_veg_M1 = os.path.join(Model1Dir, 'Lambda_veg_'+ str(i)+'.img')
rsc_M1 = os.path.join(Model1Dir, 'rsc_final_'+ str(i)+'.img')
rah_M1 = os.path.join(Model1Dir, 'rah_'+ str(i)+'.img')
u_fri_M1 = os.path.join(Model1Dir, 'u_fri_'+ str(i)+'.img')
soilm_root_M1 = os.path.join(Model1Dir, 'soilm_root_final_'+ str(i)+'.img')

if DEBUG:
    PrintCellValue(gp,eoveg_M1,'\tResult -- eoveg_M1: ')
    PrintCellValue(gp,Lambda_veg_M1,'\tResult -- Lambda_veg_M1: ')
    PrintCellValue(gp,rsc_M1,'\tResult -- rsc_M1: ')
    PrintCellValue(gp,rah_M1,'\tResult -- rah_M1: ')
    PrintCellValue(gp,u_fri_M1,'\tResult -- u_fri_M1: ')
    PrintCellValue(gp,soilm_root_M1,'\tResult -- soilm_root_M1: ')

soilm_M1=os.path.join(workDir, 'soilm_M1_'+ str(i)+'.img')
eoveg=os.path.join(workDir, 'eoveg_'+ str(i)+'.img')
Lambda_veg = os.path.join(workDir, 'Lambda_veg_'+ str(i)+'.img')                   

gp.CopyRaster(soilm_root_M1,soilm_M1)
gp.CopyRaster(eoveg_M1,eoveg)
gp.CopyRaster(Lambda_veg_M1,Lambda_veg)

Theta_sat=os.path.join(ComDir,'theta_sat_.img')
Theta_sat_eqn = 'con( '+Soil_final+' == 0 , 0.476, con( '+Soil_final+' == 1 , 0.465, con( '+Soil_final+' == 2 , 0.434, con( '+Soil_final+' == 3 , 0.439,'
Theta_sat_eqn+= 'con( '+Soil_final+' == 4 , 0.421, con( '+Soil_final+' == 5 , 0.339, con( '+Soil_final+' == 6 , 0.434, con( '+Soil_final+' == 7, 0.476,'
Theta_sat_eqn+= 'con( '+Soil_final+' == 8 , 0.387, 0.464 ) ) ) ) ) ) ) ) )'
gp.SingleOutputMapAlgebra_sa(Theta_sat_eqn,Theta_sat)

# Soil moisture at wilting point in m3/m3 
Theta_wilt=os.path.join(ComDir,'theta_wilt_.img')
Theta_wilt_eqn = 'con( '+Soil_final+' == 0 , 0.084, con( '+Soil_final+' == 1 , 0.103, con( '+Soil_final+' == 2 , 0.09, con( '+Soil_final+' == 3 , 0.11,'
Theta_wilt_eqn+= 'con( '+Soil_final+' == 4 , 0.07, con( '+Soil_final+' == 5 , 0.07, con( '+Soil_final+' == 6 , 0.1, con( '+Soil_final+' == 7, 0.12,'
Theta_wilt_eqn+= 'con( '+Soil_final+' == 8 , 0.17, 0.103 ) ) ) ) ) ) ) ) )'
gp.SingleOutputMapAlgebra_sa(Theta_wilt_eqn,Theta_wilt)

# Soil moisture at field capacity in m3/m3
Theta_ref=os.path.join(ComDir,'theta_ref_.img')
Theta_ref_eqn = 'con( '+Soil_final+' == 0 , 0.360, con( '+Soil_final+' == 1 , 0.36, con( '+Soil_final+' == 2 , 0.34, con( '+Soil_final+' == 3 , 0.329,'
Theta_ref_eqn+= 'con( '+Soil_final+' == 4 , 0.2, con( '+Soil_final+' == 5 , 0.2, con( '+Soil_final+' == 6 , 0.30, con( '+Soil_final+' == 7, 0.360,'
Theta_ref_eqn+= 'con( '+Soil_final+' == 8 , 0.36, 0.36 ) ) ) ) ) ) ) ) )'
gp.SingleOutputMapAlgebra_sa(Theta_ref_eqn,Theta_ref)

ETrF_com = os.path.join(workDir,'ETrF_com_'+ str(i)+'.img')
ETrF_com_eqn = ETrF_M+' * 1'
gp.SingleOutputMapAlgebra_sa(ETrF_com_eqn,ETrF_com)

# Calculaiton of fraction cover of soil (fc)
fc = os.path.join(workDir,'fc_'+ str(i)+'.img')
fc_eqn = '( '+NDVI+' - '+NDVImin+' ) / ( '+NDVImax+' - '+NDVImin+' )'
fc_con_eqn = 'con( '+fc_eqn+' <= '+fc_min+', '+fc_min+',  con( '+fc_eqn+' > '+fc_max+', '+fc_max+', '+fc_eqn+'))'
gp.SingleOutputMapAlgebra_sa(fc_con_eqn,fc)

LAI = os.path.join(workDir,'LAI_'+ str(i)+'.img')
LAI_eqn = 'con ( '+Landuse+' == 52, 0.7, con( '+Landuse+' == 71, 3.0, con('+LAI_M+' <= '+fc+', '+fc+','+LAI_M+' )))'
gp.SingleOutputMapAlgebra_sa(LAI_eqn,LAI)

###################################################################################################################

# Height of Canopy
hc = os.path.join(workDir, 'hc_'+ str(i)+'.img')
hc_forest = '2.5 * '+LAI
hc_sage = '0.8'
hc_other = '0.15 * '+LAI

hc_eqn = 'con( '+Landuse+' == 41, '+hc_forest+', con( '+Landuse+' == 42, '+hc_forest+', con( '+Landuse+' == 43, '+hc_forest+','
hc_eqn+= 'con( '+Landuse+' == 51, '+hc_sage+', con( '+Landuse+' == 52, '+hc_sage+', '+hc_other+')))))'
gp.SingleOutputMapAlgebra_sa(hc_eqn,hc)

# Zero plane displacement height
d= os.path.join(workDir, 'd_'+ str(i)+'.img')
d_eqn = '1.1 * '+hc+' * ln ( 1 + POW((0.2 * '+LAI+'), 0.25 ))'
d_con_eqn = 'con('+Water_flag+' == 1, 0, con( '+fc+' < '+fc_min+', 0 , '+d_eqn+'))'
gp.SingleOutputMapAlgebra_sa(d_con_eqn, d)

# Attunuation or decay coefficient
n = os.path.join(workDir, 'n_'+ str(i)+'.img')
n_eqn = '2.31 + 0.194 * '+hc
#n_eqn_con = 'con( '+hc+' <= 1, 2.5 , con( '+hc+' >= 10, 4.25, '+n_eqn+' ) )'
n_eqn_con = 'con('+Landuse+' == 52, 5, con('+hc+' <= 1, 2.5 , con( '+hc+' >= 10, 4.25, '+n_eqn+' ) ))'
gp.SingleOutputMapAlgebra_sa(n_eqn_con,n)

# Roughness length of momentum (zom)
Zom = os.path.join(workDir, 'Zom_'+ str(i)+'.img')
Zom_con_eqn =  'con('+Water_flag+' == 1,'+Zos+' , con ( '+Landuse+' == 52, 0.2, con( '+fc+' < '+fc_min+', '+Zos+' ,'+Zom_M+')))'
gp.SingleOutputMapAlgebra_sa(Zom_con_eqn,Zom)

# Integration heights 
Z1 = os.path.join(workDir,'Z1_'+ str(i)+'.img')
# Roughness lenght of heat
Z1_full = '0.1 * '+Zom
Z1_par = hc+' - '+d
Z1_int =  Z1_par+' - (( '+Z1_par+' -  '+Z1_full+') * ( '+fc+' - 0.6)) / ( 1 - 0.6)'
#Z1_con_eqn = 'con('+Water_flag+' == 1, '+Zos+', con ('+fc+' <= 0.6 , '+Z1_par+' , con ('+fc+' >= '+fc_full_veg+' , '+Z1_full+', '+Z1_int+')))'
Z1_con_eqn = 'con('+Water_flag+' == 1, '+Zos+', con('+Landuse+' == 52, '+Z1_full+', con('+Landuse+' == 71, '+Z1_full+', con ('+fc+' <= 0.6 , '+Z1_par+' , con ('+fc+' >= '+fc_full_veg+' , '+Z1_full+', '+Z1_int+')))))'
gp.SingleOutputMapAlgebra_sa(Z1_con_eqn,Z1)

###############################################################################################################################################
if DEBUG:
    PrintCellValue(gp,ETrF_com,'\tResult -- ETrF_com: ')
    PrintCellValue(gp,fc,'\tResult -- fc: ')
    PrintCellValue(gp,LAI,'\tResult -- LAI: ')
    PrintCellValue(gp,hc,'\tResult -- hc: ')
    PrintCellValue(gp,d,'\tResult -- d: ')
    PrintCellValue(gp,n,'\tResult -- n: ')
    PrintCellValue(gp,Zom,'\tResult -- Zom: ')
    PrintCellValue(gp,Z1,'\tResult -- Z1: ')
################################################################################################################################################
# Changing air temp from  to K
Tair = os.path.join(workDir, 'Tair_K_'+ str(i)+'.img')
Tair_eqn = '0.5555 * ( '+Tair_F+' - 32.0 ) + 273.0'
gp.SingleOutputMapAlgebra_sa(Tair_eqn,Tair)

Tair_oC = os.path.join(workDir, 'air_temp_oC_'+ str(i)+'.img')
gp.Plus_sa(Tair, "-273.15", Tair_oC)

# Windspeed with u and v
uz =os.path.join(workDir, 'uz_' + str(i)+'.img')
wind_eqn= 'POW (( POW ( ( '+wind_u+' ), 2) +   POW ( ( '+wind_v+' ), 2 ) ) , 0.5 )'
wind_con_eqn= 'con ( '+Water_flag+' == 1, 5.5, '+wind_eqn+')'
gp.SingleOutputMapAlgebra_sa(wind_con_eqn,uz)

#########################################################################################################################################################
# Computation of ras 
# Eddy diffusion coefficient
kh = os.path.join(workDir, 'kh_'+ str(i)+'.img')
kh_eqn = '(0.41 * 0.41 * '+uz+' * ('+hc+' - '+d+')) / (ln ((30 - ( '+d+' ) ) / '+Zom+' ))'
gp.SingleOutputMapAlgebra_sa(kh_eqn,kh)

# ras from zom to d + zom for fully vegetated surface
ras_full_ini = os.path.join(workDir, 'ras_full_ini_'+ str(i)+'.img')
ras_full_act_eqn = '('+hc+' * exp ('+n+') * ((exp ( -1 * '+ n+' * ('+Zom+' / '+hc+'))) - exp (-1 * '+n+' * ('+d+' + '+Zom+') / '+hc+')) / ('+n+' * '+kh+'))'
ras_full_ini_eqn = 'con('+ras_full_act_eqn+' <= 1, 1, '+ras_full_act_eqn+')'
gp.SingleOutputMapAlgebra_sa(ras_full_ini_eqn,ras_full_ini)

# ras from zom to d + zom for bare surface
ras_bare_ini = os.path.join(workDir, 'ras_bare_ini_'+ str(i)+'.img')
ras_bare_ini_eqn= '(ln ( 30 / '+Zos+' ) *  ln (('+d+' + '+Zom+' ) / '+Zos+'))  / (0.41 * 0.41 * '+uz+' )'
gp.SingleOutputMapAlgebra_sa(ras_bare_ini_eqn,ras_bare_ini)

# Intermediate value of ras as of parallel resistances
ras_ini = os.path.join(workDir, 'ras_ini_'+ str(i)+'.img')
ras_ini_eqn = '1 / (('+fc+' / '+ras_full_ini+') + ((1 - '+fc+') / '+ras_bare_ini+'))'
gp.SingleOutputMapAlgebra_sa(ras_ini_eqn,ras_ini)

# ras is only applicable for partial canopy
ras_neu = os.path.join(workDir, 'ras_neu_'+ str(i)+'.img')
ras_con_neu_eqn = 'con('+Water_flag+' == 1, 0, con ( '+fc+' >= '+fc_full_veg+', 0, CON ( '+ras_ini+' < 1, 1, con ('+ras_ini+' > 5000, 5000, '+ras_ini+'))))'
gp.SingleOutputMapAlgebra_sa(ras_con_neu_eqn,ras_neu)

ras = os.path.join(workDir, 'ras_'+ str(i)+'.img')
gp.CopyRaster(ras_neu, ras)

# Bulk bounary layer resistance of the vegetative elements in the canopy (rac)is also applicable for partial canopy
rac = os.path.join(workDir, 'rac_'+ str(i)+'.img')
rac_eqn = rb+' * '+fc+' / ( 2 * '+LAI+' )'
rac_con_eqn = 'CON ( '+rac_eqn+' < 1, 1, con ( '+fc+' >= '+fc_full_veg+', 0, con ('+rac_eqn+' > 5000, 5000, '+ rac_eqn+')))'
gp.SingleOutputMapAlgebra_sa(rac_con_eqn,rac)

if DEBUG:
    PrintCellValue(gp,kh,'\tResult -- kh: ')
    PrintCellValue(gp,ras_full_ini,'\tResult -- ras_full_ini: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,ras_neu,'\tResult -- ras_neu: ')
    PrintCellValue(gp,rac,'\tResult -- rac: ')
#####################################################################################################################

# Pressure in kPa
pressure=os.path.join(workDir, 'P.img')
eqn7= '101.3 * POW(((293 - '+DEM+' * 0.0065 ) / 293 ), 5.26 )'
gp.SingleOutputMapAlgebra_sa(eqn7,pressure)

# Actual vapor pressure 
ea=os.path.join(workDir, 'ea_'+ str(i)+'.img')
ea_eqn= '('+pressure+' * '+S_hum+' ) / ( 0.378 * '+S_hum+' +  0.622 )'
gp.SingleOutputMapAlgebra_sa(ea_eqn,ea)

# Psychrometric Constant, kilopascals/C from eqn 8 of FAO 56
Psyc_con=os.path.join(workDir, 'Psyc_con.img')
eqn8= '0.000665 *  '+pressure
gp.SingleOutputMapAlgebra_sa(eqn8,Psyc_con)

# bbA
bbA = os.path.join(workDir, 'bbA_'+ str(i)+'.img')
bbA_eqn = '1.0 - '+BB_Emissi
gp.SingleOutputMapAlgebra_sa(bbA_eqn,bbA)

# 
eMin = (gp.GetRasterProperties(ETins, "MINIMUM").getoutput(0))
eMax = (gp.GetRasterProperties(ETins, "MAXIMUM").getoutput(0))
eAve = (gp.GetRasterProperties(ETins, "MEAN").getoutput(0))
print 'Original ETins_hourly -- Min: ',eMin,' Max: ',eMax, ' Ave: ',eAve

# # ET hourly in mm/hr
EThourly = os.path.join(workDir, 'ETLimit_'+ str(i)+'.img')

if (low) > eMin or (high) < eMax: 
   lim_eq = 'con( '+ETins+' < '+low+', '+low+', '
   lim_eq += 'con( '+ETins+' > '+high+', '+high+', '+ETins+') )'
   gp.SingleOutputMapAlgebra_sa( lim_eq, EThourly )
else: 
   print 'Not limiting ETins_hourly. Already within limits' 
   gp.CopyRaster(ETins, EThourly)

# ET hourly in mm/sec
ETins_sec=os.path.join(workDir, 'ETins_sec_'+ str(i)+'.img')
eqn1= '0.000277778 * '+EThourly
gp.SingleOutputMapAlgebra_sa(eqn1,ETins_sec)

# Air density in kg/m3
Air_den=os.path.join(workDir, 'air_den_'+ str(i)+'.img')
eq37= pressure+' / ( 1.01 * 0.287 *  '+Tair+' )'
gp.SingleOutputMapAlgebra_sa(eq37,Air_den)

# Friction velocity (neutral condition) in m/s
u_fri=os.path.join(workDir, 'u_fri_'+ str(i)+'.img')
u_fri_neu=os.path.join(workDir, 'u_fri_neu_'+ str(i)+'.img')
eqn39= '(0.41 * '+uz+' ) / ln ((30 - ( '+d+' )) / '+Zom+' )'
gp.SingleOutputMapAlgebra_sa(eqn39,u_fri_neu)
gp.CopyRaster(u_fri_neu, u_fri)

# Aerodynamic resistance (neutral condition) in s/m
rah=os.path.join(workDir, 'rah_'+ str(i)+'.img')
rah_est=os.path.join(workDir,'rah_estimate_'+ str(i)+'.img')
eqn4= '( ln (( 30 - ( '+d+' )) / '+Zom+' ) *  ln ((30 - ( '+d+' ) ) / '+Z1+' )) / (0.41 * 0.41 * '+uz+' )'
gp.SingleOutputMapAlgebra_sa(eqn4,rah_est)

# Limiting aerodynamic resistance in s/m
rah_est_eq = 'con( '+rah_est+' < '+rahlo+', '+rahlo+', con( '+rah_est+' > '+rahhi+', '+rahhi+', '+rah_est+' ) )'
gp.SingleOutputMapAlgebra_sa(rah_est_eq, rah)

# Copying rah to rah_est
gp.CopyRaster(rah,rah_est)

######################################################################################################################################################

if DEBUG:
    PrintCellValue(gp,Tair,'\tResult -- Tair: ')
    PrintCellValue(gp,Tair_oC,'\tResult -- Tair_oC: ')
    PrintCellValue(gp,uz,'\tResult -- uz: ')
    PrintCellValue(gp,pressure,'\tResult -- pressure: ')
    PrintCellValue(gp,ea,'\tResult -- ea: ')
    PrintCellValue(gp,Psyc_con,'\tResult -- Psyc_con: ')
    PrintCellValue(gp,bbA,'\tResult -- bbA: ')
    PrintCellValue(gp,ETins_sec,'\tResult -- ETins_sec: ')
    PrintCellValue(gp,Air_den,'\tResult -- Air_den: ')
    PrintCellValue(gp,u_fri,'\tResult -- u_fri: ')
    PrintCellValue(gp,u_fri_neu,'\tResult -- u_fri_neu: ')
    PrintCellValue(gp,rah,'\tResult -- rah: ')
    PrintCellValue(gp,rah_est,'\tResult -- rah_est: ')

#####################################################################################################################################################
# Soil albedo, limiting minimum value
Albedo_soil = os.path.join(workDir, 'Albedo_soil_'+ str(i)+'.img')
Albedo_soil_eqn = Albedo_soil_min+' + (( '+Albedo_soil_max+' - '+Albedo_soil_min+' ) * (('+Theta_ref+' - '+soilm_M1+') / ('+Theta_ref+' - '+soilm_min+' )))'
Albedo_con_soil_eqn = 'con('+Water_flag+' == 1, '+Albedo+',con( '+Albedo_soil_eqn+' >= '+Albedo_soil_max+' , '+Albedo_soil_max+', con ('+Albedo_soil_eqn+' <= '+Albedo_soil_min+' , '+Albedo_soil_min+', '+Albedo_soil_eqn+' )))'

# Vegetation albdo computing from soil albedo and fc
Albedo_veg_old = os.path.join(workDir, 'Albedo_veg_old_'+ str(i)+'.img')
Albedo_veg_eqn = '( '+Albedo+' - ( 1 - '+fc+') * '+Albedo_soil+' ) / '+fc
Albedo_con_veg_eqn = 'con( '+Albedo_veg_eqn+' >= '+Albedo_veg_max+' , '+Albedo_veg_max+',  con ('+Albedo_veg_eqn+' <= '+Albedo_veg_min+' , '+Albedo_veg_min+', '+Albedo_veg_eqn+' ))'

# Limiting minimum and maximum value of vegetation albedo
Albedo_veg = os.path.join(workDir, 'Albedo_veg_'+ str(i)+'.img')
Albedo_con_new_veg_eqn = 'con( '+Albedo_veg_old+' >= '+Albedo_veg_max+', '+Albedo_veg_max+', con( '+Albedo_veg_old+' <= '+Albedo_veg_min+', '
Albedo_con_new_veg_eqn+= ''+Albedo_veg_min+' ,   '+Albedo_veg_old+'))'

gp.SingleOutputMapAlgebra_sa(Albedo_con_soil_eqn,Albedo_soil)
gp.SingleOutputMapAlgebra_sa(Albedo_con_veg_eqn,Albedo_veg_old)
gp.SingleOutputMapAlgebra_sa(Albedo_con_new_veg_eqn, Albedo_veg)

#########################################################################################################################################################
if DEBUG:
    PrintCellValue(gp,Albedo,'\tResult -- Albedo: ')
    PrintCellValue(gp,Albedo_soil,'\tResult -- Albedo_soil: ')
    PrintCellValue(gp,Albedo_veg,'\tResult -- Albedo_veg: ')
###################################################################################################
print 'VEGETATION SECTION'

#Vegetation Part
Ref_ET = '0.0002389'
ET_min = '0.000003'

ETveg_sec_pre = os.path.join(workDir,'ETveg_sec_pre_'+ str(i)+'.img')
ETveg_sec_pre_eqn  = '(1 - '+fc+')  * '+ETins_sec

# Sensible heat flux is copied to vegetation sensible heat flux for initial guess
H_Flux_rep_veg =os.path.join(workDir, 'H_Flux_rep_veg_'+ str(i)+'.img')
H_Flux_rep_veg_eqn = H_Flux+' * '+fc

# Vegeation temperature
Tc = os.path.join(workDir,'Tc_'+ str(i)+'.img')
tveg_eq = '(( '+H_Flux_rep_veg+'   * ('+rah+' + '+rac+')) / ( '+Air_den+' * '+Cp+' )) + '+Tair

# Saturation Vapor Pressure (es), kilopascals eqn 11 from FAO56
print 'Computing Saturated Vapor Pressure ... ',
eqn11_b= '0.611 * EXP((17.27 * ( '+Tc+' - 273.15)) / (('+Tc+' - 273.15) + 237.3) )'

# Latent heat of vaporization 
Lambda_eqn_veg = '(2.501 - 0.00236 * ( '+Tc+' - 273.15 )) * 1000000'

# Latent heat flux
LE_veg =os.path.join(workDir,'LE_veg_'+ str(i)+'.img')
LE_eqn_veg = Lambda_veg+' *  '+ETveg_sec_pre

# Outgoing longwave radiation
outlwr_veg = os.path.join(workDir, 'outlwr_veg_'+ str(i)+'.img')
outlwr_eq_veg = 'POW(( '+Tc+' ), 4) * '+BB_Emissi+' * '+Stefan_Boltzamn

# Net radiation
netrad_veg = os.path.join(workDir, 'netrad_veg_'+ str(i)+'.img')
netrad_eqn_veg = In_short+' - ('+Albedo_veg+' *  '+In_short+' ) + '+In_long+' - '+outlwr_veg+' - ( 1 - '+Emiss_veg+') * '+In_long

# Sensilbe heat flux
sheat_veg = os.path.join(workDir, 'sheat_veg_'+ str(i)+'.img')
sheat_eqn_veg = netrad_veg+' - '+LE_veg

sheat_veg = os.path.join(workDir, 'sheat_veg_'+ str(i)+'.img')
Tc = os.path.join(workDir,'Tc_'+ str(i)+'.img')
LE_veg =os.path.join(workDir,'LE_veg_'+ str(i)+'.img')

soilm_root_M1 = os.path.join(Model1Dir, 'soilm_root_final_'+ str(i)+'.img')
rsc_M1 = os.path.join(Model1Dir, 'rsc_final_'+ str(i)+'.img')

sheat_veg_final = os.path.join(workDir, 'sheat_veg_final_'+ str(i)+'.img')
sheat_veg_final_eqn  = '1 * '+sheat_veg

ETveg_sec = os.path.join(workDir, 'ETveg_sec_' + str(i)+'.img')
ETveg_sec_eqn = '((( '+eoveg+' -  '+ea+' ) * '+Cp+' * '+Air_den+') / (( '+rah+' + '+rsc_M1+' + '+rac+') * '+Psyc_con+')) / '+Lambda_veg
ETveg_sec_con = 'CON( '+fc+' >= '+fc_full_veg+', '+ETins_sec+' , CON ('+ETveg_sec_eqn+' > 1.15 * '+ETins_sec+', 1.15 * '+ETins_sec+', '+ETveg_sec_eqn+' ))'

ETveg_sec_ave = os.path.join(workDir,'ETveg_sec_ave_'+ str(i)+'.img')
ETveg_sec_eqn_ave= '( '+ETveg_sec_pre+' + '+ETveg_sec+' ) / 2'

# combined computation of vegetation portion

gp.SingleOutputMapAlgebra_sa(ETveg_sec_pre_eqn,ETveg_sec_pre)

gp.SingleOutputMapAlgebra_sa(H_Flux_rep_veg_eqn, H_Flux_rep_veg)
gp.SingleOutputMapAlgebra_sa(tveg_eq, Tc)

tclim_eq_veg = 'con( '+Tc+' < '+tslo+', '+tslo+', con( '+Tc+' > '+tshi+', '+tshi+', '+Tc+' ) )'
limited = os.path.join(workDir, 'limited_'+ str(i)+'.img')
errMin = (gp.GetRasterProperties(Tc, "MINIMUM").getoutput(0))
errMax = (gp.GetRasterProperties(Tc, "MAXIMUM").getoutput(0))
errAve = (gp.GetRasterProperties(Tc, "MEAN").getoutput(0))
print 'Done\n\tTc -- Min: ',errMin,' Max: ',errMax, ' Ave: ',errAve

if errMin < (float(tslo)) or errMax > (float(tshi)):
  print '\tAppling limits to Surface Temperature',
  gp.SingleOutputMapAlgebra_sa(tclim_eq_veg, limited)
  gp.CopyRaster(limited, Tc)

gp.SingleOutputMapAlgebra_sa(ETveg_sec_con, ETveg_sec)
gp.SingleOutputMapAlgebra_sa( ETveg_sec_eqn_ave,ETveg_sec_ave )

gp.SingleOutputMapAlgebra_sa(eqn11_b , eoveg)
gp.SingleOutputMapAlgebra_sa(Lambda_eqn_veg, Lambda_veg)
gp.SingleOutputMapAlgebra_sa(LE_eqn_veg , LE_veg)
gp.SingleOutputMapAlgebra_sa(outlwr_eq_veg , outlwr_veg )
gp.SingleOutputMapAlgebra_sa(netrad_eqn_veg , netrad_veg )
gp.SingleOutputMapAlgebra_sa(sheat_eqn_veg , sheat_veg )

gp.SingleOutputMapAlgebra_sa(sheat_veg_final_eqn , sheat_veg_final )

if DEBUG:
    PrintCellValue(gp,ETveg_sec,'\tResult -- ETveg_sec: ')
    
    PrintCellValue(gp,H_Flux_rep_veg,'\tResult -- H_Flux_rep_veg: ')
    PrintCellValue(gp,Tc,'\tResult -- Tc: ')
  
    PrintCellValue(gp,eoveg,'\tResult -- eoveg: ')
    PrintCellValue(gp,Lambda_veg,'\tResult -- Lambda_veg: ')
    PrintCellValue(gp,LE_veg,'\tResult -- LE_veg: ')
    PrintCellValue(gp,outlwr_veg,'\tResult -- outlwr_veg: ')
    PrintCellValue(gp,netrad_veg,'\tResult -- netrad_veg: ')
    PrintCellValue(gp,sheat_veg,'\tResult -- sheat_veg: ')
    PrintCellValue(gp,sheat_veg_final,'\tResult -- sheat_veg_final: ')
######################################################################################################
print 'SOIL SECTION'

#SOIL SECTION

# Starting value of sensible heat flux and soil ground flux is taken as METRIC value

#####################################################################################################
print ET_min

# Partition of sensible heat flux for soil
H_Flux_rep_soil =os.path.join(workDir, 'H_Flux_rep_soil_'+ str(i)+'.img')
H_Flux_rep_soil_eqn = H_Flux+' *  ( 1 - '+fc+')'

# Partition of ground heat flux for soil
G_Flux_rep_soil =os.path.join(workDir, 'G_Flux_rep_soil_'+ str(i)+'.img')
G_Flux_rep_soil_eqn = G_Flux+' * 1'

ETsoil_sec_pre = os.path.join(workDir,'ETsoil_sec_pre_'+ str(i)+'.img')
ETsoil_sec_pre_eqn  = '(1 - 0.9 * '+fc+')  * '+ETins_sec

# Surface temperature
Ts = os.path.join(workDir,'Ts_'+ str(i)+'.img')
tsurf_eq_soil = '(( '+H_Flux_rep_soil+' * ('+rah+' +  '+ras+' )) / ( '+Air_den+' * '+Cp+' )) + '+Tair

print 'Computing Saturated Vapor Pressure ... ',
eosur=os.path.join(workDir, 'eosur_'+ str(i)+'.img')
eqn11_a= '0.611 * EXP((17.27 * ( '+Ts+' - 273.15)) / (('+Ts+' - 273.15) + 237.3) )'

#saturation specific humidity at surface of water
qosur = os.path.join(workDir, 'qosur_'+ str(i)+'.img')
qosur_eqn= '(0.66 * '+eosur+') / ('+pressure+' - 0.37 * '+eosur+')'

# Latent heat of vaporization
Lambda_soil= os.path.join(workDir,'Lambda_soil_'+ str(i)+'.img')
Lambda_eqn_soil = '(2.501 - 0.00236 * ( '+Ts+' - 273.15 )) * 1000000'

# Computation of LE in W/m2
LE_soil =os.path.join(workDir,'LE_soil_'+ str(i)+'.img')
LE_eqn_soil = Lambda_soil+' *  '+ETsoil_sec_pre

# Outgoing longwave radiation in W/m2
outlwr_soil = os.path.join(workDir, 'outlwr_soil_'+ str(i)+'.img')
outlwr_eq_soil = 'POW(( '+Ts+' ), 4) * '+BB_Emissi+' * '+Stefan_Boltzamn

# Net radition in W/m2
netrad_soil = os.path.join(workDir, 'netrad_soil_'+ str(i)+'.img')
netrad_eqn_soil = In_short+' - ('+Albedo_soil+' *  '+In_short+' ) + '+In_long+' - '+outlwr_soil+' - ( 1 - '+Emiss_soil+') * '+In_long

# Sensible heat flux in W/m2
sheat_soil_1 = os.path.join(workDir, 'sheat1_soil_'+ str(i)+'.img')
sheat_eqn_1_soil = netrad_soil+' - '+G_Flux_rep_soil+' - '+LE_soil

# Ground heat flux in W/m2
gheat_soil = os.path.join(workDir,'gheat_soil_'+ str(i)+'.img')
eqn33_soil = 'Max ( 0.4 * '+sheat_soil_1+' , 0.15 * '+netrad_soil+' )'

# Ground heat flux for water, Landuse 11 in W/m2
eqngflux_soil = '( 0.9 * '+netrad_soil+' ) - 40'
congflux_soil = 'con('+Water_flag+' == 1,  '+eqngflux_soil+', '+eqn33_soil+' )'

# Ground heat flux is copied to G_Flux_ite_soil in W/m2
G_Flux_ite_soil = os.path.join(workDir, 'gheat_ite_soil_'+ str(i)+'.img')

# Sensible heat is overwrited again in sheat_soil in W/m2
sheat_soil_2 = os.path.join(workDir, 'sheat2_soil_'+ str(i)+'.img')
sheat_eqn_2_soil = netrad_soil+' - '+G_Flux_ite_soil+' - '+LE_soil

sheat_soil_final = os.path.join(workDir, 'sheat_soil_final_'+ str(i)+'.img')
sheat_soil_final_eqn  = '1 * '+sheat_soil_2

ETsoil_sec = os.path.join(workDir,'ETsoil_sec_'+ str(i)+'.img')
ET_water_eqn = '(( '+qosur+' -  '+S_hum+' )  * '+Air_den+')  / '+rah
ETsoil_sec_eqn = '( '+ETins_sec+' - ('+ETveg_sec+' * '+fc+' )) / ( 1 - '+fc+')'
ETsoil_sec_con_eqn = 'con ('+Water_flag+' == 1, '+ET_water_eqn+', con('+fc+' >= 1, '+ET_min+', con('+ETsoil_sec_eqn+' < '+ET_min+', '+ET_min+', CON ('+ETsoil_sec_eqn+' > 1.15 * '+Ref_ET+', 1.15 * '+Ref_ET+', '+ETsoil_sec_eqn+' ))))'

ETsoil_sec_ave = os.path.join(workDir,'ETsoil_sec_ave_'+ str(i)+'.img')
ETsoil_sec_eqn_ave= '( '+ETsoil_sec_pre+' + '+ETsoil_sec+' ) / 2'

# Combined computation of soil portion

gp.SingleOutputMapAlgebra_sa(H_Flux_rep_soil_eqn, H_Flux_rep_soil)
gp.SingleOutputMapAlgebra_sa(G_Flux_rep_soil_eqn, G_Flux_rep_soil)
gp.SingleOutputMapAlgebra_sa(ETsoil_sec_pre_eqn,ETsoil_sec_pre)
gp.SingleOutputMapAlgebra_sa(tsurf_eq_soil, Ts)
gp.SingleOutputMapAlgebra_sa(eqn11_a , eosur)

# Limiting the value of surface temperature, sensible heat flux and ground heat flux for soil
tslim_eq_soil = 'con( '+Ts+' < '+tslo+', '+tslo+', con( '+Ts+' > '+tshi+', '+tshi+', '+Ts+' ) )'

gp.SingleOutputMapAlgebra_sa(tsurf_eq_soil, Ts)
errMin = (gp.GetRasterProperties(Ts, "MINIMUM").getoutput(0))
errMax = (gp.GetRasterProperties(Ts, "MAXIMUM").getoutput(0))
errAve = (gp.GetRasterProperties(Ts, "MEAN").getoutput(0))
print 'Done\n\tTs -- Min: ',errMin,' Max: ',errMax, ' Ave: ',errAve

if errMin < (float(tslo)) or errMax > (float(tshi)):
  print '\tAppling limits to Surface Temperature',
  gp.SingleOutputMapAlgebra_sa(tslim_eq_soil, limited)
  gp.CopyRaster(limited, Ts)

gp.SingleOutputMapAlgebra_sa(eqn11_a , eosur)
gp.SingleOutputMapAlgebra_sa(qosur_eqn , qosur)
gp.SingleOutputMapAlgebra_sa(ETsoil_sec_con_eqn, ETsoil_sec)
gp.SingleOutputMapAlgebra_sa( ETsoil_sec_eqn_ave,ETsoil_sec_ave )
gp.SingleOutputMapAlgebra_sa(Lambda_eqn_soil,Lambda_soil)
gp.SingleOutputMapAlgebra_sa(LE_eqn_soil,LE_soil)
gp.SingleOutputMapAlgebra_sa(outlwr_eq_soil,outlwr_soil )
gp.SingleOutputMapAlgebra_sa(netrad_eqn_soil,netrad_soil )
gp.SingleOutputMapAlgebra_sa(sheat_eqn_1_soil,sheat_soil_1 )
gp.SingleOutputMapAlgebra_sa(congflux_soil, gheat_soil)
gp.CopyRaster(gheat_soil, G_Flux_ite_soil)
gp.SingleOutputMapAlgebra_sa(sheat_eqn_2_soil, sheat_soil_2)

shlim_eq_soil = 'con( '+sheat_soil_2+' < '+shlo+', '+shlo+', con( '+sheat_soil_2+' > '+shhi+', '+shhi+', '+sheat_soil_2+' ) )'

errMin = (gp.GetRasterProperties(sheat_soil_2, "MINIMUM").getoutput(0))
errMax = (gp.GetRasterProperties(sheat_soil_2, "MAXIMUM").getoutput(0))
errAve = (gp.GetRasterProperties(sheat_soil_2, "MEAN").getoutput(0))

if errMin < (float(shlo)) or errMax > (float(shhi)):
   gp.SingleOutputMapAlgebra_sa(shlim_eq_soil, limited)
   gp.CopyRaster(limited, sheat_soil_2)

gp.SingleOutputMapAlgebra_sa(sheat_soil_final_eqn, sheat_soil_final)

#############################################################################################################################################
# Saturation Vapor Pressure (es), kilopascals eqn 11 from FAO56

# Corrected value of soil surface resistance (s/m)
rss_cor = os.path.join(workDir,'rss_cor_'+ str(i)+'.img')
rss_cor_eqn = '(( '+eosur+' -  '+ea+' ) * '+Cp+' * '+Air_den+' / ( '+LE_soil+' * '+Psyc_con+' )) - '+rah+' - '+ras

# Limiting value of soil surface resistance (s/m)
rss = os.path.join(workDir, 'rss_'+ str(i)+'.img')
rss_con = 'con ( '+Water_flag+' == 1, 0, CON( '+rss_cor+' <= 35 , 35, CON( '+rss_cor+' > 5000 , 5000, '+rss_cor+' )))'

# Computing soil surface moisture in m3/m3
soilm_old = os.path.join(workDir, 'soilm_old_'+ str(i)+'.img')
#soilm_old_eqn = '(( 8.206 - ln ( '+rss+')) * '+Theta_sat+' ) / 4.225'
soilm_old_eqn = Theta_sat+' / POW((('+rss+' - 33.5) / 3.5) , (1 / 2.3))'

# Limiting value of soil surface moisture m3/m3
soilm = os.path.join(workDir, 'soilm_'+ str(i)+'.img')
soilm_con = 'con ( '+Water_flag+' == 1, '+Theta_sat+', CON( '+soilm_old+' <= 0 , 0.01, CON( '+soilm_old+' > '+Theta_ref+' , '+Theta_ref+', '+soilm_old+' )))'

gp.SingleOutputMapAlgebra_sa(rss_cor_eqn,rss_cor)
gp.SingleOutputMapAlgebra_sa(rss_con, rss)
gp.SingleOutputMapAlgebra_sa(soilm_old_eqn, soilm_old)
gp.SingleOutputMapAlgebra_sa(soilm_con, soilm )

#***********************************************************************************************************************************************
print 'Soil section results'

if DEBUG:
    PrintCellValue(gp,H_Flux,'\tResult -- H_Flux: ')
    PrintCellValue(gp,G_Flux,'\tResult -- G_Flux: ')
    PrintCellValue(gp,H_Flux_rep_soil,'\tResult -- H_Flux_rep_soil: ')
    PrintCellValue(gp,G_Flux_rep_soil,'\tResult -- G_Flux_rep_soil: ')
    PrintCellValue(gp,Ts,'\tResult -- Ts: ')
 
    PrintCellValue(gp,eosur,'\tResult -- eosur: ')
    PrintCellValue(gp,qosur,'\tResult -- qosur: ')
    PrintCellValue(gp,ETsoil_sec,'\tResult -- ETsoil_sec: ')
   
  
    PrintCellValue(gp,Lambda_soil,'\tResult -- Lambda_soil: ')
    PrintCellValue(gp,LE_soil,'\tResult -- LE_soil: ')
   
    PrintCellValue(gp,outlwr_soil,'\tResult -- outlwr_soil: ')
    PrintCellValue(gp,netrad_soil,'\tResult -- netrad_soil: ')
    PrintCellValue(gp,sheat_soil_1,'\tResult -- sheat_soil_1: ')
    PrintCellValue(gp,gheat_soil,'\tResult -- gheat_soil: ')
    PrintCellValue(gp,sheat_soil_2,'\tResult -- sheat_soil_2: ')
    PrintCellValue(gp,sheat_soil_final,'\tResult -- sheat_soil_final: ')    
    
    PrintCellValue(gp,ras,'\tResult -- ras: ')
    PrintCellValue(gp,LE_soil,'\tResult -- LE_soil: ')
    PrintCellValue(gp,rah,'\tResult -- rah: ')
    PrintCellValue(gp,ea,'\tResult -- ea: ')
    
    PrintCellValue(gp,rss_cor,'\tResult -- rss_cor: ')
    PrintCellValue(gp,rss,'\tResult -- rss: ')
    
    PrintCellValue(gp,Soil_final,'\tResult -- Soil_final2: ')
    PrintCellValue(gp,Theta_sat,'\tResult -- Theta_sat: ')
    PrintCellValue(gp,soilm_old,'\tResult -- soilm_old: ')
    PrintCellValue(gp,soilm,'\tResult -- soilm: ')
#################################################################################################################################

print 'COMBINED SECTION'

# Total sensible heat flux is computed combining soil and vegetation fraction of sensible heat flux
sheat = os.path.join(workDir, 'sheat_'+ str(i)+'.img')
sheat_eqn = sheat_veg+' * '+fc+' + (1 - '+fc+') * '+sheat_soil_2
shlim_eq = 'con('+Water_flag+' == 1, '+sheat_soil_2+', con( '+sheat_eqn+' < '+shlo+', '+shlo+', con( '+sheat_eqn+' > '+shhi+', '+shhi+', '+sheat_eqn+' ) ))'

#H_Flux_new = os.path.join(workDir,'sheat_new'+ str(i)+'.img')
# Monin Obukhov Length
L = os.path.join(workDir, 'L_'+ str(i)+'.img')
eqn50= ' - 1 * ( '+Air_den+' * '+Cp+' *  '+Tair+' * POW (( '+u_fri_neu+' ), 3 )) / ('+sheat+' *  4.02 )'
######################################################################################################################################################

## Scenario I : Surface is completely bare i.e. d = 0, zom = zos, Z1 = zos, ras = 0
## Scenario II : Surface is completely vegetated i.e. d = d, zom = zom, Z1 = Z1, ras = 0
## Scenario III : Surface is partially vegetated i.e. d = d, zom = zom, Z1 = Z1, ras = ras

# Monin Obukhov Parameter
# Correction of momentum and heat for all scenario altering d, zom, Z1
X_30m= os.path.join(workDir, 'X_30m_'+ str(i)+'.img')
eqn_X_30m='POW ( ( 1 - (( 16 * (30 - '+d+' ) ) / '+L+' )), 0.25 )'

# Psi M
psi_m_30m = os.path.join(workDir, 'psim_30m_'+ str(i)+'.img')
eqn51 = '2.0 * ln(( 1 + '+X_30m+' ) / 2.0) + '
eqn51 += 'ln (( 1 + POW (( '+X_30m+' ) , 2.0 ))  / 2.0 ) - '
eqn51 += '2.0 * ATAN( '+X_30m+' ) + 1.5708'

eqn54 = '-5 * 30 / '+L
eqn54a = '-5 * 30 / '+L
conpsi_m_30m = 'CON( '+L+' <= 0 , '+eqn51+' , '+eqn54+' )'

# Psi H
psi_h_30m = os.path.join(workDir,'psih_30m_'+ str(i)+'.img')
eqn52a = '2.0 * ln (( 1.0 + POW (( '+X_30m+' ) , 2 )) / 2.0 )'
conpsi_h_30m = 'CON( '+L+' <= 0 , '+eqn52a+' , '+eqn54a+' )'

# Corrected friction velocity (u_fri)
eqn56= '(0.41 * '+uz+' ) / ( ln ((30 - ( '+d+' )) / '+Zom+' ) - '+psi_m_30m+' )'
###########################################################################################################################################
# Correction of momentum and heat for partial canopy part of completely bare soil from zom to d+zom

X_dzom= os.path.join(workDir, 'X_dzom_'+ str(i)+'.img')
eqn_X_dzom='POW ( ( 1 - (( 16 * ('+d+' + '+Zom+' ) ) / '+L+' )), 0.25 )'

eqn52a_dzom = '2.0 * ln (( 1.0 + POW (( '+X_dzom+' ) , 2 )) / 2.0 )'
eqn54a_dzom = '-5 * ('+d+' + '+Zom+' ) / '+L

psi_h_dzom = os.path.join(workDir,'psih_dzom_'+ str(i)+'.img')
conpsi_h_dzom = 'CON( '+L+' <= 0 , '+eqn52a_dzom+' , '+eqn54a_dzom+' )'

ras_bare_ini_eqn= '((ln ( 30 / '+Zos+' ) - '+psi_m_30m+') *  (ln (('+d+' + '+Zom+' ) / '+Zos+') -  '+psi_h_dzom+')) / (0.41 * 0.41 * '+uz+' )'
###########################################################################################################################################
# Correction of heat for hc - d
# This is last part added in eqn. 57
X_hd= os.path.join(workDir, 'X_hd_'+ str(i)+'.img')
eqn_X_hd='POW ( ( 1 - (( 16 * ('+hc+' - '+d+' ) ) / '+L+' )), 0.25 )'

eqn52a_hd = '2.0 * ln (( 1.0 + POW (( '+X_hd+' ) , 2 )) / 2.0 )'
eqn54a_hd = '-5 * '+hc+' / '+L

psi_h_hd = os.path.join(workDir,'psih_hd_'+ str(i)+'.img')
conpsi_h_hd = 'con('+Water_flag+' == 1, 0, con ( '+fc+' >= '+fc_full_veg+' , 0, CON( '+L+' <= 0 , '+eqn52a_hd+' , '+eqn54a_hd+' )))'
###########################################################################################################################################
limited = os.path.join(workDir, 'limited_'+ str(i)+'.img')

# Aerodynamic resistance of water is taken as fixed value, so no convergence is required in water (rah)
eqn57 = '( ln ((30 - '+d+') / '+Z1+' ) - '+psi_h_30m+' + '+conpsi_h_hd+' ) / ( 0.41 * '+u_fri+' )'

# Limiting the value of surface temperature, sensible heat flux for vegetation
tclim_eq_veg = 'con( '+Tc+' < '+tslo+', '+tslo+', con( '+Tc+' > '+tshi+', '+tshi+', '+Tc+' ) )'
shlim_eq_veg = 'con( '+sheat_veg+' < '+shlo+', '+shlo+', con( '+sheat_veg+' > '+shhi+', '+shhi+', '+sheat_veg+' ) )'

# Limiting the value of friction velocity, aerodynamic resistance flux for soil
u_fri_lim_eq = 'con( '+u_fri+' < '+u_fri_lo+', '+u_fri_lo+', con( '+u_fri+' > '+u_fri_hi+', '+u_fri_hi+', '+u_fri+' ) )'
rahlim_eq = 'con( '+rah+' < '+rahlo+', '+rahlo+', con( '+rah+' > '+rahhi+', '+rahhi+', '+rah+' ) )'

raherror = os.path.join(workDir,'rah_error_'+ str(i)+'.img')

stat_img = os.path.join(workDir, 'status_'+ str(i)+'.img')
stat_tab = os.path.join(workDir, 'status_'+ str(i)+'.dbf')
iter_eq = 'con( '+raherror+' < '+iter_lo+', 0, con( '+raherror+' > '+iter_hi+', 0, 1 ) )'
################################################################################################################################################
# Eddy diffusivity for ras
kh = os.path.join(workDir, 'kh_'+ str(i)+'.img')
kh_eqn = '(0.41 * 0.41 * '+uz+' * ('+hc+' - '+d+')) / ((ln ((30 - ( '+d+' ) ) / '+Zom+' )) - '+psi_m_30m+') '

# ras for completely vegetated protion
ras_full_ini = os.path.join(workDir, 'ras_full_ini_'+ str(i)+'.img')
ras_full_act_eqn = '('+hc+' * exp ('+n+') * ((exp ( -1 * '+ n+' * ('+Zom+' / '+hc+'))) - exp (-1 * '+n+' * ('+d+' + '+Zom+') / '+hc+')) / ('+n+' * '+kh+'))'
ras_full_ini_eqn = 'con('+ras_full_act_eqn+' <= 1, 1, '+ras_full_act_eqn+')'

# Parallel combination of resistance
ras_ini = os.path.join(workDir, 'ras_ini_'+ str(i)+'.img')
ras_ini_eqn = '1 / (('+fc+' / '+ras_full_ini+') + ((1 - '+fc+') / '+ras_bare_ini+'))'

# No ras for complete bare surface and fully vegetated surface
ras = os.path.join(workDir, 'ras_'+ str(i)+'.img')
ras_con_eqn = 'con('+Water_flag+' == 1, 0, con ( '+fc+' >= '+fc_full_veg+', 0 , CON ( '+ras_ini+' < 1, 1, con ('+ras_ini+' > 5000, 5000, '+ras_ini+'))))'
#############################################################################################################################################
# Computing parameters in the above section
# Combined section

gp.SingleOutputMapAlgebra_sa(sheat_eqn, sheat)


sheat_rep = os.path.join(workDir, 'sheat_rep_'+ str(i)+'.img')
gp.CopyRaster(sheat,sheat_rep)
gp.SingleOutputMapAlgebra_sa(eqn50, L)
gp.SingleOutputMapAlgebra_sa(eqn_X_30m, X_30m)
gp.SingleOutputMapAlgebra_sa(conpsi_m_30m, psi_m_30m)
gp.SingleOutputMapAlgebra_sa(conpsi_h_30m, psi_h_30m)
#######################################################
gp.SingleOutputMapAlgebra_sa(eqn_X_dzom,X_dzom)
gp.SingleOutputMapAlgebra_sa(conpsi_h_dzom  , psi_h_dzom )
####################################################
gp.SingleOutputMapAlgebra_sa(eqn_X_hd,X_hd)
gp.SingleOutputMapAlgebra_sa(conpsi_h_hd,psi_h_hd)
gp.SingleOutputMapAlgebra_sa(kh_eqn,kh)
gp.SingleOutputMapAlgebra_sa(ras_bare_ini_eqn,ras_bare_ini)
gp.SingleOutputMapAlgebra_sa(ras_full_ini_eqn,ras_full_ini)
gp.SingleOutputMapAlgebra_sa(ras_ini_eqn,ras_ini)
gp.SingleOutputMapAlgebra_sa(ras_con_eqn,ras)
###########################################################
gp.SingleOutputMapAlgebra_sa(eqn56, u_fri)
u_fri_lim_eq = 'con( '+u_fri+' < '+u_fri_lo+', '+u_fri_lo+', con( '+u_fri+' > '+u_fri_hi+', '+u_fri_hi+', '+u_fri+' ) )'
errMin = (gp.GetRasterProperties(u_fri, "MINIMUM").getoutput(0))
errMax = (gp.GetRasterProperties(u_fri, "MAXIMUM").getoutput(0))
errAve = (gp.GetRasterProperties(u_fri, "MEAN").getoutput(0))

if errMin < (float(u_fri_lo)) or errMax > (float(u_fri_hi)):
   print '\tAppling limits to friction velocity',
   gp.SingleOutputMapAlgebra_sa(u_fri_lim_eq, limited)
   gp.CopyRaster(limited, u_fri)

gp.SingleOutputMapAlgebra_sa(eqn57, rah)
rahlim_eq = 'con( '+rah+' < '+rahlo+', '+rahlo+', con( '+rah+' > '+rahhi+', '+rahhi+', '+rah+' ) )'
errMin = (gp.GetRasterProperties(rah, "MINIMUM").getoutput(0))
errMax = (gp.GetRasterProperties(rah, "MAXIMUM").getoutput(0))
errAve = (gp.GetRasterProperties(rah, "MEAN").getoutput(0))

if errMin < (float(rahlo)) or errMax > (float(rahhi)):
  print '\tAppling limits to aerodynamic resistance',
  gp.SingleOutputMapAlgebra_sa(rahlim_eq, limited)
  gp.CopyRaster(limited, rah)

 
gp.Minus_sa(rah, rah_est, raherror)
gp.SingleOutputMapAlgebra_sa(iter_eq, stat_img)
gp.ZonalStatisticsAsTable_sa(stat_img, "Value", raherror, stat_tab, "DATA")

####################################################################################################################
print 'Results of combined section'

if DEBUG:
    
    PrintCellValue(gp,sheat,'\tResult -- sheat: ')
    PrintCellValue(gp,L,'\tResult -- L: ')
    PrintCellValue(gp,X_30m,'\tResult -- X_30m: ')
    PrintCellValue(gp,psi_m_30m,'\tResult -- psi_m_30m: ')
    PrintCellValue(gp,psi_h_30m,'\tResult -- psi_h_30m: ')
    PrintCellValue(gp,X_dzom,'\tResult -- X_dzom: ')
    PrintCellValue(gp,psi_h_dzom,'\tResult -- psi_h_dzom: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,X_hd,'\tResult -- X_hd: ')
    PrintCellValue(gp,psi_h_hd,'\tResult -- psi_h_hd: ')
    PrintCellValue(gp,kh,'\tResult -- kh: ')
    PrintCellValue(gp,ras_full_ini,'\tResult -- ras_full_ini: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,ras_ini,'\tResult -- ras_ini: ')
    PrintCellValue(gp,ras,'\tResult -- ras: ')
    PrintCellValue(gp,u_fri,'\tResult -- u_fri: ')
    PrintCellValue(gp,rah,'\tResult -- rah: ')

iteration = 1
Not_converged=True

while Not_converged:

  print iteration
   
   # rah is copied to rah_est
   #*******************************
  gp.CopyRaster(rah, rah_est)
  gp.CopyRaster(ETsoil_sec_ave, ETsoil_sec_pre)
  gp.CopyRaster(ETveg_sec_ave, ETveg_sec_pre)
  #*******************************

  gp.CopyRaster(soilm, soilm_M1)
  
  gp.SingleOutputMapAlgebra_sa(Albedo_con_soil_eqn,Albedo_soil)
  gp.SingleOutputMapAlgebra_sa(Albedo_con_veg_eqn,Albedo_veg_old)
  gp.SingleOutputMapAlgebra_sa(Albedo_con_new_veg_eqn, Albedo_veg)

 ###################################################################################################

  H_Flux_new_soil = os.path.join(workDir,'sheat_new_soil_'+ str(i)+'.img')
  G_Flux_new_soil = os.path.join(workDir,'gheat_new_soil_'+ str(i)+'.img')
  ux_new = os.path.join(workDir,'ux_new_'+ str(i)+'.img')
  
  H_Flux_eqn_soil= '( '+sheat_soil_2+' + '+H_Flux_rep_soil+' ) / 2'
  gp.SingleOutputMapAlgebra_sa(H_Flux_eqn_soil,H_Flux_new_soil)
  gp.CopyRaster(H_Flux_new_soil, H_Flux_rep_soil)

  G_Flux_eqn_soil= '( '+gheat_soil+' + '+G_Flux_rep_soil+' ) / 2'
  gp.SingleOutputMapAlgebra_sa(G_Flux_eqn_soil,G_Flux_new_soil)
  gp.CopyRaster(G_Flux_new_soil, G_Flux_rep_soil)

#############################################################################
  H_Flux_new_veg = os.path.join(workDir,'sheat_new_veg_'+ str(i)+'.img')

  H_Flux_eqn_veg= '( '+sheat_veg+' + '+H_Flux_rep_veg+' ) / 2'
  gp.SingleOutputMapAlgebra_sa(H_Flux_eqn_veg,H_Flux_new_veg)
  gp.CopyRaster(H_Flux_new_veg, H_Flux_rep_veg)
  #######################################################################################
  ux_eqn= '( '+u_fri+' + '+u_fri_neu+' ) / 2'
  gp.SingleOutputMapAlgebra_sa(ux_eqn,ux_new)
  gp.CopyRaster(ux_new, u_fri_neu)
##  ########################################################################################
  print 'Results of averaging parameters'
## Vegetation section
  
  gp.SingleOutputMapAlgebra_sa(ETveg_sec_con, ETveg_sec)
  gp.SingleOutputMapAlgebra_sa( ETveg_sec_eqn_ave,ETveg_sec_ave )
  gp.SingleOutputMapAlgebra_sa(tveg_eq, Tc)

  tclim_eq_veg = 'con( '+Tc+' < '+tslo+', '+tslo+', con( '+Tc+' > '+tshi+', '+tshi+', '+Tc+' ) )'
  limited = os.path.join(workDir, 'limited.img')
  errMin = (gp.GetRasterProperties(Tc, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(Tc, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(Tc, "MEAN").getoutput(0))
  
  if errMin < (float (tslo)) or errMax > (float(tshi)):
      print '\tAppling limits to Surface Temperature',
      gp.SingleOutputMapAlgebra_sa(tclim_eq_veg, limited)
      gp.CopyRaster(limited, Tc)

  gp.SingleOutputMapAlgebra_sa(eqn11_b , eoveg)
  gp.SingleOutputMapAlgebra_sa(Lambda_eqn_veg, Lambda_veg)
  gp.SingleOutputMapAlgebra_sa(outlwr_eq_veg , outlwr_veg)
  gp.SingleOutputMapAlgebra_sa(LE_eqn_veg, LE_veg)
  gp.SingleOutputMapAlgebra_sa(netrad_eqn_veg, netrad_veg)
  gp.SingleOutputMapAlgebra_sa(sheat_eqn_veg , sheat_veg )

  errMin = (gp.GetRasterProperties(sheat_veg, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(sheat_veg, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(sheat_veg, "MEAN").getoutput(0))
  if errMin < (shlo) or errMax > (shhi):
     gp.SingleOutputMapAlgebra_sa(shlim_eq_veg, limited)
     gp.CopyRaster(limited, sheat_veg)


  gp.SingleOutputMapAlgebra_sa(sheat_veg_final_eqn, sheat_veg_final)
####################################################################################################################
# Soil section
   
  
  gp.SingleOutputMapAlgebra_sa(tsurf_eq_soil, Ts)

  errMin = (gp.GetRasterProperties(Ts, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(Ts, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(Ts, "MEAN").getoutput(0))
  print 'Done\n\tTs -- Min: ',errMin,' Max: ',errMax, ' Ave: ',errAve

  if errMin < (float(tslo)) or errMax > (float(tshi)):
     print '\tAppling limits to Surface Temperature',
     gp.SingleOutputMapAlgebra_sa(tslim_eq_soil, limited)
     gp.CopyRaster(limited, Ts)

  gp.SingleOutputMapAlgebra_sa(eqn11_a , eosur)
  gp.SingleOutputMapAlgebra_sa(qosur_eqn , qosur)

  gp.SingleOutputMapAlgebra_sa(ETsoil_sec_con_eqn, ETsoil_sec)
  gp.SingleOutputMapAlgebra_sa( ETsoil_sec_eqn_ave,ETsoil_sec_ave )
  

  gp.SingleOutputMapAlgebra_sa(Lambda_eqn_soil, Lambda_soil)
  gp.SingleOutputMapAlgebra_sa(LE_eqn_soil,LE_soil)
  gp.SingleOutputMapAlgebra_sa(outlwr_eq_soil , outlwr_soil)
  gp.SingleOutputMapAlgebra_sa(netrad_eqn_soil, netrad_soil)
  gp.SingleOutputMapAlgebra_sa(sheat_eqn_1_soil, sheat_soil_1)
  gp.SingleOutputMapAlgebra_sa(congflux_soil, gheat_soil)
  #   #***************************
  gp.CopyRaster(gheat_soil, G_Flux_ite_soil)
  #***************************
  gp.SingleOutputMapAlgebra_sa(sheat_eqn_2_soil, sheat_soil_2)
  errMin = (gp.GetRasterProperties(sheat_soil_2, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(sheat_soil_2, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(sheat_soil_2, "MEAN").getoutput(0))

  if errMin < (float(shlo)) or errMax > (float(shhi)):
     gp.SingleOutputMapAlgebra_sa(shlim_eq_soil, limited)
     gp.CopyRaster(limited, sheat_soil_2)

  gp.SingleOutputMapAlgebra_sa(sheat_soil_final_eqn, sheat_soil_final)
  ##############################################################################################################################
 
  gp.SingleOutputMapAlgebra_sa(rss_cor_eqn,rss_cor)

  gp.SingleOutputMapAlgebra_sa(rss_con, rss)
  
  gp.SingleOutputMapAlgebra_sa(soilm_old_eqn, soilm_old)
  gp.SingleOutputMapAlgebra_sa(soilm_con, soilm )
#################################################################################################################################
# Combined section

  if DEBUG:
    PrintCellValue(gp,Albedo,'\tResult -- Albedo: ')
    PrintCellValue(gp,Albedo_soil,'\tResult -- Albedo_soil: ')
    PrintCellValue(gp,Albedo_veg,'\tResult -- Albedo_veg: ')
    PrintCellValue(gp,H_Flux_new_veg,'\tResult -- H_Flux_new_veg: ')
    PrintCellValue(gp,H_Flux_rep_veg,'\tResult -- H_Flux_rep_veg: ')
    PrintCellValue(gp,rah_est,'\tResult -- rah_est: ')
    PrintCellValue(gp,u_fri_neu,'\tResult -- u_fri_neu: ')
    PrintCellValue(gp,Tc,'\tResult -- Tc: ')

   
  gp.SingleOutputMapAlgebra_sa(sheat_eqn, sheat)
  
  ########################################################################################

  errMin = (gp.GetRasterProperties(sheat, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(sheat, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(sheat, "MEAN").getoutput(0))
  if errMin < (shlo) or errMax > (shhi):
     gp.SingleOutputMapAlgebra_sa(shlim_eq, limited)
     gp.CopyRaster(limited, sheat)

  gp.SingleOutputMapAlgebra_sa(eqn50, L)
   ####################################################
  gp.SingleOutputMapAlgebra_sa(eqn_X_30m, X_30m)
  gp.SingleOutputMapAlgebra_sa(conpsi_m_30m, psi_m_30m)
  gp.SingleOutputMapAlgebra_sa(conpsi_h_30m, psi_h_30m)
   #######################################################
  gp.SingleOutputMapAlgebra_sa(eqn_X_dzom,X_dzom)
  gp.SingleOutputMapAlgebra_sa(conpsi_h_dzom  , psi_h_dzom )
   ####################################################
  gp.SingleOutputMapAlgebra_sa(eqn_X_hd,X_hd)
  gp.SingleOutputMapAlgebra_sa(conpsi_h_hd,psi_h_hd)
  
  gp.SingleOutputMapAlgebra_sa(kh_eqn,kh)
  gp.SingleOutputMapAlgebra_sa(ras_bare_ini_eqn,ras_bare_ini)
  gp.SingleOutputMapAlgebra_sa(ras_full_ini_eqn,ras_full_ini)
  gp.SingleOutputMapAlgebra_sa(ras_ini_eqn,ras_ini)
  gp.SingleOutputMapAlgebra_sa(ras_con_eqn,ras)
   ###########################################################
  gp.SingleOutputMapAlgebra_sa(eqn56, u_fri)

  errMin = (gp.GetRasterProperties(u_fri, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(u_fri, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(u_fri, "MEAN").getoutput(0))
  u_fri_lim_eq = 'con( '+u_fri+' < '+u_fri_lo+', '+u_fri_lo+', con( '+u_fri+' > '+u_fri_hi+', '+u_fri_hi+', '+u_fri+' ) )'
  if errMin < (float(u_fri_lo)) or errMax > (float(u_fri_hi)):
     gp.SingleOutputMapAlgebra_sa(u_fri_lim_eq, limited)
     gp.CopyRaster(limited, u_fri)
  
  gp.SingleOutputMapAlgebra_sa(eqn57, rah)
  rahlim_eq = 'con( '+rah+' < '+rahlo+', '+rahlo+', con( '+rah+' > '+rahhi+', '+rahhi+', '+rah+' ) )'
  errMin = (gp.GetRasterProperties(rah, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(rah, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(rah, "MEAN").getoutput(0))

  if errMin < (float(rahlo)) or errMax > (float(rahhi)):
     
     print '\tAppling limits to aerodynamic resistance',
     gp.SingleOutputMapAlgebra_sa(rahlim_eq, limited)
     gp.CopyRaster(limited, rah)

  if DEBUG:

    print 'vegetation section'
   
    
    PrintCellValue(gp,H_Flux_rep_veg,'\tResult -- H_Flux_rep_veg: ')
    PrintCellValue(gp,ETveg_sec,'\tResult -- ETveg_sec: ')

    PrintCellValue(gp,Tc,'\tResult -- Tc: ')
 
    PrintCellValue(gp,eoveg,'\tResult -- eoveg: ')
    PrintCellValue(gp,Lambda_veg,'\tResult -- Lambda_veg: ')
    PrintCellValue(gp,LE_veg,'\tResult -- LE_veg: ')
    PrintCellValue(gp,outlwr_veg,'\tResult -- outlwr_veg: ')
    PrintCellValue(gp,netrad_veg,'\tResult -- netrad_veg: ')
    PrintCellValue(gp,sheat_veg,'\tResult -- sheat_veg: ')
    PrintCellValue(gp,sheat_veg_final,'\tResult -- sheat_veg_final: ')

    print 'soil section'
    PrintCellValue(gp,H_Flux_rep_soil,'\tResult -- H_Flux_rep_soil: ')
    PrintCellValue(gp,G_Flux_rep_soil,'\tResult -- G_Flux_rep_soil: ')
    PrintCellValue(gp,Ts,'\tResult -- Ts: ')
    PrintCellValue(gp,eosur,'\tResult -- eosur: ')
    PrintCellValue(gp,qosur,'\tResult -- qosur: ')
    PrintCellValue(gp,ETsoil_sec,'\tResult -- ETsoil_sec: ')
    
    PrintCellValue(gp,Lambda_soil,'\tResult -- Lambda_soil: ')
    PrintCellValue(gp,LE_soil,'\tResult -- LE_soil: ')
    PrintCellValue(gp,outlwr_soil,'\tResult -- outlwr_soil: ')
    PrintCellValue(gp,netrad_soil,'\tResult -- netrad_soil: ')
    PrintCellValue(gp,sheat_soil_1,'\tResult -- sheat_soil_1: ')
    PrintCellValue(gp,gheat_soil,'\tResult -- gheat_soil: ')
    PrintCellValue(gp,sheat_soil_2,'\tResult -- sheat_soil_2: ')
    PrintCellValue(gp,eosur,'\tResult -- eosur: ')
    PrintCellValue(gp,ras,'\tResult -- ras: ')
    PrintCellValue(gp,rah,'\tResult -- rah: ')
    PrintCellValue(gp,ea,'\tResult -- ea: ')
    PrintCellValue(gp,rss_cor,'\tResult -- rss_cor: ')
    PrintCellValue(gp,rss,'\tResult -- rss: ')
    PrintCellValue(gp,Soil_final,'\tResult -- Soil_final: ')
    PrintCellValue(gp,Theta_sat,'\tResult -- Theta_sat: ')
    PrintCellValue(gp,soilm_old,'\tResult -- soilm_old: ')
    PrintCellValue(gp,soilm,'\tResult -- soilm: ')

    print 'combined section'
    
    PrintCellValue(gp,H_Flux,'\tResult -- H_Flux: ')
    PrintCellValue(gp,G_Flux,'\tResult -- G_Flux: ')
    PrintCellValue(gp,sheat,'\tResult -- sheat: ')
    PrintCellValue(gp,sheat_rep,'\tResult -- sheat_rep: ')
    PrintCellValue(gp,L,'\tResult -- L: ')
    PrintCellValue(gp,X_30m,'\tResult -- X_30m: ')
    PrintCellValue(gp,psi_m_30m,'\tResult -- psi_m_30m: ')
    PrintCellValue(gp,psi_h_30m,'\tResult -- psi_h_30m: ')
    PrintCellValue(gp,H_Flux,'\tResult -- H_Flux: ')
    PrintCellValue(gp,G_Flux,'\tResult -- G_Flux: ')
    PrintCellValue(gp,sheat,'\tResult -- sheat: ')
    PrintCellValue(gp,X_dzom,'\tResult -- X_dzom: ')
    PrintCellValue(gp,psi_h_dzom,'\tResult -- psi_h_dzom: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,X_hd,'\tResult -- X_hd: ')
    PrintCellValue(gp,psi_h_hd,'\tResult -- psi_h_hd: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,kh,'\tResult -- kh: ')
    PrintCellValue(gp,ras_full_ini,'\tResult -- ras_full_ini: ')
    PrintCellValue(gp,ras_bare_ini,'\tResult -- ras_bare_ini: ')
    PrintCellValue(gp,ras,'\tResult -- ras: ')
    PrintCellValue(gp,u_fri,'\tResult -- u_fri: ')
    PrintCellValue(gp,rah,'\tResult -- rah: ')
#******************************************************************************
  errMin = (gp.GetRasterProperties(rah, "MINIMUM").getoutput(0))
  errMax = (gp.GetRasterProperties(rah, "MAXIMUM").getoutput(0))
  errAve = (gp.GetRasterProperties(rah, "MEAN").getoutput(0))
 
  rahlim_eq = 'con( '+rah+' < '+rahlo+', '+rahlo+', con( '+rah+' > '+rahhi+', '+rahhi+', '+rah+' ) )'

  if errMin < (float(rahlo)) or errMax > (float(rahhi)):
     gp.SingleOutputMapAlgebra_sa(rahlim_eq, limited)
     gp.CopyRaster(limited, rah)

  ####################################################################################################################################################
  gp.Minus_sa(rah, rah_est, raherror)
  gp.SingleOutputMapAlgebra_sa(iter_eq, stat_img)
  gp.ZonalStatisticsAsTable_sa(stat_img, "Value", raherror, stat_tab, "DATA")

  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  sql = '"VALUE" = 0'
  rows = gp.searchcursor(stat_tab,sql)
  row = rows.Next()
  if row:
     N_cells = row.COUNT
  else:
   
     N_cells = 0
     
  del row, rows
  sql = '"VALUE" = 1'
  rows = gp.searchcursor(stat_tab,sql)
  row = rows.Next()
  if row:
     C_cells = row.COUNT
  else:
     C_cells = 0
     
  del row, rows
  T_cells = C_cells + N_cells

  P_Convg = float(C_cells)/T_cells * 100.0000


  print 'Image Cells: ',T_cells, 'Converged Cells: ',C_cells,' Percent Converged: ',P_Convg

  print 'deleting final temporary down raster'
# Cleanup (remove) temporary files and rasters in the scratch workspace. 
  for root, dirs, files in os.walk(gp.ScratchWorkspace, topdown=False):
     
     
     for name in files:
        
        os.remove(os.path.join(root,name))
      
     for name in dirs:
        os.rmdir(os.path.join(root,name))
          
##  if((errMin < float(iter_lo) or errMax > float(iter_hi)) and P_Convg  < 99 and  iteration < 0):
  if (P_Convg  < 99 and  iteration < 15):
    
     Not_converged=True
     
     iteration += 1
     
  else:
     Not_converged=False
     print "Start of next 3 hour"
  
# Soil portion final

ETsoil_sec_final = os.path.join(workDir,'ETsoil_sec_final_'+ str(i)+'.img')
ETsoil_sec_final_eqn  = '(1 - '+fc+')  * '+ETsoil_sec

ETrF_soil = os.path.join(workDir,'ETrF_soil_'+ str(i)+'.img')
ETrF_soil_eqn= ETsoil_sec_final+' / '+Ref_ET

# Computation of soil evaporation in mm/hour
ETsoil_hour = os.path.join(workDir,'ETsoil_hour_'+ str(i)+'.img')
ETsoil_hour_eqn = '3600 * '+ETsoil_sec_final

gp.SingleOutputMapAlgebra_sa(ETsoil_sec_final_eqn,ETsoil_sec_final)
gp.SingleOutputMapAlgebra_sa( ETrF_soil_eqn,ETrF_soil )
gp.SingleOutputMapAlgebra_sa(ETsoil_hour_eqn, ETsoil_hour)

##PrintCellValue(gp,ETsoil_hour,'\tResult -- ETsoil_hour: ')
##PrintCellValue(gp,EThour_com,'\tResult -- EThour_com: ')

#_________________________________________________________________________________________________________________
# Vegetation protion
ETveg_sec_final = os.path.join(workDir,'ETveg_sec_final_'+ str(i)+'.img')
ETveg_sec_final_eqn  = fc+' * '+ETveg_sec

ETrF_veg = os.path.join(workDir,'ETrF_veg_'+ str(i)+'.img')
ETrF_veg_eqn=  ETveg_sec_final+' / '+Ref_ET

# Hourly transpiration in mm/hr
ETveg_hour = os.path.join(workDir,'ETveg_hour_'+ str(i)+'.img')
ETveg_hour_eqn = ETveg_sec_final+' * 3600'

gp.SingleOutputMapAlgebra_sa(ETveg_sec_final_eqn,ETveg_sec_final)
gp.SingleOutputMapAlgebra_sa( ETrF_veg_eqn,ETrF_veg )
gp.SingleOutputMapAlgebra_sa(ETveg_hour_eqn, ETveg_hour)

PrintCellValue(gp,ETveg_sec_final,'\tResult -- ETveg_sec_final: ')
PrintCellValue(gp,ETrF_veg,'\tResult -- ETrF_veg: ')
PrintCellValue(gp,ETveg_hour,'\tResult -- ETveg_hour: ')

#_____________________________________________________________________________________________________________________

# Combined values
# Total sensible heat flux is computed combining soil and vegetation fraction of sensible heat flux
EThour_com = os.path.join(workDir,'EThour_com_'+ str(i)+'.img')
EThour_com_1_eqn = ETveg_hour+' + '+ETsoil_hour
EThour_com_eqn = 'con('+Water_flag+' == 1, '+ETsoil_hour+', '+EThour_com_1_eqn+')'

# Combine LE in in W/m2
LE = os.path.join(workDir,'LE_'+ str(i)+'.img')
LE_1_eqn = LE_veg+' * '+fc+' + (1 - '+fc+') * '+LE_soil
LE_eqn = 'con('+Water_flag+' == 1, '+LE_soil+', '+LE_1_eqn+')'

# Combined temperature in K
T_com = os.path.join(workDir, 'T_com_'+ str(i)+'.img')
T_com_1_eqn= 'POW (( POW ( ( '+Tc+' ), 4) * '+fc+' +   POW ( ( '+Ts+' ), 4 ) * (1 - '+fc+') ) , 0.25 )'
T_com_eqn = 'con('+Water_flag+' == 1, '+Ts+', '+T_com_1_eqn+')'

# Combined net radiation in W/m2
netrad = os.path.join(workDir,'netrad_'+ str(i)+'.img')
netrad_1_eqn = netrad_veg+' * '+fc+' + (1 - '+fc+') * '+netrad_soil
netrad_eqn = 'con('+Water_flag+' == 1, '+netrad_soil+', '+netrad_1_eqn+')'

# Ground heat flux of soil poriton
gheat_com = os.path.join(workDir,'gheat_com_'+ str(i)+'.img')
gheat_com_eqn = 'con('+Water_flag+' == 1, '+gheat_soil+', '+gheat_soil+' * (1 - '+fc+'))'

gp.SingleOutputMapAlgebra_sa( EThour_com_eqn, EThour_com)
gp.SingleOutputMapAlgebra_sa(LE_eqn, LE)
gp.SingleOutputMapAlgebra_sa(netrad_eqn, netrad)
gp.SingleOutputMapAlgebra_sa(T_com_eqn, T_com)
gp.SingleOutputMapAlgebra_sa(gheat_com_eqn, gheat_com)

PrintCellValue(gp,EThour_com,'\tResult -- EThour_com: ')
PrintCellValue(gp,LE,'\tResult -- LE: ')
PrintCellValue(gp,T_com,'\tResult -- T_com: ')
PrintCellValue(gp,netrad,'\tResult -- netrad: ')
PrintCellValue(gp,gheat_com,'\tResult -- gheat_com: ')
#___________________________________________________________________________________________________________________
# Correction of water in vegetion portion

sheat_veg_wc = os.path.join(workDir,'sheat_veg_wc_'+ str(i)+'.img')
sheat_veg_wc_eqn= 'con ('+Water_flag+' == 1, '+sheat_soil_2+', '+sheat_veg+')'
gp.SingleOutputMapAlgebra_sa(sheat_veg_wc_eqn, sheat_veg_wc)

Tc_wc = os.path.join(workDir,'Tc_wc_'+ str(i)+'.img')
Tc_wc_eqn='con ('+Water_flag+' == 1, '+Ts+', '+Tc+')'
gp.SingleOutputMapAlgebra_sa(Tc_wc_eqn, Tc_wc)

LE_veg_wc = os.path.join(workDir,'LE_veg_wc_'+ str(i)+'.img')
LE_veg_wc_eqn= 'con ('+Water_flag+' == 1, '+LE_soil+', '+LE_veg+')'
gp.SingleOutputMapAlgebra_sa(LE_veg_wc_eqn, LE_veg_wc)

ETveg_hour_wc = os.path.join(workDir,'ETveg_hour_wc_'+ str(i)+'.img')
ETveg_hour_wc_eqn = 'con ('+Water_flag+' == 1, '+ETsoil_hour+', '+ETveg_hour+')'
gp.SingleOutputMapAlgebra_sa(ETveg_hour_wc_eqn, ETveg_hour_wc)

soilm_wc = os.path.join(workDir,'soilm_wc_'+ str(i)+'.img')
soilm_wc_eqn= 'con ('+Water_flag+' == 1, '+soilm_root_M1+', '+soilm+')'
gp.SingleOutputMapAlgebra_sa(soilm_wc_eqn, soilm_wc)

print 'time taken'
stop = time.time()
print 'Done ...', stop-start,'seconds'
