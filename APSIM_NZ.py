from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
from pandas import *
import datetime
from datetime import date
import numpy as np
import xml.etree.ElementTree as ET
import os
from pylab import *
import gdal
import shutil
import  matplotlib.pyplot as plt
import shapefile
import osr
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import LinearSegmentedColormap
from sklearn import linear_model, datasets
import matplotlib.colors as mcolors
from sklearn import cluster
from sklearn import preprocessing
import random
from sklearn.preprocessing import scale
params = {'mathtext.default': 'regular' }    
mpl.rcParams['axes.linewidth'] = 0.7      
plt.rcParams.update(params)
def Read():
    root_grp_MaxTemp_n = Dataset("V:\NZ\zippedDatasets\MaxTempCorr_VCSN_xaavh_1971-2000_North_Island.nc").variables['air_temperature']
    root_grp_MaxTemp_s = Dataset("V:\NZ\zippedDatasets\MaxTempCorr_VCSN_xaavh_1971-2000_South_Island.nc").variables['air_temperature']
    root_grp_MinTemp_n = Dataset("V:\NZ\zippedDatasets\MinTempCorr_VCSN_xaavh_1971-2000_North_Island.nc").variables['air_temperature']
    root_grp_MinTem_s = Dataset("V:\NZ\zippedDatasets\MinTempCorr_VCSN_xaavh_1971-2000_South_Island.nc").variables['air_temperature']
    root_grp_prep_n = Dataset("V:\NZ\zippedDatasets\TotalPrecipCorr_VCSN_xaavh_1971-2000_North_Island.nc").variables['precipitation_amount']
    root_grp_prep_s = Dataset("V:\NZ\zippedDatasets\TotalPrecipCorr_VCSN_xaavh_1971-2000_South_Island.nc").variables['precipitation_amount']
    
    Dates=date_range(date(1971,1,1), date(2000,12,31), freq='D')
    columns = ['MaxT', 'Mint','Prep']
    #north
    for row in range(root_grp_MaxTemp_n.shape[1]):
        for col in range(root_grp_MaxTemp_n.shape[2]):

            if not (root_grp_MaxTemp_n[0,row,col].data==-9999.0):
                data = np.array([root_grp_MaxTemp_n[:,row,col],root_grp_MinTemp_n[:,row,col],root_grp_prep_n[:,row,col]]).T
                df = DataFrame(data,columns=columns)
                df['Date']=Dates
                df['MaxT']=df['MaxT']-272.15
                df['Mint']=df['Mint']-272.15
                df['Date']=df['Date'].apply(lambda x:x.date())
                df[['Date','MaxT', 'Mint','Prep']].to_csv(r'V:\NZ\zippedDatasets\Finished\n_%d_%d.csv'%(row,col),index=False,float_format='%.2f')
                
                
    for row in range(root_grp_MaxTemp_s.shape[1]):
        for col in range(root_grp_MaxTemp_s.shape[2]):
            if not (root_grp_MaxTemp_s[0,row,col].data==-9999.0):
                data = np.array([root_grp_MaxTemp_s[:,row,col],root_grp_MinTem_s[:,row,col],root_grp_prep_s[:,row,col]]).T
                df = DataFrame(data,columns=columns)
                df['Date']=Dates
                df['MaxT']=df['MaxT']-272.15
                df['Mint']=df['Mint']-272.15
                df['Date']=df['Date'].apply(lambda x:x.date())

                df[['Date','MaxT', 'Mint','Prep']].to_csv(r'V:\NZ\zippedDatasets\Finished\s_%d_%d.csv'%(row,col),index=False,float_format='%.2f')  
def YearDayOfyearToDate(row):
    return (datetime.datetime(int(row['year']), 1, 1) + datetime.timedelta(int(row['day']) - 1)).date()

  
def APSIMMetToLintul():
    df=read_csv(r"\\simplace.ipf.uni-bonn.de\projects\NZ\Weather\LincolnNZ.met",delim_whitespace=True,skiprows=[0,1,2,3,4,6],header=0)
    df['Date']=df.apply(YearDayOfyearToDate, axis=1)
    df['Precipitation']=df['rain']
    df['TempMin']=df['mint']
    df['AirTemperatureMean']=df['mean']
    df['TempMax']=df['maxt']
    df['Radiation']=df['radn']*1000
    df['WindSpeed']=df['wind']
    df[['Date','Precipitation','TempMin','AirTemperatureMean','TempMax','Radiation','WindSpeed']].to_csv(r"V:\NZ\Lintul\Maize\LincolnNZ.csv",sep='\t',index=False)
def ExtractNewWeather():    
    factors=['radn','maxt','mint','rain']
    ncData={}
    
    
    Dates=date_range(date(1971,1,1), date(2000,12,31), freq='D')

    for factor in factors:
        for year in [1971,1981,1991]:
            ncData['%s_%d'%(factor,year)]=Dataset(r"V:\NZ\newWeather\%s_%d-%d.nc"%(factor,year,year+9))
    maxT=ncData['maxt_1971'].variables['air_temperature']
    for row in range(maxT.shape[1]):
        for col in range(maxT.shape[2]):
            if not (maxT[0,row,col].data==-9999.0):
                data=[]
                for factor in factors:
                    faDa=ncData['%s_%d'%(factor,1971)]
                    key=faDa.variables.keys()[3]
               
                    fData=ConverUnits(faDa.variables[key][:,row,col],factor)
                    print fData.shape
                    for year in [1981,1991]:
                        temData=ConverUnits(ncData['%s_%d'%(factor,year)].variables[key][:,row,col],factor)
                        print temData.shape
                        
                        fData=np.concatenate((fData,temData),axis=0)
                    data.append(fData)
                df = DataFrame(np.array(data).T,columns=factors)
                df['Date']=Dates
                df['year']=df['Date'].apply(lambda x:x.year)
                df['month']=df['Date'].apply(lambda x:x.month)
                
                df['day']=df['Date'].apply(lambda x:x.timetuple().tm_yday)
                
                df['avgt']=(df.maxt+df.mint)/2.0
                
                tav=df['avgt'].mean()
                amp=df[df.month==1].avgt.mean()-df[df.month==7].avgt.mean()
                PrepareAPSIMMetFromInterpolatedData(df,row,col,tav,amp)

    
def PrepareAPSIMMetFromInterpolatedData(df,row,col,tav,amp):
    head='''!title = met
[weather.met.weather]
latitude  = %.3f    (DECIMAL DEGREES)
longitude = %.3f   (DECIMAL DEGREES)
 tav =  %.2f (oC)     ! annual average ambient temperature
 amp =  %.2f (oC)     ! annual amplitude in mean monthly temperature

    year      day   radn    maxt    mint    rain
    ()        ()   (MJ/m2)  (oC)    (oC)    (mm)\n'''%(-34.375-0.05*row,166.425+0.05*col,tav,amp)
    fi=file(r'V:\NZ\APSIM\Mets\%d_%d.met'%(row,col),'w')
    fi.write(head)
    ar=np.array(df[['year','day','radn','maxt','mint','rain']])
    np.savetxt(fi,ar,fmt='%8d%8d%8.2f%8.2f%8.2f%8.2f')
    fi.close()
def ConvertToAPSIM(fN):
    df=read_csv(fN,parse_dates=[0],delimiter='\t')
    df['Year']=df['Date'].apply(lambda x:x.year)
    df['Month']=df['Date'].apply(lambda x:x.month)
    df['DOY']=df['Date'].apply(lambda x:x.timetuple().tm_yday)
    df['Radiation']=df['Radiation']/1000.0
    monthGroup=df.groupby('Month')['TempMean'].mean()
    return df[['Year','DOY','Radiation','TempMax','TempMin','Precipitation']],df['TempMean'].mean(),monthGroup[7]-monthGroup[1]
def convertMissingFilesToAPSIMMet():
    for fi in os.listdir(r'I:\MissingFiles'):
        fN=r'I:\MissingFiles\%s'%fi
        data,tav,amp=ConvertToAPSIM(fN)
        splitEles=fi.split('_')
        fiRe=int(splitEles[2][3:])
        fiCol,fiRow=splitEles[-1][1:-4].split('R')
        lon=5.8+float(fiRe)*float(fiCol)/110.567
        lat=52.5-float(fiRe)*float(fiRow)/110.567
        head='''!title = NRW
[weather.met.weather]
!station number = NRW
!station name =  %s                           
latitude  = %.2f    (DECIMAL DEGREES)
longitude = %.2f    (DECIMAL DEGREES)
tav =  %.2f (oC)     ! annual average ambient temperature
amp =  %.2f (oC)     ! annual amplitude in mean monthly temperature

year      day   radn    maxt    mint    rain
()        ()   (MJ/m2)  (oC)    (oC)    (mm)\n'''%(fi[:-4],lat,lon,tav,amp)
        apsimFi=file(r'V:\NZ\MacSur\mets\%s.met'%(fi[:-4]),'w')
        apsimFi.write(head)
        for row in data.iterrows():
            apsimFi.write('%8d%8d%8.2f%8.2f%8.2f%8.2f\n'%(row[1]['Year'],row[1]['DOY'],row[1]['Radiation'],row[1]['TempMax'],row[1]['TempMin'],row[1]['Precipitation']))
        apsimFi.close()
        
def ConverUnits(data,fac):

    if fac=='radn':
        return data*24*3600/1000000
    elif fac=='maxt' or fac=='mint':
        return data-273.15
    else:
        return data
def ProduceAPSIM(row,col):
    tree = ET.parse(r"I:\testapsim\GridCell.sim")
    for fi in os.listdir(r'V:\NZ\APSIM\Mets'):

        row,col=fi[:-4].split('_')
        
        root = tree.getroot()
        root.set('name','%s_%s'%(row,col))
        metComponent=root.findall('component')[1]

        metComponent.find('initdata/filename').text=r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\Mets\%s_%s.met'%(row,col)
        root.findall('system/component/initdata/outputfile')[0].text='%s_%s.out'%(row,col)
        tree.write(r'I:\testapsim\sims\%s_%s.sim'%(row,col))                  
def ExtractOutPuts():
    allFiles=[]
    for fi in os.listdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\simsH2'):
        if fi[-4:]=='.out':
# 
            b = os.path.getsize(r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\simsH2\%s'%fi)
#             print b
            if b>10:
                allFiles.append(fi)
    df=ReadAPSIMOutPuts_writeSpace(r"V:\NZ\APSIM\sims\%s"%allFiles[0])[['year','rotation_counter','SowingDOY','TotalBiomass','TotalYield']]
    row,col=allFiles[0][:-4].split('_')
    df['row']=row
    df['col']=col
    for fi in allFiles[1:]:
        print fi
        row,col=fi[:-4].split('_')

        tempdf=ReadAPSIMOutPuts_writeSpace(r"V:\NZ\APSIM\simsH2\%s_%s.out"%(row,col))

        tempdf['row']=row
        tempdf['col']=col        
        df=df.append(tempdf[['year','rotation_counter','SowingDOY','TotalBiomass','TotalYield','row','col']])

    df.to_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv',index=False)
            
def ReadAPSIMOutPuts(fN):
    return read_fwf(fN,widths=[15]*2+[17]+[20]+[15]*19+[16]+[15]*9,skiprows=[0,1,3],header=2,parse_dates=[0])
def ReadAPSIMOutPuts_writeSpace(fN):
    return read_csv(fN,delimiter=' ',delim_whitespace=True,parse_dates=[0],skiprows=[0,1,3])
def MapResultsYieldMean():
    df=read_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv')
    
    fac,unit='Yield',r'kg ha$^{-1}$ yr$^{-1}$'


    fig=plt.figure(figsize=(7,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=0.98, top=0.98,wspace=0.02, hspace=0.1)
    item_c=1
    f1 = fig.add_subplot(111)
    df=df[df.rotation_counter==1.0]
    dfmean=df.groupby(['row','col']).mean()
    dfmean=dfmean.reset_index()
    arr=np.zeros((dfmean.row.max()-dfmean.row.min()+1,dfmean.col.max()-dfmean.col.min()+1))
    print arr.shape
    arr[:,:]=-999
    f1.set_axis_off()
    print dfmean.row.min().astype(int),dfmean.col.min().astype(int)
    arr[dfmean.row.astype(int)-dfmean.row.min().astype(int),dfmean.col.astype(int)-dfmean.col.min().astype(int)]=dfmean.TotalYield
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=dfmean.TotalYield.min(),vmax=dfmean.TotalYield.max(), interpolation ='none')


    cax = fig.add_axes([0.2, -0.04, 0.6, 0.01])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)
    
    cbar.set_label('%s (%s)'%('Yield',unit),fontsize=8)
    #         plt.suptitle(fac,fontsize=12,x=0.7,y=1.1,transform = fig.transFigure)
    fig.savefig(r"V:\NZ\APSIM\outputs\totalYield_mean_maize2.png", dpi=300, bbox_inches='tight', pad_inches=0.0)
def MapResultsBiomassMean():
    df=read_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv')
    
    fac,unit='Yield',r'kg ha$^{-1}$ yr$^{-1}$'


    fig=plt.figure(figsize=(7,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=0.98, top=0.98,wspace=0.02, hspace=0.1)
    item_c=1
    f1 = fig.add_subplot(111)
    df=df[df.rotation_counter==1.0]
    dfmean=df.groupby(['row','col']).mean()
    dfmean=dfmean.reset_index()
    arr=np.zeros((dfmean.row.max()-dfmean.row.min()+1,dfmean.col.max()-dfmean.col.min()+1))
    print arr.shape
    arr[:,:]=-999
    f1.set_axis_off()
    print dfmean.row.min().astype(int),dfmean.col.min().astype(int)
    arr[dfmean.row.astype(int)-dfmean.row.min().astype(int),dfmean.col.astype(int)-dfmean.col.min().astype(int)]=dfmean.TotalBiomass

    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=dfmean.TotalYield.min(),vmax=dfmean.TotalYield.max(), interpolation ='none')


    cax = fig.add_axes([0.2, -0.04, 0.6, 0.01])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)
    
    cbar.set_label('%s (%s)'%('Biomass',unit),fontsize=8)
    #         plt.suptitle(fac,fontsize=12,x=0.7,y=1.1,transform = fig.transFigure)
    fig.savefig(r"V:\NZ\APSIM\outputs\totalBiomass_mean_maize2.png", dpi=300, bbox_inches='tight', pad_inches=0.0)  
def MapResultsBiomassMean_Wheat():
    df=read_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv')
    
    fac,unit='Yield',r'kg ha$^{-1}$ yr$^{-1}$'


    fig=plt.figure(figsize=(7,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=0.98, top=0.98,wspace=0.02, hspace=0.1)
    item_c=1
    f1 = fig.add_subplot(111)
    df=df[df.rotation_counter==2.0]
    dfmean=df.groupby(['row','col']).mean()
    dfmean=dfmean.reset_index()
    arr=np.zeros((dfmean.row.max()-dfmean.row.min()+1,dfmean.col.max()-dfmean.col.min()+1))
    print arr.shape
    arr[:,:]=-999
    f1.set_axis_off()
    print dfmean.row.min().astype(int),dfmean.col.min().astype(int)
    arr[dfmean.row.astype(int)-dfmean.row.min().astype(int),dfmean.col.astype(int)-dfmean.col.min().astype(int)]=dfmean.TotalBiomass

    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=dfmean.TotalBiomass.min(),vmax=dfmean.TotalBiomass.max(), interpolation ='none')


    cax = fig.add_axes([0.2, -0.04, 0.6, 0.01])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)
    
    cbar.set_label('%s (%s)'%('Biomass',unit),fontsize=8)
    #         plt.suptitle(fac,fontsize=12,x=0.7,y=1.1,transform = fig.transFigure)
    fig.savefig(r"V:\NZ\APSIM\outputs\totalBiomass_mean_wheat2.png", dpi=300, bbox_inches='tight', pad_inches=0.0)         
def MapResultsSowDOYMean_Wheat():
    df=read_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv')



    fig=plt.figure(figsize=(7,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=0.98, top=0.98,wspace=0.02, hspace=0.1)
    item_c=1
    f1 = fig.add_subplot(111)
    df=df[df.rotation_counter==2.0]
    dfmean=df.groupby(['row','col']).mean()
    dfmean=dfmean.reset_index()
    arr=np.zeros((dfmean.row.max()-dfmean.row.min()+1,dfmean.col.max()-dfmean.col.min()+1))
    print arr.shape
    arr[:,:]=-999
    f1.set_axis_off()
    print dfmean.row.min().astype(int),dfmean.col.min().astype(int)
    arr[dfmean.row.astype(int)-dfmean.row.min().astype(int),dfmean.col.astype(int)-dfmean.col.min().astype(int)]=dfmean.SowingDOY


    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=dfmean.SowingDOY.min(),vmax=dfmean.SowingDOY.max(), interpolation ='none')


    cax = fig.add_axes([0.2, -0.04, 0.6, 0.01])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)
    
    cbar.set_label('%s (Day of Year)'%('Sowing DOY'),fontsize=8)
    #         plt.suptitle(fac,fontsize=12,x=0.7,y=1.1,transform = fig.transFigure)
    fig.savefig(r"V:\NZ\APSIM\outputs\SowingDOY_mean_wheat2.png", dpi=300, bbox_inches='tight', pad_inches=0.0)  
def MapResultsSowDOYMean_maize():
    df=read_csv(r'V:\NZ\APSIM\outputs\CollectedOutput2.csv')
    fig=plt.figure(figsize=(7,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=0.98, top=0.98,wspace=0.02, hspace=0.1)
    item_c=1
    f1 = fig.add_subplot(111)
    df=df[df.rotation_counter==1]
    dfmean=df.groupby(['row','col']).mean()
    dfmean=dfmean.reset_index()
    arr=np.zeros((dfmean.row.max()-dfmean.row.min()+1,dfmean.col.max()-dfmean.col.min()+1))
    print arr.shape
    arr[:,:]=-999
    f1.set_axis_off()
    print dfmean.row.min().astype(int),dfmean.col.min().astype(int)
    arr[dfmean.row.astype(int)-dfmean.row.min().astype(int),dfmean.col.astype(int)-dfmean.col.min().astype(int)]=dfmean.SowingDOY


    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=dfmean.SowingDOY.min(),vmax=dfmean.SowingDOY.max(), interpolation ='none')


    cax = fig.add_axes([0.2, -0.04, 0.6, 0.01])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)
    
    cbar.set_label('%s (Day of Year)'%('Sowing DOY'),fontsize=8)
    #         plt.suptitle(fac,fontsize=12,x=0.7,y=1.1,transform = fig.transFigure)
    fig.savefig(r"V:\NZ\APSIM\outputs\SowingDOY_mean_maize2.png", dpi=300, bbox_inches='tight', pad_inches=0.0) 
def GetCol(field):
    splited=field[:-4].split('_')
    return int(splited[-1][1:].split('R')[0][:])
def GetRow(field):
    splited=field[:-4].split('_')
    return int(splited[-1].split('R')[-1][:])
def GetResolution(field):
    splited=field.split('_')
    return int(splited[2][3:])
def ProduceSimForMACSUR():
    sims=[fi for fi in os.listdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\simMACSUR') if fi[-4:]=='.sim']
    metFis=[fi[:-4]+'.met' for fi in os.listdir(r'I:\MissingFiles')]
    for sim in sims:
        print sim
        tree = ET.parse(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\simMACSUR\%s'%sim)
#         shutil.rmtree(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])
#         os.mkdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])

        for met in metFis:
            row,col,res=GetRow(met),GetCol(met),GetResolution(met)
            
            root = tree.getroot()
            root.set('name','%d_%d_%d'%(res,row,col))
            metComponent=root.findall('component')[1]
     
            metComponent.find('initdata/filename').text=r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\Mets\%s'%(met)
            root.findall('system/component/initdata/outputfile')[0].text='%d_%d_%d.out'%(res,row,col)
            tree.write(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s\%d_%d_%d.sim'%(sim[:-4],res,row,col))  
def ExtractOutPutsMACSUR():
    for fold in ['M_PLW','M_PLP','M_PLN','W_PLP','W_PLW','W_PLN']:
        allFiles_res={10:[],25:[],50:[],100:[]}
        folder=r"V:\NZ\MacSur\%s"%fold
        print folder
        secna=folder.split('\\')[-1]
        print fold
        for fi in os.listdir(folder):
            print fi
            if fi[-4:]=='.out':
                res,row,col=fi[:-4].split('_')
                res,row,col=int(res),int(row),int(col)
                allFiles_res[res].append(fi)
        for res in [10,25,50,100]:
            print res
            allFiles=allFiles_res[res]
            print allFiles
            df=ReadAPSIMOutPuts_writeSpace(r"%s\%s"%(folder,allFiles[0]))
            
            res,row,col=allFiles[0][:-4].split('_')
            df['row']=row
            df['col']=col
            df['res']=res
            for fi in allFiles[1:]:
                print fi
                res,row,col=fi[:-4].split('_')
        
                tempdf=ReadAPSIMOutPuts_writeSpace(r"%s\%s"%(folder,fi))
        
                tempdf['row']=row
                tempdf['col']=col   
                tempdf['res']=res     
                df=df.append(tempdf)
            
            df[['row','col','res']+df.columns[1:-3].tolist()].to_csv(r'V:\NZ\MacSur\outputs\%s_%s.csv'%(secna,res),index=False) 
def extracAPSIMMacsur_1x1():
    for crop in ['Wheat','Maize']:
        print crop
        folder=r"V:\NZ\MacSur\1x1%s"%crop
        for mana in ['PLP','PLW','PLN']:
            print mana
            allFis=[fi for fi in os.listdir(folder) if fi[-4:]=='.out' and fi.find(mana)>-1]
            cropab,secna,res,row,col=allFis[0][:-4].split('_')
            
            df=ReadAPSIMOutPuts_writeSpace(folder+'\\%s'%allFis[0])
            df['row']=row
            df['col']=col
            df['res']=res
            for fi in allFis[1:]:
                print fi
                cropab,secna,res,row,col=fi[:-4].split('_')
                fipath=folder+'\\%s'%fi
                tempdf=ReadAPSIMOutPuts_writeSpace(fipath)
                tempdf['row']=row
                tempdf['col']=col
                tempdf['res']=res
                df=df.append(tempdf)
            df[['row','col','res']+df.columns[1:-3].tolist()].to_csv(r'V:\NZ\MacSur\outputs\%s_%s_%s.csv'%(crop[0],mana,res),index=False)  
def extracAPSIMMacsur_1x1_maize():
    for crop in ['Wheat','Maize'][1:]:
        print crop
        folder=r"V:\NZ\MacSur\1x1%s"%crop
        for mana in ['PLP','PLW','PLN']:
            print mana
            allFis=[fi for fi in os.listdir(folder) if fi[-4:]=='.out' and fi.find(mana)>-1]
            cropab,secna,res,row,col=allFis[0][:-4].split('_')
            
            df=ReadAPSIMOutPuts_writeSpace(folder+'\\%s'%allFis[0])
            df['row']=row
            df['col']=col
            df['res']=res
            for fi in allFis[1:]:
                print fi
                cropab,secna,res,row,col=fi[:-4].split('_')
                fipath=folder+'\\%s'%fi
                tempdf=ReadAPSIMOutPuts_writeSpace(fipath)
                tempdf['row']=row
                tempdf['col']=col
                tempdf['res']=res
                df=df.append(tempdf)
            df[['row','col','res']+df.columns[1:-3].tolist()].to_csv(r'V:\NZ\MacSur\outputs\%s_%s_%s.csv'%(crop[0],mana,res),index=False)  
def rerange10_100():
    for fi in os.listdir(r'V:\NZ\MacSur\outputs - Copy'):
        res=fi[:-4].split('_')[-1]
        df=read_csv(r'V:\NZ\MacSur\outputs - Copy\%s'%fi)
        df['res']=res
        print df.columns
        df[['row','col','res']+df.columns[0:-3].tolist()].to_csv(r'V:\NZ\MacSur\outputs\%s'%(fi),index=False) 
def MapAll():
    folder=r'V:\NZ\MacSur\outputs'
    mapfolder=r'V:\NZ\MacSur\maps'
    files=os.listdir(folder)
    for fi in files:

        df=read_csv(r'%s\%s'%(folder,fi))
        for field in df.columns:
            if field!='row' and field!='col' and field!='year' :
                map(df,field,mapfolder,fi)
def map(df,field,mapfolder,fi):
    df=df.groupby(['row','col']).mean().reset_index()
    fig=plt.figure(figsize=(4,4))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
     
     
    f1 = fig.add_subplot(1,1,1)
    arr=np.zeros((df.row.max()-df.row.min()+1,df.col.max()+1))
    arr[:,:]=-999
    f1.set_axis_off()            
    arr[df.row-df.row.min(),df.col]=df[field].round(1)
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)

    image=f1.imshow(ArrayMasked,cmap=cm.RdYlGn,vmin=df[field].min(),vmax=df[field].max(),interpolation ='none')
    
    cax = fig.add_axes([0.1, -0.08, 0.7, 0.02])#
    cbar=fig.colorbar(image, cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5,size=5)
    cbar.set_label('%s'%(field),fontsize=8) 
    
    fig.savefig(r"%s\%s_%s.png"%(mapfolder,fi[:-4],field), dpi=300, bbox_inches='tight', pad_inches=0.0) 
    close()
def ClearResults():
    for fold in ['M_PLW','M_PLP','M_PLN','W_PLP','W_PLW','W_PLN']:
        folder=r"V:\NZ\MacSur\%s"%fold
        for fi in os.listdir(folder):
            if fi[-4:]=='.out':                   
                os.remove(folder+'\\'+fi) 
def ProduceSimForMACSUR_1x1():
    sims=[fi for fi in os.listdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\simMACSUR') if fi[-4:]=='.sim']
#     metFis=[fi for fi in os.listdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\Mets') if fi.find('_RES1_')>-1]
    metFis=[fi[:-4]+'.met' for fi in os.listdir(r'I:\MissingFiles')]#only the left ones
    for sim in sims:
        if sim.find('M_')>-1:
            print sim
            tree = ET.parse(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\simMACSUR\%s'%sim)
    #         shutil.rmtree(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])
    #         os.mkdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])
    
            for met in metFis:
                row,col,res=GetRow(met),GetCol(met),GetResolution(met)
                
                root = tree.getroot()
                root.set('name','%d_%d_%d'%(res,row,col))
                metComponent=root.findall('component')[1]
         
                metComponent.find('initdata/filename').text=r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\Mets\%s'%(met)
                root.findall('system/component/initdata/outputfile')[0].text='%s_%d_%d_%d.out'%(sim[:-4],res,row,col)
                tree.write(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\1x1Maize\%s_%d_%d_%d.sim'%(sim[:-4],res,row,col))            
        elif sim.find('W_')>-1:
            print sim
            tree = ET.parse(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\simMACSUR\%s'%sim)
    #         shutil.rmtree(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])
    #         os.mkdir(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\%s'%sim[:-4])
    
            for met in metFis:
                row,col,res=GetRow(met),GetCol(met),GetResolution(met)
                
                root = tree.getroot()
                root.set('name','%d_%d_%d'%(res,row,col))
                metComponent=root.findall('component')[1]
         
                metComponent.find('initdata/filename').text=r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\Mets\%s'%(met)
                root.findall('system/component/initdata/outputfile')[0].text='%s_%d_%d_%d.out'%(sim[:-4],res,row,col)
                tree.write(r'\\simplace.ipf.uni-bonn.de\projects\NZ\MacSur\1x1Wheat\%s_%d_%d_%d.sim'%(sim[:-4],res,row,col))                         
def ChangeDbf():

    e = shapefile.Editor(shapefile=r"V:\NZ\GIS_layers\LandUse\nzlri-land-use-capability")
    for s in e.records:

        if s[10]=='1' or s[10]=='2' or s[10]=='3' or s[10]=='4' or s[10]=='5':
            s[16]=1
        else:
            s[16]=0
    e.save(r"V:\NZ\GIS_layers\LandUse\nzlri-land-use-capability")
def creatLookUpTable():
    landuse=gdal.Open(r"V:\NZ\zoneRaster\landuse").ReadAsArray()
    pioneer=gdal.Open(r"V:\NZ\zoneRaster\zonePioneer.img").ReadAsArray()
    cliamte=Dataset(r"V:\NZ\newWeather\maxt_1971-1980.nc")
    mask=cliamte.variables['air_temperature'][0,:,:].mask

    data={}

    fi=file(r'v:\nz\row_col_zone.csv','w')
    fi.write('row,col,zone\n')
    for row in range(260):
        for col in range(243):


            if (not mask[row,col]) and landuse[row,col]==1 and pioneer[row,col]<100:
                fi.write('%d,%d,%d\n'%(row,col,pioneer[row,col]))

    fi.close()
def SimsForZoneSowring():
    tree = ET.parse(r"\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\GridCell3.sim")
    root = tree.getroot()
    df_row_col_zone=read_csv(r"\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\row_col_zone.csv")
    df_zone_cul=read_csv(r"\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\zone_cultivar.csv")
    
    
    for row,col,zone in np.array(df_row_col_zone):
        if row>=120:
            for sowdate,fir,sec in zip(range(1,6),['1-Sep','1-Oct','1-Nov','1-Dec','29-Dec'],['2-Sep','2-Oct','2-Nov','2-Dec','30-Dec']):
                for hb in range(1,4):
                    if not os.path.exists(r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\120_240\%d_%d_%d_%d_%d.sim'%(row,col,zone,sowdate,hb)):
                        try:
                            p1,p2= np.array(df_zone_cul[(df_zone_cul.Region==zone)&(df_zone_cul.SowDate==sowdate)&(df_zone_cul.Hybrid==hb)][['P1','P2']])[0]
            
                            root.set('name','%d_%d_%d_%d_%d'%(row,col,zone,sowdate,hb))
                            root.find(".//*[@name='met']/initdata/filename").text=r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\Mets\%d_%d.met'%(row,col)
                            root.find(".//*[@name='MainCrop']/initdata/ui/StartSow").text=fir
                            root.find(".//*[@name='MainCrop']/initdata/ui/EndSow").text=sec
                            root.find(".//*[@name='maize']/initdata/H_1/tt_emerg_to_endjuv").text='%d'%p1
                            root.find(".//*[@name='maize']/initdata/H_1/tt_flower_to_maturity").text='%d'%p2
                            root.findall('system/component/initdata/outputfile')[0].text='%d_%d_%d_%d_%d.out'%(row,col,zone,sowdate,hb)
                            
                            tree.write(r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\120_240\%d_%d_%d_%d_%d.sim'%(row,col,zone,sowdate,hb))  
                        except:
                            print r'\\simplace.ipf.uni-bonn.de\projects\NZ\APSIM\120_240\%d_%d_%d_%d_%d.sim'%(row,col,zone,sowdate,hb)
def Move180Then():
    fis=os.listdir(r'V:\NZ\APSIM\120_240')
    for fi in fis:
        if fi[-4:]=='.sim':
            if  not os.path.exists(r'V:\NZ\APSIM\120_240\%s.out'%(fi[:-4])):
                if int(fi.split('_')[0])>200:
                    shutil.move(r'V:\NZ\APSIM\120_240\%s'%fi, r'V:\NZ\APSIM\180_240') 
def moveBack():            
    fis=os.listdir(r'V:\NZ\APSIM\180_240')
    for fi in fis:
        if int(fi.split('_')[0])<200:
            shutil.move(r'V:\NZ\APSIM\180_240\%s'%fi, r'V:\NZ\APSIM\120_240')
def calcuateMeanRadAndT():
    luma=gdal.Open(r"V:\NZ\GIS_layers\LandUse\LandUseFilter.tif").ReadAsArray()
    variables=['MaxTempCorr','MinTempCorr','SurfRad','TotalPrecipCorr'][-1:]
    vaNs=['air_temperature','air_temperature','downwelling_shortwave_flux_in_air','precipitation_amount'][-1:]
    for variable,vaN in zip(variables,vaNs):
        amean=np.zeros((3,260,243))
        for count,year in enumerate(range(1971,2001,10)):
            print r"V:\NZ\weatherHistorical\ERA-40_revised(Nov2014)\%s_VCSN_xaavh_%d-%d.nc"%(variable,year,year+9)
            ds=Dataset(r"V:\NZ\weatherHistorical\ERA-40_revised(Nov2014)\%s_VCSN_xaavh_%d-%d.nc"%(variable,year,year+9))
            datayr=ds.variables[vaN][:,:,:]
            dmean= np.mean(datayr,axis=0)
            amean[count,:,:]=dmean
        tmean=np.mean(amean,axis=0)
        imshow(tmean);colorbar();show()
        if variable in variables[:2]:
            tmean=np.ma.masked_where((luma<3)|(tmean==0),tmean)
        else:
            tmean=np.ma.masked_where(luma<3,tmean)
        tmean=convertU(tmean,variable)
        tmean.dump(r'V:\NZ\weatherHistorical\summerized\%s.npy'%(variable))
def convertU(data,fac):
        if fac=='SurfRad':
            return data*24*3600/1000000
        elif fac=='MaxTempCorr' or fac=='MinTempCorr':
            return data-273.15
        elif fac=='TotalPrecipCorr':
            return data
def averageTmaxTmin():
    tmin=np.load("V:\NZ\weatherHistorical\summerized\MinTempCorr.npy")
    tmax=np.load("V:\NZ\weatherHistorical\summerized\MaxTempCorr.npy")
    T=(tmin+tmax)*0.5
    imshow(T);show()
    T.dump(r"V:\NZ\weatherHistorical\summerized\air_temperature.npy")

def Makexxyy():
    ds=gdal.Open(r"V:\NZ\GIS_layers\LandUse\FilterLU.tif")
    band=ds.GetRasterBand(1)
    no_data=band.GetNoDataValue()
    Arr=ds.ReadAsArray()

    imshow(Arr);show()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    
    xres = gt[1]
    yres = gt[5]
    print 2
    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
    ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
    ymax = gt[3] - yres * 0.5
     
    print xmin,xmax,xres
     
    ds = None
 
    # create a grid of xy coordinates in the original projection
    xy_source = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
    print xy_source.shape
    xx,yy=xy_source
    print xx.shape
    extent=[xmin,ymin,xmax,ymax]
    baseName=r"V:\NZ\GIS_layers\NZL_adm\NZL_adm2"
    m = Basemap(projection='laea',\
        llcrnrlat=extent[1],urcrnrlat=extent[3],\
        llcrnrlon=extent[0],urcrnrlon=extent[2],lat_0=(extent[1]+extent[3])/2,lon_0=(extent[0]+extent[2])/2,\
        resolution='l')
     
    inproj = osr.SpatialReference()
    inproj.ImportFromWkt(proj)
    print xy_source
     
    outproj = osr.SpatialReference()
    outproj.ImportFromProj4(m.proj4string)   
    print outproj
    xx, yy = convertXY(xy_source, inproj, outproj)

    np.save(r'V:\NZ\GIS_layers\xx.npy',xx)
    np.save(r'V:\NZ\GIS_layers\yy.npy',yy)    
def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the Climate object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy
def ClusterTandR():
    
    from sklearn import metrics
    dat=np.load(r"V:\NZ\weatherHistorical\summerized\air_temperature.npy")
    dat[dat.mask]=-999
    
    dar=np.load(r"V:\NZ\weatherHistorical\summerized\SurfRad.npy")
    dar*=365.25
    dar[dar.mask]=-999

    darain=np.load(r"V:\NZ\weatherHistorical\summerized\TotalPrecipCorr.npy")
    darain*=365.25
    darain[darain.mask]=-999
    
    ind=np.indices(dat.shape)
    valuesT=np.dstack((ind[0,...].flatten(),ind[1,...].flatten(),dat.flatten()))
    dfT=DataFrame(valuesT[0,...],columns=['row','col','TMean'])
    dfT=dfT[dfT.TMean!=-999]    
    valuesR=np.dstack((ind[0,...].flatten(),ind[1,...].flatten(),dar.flatten()))
    dfR=DataFrame(valuesR[0,...],columns=['row','col','Radn'])
    dfR=dfR[dfR.Radn!=-999] 
    
    valuesRain=np.dstack((ind[0,...].flatten(),ind[1,...].flatten(),darain.flatten()))
    dfRain=DataFrame(valuesRain[0,...],columns=['row','col','Rain'])

    dfRain=dfRain[dfRain.Rain!=-999]     
    
    df=dfT.set_index(['row','col']).join(dfR.set_index(['row','col'])).join(dfRain.set_index(['row','col'])).reset_index().dropna()
    print df
    k_means = cluster.KMeans(n_clusters=3, n_init=1)
    data= np.array(df['TMean'].reshape(-1, 1))
    k_means.fit(data)
    values = list(k_means.cluster_centers_.squeeze())
    labels = k_means.labels_
    orderedLabel= [sorted(values).index(i) for i in values] 
    labeldic=dict(zip(range(3),orderedLabel))
    df['TC']=labels
    df.TC=df.TC.apply(lambda x:labeldic[x])
    
    k_means.fit(df['Radn'].reshape(-1, 1))
    values = list(k_means.cluster_centers_.squeeze())
    labels = k_means.labels_
    orderedLabel= [sorted(values).index(i) for i in values] 
    labeldic=dict(zip(range(3),orderedLabel))
    df['RC']=labels
    df.RC=df.RC.apply(lambda x:labeldic[x])

    df['Combine']=''
    df.loc[(df.TC==0)&(df.RC==0),'Combine']='LTLR'
    df.loc[(df.TC==0)&(df.RC==1),'Combine']='LTMR'
    df.loc[(df.TC==0)&(df.RC==2),'Combine']='LTHR'
    df.loc[(df.TC==1)&(df.RC==0),'Combine']='MTLR'
    df.loc[(df.TC==1)&(df.RC==1),'Combine']='MTMR'
    df.loc[(df.TC==1)&(df.RC==2),'Combine']='MTHR'
    df.loc[(df.TC==2)&(df.RC==0),'Combine']='HTLR'
    df.loc[(df.TC==2)&(df.RC==1),'Combine']='HTMR'
    df.loc[(df.TC==2)&(df.RC==2),'Combine']='HTHR'
    mapDic={'LTLR':0,'LTMR':1,'LTHR':2,'MTLR':3,'MTMR':4,'MTHR':5,'HTLR':6,'HTMR':7,'HTHR':8}
    df['CN']=df.Combine.apply(lambda x:mapDic[x])
    df.to_csv(r'V:\NZ\weatherHistorical\summerized\Cluster.csv',index=False)
#     arr=np.zeros((260,243))
#     arr[:,:]=-1
#     arr[df.row,df.col]=df.CN
#     arr=np.ma.masked_where(arr<0,arr)
#     imshow(arr,cmap=plt.cm.get_cmap('RdBu_r',9));colorbar();show()
def APSIM_Cluster_file():
    dfc=read_csv(r'V:\NZ\weatherHistorical\summerized\Cluster.csv')
    dfb=read_csv(r"V:\NZ\sowByHybridResults_2015\APSIM_Res.csv")
    f = {'TotalBiomass':['std','mean'],'HarvestIndex':['std','mean']}
    dfb=dfb.groupby(['row','col','thisHyb','thisSow','thisLong','thisLat']).agg(f)
    dfb.columns = ['_'.join(col).strip() for col in dfb.columns.values]
    dfb=dfb.reset_index()
    print dfb
    dfb.row-=1
    dfb.col-=1#r start counting from 1 while python from 0
    df=dfb.set_index(['row','col']).join(dfc.set_index(['row','col'])).dropna()
    df.to_csv(r'V:\NZ\sowByHybridResults_2015\RES_CLuster.csv')
def SD_MeanacrossSH():
    dfb=read_csv(r"V:\NZ\sowByHybridResults_2015\APSIM_Res.csv")
    dfb=dfb.groupby(['row','col','thisHyb','thisSow','thisLong','thisLat']).mean().reset_index()
    f = {'TotalBiomass':['std','mean'],'HarvestIndex':['std','mean']}
    dfb=dfb.groupby(['row','col']).agg(f)
    dfb.columns = ['_'.join(col).strip() for col in dfb.columns.values]
    dfb=dfb.reset_index()
    dfb['TB_cv']=dfb.TotalBiomass_std/dfb.TotalBiomass_mean*100
    dfb['HI_cv']=dfb.HarvestIndex_std/dfb.HarvestIndex_mean*100
    dfb.row-=1
    dfb.col-=1#r start counting from 1 while python from 0
    dfb.to_csv(r'V:\NZ\sowByHybridResults_2015\CVAcrossSowHy.csv',index=False)
def CountourAll_meanBiomass():
    df=read_csv(r'V:\NZ\sowByHybridResults_2015\RES_CLuster.csv')
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    df=df.groupby(['thisHyb','thisSow','TC','RC','CN']).mean().reset_index()
    fig=plt.figure(figsize=(7.5,7.5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.07, hspace=0.05)
    item_c=1

        
    yi, xi = np.mgrid[slice(1, 5 + 1, 1),
                            slice(1, 5 + 1, 1)]
    for TC,tcN in zip(range(3),['LT','MT','HT']):
        for RC,rcN in zip(range(3),['LR','MR','HR']):
            f1=fig.add_subplot(3,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
            f1.set_xticks(range(1,6))
            f1.set_yticks(range(1,6))
#             f1.set_ylim((0,6))
            if TC==0:
                f1.set_title(rcN,fontsize=10)
            if RC==0:
                f1.set_ylabel(tcN,fontsize=10)
            if TC<2:
                f1.set_xticklabels(())
            else:
                f1.set_xticklabels(('s1','s2','s3','s4','s5'))
            if RC>0:
                f1.set_yticklabels(())
            else:
                f1.set_yticklabels(('h1','h2','h3','h4','h5'))
            dftr=df[(df.TC==TC)&(df.RC==RC)]
            dftr=dftr[['HB','Sow','TotalBiomass_mean']]
            dftr=dftr.pivot(values='TotalBiomass_mean',columns='Sow',index='HB')
#             dftr=dftr.iloc[::-1]
            zi=np.array(dftr)/1000
            CS = f1.contourf(xi,yi,zi,20,cmap=plt.cm.RdYlGn,alpha=1,vmin=10,vmax=26,levels=np.arange(10,27,1),extend='both')
            con=f1.contour(CS,levels=np.arange(10,27,2),colors='k')
            plt.clabel(con, levels=np.arange(10,27,2), inline=True,fmt='%1.1f',fontsize=10,colors='black')
    cax = fig.add_axes([0,-0.07,1,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal')
    cbar.ax.tick_params(labelsize=10,size=11)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Total\ biomass\ (t\ ha^{-1})$',fontsize=12)
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\CountourAll_meanBiomass.png", dpi=300, bbox_inches='tight', pad_inches=0.05) 
def CountourAll_meanHI():
    df=read_csv(r'V:\NZ\sowByHybridResults_2015\RES_CLuster.csv')
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    df=df.groupby(['thisHyb','thisSow','TC','RC','CN']).mean().reset_index()
    fig=plt.figure(figsize=(7.5,7.5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.07, hspace=0.05)
    item_c=1

        
    yi, xi = np.mgrid[slice(1, 5 + 1, 1),
                            slice(1, 5 + 1, 1)]
    for TC,tcN in zip(range(3),['LT','MT','HT']):
        for RC,rcN in zip(range(3),['LR','MR','HR']):
            f1=fig.add_subplot(3,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
            f1.set_xticks(range(1,6))
            f1.set_yticks(range(1,6))
#             f1.set_ylim((0,6))
            if TC==0:
                f1.set_title(rcN,fontsize=10)
            if RC==0:
                f1.set_ylabel(tcN,fontsize=10)
            if TC<2:
                f1.set_xticklabels(())
            else:
                f1.set_xticklabels(('s1','s2','s3','s4','s5'))
            if RC>0:
                f1.set_yticklabels(())
            else:
                f1.set_yticklabels(('h1','h2','h3','h4','h5'))
            dftr=df[(df.TC==TC)&(df.RC==RC)]
            dftr=dftr[['HB','Sow','HarvestIndex_mean']]
            dftr=dftr.pivot(values='HarvestIndex_mean',columns='Sow',index='HB')
#             dftr=dftr.iloc[::-1]
            zi=np.array(dftr)
            CS = f1.contourf(xi,yi,zi,20,cmap=plt.cm.RdYlGn,alpha=1,vmin=0.1,vmax=0.6,levels=np.arange(0.1,0.6,0.05),extend='both')
            con=f1.contour(CS,levels=np.arange(0.1,0.6,0.1),colors='k')
            plt.clabel(con, levels=np.arange(0.1,0.6,0.1), inline=True,fmt='%1.1f',fontsize=10,colors='black')
    cax = fig.add_axes([0,-0.07,1,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal')
    cbar.ax.tick_params(labelsize=10,size=11)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Harvest\ index$',fontsize=12)
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\CountourAll_meanHI.png", dpi=300, bbox_inches='tight', pad_inches=0.05) 
def CountourAll_stdBiomass():
    df=read_csv(r'V:\NZ\sowByHybridResults_2015\RES_CLuster.csv')
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    df=df.groupby(['thisHyb','thisSow','TC','RC','CN']).mean().reset_index()
    fig=plt.figure(figsize=(7.5,7.5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.07, hspace=0.05)
    item_c=1

        
    yi, xi = np.mgrid[slice(1, 5 + 1, 1),
                            slice(1, 5 + 1, 1)]
    for TC,tcN in zip(range(3),['LT','MT','HT']):
        for RC,rcN in zip(range(3),['LR','MR','HR']):
            f1=fig.add_subplot(3,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
            f1.set_xticks(range(1,6))
            f1.set_yticks(range(1,6))
#             f1.set_ylim((0,6))
            if TC==0:
                f1.set_title(rcN,fontsize=10)
            if RC==0:
                f1.set_ylabel(tcN,fontsize=10)
            if TC<2:
                f1.set_xticklabels(())
            else:
                f1.set_xticklabels(('s1','s2','s3','s4','s5'))
            if RC>0:
                f1.set_yticklabels(())
            else:
                f1.set_yticklabels(('h1','h2','h3','h4','h5'))
            dftr=df[(df.TC==TC)&(df.RC==RC)]
            dftr=dftr[['HB','Sow','TotalBiomass_std']]
            dftr=dftr.pivot(values='TotalBiomass_std',columns='Sow',index='HB')
#             dftr=dftr.iloc[::-1]
            zi=np.array(dftr)/1000
            CS = f1.contourf(xi,yi,zi,20,cmap=plt.cm.RdYlGn_r,alpha=1,vmin=1,vmax=3,levels=np.arange(1,3.1,0.2),extend='both')
            con=f1.contour(CS,levels=np.arange(1,3.1,0.4),colors='k')
            plt.clabel(con, levels=np.arange(1,3.1,0.4), inline=True,fmt='%1.1f',fontsize=10,colors='black')
    cax = fig.add_axes([0,-0.07,1,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal')
    cbar.ax.tick_params(labelsize=10,size=11)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Standard\ deviation\ of\ total\ biomass\ (t\ ha^{-1})$',fontsize=12)
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\CountourAll_stdBiomass.png", dpi=300, bbox_inches='tight', pad_inches=0.05) 
def GetBase():
    xx, yy = np.load(r"V:\NZ\GIS_layers\xx.npy"),np.load(r"V:\NZ\GIS_layers\yy.npy")

    ds=gdal.Open(r"V:\NZ\GIS_layers\LandUse\FilterLU.tif")
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    
    xres = gt[1]
    yres = gt[5]
    band=ds.GetRasterBand(1)
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
    ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
    ymax = gt[3] - yres * 0.5    
    no_data=band.GetNoDataValue()
    cliArr=ds.ReadAsArray()
    extent=[-10.662957,34.803539,67,71.185448]
    baseName=r"V:\NZ\GIS_layers\NZL_adm\NZL_adm1"
    extent=[xmin,ymin,xmax,ymax]

    m = Basemap(projection='tmerc',\
        llcrnrlat=extent[1],urcrnrlat=extent[3],\
        llcrnrlon=extent[0],urcrnrlon=extent[2],lat_0=(extent[1]+extent[3])/2,lon_0=(extent[0]+extent[2])/2,\
        resolution='l')
    ind=np.indices(cliArr.shape)
    values=np.dstack((ind[0,...].flatten(),ind[1,...].flatten(),cliArr.flatten()))
    df_zone=DataFrame(values[0,...],columns=['row','col','Value'])
    df_zone=df_zone[df_zone.Value!=no_data]    
    df_zone=df_zone[df_zone.Value==3] 
    return df_zone,cliArr,m,xx,yy,baseName            
def mapclusters():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\weatherHistorical\summerized\Cluster.csv")
    fig=plt.figure(figsize=(7.5,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    item_c=1
    mapDic={'LTLR':0,'LTMR':1,'LTHR':2,'MTLR':3,'MTMR':4,'MTHR':5,'HTLR':6,'HTMR':7,'HTHR':8}
    mapDic=dict (zip(mapDic.values(),mapDic.keys()))
    for clu,level,fac,tit,Types,unit in zip(['TC','RC','CN'],[3,3,9],[['TMean'],['Radn'],['TMean','Radn']],['Temperature','Global radiation','Temperature and \nglobal radiation'],[['LT','MT','HT'],['LR','MR','HR'],[mapDic[k] for k in range(9)]],[r'$^\circ$C',r'MJ m$^{-2}$ yr$^{-1}$','']):
        print clu
        f1=fig.add_subplot(1,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
        arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
        arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc[clu]
        ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
        m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r')  #,norm=MidPointNorm()
        m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        f1.text(0.0,0.8,tit,fontsize=8,transform=f1.transAxes)
        pats=[]
        for i,Type in enumerate(Types):
            red_patch = mpatches.Patch(color=cm.get_cmap('RdYlGn_r', len(Types))(i), label=u'%s'%Type.decode('utf-8', "replace"))
            pats.append(red_patch)
        plt.legend(handles=pats,ncol=1,fontsize=6.5,bbox_to_anchor=(0.27,0.38,0,0.37),markerscale=0.3,frameon=True) 

        if clu=='TC':
            f1pos=f1.get_position()
            f2 = fig.add_axes([0.22,0.23,0.11,0.18]);f2.tick_params(labelsize=6,size=2);#
            data=[]
            for co in range(level):
                data.append(df.loc[df[clu]==co,fac].values)
            f2.boxplot(data)
            f2.set_xticklabels(Types)
            f2.set_ylabel(r'%s (%s)'%(tit,unit),fontsize=6) 
        elif clu=='RC':
            f1pos=f1.get_position()
            f2 = fig.add_axes([0.56,0.23,0.11,0.16]);f2.tick_params(labelsize=6,size=2);#
            data=[]
            for co in range(level):
                data.append(df.loc[df[clu]==co,fac].values)
            f2.boxplot(data)
            f2.set_xticklabels(Types)
            f2.set_ylabel(r'%s (%s)'%(tit,unit),fontsize=6) 
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\Cluster.png", dpi=300, bbox_inches='tight', pad_inches=0.0)
def mapTRRain():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\weatherHistorical\summerized\Cluster.csv")
    fig=plt.figure(figsize=(7.5,4))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    item_c=1

    for clu,tit,unit,vmin,vmax,cmap in zip(['TMean','Radn','Rain'],['Temperature','Global radiation','Precipitation'],[r'$^\circ$C',r'MJ m$^{-2}$ yr$^{-1}$',r'mm yr$^{-1}$'],[8,3200,400],[16,6000,1200],['RdYlGn_r','RdYlGn_r','RdYlGn']):
        print clu
        f1=fig.add_subplot(1,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
        arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
        arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc[clu]
        ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
        print xx
        image=m.pcolor(xx,yy,ArrayMasked.T,alpha=1,cmap=cmap,vmin=vmin,vmax=vmax)  #,norm=MidPointNorm()
        m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        f1pos=f1.get_position()
        cax = fig.add_axes([f1pos.x0+0.01,f1pos.y0+0.37,0.01,(f1pos.y1-f1pos.y0)*0.55])#
        cbar=colorbar(image,cax,orientation='vertical',extend='both')
        cbar.ax.tick_params(labelsize=7,size=5)
        cbar.locator = MaxNLocator( nbins = 4)
     
        cbar.ax.get_xaxis().labelpad = -40
        cbar.set_label(r'%s (%s)'%(tit,unit),fontsize=12)       

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\TRRain.png", dpi=300, bbox_inches='tight', pad_inches=0.0)
def mapTRRain_Contour():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\weatherHistorical\summerized\Cluster.csv")
    fig=plt.figure(figsize=(7.5,4))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    item_c=1
    import scipy.ndimage as ndimage
    for clu,tit,unit,vmin,vmax,level,cmap in zip(['TMean','Radn','Rain'],['Temperature','Global radiation','Precipitation'],[r'$^\circ$C',r'MJ m$^{-2}$ yr$^{-1}$',r'mm yr$^{-1}$'],[8,3200,400],[16,6000,1200],[0.5,300,50],['RdYlGn_r','RdYlGn_r','RdYlGn']):
        print clu
        f1=fig.add_subplot(1,3,item_c);item_c+=1;f1.tick_params(labelsize=8,size=3);
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
        arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
        arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc[clu]
        
        ArrayMasked = np.ma.masked_where(arr == -999.0, arr)

        image=m.pcolor(xx.T,yy.T,ArrayMasked,cmap=cm.get_cmap(cmap, 20)  ,vmin=vmin,vmax=vmax,latlon=True)  #,norm=MidPointNorm()
        m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)


        f1pos=f1.get_position()
        cax = fig.add_axes([f1pos.x0+0.01,f1pos.y0+0.37,0.01,(f1pos.y1-f1pos.y0)*0.55])#
        cbar=colorbar(image,cax,orientation='vertical',extend='both')
        cbar.ax.tick_params(labelsize=7,size=5)
        cbar.locator = MaxNLocator( nbins = 4)
     
        cbar.ax.get_xaxis().labelpad = -40
        cbar.set_label(r'%s (%s)'%(tit,unit),fontsize=12)       

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\TRRain_T.png", dpi=500, bbox_inches='tight', pad_inches=0.0)
def mapTheBestCombination():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_mean==x.TotalBiomass_mean.max()])
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.5,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,1,1)
    print df
    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.04+0.061*(hb-1),0.65,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.37,0.38,0,0.37),markerscale=0.3,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\Best_combine.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapTheBestCombination_std_biomass():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_std==x.TotalBiomass_std.max()])
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.5,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,1,1)
    print df
    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.04+0.061*(hb-1),0.65,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.37,0.38,0,0.37),markerscale=0.3,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\mapTheBestCombination_std_biomass_combine.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapTheBestCombination_std_biomass_min():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_std==x.TotalBiomass_std.min()])
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.5,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,1,1)
    print df
    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.04+0.061*(hb-1),0.65,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.37,0.38,0,0.37),markerscale=0.3,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\mapTheBestCombination_std_biomass_min.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapTheBestCombination_std():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    dfh=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_mean==x.TotalBiomass_mean.max()])
    df=dfh.set_index(['row','col','thisHyb','thisSow']).join(df.set_index(['row','col','thisHyb','thisSow']),rsuffix='_or').dropna().reset_index()
    df= df[['row','col','thisHyb','thisSow','TotalBiomass_std_or']]
    df.row=df.row.astype(int)
    df.col=df.col.astype(int)
    fig=plt.figure(figsize=(6,7.5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,1,1)

    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['TotalBiomass_std_or']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)/1000
#     imshow(ArrayMasked);show()
    print np.min(ArrayMasked),np.max(ArrayMasked)
    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=1,vmax=2) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.1,-0.07,0.8,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal')
    cbar.ax.tick_params(labelsize=10,size=10)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Standard\ deviation\ of\ total\ biomass\ (t\ ha^{-1})$',fontsize=12)
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\Best_combine_std.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapBestCombinationAndSTD():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_mean==x.TotalBiomass_mean.max()])
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.8,5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,3,1)

    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    f1.text(0.8,0.03,'(a)', fontsize=8,transform=f1.transAxes)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.08+0.095*(hb-1),0.63,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.6,0.42,0,0.37),markerscale=0.1,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    f1=fig.add_subplot(1,3,2)
    f1.text(0.8,0.03,'(b)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['TotalBiomass_mean']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)/1000
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn',vmin=18,vmax=30) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.34,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Total\ biomass\ (t\ ha^{-1})$',fontsize=12)    
    
    f1=fig.add_subplot(1,3,3)
    f1.text(0.8,0.03,'(c)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['TotalBiomass_std']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)/1000
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=1,vmax=2) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.675,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$SD\ of\ total\ biomass\ (t\ ha^{-1})$',fontsize=12)    
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\Best_combine_and_STD.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapBestCombinationAndCV_Biomass():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df['TotalBiomass_cv']=df['TotalBiomass_std']/df['TotalBiomass_mean']*100
#     df['HarvestIndex_cv']=df['HarvestIndex_std']/df['HarvestIndex_mean']*100
    df=df.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_mean==x.TotalBiomass_mean.max()])
    print df.TotalBiomass_cv.max()
    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.8,5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,3,1)

    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    f1.text(0.8,0.03,'(a)', fontsize=8,transform=f1.transAxes)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.08+0.095*(hb-1),0.63,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.6,0.42,0,0.37),markerscale=0.1,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    f1=fig.add_subplot(1,3,2)
    f1.text(0.8,0.03,'(b)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['TotalBiomass_mean']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)/1000
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn',vmin=18,vmax=30) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.34,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$Total\ biomass\ (t\ ha^{-1})$',fontsize=12)    
    
    f1=fig.add_subplot(1,3,3)
    f1.text(0.8,0.03,'(c)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['TotalBiomass_cv']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=5,vmax=10) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.675,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$CV(\%)$',fontsize=12)    
    fig.savefig(r"V:\NZ\sowByHybridResults_2015\mapBestCombinationAndCV.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapBestCombinationAndSTD_HI():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df=df.groupby(['row','col']).apply(lambda x:x[x.HarvestIndex_mean==x.HarvestIndex_mean.max()])

    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.8,5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,3,1)

    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    f1.text(0.8,0.03,'(a)', fontsize=8,transform=f1.transAxes)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.08+0.095*(hb-1),0.63,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.6,0.42,0,0.37),markerscale=0.1,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    f1=fig.add_subplot(1,3,2)
    f1.text(0.8,0.03,'(b)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['HarvestIndex_mean']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn',vmin=0.4,vmax=0.55) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.34,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$HI\ (t\ t^{-1})$',fontsize=12)    
    
    f1=fig.add_subplot(1,3,3)
    f1.text(0.8,0.03,'(c)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['HarvestIndex_std']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=0,vmax=0.1) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.675,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$SD\ of\ HI\ (t\ t^{-1})$',fontsize=12)    

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\Best_combine_and_STD_HI.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapBestCombinationAndCV_HI():
    import matplotlib.patches as mpatches
    df=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")
    df['HarvestIndex_cv']=df['HarvestIndex_std']/df['HarvestIndex_mean']*100
    df=df.groupby(['row','col']).apply(lambda x:x[x.HarvestIndex_mean==x.HarvestIndex_mean.max()])

    df['HB']=df.thisHyb.apply(lambda x:int(x[1:]))
    df['Sow']=df.thisSow.apply(lambda x:int(x[1:]))
    fig=plt.figure(figsize=(7.8,5))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)
    f1=fig.add_subplot(1,3,1)


    pats=[]
    mycolors=['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000']
    from matplotlib import colors
    cmap = colors.ListedColormap(['#80CA33', '#308DBF','#5D30BF','#E0A738','#FF0000'])
    bounds=range(1,6)
    f1.text(0.8,0.03,'(a)', fontsize=8,transform=f1.transAxes)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    for hb in range(1,6):
        print hb
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfs=df[df.HB==hb]
        if dfs.Sow.count()>0:
            dfzc=df_zone.set_index(['row','col']).join(dfs.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
            arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
            arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Sow']
            ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
            m.pcolormesh(xx,yy,ArrayMasked.T,alpha=hb*0.2,cmap=cmap, norm=norm,vmin=1,vmax=5) 
            m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
        labels=[]
        for i in range(1,6):
            if hb==5:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2,label='s%d'%i)
            else:
                red_patch = mpatches.Patch(color=mycolors[i-1],alpha=hb*0.2+0.11)
            pats.append(red_patch)
        f1.text(0.08+0.095*(hb-1),0.63,'h%d'%hb,fontsize=8,transform=f1.transAxes)
    plt.legend(handles=pats,ncol=5,fontsize=6.5,bbox_to_anchor=(0.6,0.42,0,0.37),markerscale=0.1,frameon=False,columnspacing=0,labelspacing=0,borderaxespad=0.1) 
    
    f1=fig.add_subplot(1,3,2)
    f1.text(0.8,0.03,'(b)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['HarvestIndex_mean']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn',vmin=0.4,vmax=0.5) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.34,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$HI\ (t\ t^{-1})$',fontsize=12)    
    
    f1=fig.add_subplot(1,3,3)
    f1.text(0.8,0.03,'(c)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['HarvestIndex_cv']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=5,vmax=10) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.675,0.4,0.01,0.45])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$CV(\%)$',fontsize=12)      

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\mapBestCombinationAndCV_HI.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def mapBestCombinationAndSTD_HI_Delta():
    import matplotlib.patches as mpatches
    dfa=read_csv(r"V:\NZ\sowByHybridResults_2015\RES_CLuster.csv")

    dfa['HB']=dfa.thisHyb.apply(lambda x:int(x[1:]))
    dfa['Sow']=dfa.thisSow.apply(lambda x:int(x[1:]))
    dfb=dfa.groupby(['row','col']).apply(lambda x:x[x.TotalBiomass_mean==x.TotalBiomass_mean.max()])
    dfhi=dfa.groupby(['row','col']).apply(lambda x:x[x.HarvestIndex_mean==x.HarvestIndex_mean.max()])
    df33=dfa[(dfa.HB==3)&(dfa.Sow==3)]
    
    df=dfb.set_index(['row','col']).join(dfhi.set_index(['row','col']),rsuffix='_hi').join(df33.set_index(['row','col']),rsuffix='_33').reset_index()

    df['Delta_TB']=df.TotalBiomass_mean-df.TotalBiomass_mean_33
    df['Delta_HI']=df.HarvestIndex_mean_hi-df.HarvestIndex_mean_33
    df[['row','col','Delta_TB','Delta_HI']].to_csv(r'Q:\NZ\sowByHybridResults_2015\Delta.csv',index=False)
    print df.Delta_TB
    fig=plt.figure(figsize=(7.8,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)

    
    f1=fig.add_subplot(1,2,1)
    f1.text(0.8,0.03,'(a)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()

    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Delta_TB']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)/1000
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=2,vmax=10) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.02,0.35,0.01,0.6])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='max')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$\Delta_{total\ biomass} (t\ ha^{-1})$',fontsize=12)    
    
    f1=fig.add_subplot(1,2,2)
    f1.text(0.8,0.03,'(b)', fontsize=8,transform=f1.transAxes)
    df_zone,cliArr,m,xx,yy,baseName=GetBase()
    dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
    print dfzc.columns
    arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
    arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['Delta_HI']
    ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
#     imshow(ArrayMasked);show()

    CS=m.pcolormesh(xx,yy,ArrayMasked.T,alpha=1,cmap='RdYlGn_r',vmin=0,vmax=0.05) 
    m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)

    cax = fig.add_axes([0.52,0.35,0.01,0.6])#
    cbar=colorbar(CS,cax,orientation='vertical',extend='max')
    cbar.ax.tick_params(labelsize=7,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$\Delta _{HI}(t\ t^{-1})$',fontsize=12)    

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\mapBestCombinationAndSTD_HI_Delta.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def MapCVBiomassHIAcrossSowHy():
    df=read_csv(r"Q:\NZ\sowByHybridResults_2015\CVAcrossSowHy.csv")
    fig=plt.figure(figsize=(7.8,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)

    for va,tag,item,VN in zip(['TB','HI'],['(a)','(b)'],[1,2],['Total biomass','HI']):
        f1=fig.add_subplot(1,2,item)
        f1.text(0.8,0.03,tag, fontsize=8,transform=f1.transAxes)
        f1.text(0.03,0.8,VN, fontsize=10,transform=f1.transAxes)
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
        arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
        arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['%s_cv'%va]
        ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
        CS=m.pcolormesh(xx,yy,ArrayMasked.T,cmap='RdYlGn_r',vmin=0,vmax=50) 
        m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
    cax = fig.add_axes([0.1,0.,0.8,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal',extend='max')
    cbar.ax.tick_params(labelsize=8,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$CV (\%)$',fontsize=12)    

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\CV_B_HI.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
def MapCVBiomassHIAcrossSowHy():
    df=read_csv(r"Q:\NZ\sowByHybridResults_2015\CVAcrossSowHy.csv")
    fig=plt.figure(figsize=(7.8,6))
    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1,wspace=0.0, hspace=0.0)

    for va,tag,item,VN in zip(['TB','HI'],['(a)','(b)'],[1,2],['Total biomass','HI']):
        f1=fig.add_subplot(1,2,item)
        f1.text(0.8,0.03,tag, fontsize=8,transform=f1.transAxes)
        f1.text(0.03,0.8,VN, fontsize=10,transform=f1.transAxes)
        df_zone,cliArr,m,xx,yy,baseName=GetBase()
        dfzc=df_zone.set_index(['row','col']).join(df.set_index(['row','col']),rsuffix='_r').reset_index().dropna()
        arr=np.zeros(cliArr.shape);arr[:,:]=-999.0
        arr[dfzc.row.astype(int),dfzc.col.astype(int)]=dfzc['%s_cv'%va]
        ArrayMasked = np.ma.masked_where(arr == -999.0, arr)
        CS=m.pcolormesh(xx,yy,ArrayMasked.T,cmap='RdYlGn_r',vmin=0,vmax=50) 
        m.readshapefile(baseName,'base',linewidth=0.1,drawbounds=True)
    cax = fig.add_axes([0.1,0.,0.8,0.02])#
    cbar=colorbar(CS,cax,orientation='horizontal',extend='max')
    cbar.ax.tick_params(labelsize=8,size=6)
    cbar.locator = MaxNLocator( nbins = 4)
    cbar.set_label(r'$CV (\%)$',fontsize=12)    

    fig.savefig(r"V:\NZ\sowByHybridResults_2015\CV_B_HI.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
if __name__ == '__main__':
    SD_MeanacrossSH()
    MapCVBiomassHIAcrossSowHy()


        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
