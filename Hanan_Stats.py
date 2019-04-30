import geopandas as gpd
import pandas as pd
import glob
import os
from rasterstats import zonal_stats

rsterDir=r'C:\Students\Hanan\Thesis_work\Data\Processing\Lebanon\LB_Planet\LB_Rasters'
mydir=r'C:\Students\Hanan\Thesis_work\Data\Processing\Lebanon\LB_Planet\LB_Rasters'
shpDir=r'C:\Students\Hanan\Thesis_work\Data\Processing\Lebanon\LB_Planet\shp'
procDir=r'C:\Students\Hanan\Thesis_work\Data\Processing\Idaho\ID_Planet\GridStats'

os.chdir(shpDir)
for shp in glob.glob('*.shp'):
    print shp
    shpNmeList=shp.split('.')
    shpNme=shpNmeList[0]
    rasterList=[]
    gdf=gpd.read_file(shp)
    
    gdf['FID']=None
    counterShp=0
    while counterShp<len(gdf):
        gdf.loc[counterShp,'FID']=shpNme+'_'+str(counterShp)
        counterShp+=1
    os.chdir(rsterDir)
##    for rast in glob.glob('*'):
##        rasterList.append(rast)
    statOpsList=['mean','sum']
    for statVar in statOpsList:
        print statVar
        cols=[]
        cols.append('FID')
        for rster in rasterList:    
            rster=rster.split('.')
            rster=rster[0]          
            cols.append(rster)
            
        df=pd.DataFrame(columns=cols)
        counter=0
        while counter<len(gdf):
            df.loc[counter,'FID']=gdf.loc[counter,'FID']
            counter+=1
        for rster in rasterList:
            os.chdir(rsterDir)
            stats=zonal_stats(gdf,rster,stats=statVar)
            rster=rster.split('.')
            rster=rster[0]
            print rster

            counter1=0
            while counter1<len(stats):
                stat=stats[counter1]
                fld=df.loc[counter1,'FID']
                counter2=0
                while counter2<len(df):
                    if df.loc[counter2,'FID']==fld:
                        df.loc[counter2,rster]=stat[statVar]
                    counter2+=1
                counter1+=1
        os.chdir(mydir)
        df.to_csv(shpNme+'_'+statVar+'.csv')    
        if statVar=='mean':
            for col in cols:
                if col!='FID':
                    df[col]=df[col].astype(float)
           
            gdfFinal=gdf.join(df.set_index('FID'),on='FID')
            print list(gdfFinal)
                      
            gdfFinal.to_file(shpNme+'_stats.shp',driver='ESRI Shapefile')
                        
        
    os.chdir(shpDir)
        
