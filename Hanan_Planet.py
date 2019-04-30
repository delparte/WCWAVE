# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:33:11 2018

@author: abouhana
"""

import os
import arcpy
import re

from arcpy import env
from arcpy.sa import *

#initial directory with all the subfolders and raw satellite imagery
mydir = r"C:\Students\Hanan\Thesis_work\Data\4Band\Lebanon"
#directory to store the individual corrected bands and indices 
rasters = r"C:\Students\Hanan\Thesis_work\Data\Processing\Lebanon\LB_Planet\LB_Rasters"

env.workspace =  mydir
env.overwriteOutput = True

#lists for images & corresponding xml files 
P_4Band = []
P_Metadata = []

#iterate through folders and files to find the tif & xml
for root,dirs,files in os.walk(mydir):
    for name in files:
        if name.endswith("MS.tif"):
            #taking the full path with the extension to be able to access it later
            mypath = root+"\\"+name
            P_4Band.append(mypath)
            print(name)
            
        elif name.endswith("metadata.xml"):
            mypath = root+"\\"+name
            P_Metadata.append(mypath)
            print(name)            

#progress report purposes
print(P_4Band)
print("---------------------------")

print(P_Metadata)
print("---------------------------")

#list element index counter
j=0
temp1 = 1

#reading through xml files 
for entry in P_Metadata:

    #read xml file & extract correction coefficients
    print (entry)
    doc = open(entry,"r")
    i = 0
    while i<208:
        var = doc.readline()
        if i==183:        
            temp = var            
            a_2=re.sub("<ps:reflectanceCoefficient>", "", temp)
            a_2=re.sub("</ps:reflectanceCoefficient>","",a_2)
            print (i)
            print (a_2)
            
        elif i==192:
            temp = var
            a_3=re.sub("<ps:reflectanceCoefficient>", "", temp)
            a_3=re.sub("</ps:reflectanceCoefficient>","",a_3)
            print (i)
            print (a_3)
            
        elif i==201:
            temp = var
            a_4=re.sub("<ps:reflectanceCoefficient>", "", temp)
            a_4=re.sub("</ps:reflectanceCoefficient>","",a_4)
            print (i)
            print (a_4)          
       
        i+=1
    
    print("---------------------------")
            
    print(j)
    
    print("---------------------------")
    print("---------------------------")

    #to see which file its at
    print(entry)

    #rename date using format _DDMM      
    rename = entry[-43:-39]       
    
    a = rename[2:]
    b = rename[:2]

    rename = a+b

    print(rename)

    #pulll out the image file that corresponds to the xml file      
    x = P_4Band[j]

    #to verify that it has the matching set 
    print(x)

    arcpy.env.workspace = entry
    
    #first run 
    if j==0:
        
        #do atmoshperic correction & indices' calculations here
        
        Green = Raster(os.path.join(x,"band_2"))
        Green.save("B02_"+rename+"_old.tif")
        Green_correct = Green*float(a_2)
        env.workspace = rasters
        Green_correct.save("B02_"+rename+".tif")  
        print("Green band corrected & saved")
        env.workspace = entry
            
            
        Red = Raster(os.path.join(x,"band_3"))
        Red.save("B03_"+rename+"_old.tif")
        Red_correct = Red*float(a_3)
        env.workspace = rasters
        Red_correct.save("B03_"+rename+".tif")
        print("Red band corrected & saved")
        env.workspace = entry
            
            
        NIR = Raster(os.path.join(x,"band_4"))
        NIR.save("B04_"+rename+"_old.tif")    
        NIR_correct = NIR*float(a_4)
        env.workspace = rasters
        NIR_correct.save("B04_"+rename+".tif")      
        print("NIR band corrected & saved")    
        print("---------------------------")
        
        
        Green = Green_correct
        Red = Red_correct
        NIR = NIR_correct
        
        
        NDVI = ((NIR-Red)/(NIR+Red))
        NDVI.save("NDVI_"+rename+".tif")
        print("NDVI calculated")


        GNDVI = ((NIR-Green)/(NIR+Green))
        GNDVI.save("GNDVI_"+rename+".tif") 
        print("GNDVI calculated")
        

        SAVI = (((NIR-Red)/(NIR+Red+0.5))*1.5)
        SAVI.save("SAVI_"+rename+".tif") 
        print("SAVI calculated")
        
        
        MSAVI2 = ((2*NIR)+1-SquareRoot(Square((2*NIR)+1)-8*(NIR-Red)))/2
        MSAVI2.save("MSAVI_"+rename+".tif")    
        print("MSAVI2 calculated")
        
        
        NDWI = ((Green-NIR)/(Green+NIR))
        NDWI.save("NDWI_"+rename+".tif")
        print("NDWI calculated")
        
        print("---------------------------")
        print("---------------------------")

        env.workspace = entry
   
    if j!=0:
                     
        #read the date of the entry & the one before it to see if same date but different scene
        c1 = P_Metadata[j-1]
        c2 = P_Metadata[j]

        rename1 = c1[-43:-39]
        rename2 = c2[-43:-39]

        a1 = rename1[2:]
        a2 = rename2[2:]

        b1 = rename1[:2]
        b2 = rename2[:2]

        rename1 = a1+b1
        rename2 = a2+b2
        
        if rename2 != rename1:
            
            Green = Raster(os.path.join(x,"band_2"))
            Green.save("B02_"+rename+"_old.tif")
            Green_correct = Green*float(a_2)
            env.workspace = rasters
            Green_correct.save("B02_"+rename+".tif")  
            print("Green band corrected & saved")
            env.workspace = entry
                
                
            Red = Raster(os.path.join(x,"band_3"))
            Red.save("B03_"+rename+"_old.tif")
            Red_correct = Red*float(a_3)
            env.workspace = rasters
            Red_correct.save("B03_"+rename+".tif")
            print("Red band corrected & saved")
            env.workspace = entry
                
                
            NIR = Raster(os.path.join(x,"band_4"))
            NIR.save("B04_"+rename+"_old.tif")    
            NIR_correct = NIR*float(a_4)
            env.workspace = rasters
            NIR_correct.save("B04_"+rename+".tif")      
            print("NIR band corrected & saved")    
            print("---------------------------")
            
            
            Green = Green_correct
            Red = Red_correct
            NIR = NIR_correct
            
            
            NDVI = ((NIR-Red)/(NIR+Red))
            NDVI.save("NDVI_"+rename+".tif")
            print("NDVI calculated")


            GNDVI = ((NIR-Green)/(NIR+Green))
            GNDVI.save("GNDVI_"+rename+".tif") 
            print("GNDVI calculated")
            

            SAVI = (((NIR-Red)/(NIR+Red+0.5))*1.5)
            SAVI.save("SAVI_"+rename+".tif") 
            print("SAVI calculated")
            
            
            MSAVI2 = ((2*NIR)+1-SquareRoot(Square((2*NIR)+1)-8*(NIR-Red)))/2
            MSAVI2.save("MSAVI_"+rename+".tif")    
            print("MSAVI2 calculated")
            
            
            NDWI = ((Green-NIR)/(Green+NIR))
            NDWI.save("NDWI_"+rename+".tif")
            print("NDWI calculated")
            
            print("---------------------------")
            print("---------------------------")

            temp1 = 1
            env.workspace = entry

        elif rename2 == rename1:

            #if same date but different scene use different naming system using additional parameter counter
                       
            Green = Raster(os.path.join(x,"band_2"))
            Green.save("B02_"+rename+"_old.tif")
            Green_correct = Green*float(a_2)
            env.workspace = rasters
            Green_correct.save("B02_"+rename+"_"+str(temp1)+".tif")  
            print("Green band corrected & saved")
            env.workspace = entry
                
                
            Red = Raster(os.path.join(x,"band_3"))
            Red.save("B03_"+rename+"_old.tif")
            Red_correct = Red*float(a_3)
            env.workspace = rasters
            Red_correct.save("B03_"+rename+"_"+str(temp1)+".tif")
            print("Red band corrected & saved")
            env.workspace = entry
                
                
            NIR = Raster(os.path.join(x,"band_4"))
            NIR.save("B04_"+rename+"_old.tif")    
            NIR_correct = NIR*float(a_4)
            env.workspace = rasters
            NIR_correct.save("B04_"+rename+"_"+str(temp1)+".tif")      
            print("NIR band corrected & saved")    
            print("---------------------------")
            
            
            Green = Green_correct
            Red = Red_correct
            NIR = NIR_correct
            
            
            NDVI = ((NIR-Red)/(NIR+Red))
            NDVI.save("NDVI_"+rename+"_"+str(temp1)+".tif")
            print("NDVI calculated")


            GNDVI = ((NIR-Green)/(NIR+Green))
            GNDVI.save("GNDVI_"+rename+"_"+str(temp1)+".tif") 
            print("GNDVI calculated")
            

            SAVI = (((NIR-Red)/(NIR+Red+0.5))*1.5)
            SAVI.save("SAVI_"+rename+"_"+str(temp1)+".tif") 
            print("SAVI calculated")
            
            
            MSAVI2 = ((2*NIR)+1-SquareRoot(Square((2*NIR)+1)-8*(NIR-Red)))/2
            MSAVI2.save("MSAVI_"+rename+"_"+str(temp1)+".tif")    
            print("MSAVI2 calculated")
            
            
            NDWI = ((Green-NIR)/(Green+NIR))
            NDWI.save("NDWI_"+rename+"_"+str(temp1)+".tif")
            print("NDWI calculated")
            
            print("---------------------------")
            print("---------------------------")

            temp1+=1
            env.workspace = entry
    
    j+=1
    
