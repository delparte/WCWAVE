'''
Copyright (c) 2015 Wesley Joel Johansen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os
import sys
import traceback
from numpy import *
import mysql.connector
from mysql.connector import errorcode
import datetime
print("Start: " + str(datetime.datetime.now()))

#Get path to .dat files
pathDataFolder = r'C:\AA_Thesis_Project\ZZ_MySQL_Work\FTP_Data\snow_depth'
#List files in folder
lsFiles = os.listdir(pathDataFolder)

#connect to MySQL database
try:
  cnx = mysql.connector.connect(user='root', password='',
                                host='localhost',
                                database='rc_data')
except mysql.connector.Error as err:

    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        print("Something is wrong with your user name or password")
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        print("Database does not exist")
    else:
        print(err)
else:
    print("Connection successful")

    #Iterate over .dat files
    for file in lsFiles:
        #Check if file ends with ".dat"
        if not file[-4:] == ".dat":
            pass
        else:
            lsHeader = []
            boolHeader = True
            boolFirstLine = True
            #open file and read one line at a time
            fullPath = os.path.join(pathDataFolder, file)
            for line in open(fullPath,'r'):
                sLine = line.rstrip()
                if sLine[0] == "#" and boolFirstLine == True:
                        sSiteKey = sLine.split(" ")[-1]
                        boolFirstLine = False
                elif sLine[0] != "#" and boolHeader == True:
                    sHeader = sLine
                    #print(sHeader)
                    lsHeader = sHeader.split(",")
                    lsHeader[0] = "date_time"
                    lsHeader.insert(0,'Site_Key')
                    boolHeader = False
                elif sLine[0] != "#":
                    try:
                        sLine = sSiteKey + "," + sLine
                        lsLine = sLine.split(",")
                        #Format date_time
                        try:
                            dtTemp = datetime.datetime.strptime(lsLine[1], "%Y-%m-%d %H:%M")
                            lsLine[1] = dtTemp.strftime("%Y-%m-%d %H:%M:00")
                        except:
                            pass
                        cur = cnx.cursor()
                        query = 'INSERT INTO snow_depth ({0}) VALUES ({1})'
                        sTemp = "%s,"
                        sTemp2 = sTemp * len(lsLine)
                        sTemp3 = sTemp2[:-1]
                        query = query.format(','.join(lsHeader), sTemp3)
                        if len(sLine) > 10:
                            cur.execute(query, lsLine)

                    except:
                        print("### Did not execute " + sLine)
        cnx.commit()
    cnx.close()
print("Finish: " + str(datetime.datetime.now()))
