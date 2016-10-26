import sqlite3 as sql
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
bad_station = []

conn = sql.connect("{0}/rc_data.db".format(dir_path))
conn.execute('''CREATE TABLE IF NOT EXISTS 'weather' ('Site_Key' TEXT NOT NULL,
'date_time' TEXT NOT NULL,
'wy' INTEGER NOT NULL DEFAULT '1900',
'wd' INTEGER NOT NULL DEFAULT '1',
'year' INTEGER NOT NULL DEFAULT '1900',
'month' INTEGER NOT NULL DEFAULT '1',
'day' INTEGER NOT NULL DEFAULT '1',
'hour' INTEGER NOT NULL DEFAULT '1',
'minute' INTEGER NOT NULL DEFAULT '0',
'tmp3' REAL NOT NULL DEFAULT '-999.0',
'hum3' REAL NOT NULL DEFAULT '-999.0',
'vap3' REAL NOT NULL DEFAULT '-999.0',
'dpt3' REAL NOT NULL DEFAULT '-999.0',
'sol' REAL NOT NULL DEFAULT '-999.0',
'wnd3sa' REAL NOT NULL DEFAULT '-999.0',
'wnd3d' REAL NOT NULL DEFAULT '-999.0',
PRIMARY KEY('site_key','date_time'))''')

conn.execute('''CREATE TABLE IF NOT EXISTS 'soil_temperature'('Site_Key' TEXT NOT NULL,
'date_time' TEXT NOT NULL,
'stm002_5' REAL NOT NULL DEFAULT '-999.0',
'stm005' REAL NOT NULL DEFAULT '-999.0',
'stm010' REAL NOT NULL DEFAULT '-999.0',
'stm015' REAL NOT NULL DEFAULT '-999.0',
'stm020' REAL NOT NULL DEFAULT '-999.0',
'stm030' REAL NOT NULL DEFAULT '-999.0',
'stm040' REAL NOT NULL DEFAULT '-999.0',
'stm050' REAL NOT NULL DEFAULT '-999.0',
'stm055' REAL NOT NULL DEFAULT '-999.0',
'stm060' REAL NOT NULL DEFAULT '-999.0',
'stm070' REAL NOT NULL DEFAULT '-999.0',
'stm090' REAL NOT NULL DEFAULT '-999.0',
'stm120' REAL NOT NULL DEFAULT '-999.0',
'stm180'  REAL NOT NULL DEFAULT '-999.0',
PRIMARY KEY('site_key','date_time'))''')

conn.execute('''CREATE TABLE IF NOT EXISTS 'snow_depth' ('Site_Key' TEXT NOT NULL,
'date_time' TEXT NOT NULL,
'wy' INTEGER NOT NULL DEFAULT '1900',
'wd' INTEGER NULL DEFAULT '1',
'year' INTEGER NOT NULL DEFAULT '1900',
'month' INTEGER NOT NULL DEFAULT '1',
'day' INTEGER NOT NULL DEFAULT '1',
'hour' INTEGER NOT NULL DEFAULT '1',
'minute' INTEGER NOT NULL DEFAULT '0',
'snowdepth' REAL NOT NULL DEFAULT '-999.0',
PRIMARY KEY('site_key','date_time'))''')

conn.execute('''CREATE TABLE IF NOT EXISTS 'precipitation' ('Site_Key' TEXT NOT NULL,
'date_time' TEXT NOT NULL,
'ppts' REAL NOT NULL DEFAULT '-999.0',
'pptu' REAL NOT NULL DEFAULT '-999.0',
'ppta' REAL NOT NULL DEFAULT '-999.0',
PRIMARY KEY('site_key','date_time'))''')

def main():
    path_data_folder_climate = "C:/ReynoldsCreek/Watersheds/RC/ftp_data/climate/hourly-model-reference"
    print("weather")
    insertData("weather", path_data_folder_climate)
    path_data_folder_snow_depth = "C:/ReynoldsCreek/Watersheds/RC/ftp_data/snow-depth"
    print("snow depth")
    insertData("snow_depth", path_data_folder_snow_depth)
    path_data_folder_precipitation = "C:/ReynoldsCreek/Watersheds/RC/ftp_data/precip/hourly"
    print("Precipitation")
    insertData("precipitation", path_data_folder_precipitation)
    path_data_folder_soil_temp = "C:/ReynoldsCreek/Watersheds/RC/ftp_data/soil_temp/hourly"
    print("soil temperature")
    insertData("soil_temperature", path_data_folder_soil_temp)

def insertData(table, path_dir):
    ls_files = os.listdir(path_dir)

    for file in ls_files:
        if file[-4:] == ".dat":
            file_location = '{0}/{1}'.format(path_dir,file)
            with open(file_location, 'rb') as f:
                ls_header = []
                for line in f:
                    line = line.rstrip()
                    ls_data = []
                    if line[0] == '#':
                        line = line.split()
                        if "Site:" in line or "Site" and "Key" in line:
                            ##print(line)
                            site_key = line[-1]
                    else:
                        line = line.split(',')
                    if line[0]== 'datetime':
                        ls_header = line
                        ls_header[0] = 'date_time'
                        ls_header.insert(0, 'Site_Key')
                    elif line[0][0] != "#" and len(line) > 1:
                        ls_data = line
                        ls_data.insert(0, site_key)
                        placeholders = ['?'] * len(ls_data)
                        placeholders = tuple(placeholders)
                        ph = "("
                        for i in placeholders:
                            ph = ph + i +','
                        ph = ph[:-1] + ')'
                        query = '''INSERT INTO {0}{1} VALUES{2}'''.format(table,tuple(ls_header),ph)
                        try:
                            conn.execute(query,ls_data)
                        except sql.IntegrityError:
                            print(ls_data)
                            print("Not added. Likely a duplicate")
                        except sql.OperationalError as e:
                            if ls_data[0] not in bad_station:
                                print("Error: {0}, {1},\n {2}".format(ls_data[0],table, e))
                                bad_station.append(ls_data[0])
                conn.commit()

if __name__ == "__main__":
    main()
    conn.close()


