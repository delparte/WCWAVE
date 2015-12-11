import sqlite3 as sql

param = 'soil_temp'

data = open('demo_data_{0}.csv'.format(param), 'r')

c = sql.connect('demo.db')
if param == 'snow_depth':
    count = 0
    for row in data:
        if count == 0:
                count = 1
        else:
            d = row.split(',')
            query = 'insert into snow_depth(Site_Key, date_time, WY, Year, Month, Day, Hour, Minute, zs ) values(?,?,?,?,?,?,?,?,?)'
            d[8] = d[8].rstrip('\n')
            c.execute(query, (d[0], d[1], d[2], d[3], d[4], d[5],d[6],d[7],d[8]))
            c.commit()
            print d

if param == 'weather':
    count = 0
    for row in data:
        if count == 0:
            count = 1
        else:
            d = row.split(',')
            query = 'insert into weather(Site_Key, date_time, wy, year, month, day, hour, minute, ta, rh, ea, td, si, ws, wd ) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
            d[14] = d[14].rstrip('\n')
            c.execute(query, (d[0], d[1], d[2], d[3], d[4], d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14]))
            c.commit()
            print d

if param == 'soil_temp':
    count = 0
    for row in data:
        if count == 0:
            count = 1
        else:
            d = row.split(',')
            query = 'insert into soil_temperature(Site_Key, date_time, stm005) values(?,?,?)'
            c.execute(query, (d[0], d[1], d[3]))
            c.commit()
            print d

if param == 'precipt':
    count = 0
    for row in data:
        if count == 0:
            count = 1
        else:
            d = row.split(',')
            query = 'insert into precipitation(Site_Key, date_time, ppts, pptu, ppta) values(?,?,?,?,?)'
            d[4] = d[4].rstrip('\n')
            c.execute(query, (d[0], d[1], d[2], d[3], d[4]))
            c.commit()
            print d
