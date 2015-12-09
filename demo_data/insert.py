import sqlite3 as sql

data = open('demo_data_{0}.csv'.format('snow_depth'), 'r')

c = sql.connect('demo.db')
count = 0
for row in data:
	if count == 0:
		count = 1
	else:
		d = row.split(',')
		query = 'insert into snow_depth(Site_Key, date_time, WY, Year, Month, Day, Hour, Minute, zs ) values(?,?,?,?,?,?,?,?,?)'
		c.execute(query, (d[0], d[1], d[2], d[3], d[4], d[5],d[6],d[7],d[8]))
		c.commit()
		print d
