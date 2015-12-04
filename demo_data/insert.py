import sqlite3 as sql

data = open('demo_data_{0}.csv'.format('soil_temp'), 'r')

c = sql.connect('demo.db')
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