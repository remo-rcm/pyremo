

import pandas as pd

remo_table = 'domains.csv'

cordex = 'cordex.csv'


remo = pd.read_csv(remo_table)
print(remo)
cordex = pd.read_csv(cordex, index_col='short_name')
print(cordex)


remo.rename(columns={'Short Name':'short_name', 'West (Cell Center)':'ll_lon', 'South (Cell Center)':'ll_lat', 
    'Nr. of Long. Cells (X-Size)':'nlon','Nr. of Lat. Cells (Y-Size)':'nlat', 
    'Long. North Pole.':'pollon', 'Lat. North Pole':'pollat', 'Long. Res.':'dlon', 'Lat. Res.':'dlat'}, inplace=True)

print(remo)
#print(remo['Short Name'])
remo.to_csv('remo-domains.csv')


