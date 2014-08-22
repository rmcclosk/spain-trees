#!/usr/bin/env python

"""
Parse Spain data
"""
import csv
import sys

reader = csv.DictReader(open('00_raw.csv'))

data = {}

for row in reader:
    day = int(row["Time point"].split()[-1])
    count = int(row["Count"])
    try:
        fpr = float(row['G2P_FPR%'])
    except ValueError:
        continue
    patid = row["Patient"]

    if patid not in data:
        data.update({patid: {}})
    
    if day not in data[patid]:
        data[patid].update({day: dict(zip(range(5), [0]*5))})
    
    for i, cutoff in enumerate([2, 3.5, 10, 20, 100]):
        if fpr < cutoff:
            data[patid][day][i] += count
            break

print('patid,day,FPR.lt.2,FPR.lt.3.5,FPR.lt.10,FPR.lt.20,FPR.lt.100')

for patid in data.keys():
    for day in data[patid].keys():
        bin_values = [str(x) for x in data[patid][day].values()]
        print('{},{},{}'.format(patid, day, ",".join(bin_values)))
