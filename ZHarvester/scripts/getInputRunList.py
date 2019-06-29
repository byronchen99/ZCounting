#This script scan input byls csv file and generate a run list to be put as parameter in batch setup 
import csv
import os,sys

runlist=[]
previousRun=299918

with open("../2017LumiByLS_hfet_trig_PU.csv", "rb") as input_file:
    reader = csv.reader(input_file, delimiter=",")
    for row in reader:
        if row[0] == '#Data tag : v1 ':
            continue
        if row[0] == '#run:fill':
            continue 
	if len(row) != 8:
            continue;
        currentRun = int(row[0].split(':')[0])
        if currentRun > previousRun:
            runlist.append(currentRun)
            previousRun = currentRun

with open("runlist.csv", "w+") as to_file:
    writer = csv.writer(to_file, delimiter=",")
    writer.writerow(runlist)
