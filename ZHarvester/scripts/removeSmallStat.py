#The script does a few things:
#1) remove small stat. points and write into a new csv file
#2) save luminosity values for all points (both passing and failing the lulminosity cut) and sort them and write into a csv files in case to do further checks
#3) sum over stat. amoung all passing points and all failing points to know how much we are going to remove

import csv
import operator

lumi_pass = [] #Buffer list 
lumi_all  = [] #Buffer list 

with open("ZCounting2017Total.csv", "rb") as input_file:
    reader = csv.reader(input_file, delimiter=",")
    for row in reader:
        lumi_all.append(row[5])
        if float(row[5]) > 6.0: #Here is the luminosity cut to remove small stat. measurements
            lumi_pass.append(row)

with open("ZCounting2017Total_RemoveSmallStat.csv", "w+") as to_file:
    writer = csv.writer(to_file, delimiter=",")
    for new_row in lumi_pass:
       writer.writerow(new_row)

lumi_all.sort(key=float)

totPass = 0
totFail = 0
tot = 0

with open("checkstat.csv", "w+") as stat_file:
    writer2 = csv.writer(stat_file)
    for lumi in lumi_all:
        writer2.writerow([lumi])

        tot += float(lumi)
        if float(lumi) < 6.0: #Here is the luminosity cut to remove small stat. measurements
            totFail += float(lumi)
        else:
            totPass += float(lumi)

print "Total RecLumi: "+str(tot)
print "RecLumi > LumiCut: "+str(totPass)
print "RecLumi < LumiCut: "+str(totFail)
