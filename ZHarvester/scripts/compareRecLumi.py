#This script compares triggerpath dependent deadtime and general deadtime, and print out different points and write into a csv file
import csv
import os,sys

f1 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_trig.csv', 'r')
f2 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_notrig.csv', 'r')
f3 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_diff.csv', 'w')

c1 = csv.reader(f1, delimiter=",")
c2 = csv.reader(f2, delimiter=",")
c3 = csv.writer(f3, delimiter=",")

triggerlist = list(c1)
diff_row = []

for notrigger_row in c2:
    if notrigger_row[0] == '#Data tag : v1 ':
        continue;
    if notrigger_row[0] == '#run:fill':
        continue; 
    found = False
    for trigger_row in triggerlist:
        if trigger_row[0] == notrigger_row[0] and trigger_row[1] == notrigger_row[1]:
            if float(notrigger_row[6]) != float(trigger_row[5]):
                print "luminosity difference spotted:"
                print trigger_row[0],trigger_row[1]
                print float(notrigger_row[6]) - float(trigger_row[5])

                diff_row = trigger_row
                diff_row.append(",")
                diff_row.append(float(notrigger_row[6]) - float(trigger_row[5]))
                c3.writerow(diff_row)
            found = True
            break
    if found == False:
        print "not in trigger csv"
        print notrigger_row
 
        diff_row = notrigger_row
        c3.writerow(diff_row)


f1.close()
f2.close()
f3.close()
