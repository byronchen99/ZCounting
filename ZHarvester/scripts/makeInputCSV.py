#This script started from input byls csv files with triggerpath dependent deadtime, and insert average pile-up info from the other byls csv files
#One need to add by hand two head lines into the new csv file
#Data tag : v1 , Norm tag: None
#run:fill,ls,time,trigger,delivered(/pb),recorded(/pb),source,avgpu

import csv
import os,sys

#f1 = file('/eos/cms/store/group/comm_luminosity/ZCounting/2017LumiByLS_trig.csv', 'r')
#f2 = file('/eos/cms/store/group/comm_luminosity/ZCounting/2017LumiByLS_notrig_PU.csv', 'r')
#f3 = file('/eos/cms/store/group/comm_luminosity/ZCounting/2017LumiByLS_trig_PU.csv', 'w')
f1 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_trig.csv', 'r')
f2 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_notrig.csv', 'r')
f3 = file('/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_trig_PU.csv', 'w')

c1 = csv.reader(f1, delimiter=",")
c2 = csv.reader(f2, delimiter=",")
c3 = csv.writer(f3, delimiter=",")

notriggerlist = list(c2)

for trigger_row in c1:
    if trigger_row[0] == '#Data tag : v1 ':
        c3.writerow(trigger_row)
        continue;
    if trigger_row[0] == '#run:fill':
        results_row = trigger_row
        results_row.append(",avgpu")
        c3.writerow(trigger_row)
        continue; 
    found = False
    for notrigger_row in notriggerlist:
        results_row = trigger_row
        if notrigger_row[0] == trigger_row[0] and notrigger_row[1] == trigger_row[1]:
            results_row.append(notrigger_row[7])
            found = True
            break
    c3.writerow(results_row)

f1.close()
f2.close()
f3.close()
