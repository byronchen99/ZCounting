#This script write a fill list into a new csv file, which can be used as the input for plotting marco
#plot_Zrate_Data.py -f "filllist"

import csv

fill_list = [] #Buffer list 
prevous_fill = 0

with open("ZCounting2017Total_RemoveSmallStat.csv", "rb") as input_file:
    reader = csv.reader(input_file, delimiter=",")
    for row in reader:
        current_fill = int(row[0])
        if current_fill > prevous_fill:
            fill_list.append(current_fill) 
            prevous_fill = current_fill
            #print current_fill
        else:
            continue

with open("filllist.csv", "w+") as to_file:
    writer = csv.writer(to_file, delimiter=",")
    writer.writerow(fill_list)
