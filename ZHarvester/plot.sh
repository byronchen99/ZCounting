
PROJECTDIR=$1
REFBYLS=$2
SIGMA=$3

python Plotting/plot_ZEfficiency_DataCMS.py -c $PROJECTDIR/csvFiles/Mergedcsvfile_perMeasurement.csv -s $PROJECTDIR/ZEfficiency

python Plotting/plot_Zrate_DataCMS.py -c $PROJECTDIR/csvFiles/Mergedcsvfile_perMeasurement.csv -s $PROJECTDIR/Zrate -r $REFBYLS

python Plotting/makeZvsL.py -c $PROJECTDIR/csvFiles/Mergedcsvfile_perMeasurement.csv -s $PROJECTDIR/PlotsZvsL --xsec $SIGMA
