from python.logging import child_logger
log = child_logger(__name__)

from datetime import datetime
import pandas as pd
import numpy as np
import uncertainties as unc
import glob
import uproot 
import ROOT
import boost_histogram as bh
from hist import Hist
import json, os

import pdb

# ------------------------------------------------------------------------------
def load_input_csv(byLS_data):
    """
    load input by lumisection CSV file, convert it and return the data 
    
    Parameters
    ----------
    byLS_data : str
        path to the file by lumisection CSV file
    """
    
    byLS_file = open(str(byLS_data))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_data, sep=',', low_memory=False,
        skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))
        
    log.info("formatting csv file...")    # formatting the csv
    byLS_data[['run', 'fill']] = byLS_data['#run:fill'].str.split(':', expand=True).apply(pd.to_numeric)
    byLS_data['ls'] = byLS_data['ls'].str.split(':', expand=True)[0].apply(pd.to_numeric)   

    if 'delivered(/ub)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/ub)'] = byLS_data['delivered(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data['recorded(/ub)'] = byLS_data['recorded(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/ub)': 'delivered(/pb)', 'recorded(/ub)': 'recorded(/pb)'})
    elif 'delivered(/fb)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/fb)'] = byLS_data['delivered(/fb)'].apply(lambda x: x * 1000.)
        byLS_data['recorded(/fb)'] = byLS_data['recorded(/fb)'].apply(lambda x: x * 1000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/fb)': 'delivered(/pb)', 'recorded(/fb)': 'recorded(/pb)'})

    # if there are multiple entries of the same ls (for example from different triggers), 
    #   only keep the one with the highest luminosity.
    byLS_data = byLS_data.sort_values(['fill', 'run', 'ls', 'delivered(/pb)', 'recorded(/pb)'])
    byLS_data = byLS_data.drop_duplicates(['fill', 'run', 'ls'])
    
    return byLS_data

# ------------------------------------------------------------------------------
def to_DateTime(time, string_format = "yy/mm/dd"):

    # converts time string to python datetime
    time = time.split(" ")

    if string_format == "yy/mm/dd":     # format agreed upon ATLAS and CMS
        year, month, day = [int(x) for x in time[0].split("/")]
        year += 2000
        hour, min, sec   = [int(x) for x in time[1].split(":")]
    elif string_format == "mm/dd/yy":   # brilcalc format
        month, day, year = [int(x) for x in time[0].split("/")]
        year += 2000
        hour, min, sec   = [int(x) for x in time[1].split(":")]
    else:   # default format
        year, month, day = [int(x) for x in time[0].split("-")]
        hour, min, sec   = [int(float(x)) for x in time[1].split(":")]
    
    return datetime(year, month, day, hour, min, sec)

# ------------------------------------------------------------------------------
def getFileName(directory, run):
    """
    Get root file with histograms
    
    Parameters
    ----------
    directory : str
        directory where the root files are stored
    run : integer
        run number for the rootfile to load
    """

    # check if run was processed already
    eosFileList = glob.glob(directory + '/*' + str(run) + '*.root')

    # look one level deeper
    if len(eosFileList) == 0:
        eosFileList = glob.glob(f'{directory}/000{str(run)[:-2]}xx/*{run}*.root')
    if len(eosFileList) == 0:
        eosFileList = glob.glob(f'{directory}/Muon/000{str(run)[:-2]}xx/*{run}*.root')
    if len(eosFileList) == 0:
        eosFileList = glob.glob(f'{directory}/SingleMuon/000{str(run)[:-2]}xx/*{run}*.root')

    if not len(eosFileList) > 0:
        log.warning(f"The file does not (yet) exist for run: {run}")
        return None
        
    elif len(eosFileList) > 1:
        log.warning(f"Multiple files found for run: {run}")
        return None
    else:
        return eosFileList[0]

# ------------------------------------------------------------------------------
def load_histogram(
    name,
    fileName,
    lumisections=[0,],
    run=0, 
    prefix="", suffix="", 
    MassBin=50, MassMin=66, MassMax=116, 
    pileup=False
):
    """
    load 2D histograms, project the specified lumisections to the mass axis: 
    - if pileup=False with (mass, lumisections); rebinning possible
    - if pileup=True, (PU, lumisections)
    
    Parameters
    ----------
    name : str
        name of the histogram to load
    fileName : str
        file where the histogram is stored
    lumisections : list
        list of lumisections that are taken from the histogram
    run : integer
        run number 
    prefix/suffix : str
        prefix for naming         
    MassBin/MassMin/MassMax : int
        For rebinning, Number of bins / Lower bound / Upper bound 
    pileup : Boolean
        If the pileup histogram is to be returned
    """

    with uproot.open("{0}:{1}{2}".format(fileName, prefix, name)) as th_:

        h_ = th_.to_numpy()
        
        # x-axis is mass or primary vertices
        bins_x = h_[2]

        # select lumisections
        h_ = np.sum(h_[0][lumisections], axis=0)

        if sum(h_) <= 0:
            log.warning(f"No entries found in histogram {name}")

        if pileup:
            return h_
        
        bins_new = np.linspace(MassMin, MassMax, MassBin+1)

        # no rebinning if bins are already the same
        if bins_new == bins_x or (len(bins_new)==len(bins_x) and (bins_new == bins_x).all()):
            return h_
        
        binWidth = bins_x[1:]-bins_x[:-1]
        binWidthNew = (MassMax - MassMin)/MassBin

        if (binWidth[0] == binWidth).all() == False:
            log.error("variable bin sizes, no rebinning possible!")
            return None

        if binWidth[0] > binWidthNew:
            log.error("bin width is smaller than original binning, no rebinning possible!")
            return None

        # same bin widths or larger, select different range
        bin_lo = (bins_new[0] - bins_x[0])/binWidth[0]
        bin_hi = (bins_new[-1] - bins_x[-1])/binWidth[0]

        if not bin_lo.is_integer() or not bin_hi.is_integer() or bin_lo<0 or bin_hi>0:
            log.error("can not rebin to new range!")
            return None

        h_ = h_[int(bin_lo):int(bin_hi)]

        if binWidth[0] < binWidthNew:
            f = binWidthNew/binWidth[0]
            l = len(h_)/f

            if not f.is_integer() or not l.is_integer():
                log.error("can not rebin to new width!")
                return None           

            h_ = h_.reshape(int(l), int(f)).sum(axis=1)
        
        return h_

# ------------------------------------------------------------------------------
def np_to_hist(y, nBins=None, binLow=None, binHigh=None, histtype="TH1", name=""):
    if histtype == "np":
        return y

    if histtype == "boost":
        hist = bh.Histogram(bh.axis.Regular(bins=nBins, start=binLow, stop=binHigh))
        hist[:] = y
        return hist

    if histtype == "hist":
        hist = Hist.new.Regular(nBins, binLow, binHigh, name="x").Double()
        hist[:] = y
        return hist

    if histtype == "TH1":
        hist = ROOT.TH1D(name, name, nBins, binLow, binHigh)
        for i, v in enumerate(y):
            hist.SetBinContent(i, v)
        return hist


# ------------------------------------------------------------------------------
def get_ls_for_next_measurement(
    lumisections, luminosities=None, zcounts=None, 
    lumiPerMeasurement=20, lsPerMeasurement=100, 
    threshold = 0.02
):
    """
    generator that takes the set of lumisections that are process 
    and yields slizes of lists of these lumisections 
    that should be used in the next measurement. 
    It can be run in two modes: 
    - If the parameter `luminosities` is not specified, the number of lsPerMeasurement is taken as criteria
    - If the parameter `luminosities` is specified, the amount of lumiPerMeasurementis taken as criteria
    
    Parameters
    ----------
    lumisections : list
        The list of lumisections that are going to be processed
    luminosities : list, optional
        The list of luminosity in \pb for each lumisection. 
    zcounts : list, optional
        The list of z boson counts for each lumisection.     
    lumiPerMeasurement : float, optional
        The amount of luminosity in \pb required for a measurement
    lsPerMeasurement : int, optional
        The number of lumisections required for a measurement
    threshold : float, optional
        If the luminosity in one lumisection is above this value and the number z counts is zero, the ls is skipped 
    """
    
    if luminosities:
        # make measurement based on number of lumisections
        if len(lumisections) != len(luminosities):
            log.error("Same length of lumisections and luminosities is required!")
        
    while len(lumisections) > 0:
        
        # merge data to one measuement if remaining luminosity is too less for two measuements
        if luminosities:
            mergeMeasurements_ = sum(luminosities) < 1.5 * lumiPerMeasurement
        else:
            mergeMeasurements_ = len(lumisections) < 1.5 * lsPerMeasurement        
        
        recLumi_ = 0
        # produce list_good_ls_ with lumisections that are used for one measurement
        list_good_ls_ = []
        while len(lumisections) > 0:
            
            if luminosities and zcounts:
                # consider lumisections where we would expect to have at least any z count 
                #   (for lumi > 0.01 /pb we expect 0.01*500 = 5 Z bosons, the probability to find 0 is < 1%)
                #   (for lumi > 0.02 /pb we expect 0.02*500 = 10 Z bosons, the probability to find 0 is < 0.01%)
                # sort out lumisections without any Z candidate (maybe trigger was off)
                if luminosities[0] > threshold and zcounts[0] == 0:
                    log.warning("Zero Z boson candidates found {0}/pb while we would expect {1} -> skip lumi section {2}".format(luminosities[0], luminosities[0]*500, lumisections[0]))
                    del lumisections[0]
                    del luminosities[0]
                    del zcounts[0]
                    continue

            list_good_ls_.append(lumisections[0])
            del lumisections[0]
            
            if zcounts:
                del zcounts[0]
                
            if luminosities:
                recLumi_ += luminosities[0]
                del luminosities[0]

                
                # if we have collected enough luminosity
                if not mergeMeasurements_ and recLumi_ >= lumiPerMeasurement:
                    break
            else:
                # if we have collected enough luminosity sections
                if not mergeMeasurements_ and len(list_good_ls_) >= lsPerMeasurement:
                    break

        log.debug("Selected lumi section list {0}".format(list_good_ls_))

        yield list_good_ls_

# ------------------------------------------------------------------------------
def getCorrelation(hPV_data, correlationsFileName, which, region="I"):
    """
    calculate the correlation factor from a resource file by folding the PV histogram with the correlation histogram  

    Parameters
    ----------
    hPV_data : TH1
        1D histogram with number of primary vertices in data
    correlationsFileName : str
        path to the file with correlation factors as function of number of primary vertices
    which: str
        which correlation factor to be read, can be: "IO", "ID", "Acceptance"
    region: str
        eta region for the correlation factor, can be: "I", "BB", "BE", "EE"
    """
    ROOT.gROOT.SetBatch(True) # disable root prompts

    tfile = ROOT.TFile(correlationsFileName,"READ")

    # normalize pileup histogram
    hPV_data.Scale(1./hPV_data.Integral())
    avgPV = hPV_data.GetMean()
    hC = tfile.Get("C{0}_{1}".format(which, region))

    hC.Multiply(hPV_data)

    corr = hC.Integral()

    tfile.Close()
    log.info("Correlation coefficienct c{0}_{1} = {2}".format(which, region, corr))
    log.info("For average primary vertices of <pv> = {0}".format(avgPV))
    return corr

# ------------------------------------------------------------------------------
def writePerRunCSV(results, outCSVDir, run):
    ## Write per measurement csv file - one per run
    log.info("Writing per Run CSV file")
    results = pd.concat([pd.DataFrame([result]) for result in results], ignore_index=True, sort=False)

    with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
        results.to_csv(file, index=False)

# ------------------------------------------------------------------------------
def writeSummaryCSV(outCSVDir, outName="Mergedcsvfile", writeByLS=False, keys=None):
    """
    Collect "by LS" (and "by measurement") csv files and write one big csv file
    
    Parameters
    ----------
    outCSVDir : str
        Directory where the csv files of each run are located
    writeByLS : boolean, optional
        Boolean if a summary should be created for csv files that contain the Z rates for each lumisection
    """
    
    log.info("Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True, sort=False)
    
    if keys:
        df_merged = df_merged[keys]
    
    with open(outCSVDir + '/' + outName + '_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    if writeByLS:
        log.info("Writing overall CSV file per LS")
        rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
        csvList = []
        # add measurement label to the csv list
        for m in rateFileList:
            csv = pd.read_csv(m)
            measurement = int(m.split("_")[-1][:-4])
            fill = int(m.split("_")[-2].split("csvfile")[-1])
            csv["measurement"] = measurement
            csvList.append(csv)
        df_merged = pd.concat(csvList, ignore_index=True, sort=False)

        # df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

        with open(outCSVDir + '/' + outName + '_perLS.csv', 'w') as file:
            df_merged.to_csv(file, index=False)

# ------------------------------------------------------------------------------
def open_workspace(directory, filename, m):
    # loading root workspace
    file_result = f"{directory}/{filename}_{m}.root"
    
    if not os.path.isfile(file_result):
        log.warning("No result for `{0}`".format(file_result))
        return None, None
    
    f = ROOT.TFile(file_result,"READ")
    w = f.Get("workspace")
    return f, w

def open_workspace_yield(directory, filename, m):
    # reading fit result from root workspace

    f, w = open_workspace(directory, "workspace_yield_"+filename, m)
    
    Nsig = w.var("Nsig").getVal()
    chi2 = w.arg("chi2").getVal()
    
    f.Close()
    f.Delete()
    w.Delete()
    return Nsig, chi2

def unorm(value):
    if value > 0:
        return unc.ufloat(value, np.sqrt(value))
    else: 
        return unc.ufloat(value, value)

# ------------------------------------------------------------------------------
def getEra(run):        
    """
    return the era for a given run number
    
    Parameters
    ----------
    run : int
        Drun number
    """

    if run <= 271658:
        return "Run2016A"
    if run <= 275376:
        return "Run2016B"        
    if run <= 276283:
        return "Run2016C"
    if run <= 276811:
        return "Run2016D"
    if run <= 277420:
        return "Run2016E"
    if run <= 278808:
        return "Run2016F"        
    if run <= 280385:
        return "Run2016G"
    if run <= 284044:
        return "Run2016H"        
    if run <= 297019:
        return "Run2017A"
    if run <= 299329:
        return "Run2017B"        
    if run <= 302029:
        return "Run2017C"
    if run <= 303434:
        return "Run2017D"
    if run <= 304826:
        return "Run2017E"
    if run <= 306462:
        return "Run2017F"        
    if run <= 306826:
        return "Run2017G"
    if run <= 307082:
        return "Run2017H"           
    if run <= 316995:
        return "Run2018A"
    if run <= 319312:
        return "Run2018B"        
    if run <= 320393:
        return "Run2018C"
    if run <= 325273:
        return "Run2018D"     
    if run <= 355769:
        return "Run2022B"   
    if run <= 357482:
        return "Run2022C"
    if run <= 357900:
        return "Run2022D"

    return "Run2022" #default 

# ------------------------------------------------------------------------------
def chart_to_js(_chart, _name, _data="data.csv"):
    
    chart_json = _chart.to_dict()
    # delete in line data
    del chart_json["datasets"]
    # add external data resource
    chart_json["data"] = {"url": _data}
    with open(_name+'.js','a') as file:
        file.write("""(function(vegaEmbed) {
    var spec = """)
        json.dump(chart_json, file, indent=4)
        file.write(""";
    var embedOpt = {\"renderer\": \"svg\", \"mode\": \"vega-lite\"};
    function showError(el, error){
        el.innerHTML = ('<div class=\"error\" style=\"color:red;\">'
        + '<p>JavaScript Error: ' + error.message + '</p>'
        + \"<p>This usually means there's a typo in your chart specification. \"
        + \"See the javascript console for the full traceback.</p>\"
        + '</div>');
        throw error;
    }""")
        file.write("const el = document.getElementById('vis');\n")
        file.write("vegaEmbed(\"#vis\", spec, embedOpt).catch(error => showError(el, error));\n")
        file.write("})(vegaEmbed);\n")

def chart_to_html(_chart, _name, **kwargs):

    with open(_name+'.js','w') as file:
        file.write("""
<!DOCTYPE html>
<html>
<head>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega@5"></script>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega-lite@4.8.1"></script>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega-embed@6"></script>
</head>
<body>
  <div id="vis"></div>
  <script type="text/javascript">
    """)
    chart_to_js(_chart, _name, **kwargs)
    with open(_name+'.js','a') as file:
        file.write("""
  </script>
</body>
</html>
""")

    os.system("mv {0}.js {0}.html".format(_name))
