
# ------------------------------------------------------------------------------
## load input by lumisection CSV file, convert it and return the data 
def load_input_csv(byLS_data):

    import pandas as pd
    
    byLS_file = open(str(byLS_data))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_data, sep=',', low_memory=False,
        skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))
        
    print(" formatting csv file...")    # formatting the csv
    byLS_data[['run', 'fill']] = byLS_data['#run:fill'].str.split(':', expand=True).apply(pd.to_numeric)
    byLS_data['ls'] = byLS_data['ls'].str.split(':', expand=True)[0].apply(pd.to_numeric)
    byLS_data = byLS_data.drop(['#run:fill', 'hltpath', 'source'], axis=1)

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
def load_histogram(
    name,      # name of the histogram to load
    fileName,  # file where the histogram is stored
    lumisections=[0,], # list of lumisections
    run=0, 
    prefix="", suffix="", 
    MassBin=50, MassMin=66, MassMax=116, 
    pileup=False
):
    import ROOT 
    
    file_ = ROOT.TFile(fileName)
    
    h_X_ls = file_.Get("{0}{1}".format(prefix, name)).Clone("{0}{1}{2}".format(prefix, name, suffix))
    h_X_ls.SetDirectory(0)
    h_X = h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}".format(name, run, suffix), lumisections[0], lumisections[0], "e")
    
    for ls in lumisections[1:]:
        h_X.Add(h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}_{3}".format(name, run, ls, suffix), ls, ls, "e"))
    
    if pileup:
        h_X.SetDirectory(0)
        return h_X
        
    # create new histogram in correct bin range
    hNew = ROOT.TH1D("h_mass_{0}_{1}_{2}_{3}".format(
        name, run, lumisections[0], suffix), "",MassBin, MassMin, MassMax)
    
    for ibin in range(0, h_X.GetNbinsX()+1):
        binCenter = h_X.GetBinCenter(ibin)
        if binCenter < MassMin:
            continue
        elif binCenter > MassMax:
            break
        else:
            content = h_X.GetBinContent(ibin)
            newBin = hNew.FindBin(binCenter)
            hNew.SetBinContent(newBin, hNew.GetBinContent(newBin) + content)
        
    hNew.SetDirectory(0)

    return hNew

# ------------------------------------------------------------------------------
## Collect "by LS" and "by measurement" csv files and write big csv file
def writeSummaryCSV(outCSVDir, writeByLS=True):
    import logging as log
    import pandas as pd
    import glob
    
    print(" ===Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True, sort=False)

    with open(outCSVDir + '/Mergedcsvfile_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    if writeByLS:
        print(" ===Writing overall CSV file per LS")
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

        with open(outCSVDir + '/Mergedcsvfile_perLS.csv', 'w') as file:
            df_merged.to_csv(file, index=False)

# ------------------------------------------------------------------------------
def chart_to_js(_chart, _name, _data="data.csv"):
    import json, os
    
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
