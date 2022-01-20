
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
