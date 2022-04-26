Z Counting
==========
The project is divided into 4 packages.

ZCountAnalyze
----------------
The `ZCountAnalyze` package is used to produce ntuples from cmssw files in MiniAOD or AOD format
and is run over MC or data.
These ntuples contain truth level information and can be used to analyze the method based on simulation, or re-analyze the data in a more flexible way.
It contains scripts in the `ZCountAnalyze/python/Analyze` folder:
* to investigate certain quantities of the method based on truth level information.
* to produce a `.json` file with the correction of the tag and probe muon efficiency based on truth level.
* to produce MC templates which can then be used in the fitting procedure of the tag and probe measurement (with the `ZHarvester`).
* to investigate ntuples (w/o truth level) and compare it with MC.

ZHarvester
-------------
This is the "core" package of the project. 
It can be used with the histograms from the `ZCounting DQM module` or the ones produced from the ntuples obtained from the `ZCountAnalyze` package.
* Here, the measurement of the tag-and-probe muon efficiencies and calculation of e.g. the Z Rate can be done.
    * based on RooFit
    * HTCondor job submission is supported
* Several scipts are available to provide additional inputs and produce plots

ZUtils
-------------
Contains code that is used by multiple packages in this project in common.

BLUE
-------------
Combinations of Z counts and reference luminosity with the BLUE tool
