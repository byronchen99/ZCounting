Z Counting
==========
The project is divided into 4 packages.

TnPPairTreeProducer
----------------------
The `TnPPairTreeProducer` packege is used to produce tag-and-probe trees from cmssw files in AOD format and is mainly run over data.
The production is done via `crab3`. Scripts are available in `TnPPairTreeProducer/production/`.
It replaces the part of the `ZCounting DQM module` in cmssw (https://github.com/cms-sw/cmssw/tree/master/DQMOffline/Lumi) for offline analyses.
Instead of producing histograms (like the `ZCounting DQM module` does), the output is stored in root trees which makes it more flexible.
It contains scripts in the `TnPPairTreeProducer/python/Analyze` folder:
* with the `histProducer.py` the tag-and-probe trees can be filled into histograms which have the same format as the ones from the `ZCounting DQM module`.
* with the `produceBkgTemplate.py` a kernel density estimation is performed on the tag-and-probe trees and stored as a histogram.
this can be used on a tag-and-probe same sign selection to generate data driven QCD background templates.

ZCountAnalyze
----------------
The `ZCountAnalyze` package is used to produce muon pair trees from cmssw files in MiniAOD format
and is usually run over MC only.
These muon pair trees contain truth level information and can be used to analyze the method based on simulation.
It contains scripts in the `ZCountAnalyze/python/Analyze` folder:
* to investigate certain quantities of the method based on truth level information.
* to produce a `.json` file with the correction of the tag and probe muon efficiency based on truth level.
* to produce MC templates which can then be used in the fitting procedure of the tag and probe measurement (with the `ZHarvester`).
* to investigate tag-and-probe trees (w/o truth level) and compare it with MC.

ZHarvester
-------------
This is the "core" package of the project. 
It can be used with the histograms from the `ZCounting DQM module` or the ones produced in the `TnPPairTreeProducer` package.
* Here, the measurement of the tag-and-probe muon efficiencies and calculation of e.g. the Z Rate can be done.
    * based on RooFit
    * HTCondor job submission is supported
* Several scipts are available to provide additional inputs and produce plots

ZUtils
-------------
Contains code that is used by multiple packages in this project in common.
