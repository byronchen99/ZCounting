import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCountingAOD',
    pileupSummaryInfoCollection=cms.InputTag("addPileupInfo"),

    era = cms.string("None"),
    # Primary vertex - source
    edmPVName = cms.InputTag("offlinePrimaryVertices"),
    
    # Muons - source
    reco_muons = cms.InputTag("muons","","RECO"),

    roccorFile = cms.string(""),   

    # Standalone muon tracks - source
    reco_standalones = cms.InputTag("standAloneMuons","","RECO"),
    reco_standalonesUpdated = cms.InputTag("standAloneMuons","UpdatedAtVtx","RECO"),

    # Tracks - source
    reco_tracks = cms.InputTag("generalTracks","","RECO"),

    # Electrons - source
    reco_electrons = cms.InputTag("gedGsfElectrons","","RECO"),
    reco_superclusters = cms.InputTag("particleFlowEGamma","","RECO"),

    # Trigger - sources
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
    TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
    
    emulateTrigger = cms.untracked.bool(False),

    # Muon - triggers
    # 2^{x}
    muon_trigger_patterns=cms.vstring(
                                # x =
        "HLT_IsoMu24_v*",       # 0
        "HLT_IsoMu27_v*",       # 1
        "HLT_IsoTkMu24_v*",     # 2
        "HLT_IsoTkMu27_v*",     # 3
        "HLT_Mu17_v*",          # 4
        "HLT_Mu19_v*",          # 5
        "HLT_Mu20_v*",          # 6
        "HLT_Mu27_v*",          # 8
        "HLT_L1SingleMu18_v*",  # 9
        "HLT_L1SingleMu25_v*",  # 10
        # "HLT_IsoMu30_v*",       # 6
        # "HLT_IsoTkMu30_v*",     # 7
        # "HLT_Mu17_v*",          # 8
        # "HLT_TkMu17_v*",        # 9
        # "HLT_Mu19_v*",          # 10
        # "HLT_TkMu19_v*",        # 11
        # "HLT_Mu20_v*",          # 12
        # "HLT_TkMu20_v*",        # 13
        # "HLT_Mu27_v*",          # 14
        # "HLT_TkMu27_v*",        # 15
        # "HLT_IsoMu22_eta2p1",   # 16
        # "HLT_IsoTkMu22_eta2p1", # 17
        # "HLT_IsoMu24_eta2p1",   # 18
        # "HLT_IsoTkMu24_eta2p1", # 19
        # "HLT_Mu50_v*",          # 20
        # "HLT_TkMu50_v*"         # 21
    ),
    muon_trigger_DRMAX=cms.untracked.double(0.1),

    # MET - triggers
    # 2^{x}
    met_trigger_patterns = cms.vstring(
                                                            # x =
        "HLT_MET200_v*",                                    #0
        "HLT_MET250_v*",                                    #1
        "HLT_MET300_v*",                                    #2
        "HLT_PFMET90_PFMHT90_IDTight_v*",                   #3
        "HLT_PFMET100_PFMHT100_IDTight_v*",                 #4
        "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v*",          #5
        "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v*", #6
        "HLT_PFMET110_PFMHT110_IDTight_v*",                 #7
        "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v*", #8
        "HLT_PFMET120_BTagCSV_p067_v*",                     #9
        "HLT_PFMET120_PFMHT120_IDTight_v*",                 #10
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*",          #11
        "HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v*", #12
        "HLT_PFMET130_PFMHT130_IDTight_v*",                 #13
        "HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v*", #14
        "HLT_PFMET140_PFMHT140_IDTight_v*",                 #15
        "HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v*", #16
        "HLT_PFMET170_NotCleaned_v*",                       #17
        "HLT_PFMET170_HBHECleaned_v*",                      #18
        "HLT_PFMET170_JetIdCleaned_v*",                     #19
        "HLT_PFMET170_NoiseCleaned_v*",                     #20
        "HLT_PFMET200_NotCleaned_v*",                       #21
        "HLT_PFMET200_HBHECleaned_v*",                      #22
        "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",             #23
        "HLT_PFMET250_HBHECleaned_v*",                      #24
        "HLT_PFMET300_v*",                                  #25
        "HLT_PFMET300_HBHECleaned_v*",                      #26
        "HLT_PFMET400_v*",                                  #27
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",         #28
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",  #29
        "HLT_PFMETTypeOne110_PFMHT110_IDTight_v*"           #30
    ),

    # Electron - triggers
    # 2^{x}
    electron_trigger_patterns=cms.vstring(
                                # x =
        "HLT_Ele27_WPTight_Gsf_v*",      # 0   lowest unprescaled in 2016
        "HLT_Ele32_WPTight_Gsf_v*"       # 1   exists in 2017 (But must be emulated https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#Emulation_of_HLT_Ele32_WPTight_G) and 2018

    ),

    # Electrons - 
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt"),
    beamspotName = cms.InputTag('offlineBeamSpot'),
    conversionsName = cms.InputTag('conversions'),
    rhoName = cms.InputTag('fixedGridRhoFastjetAll'),

    #Flags
    isData = cms.untracked.bool(False),
    hasGenZ = cms.untracked.bool(False),
    hasGenTt = cms.untracked.bool(False),

    store_muons = cms.untracked.bool(True),
    store_electrons = cms.untracked.bool(True),    

    # For MC only
    genParticles         = cms.InputTag("genParticles"),
    genEventInfo         = cms.InputTag('generator'),
    lheEventInfo         = cms.InputTag('externalLHEProducer'),
    lheRunInfo           = cms.InputTag('externalLHEProducer'),
    genZLeptonCollection = cms.InputTag(''),
    genTtCollection = cms.InputTag(''),
    
    particleLevelLeptonCollection = cms.InputTag("particleLevel:leptons"),
    printLHE = cms.bool(False),
    
    # settings of UL aMC@NLO
    # (muRmuF up, muRmuF down, muR up, muR down, muF up, muF down)
    genWeights = cms.vint32(1005, 1009, 1004, 1007, 1002, 1003),
    pdfWeights = cms.vint32(1214, 1316)

)
