import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCounting',
    pileupSummaryInfoCollection=cms.InputTag("slimmedAddPileupInfo"),

    era = cms.string("None"),
    # Primary vertex - source
    edmPVName = cms.InputTag("offlineSlimmedPrimaryVertices"),

    # Trigger - sources
    trigger_bits = cms.InputTag("TriggerResults","","HLT"),
    trigger_objects = cms.InputTag("slimmedPatTrigger"),

    # MET - sources
    pat_met_PF = cms.InputTag("slimmedMETs"),
    pat_met_puppi = cms.InputTag("slimmedMETsPuppi"),

    # Muons - source
    pat_muons = cms.InputTag("slimmedMuons","","PAT"),

    # Muon - triggers
    # 2^{x}
    muon_trigger_patterns=cms.vstring(
                                # x =
        "HLT_L1SingleMu18_v*",  # 0
        "HLT_L1SingleMu25_v*",  # 1
        "HLT_IsoMu24_v*",       # 2
        "HLT_IsoTkMu24_v*",     # 3
        "HLT_IsoMu27_v*",       # 4
        "HLT_IsoTkMu27_v*",     # 5
        "HLT_IsoMu30_v*",       # 6
        "HLT_IsoTkMu30_v*",     # 7
        "HLT_Mu17_v*",          # 8
        "HLT_TkMu17_v*",        # 9
        "HLT_Mu19_v*",          # 10
        "HLT_TkMu19_v*",        # 11
        "HLT_Mu20_v*",          # 12
        "HLT_TkMu20_v*",        # 13
        "HLT_Mu27_v*",          # 14
        "HLT_TkMu27_v*",        # 15
        "HLT_IsoMu22_eta2p1",   # 16
        "HLT_IsoTkMu22_eta2p1", # 17
        "HLT_IsoMu24_eta2p1",   # 18
        "HLT_IsoTkMu24_eta2p1", # 19
        "HLT_Mu50_v*",          # 20
        "HLT_TkMu50_v*"         # 21
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
    # more muon triggers
    met_trigger_patterns_ext = cms.vstring(
                                                            # x =
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_v*",          #0
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v*",   #1
        "HLT_PFMETTypeOne130_PFMHT130_IDTight_v*",          #2
        "HLT_PFMETTypeOne140_PFMHT140_IDTight_v*",          #3
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*",      #4
        "HLT_TripleJet110_35_35_Mjj650_PFMET110_v*",        #5
        "HLT_TripleJet110_35_35_Mjj650_PFMET120_v*",        #6
        "HLT_TripleJet110_35_35_Mjj650_PFMET130_v*",        #7
        "HLT_CaloMET70_HBHECleaned_v*",                     #8
        "HLT_CaloMET80_NotCleaned_v*",                      #9
        "HLT_CaloMET80_HBHECleaned_v*",                     #10
        "HLT_CaloMHT90_v*",                                 #11
        "HLT_CaloMET90_NotCleaned_v*",                      #12
        "HLT_CaloMET90_HBHECleaned_v*",                     #13
        "HLT_CaloMET100_NotCleaned_v*",                     #14
        "HLT_CaloMET100_HBHECleaned_v*",                    #15
        "HLT_CaloMET110_NotCleaned_v*",                     #16
        "HLT_CaloMET250_NotCleaned_v*",                     #17
        "HLT_CaloMET250_HBHECleaned_v*",                    #18
        "HLT_CaloMET300_HBHECleaned_v*",                    #19
        "HLT_CaloMET350_HBHECleaned_v*",                    #20
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_v*",     #21
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v*",    #22
        "HLT_DiJet110_35_Mjj650_PFMET110_v*",                           #23
        "HLT_DiJet110_35_Mjj650_PFMET120_v*",                           #24
        "HLT_DiJet110_35_Mjj650_PFMET130_v*",                           #25
        "HLT_L1ETMHadSeeds_v*",                                         #26
        "HLT_DiCentralPFJet55_PFMET110_v*",                             #27
        "HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v*",         #28
        "HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v*",         #29
        "HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v*",         #30
        "HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v*",         #31
        "HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v*"          #32
    ),


    #Flags
    isData = cms.untracked.bool(False),
    hasGenZ = cms.untracked.bool(False),
    hasGenTt = cms.untracked.bool(False),

    # For MC only
    genParticles         = cms.InputTag("prunedGenParticles"),
    genEventInfo         = cms.InputTag('generator'),
    lheEventInfo         = cms.InputTag('externalLHEProducer'),
    genZLeptonCollection = cms.InputTag(''),
    genTtCollection = cms.InputTag(''),

    # settings of UL aMC@NLO
    genWeights = cms.vint32(1005, 1009, 1004, 1007, 1002, 1003),
    pdfWeights = cms.vint32(1214, 1316)

)
