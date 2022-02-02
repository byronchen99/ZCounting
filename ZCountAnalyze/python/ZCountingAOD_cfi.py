import FWCore.ParameterSet.Config as cms

zcounting = cms.EDAnalyzer('ZCountingAOD',
    pileupSummaryInfoCollection=cms.InputTag("addPileupInfo"),

    era = cms.string("None"),
    # Primary vertex - source
    edmPVName = cms.InputTag("offlinePrimaryVertices"),
    
    # Muons - source
    reco_muons = cms.InputTag("muons","","RECO"),

    # Trigger - sources
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT'),
    TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),

    # Muon - triggers
    # 2^{x}
    muon_trigger_patterns=cms.vstring(
                                # x =
        # "HLT_L1SingleMu18_v*",  # 0
        # "HLT_L1SingleMu25_v*",  # 1
        "HLT_IsoMu24_v*",       # 2
        # "HLT_IsoTkMu24_v*",     # 3
        "HLT_IsoMu27_v*",       # 4
        # "HLT_IsoTkMu27_v*",     # 5
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


    #Flags
    isData = cms.untracked.bool(False),
    hasGenZ = cms.untracked.bool(False),
    hasGenTt = cms.untracked.bool(False),

    # For MC only
    genParticles         = cms.InputTag("genParticles"),
    genEventInfo         = cms.InputTag('generator'),
    lheEventInfo         = cms.InputTag('externalLHEProducer'),
    genZLeptonCollection = cms.InputTag(''),
    genTtCollection = cms.InputTag(''),

    # settings of UL aMC@NLO
    genWeights = cms.vint32(1005, 1009, 1004, 1007, 1002, 1003),
    pdfWeights = cms.vint32(1214, 1316)

)
