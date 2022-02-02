// system include files
#include <memory>

// CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "ZCounting/ZUtils/interface/triggertool.h"
#include "ZCounting/ZUtils/interface/Helper.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

//
// class declaration
//


class TnPPairTreeProducer :
    public edm::EDAnalyzer
{
public:
    TnPPairTreeProducer(const edm::ParameterSet&);
    ~TnPPairTreeProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum MuonIDTypes { NoneID, LooseID, MediumID, TightID, CustomID};
    enum MuonIsoTypes { NoneIso, TrackerIso, PFIso };

private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool isMuonTrigger();
    bool isMuonTriggerObj(const double eta, const double phi);
    bool isMuonTriggerObjEmulated(const double eta, const double phi, const long long unsigned eventNumber);
    bool passMuonIso(const reco::Muon&);
    bool passMuonID(const reco::Muon&, const reco::Vertex&);

    void clearVariables();
    void clearTagVariables();
    void clearProbeVariables();

    // ----------member data ---------------------------

    const double MUON_MASS = 0.105658369;

    edm::Service<TFileService> fs;
    TTree *tree_;

    // --- input
    double MassMin_;
    double MassMax_;

    // trigger
    const edm::InputTag triggerResultsInputTag_;
    
    HLTConfigProvider hltConfigProvider_;
    triggertool *triggers;

    bool emulateTrigger_;
    TFile *_fileHLTEmulation = 0;

    // Tracks
    std::string fTrackName;
    edm::EDGetTokenT<std::vector<reco::Track>> fTrackName_token;

    // Muons
    std::string fMuonName;
    edm::EDGetTokenT<std::vector<reco::Muon>> fMuonName_token;
    std::vector<std::string> fMuonHLTNames;
    double fMuonHLTDRMAX;

    std::string IDTypestr_;
    std::string IsoTypestr_;
    MuonIDTypes IDType_{NoneID};
    MuonIsoTypes IsoType_{NoneIso};
    double IsoCut_;
    double DxyCut_;
    double DzCut_;

    double PtCutL1_;
    double EtaCutL1_;
    double PtCutL2_;
    double EtaCutL2_;

    // PV
    std::string fPVName;
    edm::EDGetTokenT<std::vector<reco::Vertex>> fPVName_token;
    double VtxNTracksFitCut_;
    double VtxNdofCut_;
    double VtxAbsZCut_;
    double VtxRhoCut_;

    // --- output
    int nPV_;
    unsigned run_;
    unsigned ls_;
    long long unsigned eventNumber_;

    float pt1_;
    float eta1_;
    float phi1_;
    float q1_;
    float pfIso1_;
    float tkIso1_;
    float dxy1_;        // transvers distance of muon from nominal interaction point
    float dz1_;         // longitudinal distance of muon from nominal interaction point    
    bool is1IsoMu27_;

    float pt2_;
    float eta2_;
    float phi2_;
    float q2_;
    float pfIso2_;
    float tkIso2_;
    float dxy2_;        // transvers distance of muon from nominal interaction point
    float dz2_;         // longitudinal distance of muon from nominal interaction point  
    float nTrackerLayers2_;
    float nValidPixelHits2_;
    int trackAlgo2_;    // algorithm used in the track reconstruction, to identify muon seeded track
    bool is2IsoMu27_;


    bool is2HLT_;
    bool isSel_;
    bool isGlo_;
    bool isSta_;
    bool isTrk_;

    float dilepMass_;
    float delR_;

    float drhoPV_;      // transvers distance of PV from nominal interaction point
    float dzPV_;       // longitudinal distance of PV from nominal interaction point

    // --- things for trigger emulation
    // filter tag for HLT_IsoMu24_v11:

};


//define this as a plug-in
DEFINE_FWK_MODULE(TnPPairTreeProducer);
