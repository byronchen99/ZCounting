
// system include files
#include <memory>

// CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "ZCounting/TnPPairTreeProducer/interface/triggertool.h"
#include "ZCounting/TnPPairTreeProducer/interface/MuonDefinitions.h"

// ROOT includes
#include "TTree.h"


//
// class declaration
//

class MuonTreeProducer :
    public edm::EDAnalyzer,
    public MuonDefinitions
{
public:
    MuonTreeProducer(const edm::ParameterSet&);
    ~MuonTreeProducer() {};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    bool isMuonTrigger();
    bool isMuonTriggerObj(const double eta, const double phi);

    void clearVariables();
    void clearMuonVariables();

    // ----------member data ---------------------------
    edm::Service<TFileService> fs;
    TTree *tree_;

    HLTConfigProvider hltConfigProvider_;
    triggertool *triggers;

    // --- input

    // trigger
    const edm::InputTag triggerResultsInputTag_;

    // Tracks
    const std::string fTrackName;
    const edm::EDGetTokenT<std::vector<reco::Track>> fTrackName_token;

    // Muons
    const std::string fMuonName;
    edm::EDGetTokenT<std::vector<reco::Muon>> fMuonName_token;
    const std::vector<std::string> fMuonHLTNames;

    const double PtCutL_;
    const double EtaCutL_;

    // PV
    const std::string fPVName;
    edm::EDGetTokenT<std::vector<reco::Vertex>> fPVName_token;
    const double VtxNTracksFitCut_;
    const double VtxNdofCut_;
    const double VtxAbsZCut_;
    const double VtxRhoCut_;

    // --- output
    int nPV_;
    unsigned run_;
    unsigned ls_;
    unsigned eventNumber_;

    float pt_;
    float eta_;
    float phi_;
    float q_;
    float pfIso_;
    float tkIso_;
    float dxy_;
    float dz_;

    bool is_HLT_Mu17_;
    bool is_HLT_IsoMu27_;

};


//define this as a plug-in
DEFINE_FWK_MODULE(MuonTreeProducer);
