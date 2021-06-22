// -*- C++ -*-
//
// Package:    ZCounting/ZCountAnalyze
// Class:      ZCounting
//
/**\class ZCounting ZCounting.cc ZCounting/ZCountAnalyze/plugins/ZCounting.cc

 Description:
     Purpose is to derive MC corrections of Z counting efficiency in respect of true (MC) Z efficiency

     Select events which have:
        - a decay of Z to mu mu
        - available generator level information of both muons
     Store:
        - information of the gen muons
        - information of the corresponding reco muons (if given)



 Implementation:

*/
//
// Original Author:  David Walter
//         Created:  Sat, 29 Jun 2019 13:26:06 GMT
//
//


// system include files
#include <memory>

// CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/RegexMatch.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

// ROOT includes
#include "TTree.h"
#include "Math/Vector4D.h"
#include <TLorentzVector.h>

// typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LV
//
// class declaration
//


class ZCounting : //public edm::one::EDAnalyzer<edm::one::SharedResources>
    public edm::EDAnalyzer
{
public:
    explicit ZCounting(const edm::ParameterSet&);
    ~ZCounting();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void clearVariables();

    std::string get_triggerPath(std::string, const edm::TriggerNames&);
    std::vector<pat::TriggerObjectStandAlone> get_muonTriggerObjects(const std::vector<pat::TriggerObjectStandAlone> &,
                                                                     const edm::TriggerNames&,
                                                                     const std::string&
                                                                     );
    bool isTriggerObject(const std::vector<pat::TriggerObjectStandAlone>&, const pat::Muon&);
    bool isGoodPV(const reco::Vertex&);
    bool isValidTrack(const reco::Track&);
    bool customIsTightMuon(const pat::Muon&);
    double pfIso(const pat::Muon&);
    double tkIso(const pat::Muon&);
    int getMuonID(const pat::Muon&, const reco::Vertex&);
    double dxy(const pat::Muon&, const reco::Vertex&);
    double dz(const pat::Muon&, const reco::Vertex&);
    double pointsDistance(const reco::Candidate::Point&, const reco::Candidate::Point&);
    bool isPVClosestVertex(const std::vector<reco::Vertex>&, const pat::Muon&);

    bool isTau(const reco::GenParticle* lepton)const;
    const reco::GenParticle* tauDaughter(const reco::GenParticle* tau)const;


    // ----------member data ---------------------------

    edm::Service<TFileService> fs;
    TTree *tree_;

    const double MUON_MASS = 0.105658369;

    TLorentzVector vMuon;
    TLorentzVector vAntiMuon;

    // --- input
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
    edm::EDGetTokenT<pat::MuonCollection> muonCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > pvCollection_;
    edm::EDGetTokenT<std::vector<GenZDecayProperties> > genZCollection_;
    edm::EDGetTokenT<TtGenEvent> genTtEvent_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoCollection_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;

    edm::EDGetTokenT< double > prefweightECAL_token;
    edm::EDGetTokenT< double > prefweightupECAL_token;
    edm::EDGetTokenT< double > prefweightdownECAL_token;

    edm::EDGetTokenT< double > prefweightMuon_token;
    edm::EDGetTokenT< double > prefweightupMuon_token;
    edm::EDGetTokenT< double > prefweightdownMuon_token;
    edm::EDGetTokenT< double > prefweightupSystMuon_token;
    edm::EDGetTokenT< double > prefweightdownSystMuon_token;
    edm::EDGetTokenT< double > prefweightupStatMuon_token;
    edm::EDGetTokenT< double > prefweightdownStatMuon_token;

    edm::EDGetTokenT< double > prefweightMuon2016H_token;
    edm::EDGetTokenT< double > prefweightupMuon2016H_token;
    edm::EDGetTokenT< double > prefweightdownMuon2016H_token;
    edm::EDGetTokenT< double > prefweightupSystMuon2016H_token;
    edm::EDGetTokenT< double > prefweightdownSystMuon2016H_token;
    edm::EDGetTokenT< double > prefweightupStatMuon2016H_token;
    edm::EDGetTokenT< double > prefweightdownStatMuon2016H_token;

    std::string era_;

    bool hltChanged_;
    std::vector<std::string> muonTriggerPatterns_;
    std::vector<std::string> muonTriggerPaths_;
    // max dR matching between muon and hlt object
    double DRMAX;

    // flags
    bool hasGenZ_;
    bool hasGenTt_;

    // primary vertex cuts
    double VtxNTracksFitCut_;
    double VtxNdofCut_;
    double VtxAbsZCut_;
    double VtxRhoCut_;

    // --- output

    // ... for event info
    unsigned int runNumber_;
    unsigned int lumiBlock_;
    unsigned int eventNumber_;

    float eventweight_;

    float prefiringweightECAL_;
    float prefiringweightECALup_;
    float prefiringweightECALdown_;
    float prefiringweightMuon_;
    float prefiringweightMuonup_;
    float prefiringweightMuondown_;
    float prefiringweightMuonupSyst_;
    float prefiringweightMuondownSyst_;
    float prefiringweightMuonupStat_;
    float prefiringweightMuondownStat_;

    float prefiringweightMuon2016H_;
    float prefiringweightMuonup2016H_;
    float prefiringweightMuondown2016H_;
    float prefiringweightMuonupSyst2016H_;
    float prefiringweightMuondownSyst2016H_;
    float prefiringweightMuonupStat2016H_;
    float prefiringweightMuondownStat2016H_;

    int nPU_;
    int nPV_;
    int decayMode_;
    float z_genMass_;
    float z_recoMass_;

    int muon_genRecoMatches_;
    int muon_genRecoObj_;
    float muon_genPt_;
    float muon_genEta_;
    float muon_genPhi_;
    float muon_genVtxToPV_;

    int antiMuon_genRecoMatches_;
    int antiMuon_genRecoObj_;
    float antiMuon_genPt_;
    float antiMuon_genEta_;
    float antiMuon_genPhi_;
    float antiMuon_genVtxToPV_;

    unsigned int nMuon_;
    std::vector<int> muon_charge_;
    std::vector<float> muon_tkRelIso_;
    std::vector<float> muon_pfRelIso04_all_;
    std::vector<int> muon_ID_;
    std::vector<int> muon_genPartFlav_;
    std::vector<int> muon_triggerBits_;
    std::vector<float> muon_dxy_;
    std::vector<float> muon_dz_;
    std::vector<float> muon_pt_;
    std::vector<float> muon_eta_;
    std::vector<float> muon_phi_;
    std::vector<bool> muon_isFromPV_;
    std::vector<bool> muon_isMedium_;
    std::vector<bool> muon_isStandalone_;
    std::vector<bool> muon_isTracker_;
    std::vector<bool> muon_isGlobal_;
    std::vector<bool> muon_isPFCand_;
    std::vector<float> muon_normChi2_;
    std::vector<int> muon_nTrackHits_;
    std::vector<int> muon_nStations_;
    std::vector<float> muon_SegmentCompatibility_;
    std::vector<float> muon_chi2LocalPosition_;
    std::vector<float> muon_trkKink_;
    std::vector<int> muon_nPixelHits_;
    std::vector<int> muon_nTrackerLayers_;
    std::vector<float> muon_validFraction_;

};

//
// constructors and destructor
//
ZCounting::ZCounting(const edm::ParameterSet& iConfig):
    triggerBits_     (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigger_bits"))),
    triggerObjects_  (consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("trigger_objects"))),
    muonCollection_  (consumes<pat::MuonCollection> (iConfig.getParameter<edm::InputTag>("pat_muons"))),
    pvCollection_    (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("edmPVName"))),
    genZCollection_  (consumes<std::vector<GenZDecayProperties> > (iConfig.getParameter<edm::InputTag>("genZLeptonCollection"))),
    genTtEvent_  (consumes<TtGenEvent> (iConfig.getParameter<edm::InputTag>("genTtCollection"))),
    pileupInfoCollection_  (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfoCollection"))),
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")))
{
    LogDebug("ZCounting")<<"ZCounting(...)";

    era_ = iConfig.getParameter<std::string>("era");

    hasGenZ_ = iConfig.getUntrackedParameter<bool>("hasGenZ");
    hasGenTt_ = iConfig.getUntrackedParameter<bool>("hasGenTt");

    muonTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("muon_trigger_patterns");
    DRMAX = iConfig.getUntrackedParameter<double>("muon_trigger_DRMAX");

    VtxNTracksFitCut_ = iConfig.getUntrackedParameter<double>("VtxNTracksFitMin");
    VtxNdofCut_       = iConfig.getUntrackedParameter<double>("VtxNdofMin");
    VtxAbsZCut_       = iConfig.getUntrackedParameter<double>("VtxAbsZMax");
    VtxRhoCut_        = iConfig.getUntrackedParameter<double>("VtxRhoMax");
    hltChanged_ = true;

    prefweightECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECAL"));
    prefweightupECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECALUp"));
    prefweightdownECAL_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECALDown"));

    prefweightMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuon"));
    prefweightupMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonUp"));
    prefweightdownMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonDown"));
    prefweightupSystMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonSystUp"));
    prefweightdownSystMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonSystDown"));
    prefweightupStatMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonStatUp"));
    prefweightdownStatMuon_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuonStatDown"));

    if(era_ == "2016postVFP"){
        std::cout<<"2016 post VFP era"<<std::endl;

        prefweightMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuon"));
        prefweightupMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonUp"));
        prefweightdownMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonDown"));
        prefweightupSystMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonSystUp"));
        prefweightdownSystMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonSystDown"));
        prefweightupStatMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonStatUp"));
        prefweightdownStatMuon2016H_token = consumes< double >(edm::InputTag("prefiringweight2016H:nonPrefiringProbMuonStatDown"));
    }

}


ZCounting::~ZCounting()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZCounting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    LogDebug("ZCounting::analyze");

    this->clearVariables();

    // Event info
    runNumber_ = iEvent.id().run();
    lumiBlock_ = iEvent.id().luminosityBlock();
    eventNumber_ = iEvent.id().event();


    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
    edm::Handle<pat::MuonCollection> muonCollection;
    edm::Handle<std::vector<reco::Vertex> > pvCollection;
    edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoCollection;

    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(muonCollection_, muonCollection);
    iEvent.getByToken(pvCollection_, pvCollection);
    iEvent.getByToken(pileupInfoCollection_, pileupInfoCollection);

    // --- >>> prefiring
    edm::Handle< double > theprefweightECAL;
    edm::Handle< double > theprefweightupECAL;
    edm::Handle< double > theprefweightdownECAL;
    edm::Handle< double > theprefweightMuon;
    edm::Handle< double > theprefweightupMuon;
    edm::Handle< double > theprefweightdownMuon;
    edm::Handle< double > theprefweightupsystMuon;
    edm::Handle< double > theprefweightdownsystMuon;
    edm::Handle< double > theprefweightupstatMuon;
    edm::Handle< double > theprefweightdownstatMuon;

    iEvent.getByToken(prefweightECAL_token, theprefweightECAL ) ;
    iEvent.getByToken(prefweightupECAL_token, theprefweightupECAL ) ;
    iEvent.getByToken(prefweightdownECAL_token, theprefweightdownECAL ) ;
    iEvent.getByToken(prefweightMuon_token, theprefweightMuon ) ;
    iEvent.getByToken(prefweightupMuon_token, theprefweightupMuon ) ;
    iEvent.getByToken(prefweightdownMuon_token, theprefweightdownMuon ) ;
    iEvent.getByToken(prefweightupSystMuon_token, theprefweightupsystMuon ) ;
    iEvent.getByToken(prefweightdownSystMuon_token, theprefweightdownsystMuon ) ;
    iEvent.getByToken(prefweightupStatMuon_token, theprefweightupstatMuon ) ;
    iEvent.getByToken(prefweightdownStatMuon_token, theprefweightdownstatMuon ) ;

    prefiringweightECAL_ =(*theprefweightECAL);
    prefiringweightECALup_ =(*theprefweightupECAL);
    prefiringweightECALdown_ =(*theprefweightdownECAL);
    prefiringweightMuon_ =(*theprefweightMuon);
    prefiringweightMuonup_ =(*theprefweightupMuon);
    prefiringweightMuondown_ =(*theprefweightdownMuon);
    prefiringweightMuonupSyst_ =(*theprefweightupMuon);
    prefiringweightMuondownSyst_ =(*theprefweightdownMuon);
    prefiringweightMuonupStat_ =(*theprefweightupMuon);
    prefiringweightMuondownStat_ =(*theprefweightdownMuon);

    if(era_ == "2016postVFP"){
        edm::Handle< double > theprefweightMuon2016H;
        edm::Handle< double > theprefweightupMuon2016H;
        edm::Handle< double > theprefweightdownMuon2016H;
        edm::Handle< double > theprefweightupsystMuon2016H;
        edm::Handle< double > theprefweightdownsystMuon2016H;
        edm::Handle< double > theprefweightupstatMuon2016H;
        edm::Handle< double > theprefweightdownstatMuon2016H;

        iEvent.getByToken(prefweightMuon2016H_token, theprefweightMuon2016H);
        iEvent.getByToken(prefweightupMuon2016H_token, theprefweightupMuon2016H);
        iEvent.getByToken(prefweightdownMuon2016H_token, theprefweightdownMuon2016H);
        iEvent.getByToken(prefweightupSystMuon2016H_token, theprefweightupsystMuon2016H);
        iEvent.getByToken(prefweightdownSystMuon2016H_token, theprefweightdownsystMuon2016H);
        iEvent.getByToken(prefweightupStatMuon2016H_token, theprefweightupstatMuon2016H);
        iEvent.getByToken(prefweightdownStatMuon2016H_token, theprefweightdownstatMuon2016H);

        prefiringweightMuon2016H_ =(*theprefweightMuon2016H);
        prefiringweightMuonup2016H_ =(*theprefweightupMuon2016H);
        prefiringweightMuondown2016H_ =(*theprefweightdownMuon2016H);
        prefiringweightMuonupSyst2016H_ =(*theprefweightupsystMuon2016H);
        prefiringweightMuondownSyst2016H_ =(*theprefweightdownsystMuon2016H);
        prefiringweightMuonupStat2016H_ =(*theprefweightupstatMuon2016H);
        prefiringweightMuondownStat2016H_ =(*theprefweightdownstatMuon2016H);
    }

    // <<< ---

    const reco::GenParticle* genLepton = 0;
    const reco::GenParticle* genAntiLepton = 0;
    if(!iEvent.isRealData()){

        edm::Handle<GenEventInfoProduct> evt_info;
        iEvent.getByToken(genEventInfo_, evt_info);

        eventweight_ = evt_info->weight();

        if(hasGenZ_){
            edm::Handle<std::vector<GenZDecayProperties> > genZCollection;
            iEvent.getByToken(genZCollection_, genZCollection);
            decayMode_ = (genZCollection->size() != 0 ? genZCollection->at(0).decayMode() : 0);
            genLepton = genZCollection->at(0).stableLepton();
            genAntiLepton = genZCollection->at(0).stableAntiLepton();

        }
        if(hasGenTt_){
            LogDebug("ZCounting::analyze")<<"hasGenZ";

            edm::Handle<TtGenEvent> genTtEvent;
            iEvent.getByToken(genTtEvent_, genTtEvent);
            if(genTtEvent->lepton()){
                LogDebug("ZCounting::analyze")<<"hasLepton";
                genLepton = genTtEvent->lepton();
                if(this->isTau(genLepton)){
                    decayMode_ += 150000;
                    genLepton = this->tauDaughter(genTtEvent->lepton());
                }

            }
            if(genTtEvent->leptonBar()){
                LogDebug("ZCounting::analyze")<<"hasAntiLepton";
                genAntiLepton = genTtEvent->leptonBar();
                if(this->isTau(genAntiLepton)){
                    decayMode_ += 150000;
                    genAntiLepton = this->tauDaughter(genTtEvent->leptonBar());
                }
            }
            LogDebug("ZCounting::analyze")<<"setDecayMode";
            if(genLepton)
                decayMode_ += 100*std::abs(genLepton->pdgId());
            if(genAntiLepton)
                decayMode_ += std::abs(genAntiLepton->pdgId());
        }
    }

    // --- get muon trigger names
    LogDebug("ZCounting") << "get trigger names";
    const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
    if(hltChanged_ == true){
        muonTriggerPaths_.clear();
        for(const std::string pattern: muonTriggerPatterns_){
            muonTriggerPaths_.push_back(get_triggerPath(pattern, triggerNames));
        }
        hltChanged_=false;
    }

    // PV selection
    const std::vector<reco::Vertex> *pvCol = pvCollection.product();
    reco::Vertex pv = *pvCol->begin();
    nPV_ = 0;
    for (auto const& itVtx : *pvCol) {
        if(!isGoodPV(itVtx))
            continue;
        if (nPV_ == 0)
            pv = itVtx;
        nPV_++;
    }

    if(genLepton){
        muon_genPt_ = genLepton->pt();
        muon_genEta_ = genLepton->eta();
        muon_genPhi_ = genLepton->phi();
        if(nPV_){
            muon_genVtxToPV_ = pointsDistance(pv.position(), genLepton->vertex());
        }
    }
    if(genAntiLepton){
        antiMuon_genPt_ = genAntiLepton->pt();
        antiMuon_genEta_ = genAntiLepton->eta();
        antiMuon_genPhi_ = genAntiLepton->phi();
        if(nPV_){
            antiMuon_genVtxToPV_ = pointsDistance(pv.position(), genAntiLepton->vertex());
        }
    }
    if(genLepton && genAntiLepton){
        z_genMass_ = (genLepton->p4() + genAntiLepton->p4()).M();
    }

    // --- store all reco muons
    LogDebug("ZCounting") << "find reco muons";
    for (pat::Muon mu : *muonCollection){
        if(std::abs(mu.eta()) > 2.4) continue;
        if(mu.pt() < 15) continue;

        muon_charge_.push_back(mu.charge());
        muon_pt_.push_back(mu.pt());
        muon_eta_.push_back(mu.eta());
        muon_phi_.push_back(mu.phi());
        muon_ID_.push_back(getMuonID(mu, pv));

        muon_isMedium_.push_back(mu.isMediumMuon());
        muon_isStandalone_.push_back(mu.isStandAloneMuon());
        muon_isTracker_.push_back(mu.isTrackerMuon());
        muon_isGlobal_.push_back(mu.isGlobalMuon());
        muon_isPFCand_.push_back(mu.isPFMuon());
        muon_nStations_.push_back(mu.numberOfMatchedStations());
        muon_SegmentCompatibility_.push_back(mu.segmentCompatibility());
        muon_chi2LocalPosition_.push_back(mu.combinedQuality().chi2LocalPosition);
        muon_trkKink_.push_back(mu.combinedQuality().trkKink);

        if(mu.isTrackerMuon()){
            muon_nPixelHits_.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
            muon_nTrackerLayers_.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
            muon_validFraction_.push_back(mu.innerTrack()->validFraction());
        }
        else{
            muon_nPixelHits_.push_back(-1);
            muon_nTrackerLayers_.push_back(-1);
            muon_validFraction_.push_back(-1);
        }

        if(mu.isGlobalMuon()){
            muon_normChi2_.push_back(mu.globalTrack()->normalizedChi2());
            muon_nTrackHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
        }
        else{
            muon_normChi2_.push_back(-1);
            muon_nTrackHits_.push_back(-1);
        }

        muon_tkRelIso_.push_back(tkIso(mu));
        muon_pfRelIso04_all_.push_back(pfIso(mu));
        muon_dxy_.push_back(dxy(mu, pv));
        muon_dz_.push_back(dz(mu, pv));
        muon_isFromPV_.push_back(isPVClosestVertex(*pvCol, mu));

        int bits_ = 0;
        for(unsigned j = 0; j < muonTriggerPaths_.size(); ++j){
            std::vector<pat::TriggerObjectStandAlone> tObjCol = get_muonTriggerObjects(*triggerObjects, triggerNames, muonTriggerPaths_.at(j));
            bits_ += std::pow(2,j) * isTriggerObject(tObjCol, mu);
        }
        muon_triggerBits_.push_back(bits_);

        if(!iEvent.isRealData()){
            muon_genPartFlav_.push_back(mu.simFlavour());

            if(genLepton && mu.pdgId() == 13 && reco::deltaR(mu.eta(), mu.phi(), genLepton->eta(), genLepton->phi()) < 0.03){
                muon_genRecoMatches_++;
                if(muon_genRecoObj_ == -1
                    || std::abs(mu.pt() - muon_genPt_) < std::abs(muon_pt_[muon_genRecoObj_] - muon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    muon_genRecoObj_ = nMuon_;
                }
            }

            if(genAntiLepton && mu.pdgId() == -13 && reco::deltaR(mu.eta(), mu.phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                antiMuon_genRecoMatches_++;
                if(antiMuon_genRecoObj_ == -1
                    || std::abs(mu.pt() - antiMuon_genPt_) < std::abs(muon_pt_[antiMuon_genRecoObj_] - antiMuon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    antiMuon_genRecoObj_ = nMuon_;
                }
            }
        }

        nMuon_++;

    }

    if(antiMuon_genRecoObj_ != -1 && muon_genRecoObj_ != -1){
        vMuon.SetPtEtaPhiM(muon_pt_[muon_genRecoObj_], muon_eta_[muon_genRecoObj_], muon_phi_[muon_genRecoObj_], MUON_MASS);
        vAntiMuon.SetPtEtaPhiM(muon_pt_[antiMuon_genRecoObj_], muon_eta_[antiMuon_genRecoObj_], muon_phi_[antiMuon_genRecoObj_], MUON_MASS);

        z_recoMass_ = (vMuon + vAntiMuon).M();
    }

    for(std::vector<PileupSummaryInfo>::const_iterator puI = pileupInfoCollection->begin(); puI != pileupInfoCollection->end(); ++puI)
    {
        if(puI->getBunchCrossing() == 0){
            nPU_ = puI->getTrueNumInteractions();
        }
    }

    tree_->Fill();
}

// ------------ method called once each run just before starting event loop  ------------
void ZCounting::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    // edm::LogVerbatim("ZCounting") << "now at "<<iRun.id();
    std::cout<< "new run"<<std::endl;
    hltChanged_ = true;
}

// ------------ method called once each job just before starting event loop  ------------
void
ZCounting::beginJob()
{
    LogDebug("ZCounting")<<"beginJob()";

    if( !fs ){
        edm::LogError("ZCounting") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    // Event info
    tree_->Branch("runNumber", &runNumber_, "runNumber/i");
    tree_->Branch("lumiBlock", &lumiBlock_,"lumiBlock/i");
    tree_->Branch("eventNumber", &eventNumber_, "eventNumber/i");

    tree_->Branch("nPV", &nPV_,"nPV_/i");
    tree_->Branch("nPU", &nPU_,"nPU_/i");

    // weights
    // ME weight
    tree_->Branch("eventweight", &eventweight_, "eventweight_/f");
    // prefire weights
    tree_->Branch("prefiringweightECAL",         &prefiringweightECAL_,     "prefiringweightECAL_/f");
    tree_->Branch("prefiringweightECALup",       &prefiringweightECALup_,   "prefiringweightECALup_/f");
    tree_->Branch("prefiringweightECALdown",     &prefiringweightECALdown_, "prefiringweightECALdown_/f");
    tree_->Branch("prefiringweightMuon",         &prefiringweightMuon_,         "prefiringweightMuon_/f");
    tree_->Branch("prefiringweightMuonup",       &prefiringweightMuonup_,       "prefiringweightMuonup_/f");
    tree_->Branch("prefiringweightMuondown",     &prefiringweightMuondown_,     "prefiringweightMuondown_/f");
    tree_->Branch("prefiringweightMuonupSyst",   &prefiringweightMuonupSyst_,   "prefiringweightMuonupSyst_/f");
    tree_->Branch("prefiringweightMuondownSyst", &prefiringweightMuondownSyst_, "prefiringweightMuondownSyst_/f");
    tree_->Branch("prefiringweightMuonupStat",   &prefiringweightMuonupStat_,   "prefiringweightMuonupStat_/f");
    tree_->Branch("prefiringweightMuondownStat", &prefiringweightMuondownStat_, "prefiringweightMuondownStat_/f");

    if(era_ == "2016postVFP"){
        tree_->Branch("prefiringweightMuon2016H",         &prefiringweightMuon2016H_,         "prefiringweightMuon2016H_/f");
        tree_->Branch("prefiringweightMuonup2016H",       &prefiringweightMuonup2016H_,       "prefiringweightMuonup2016H_/f");
        tree_->Branch("prefiringweightMuondown2016H",     &prefiringweightMuondown2016H_,     "prefiringweightMuondown2016H_/f");
        tree_->Branch("prefiringweightMuonupSyst2016H",   &prefiringweightMuonupSyst2016H_,   "prefiringweightMuonupSyst2016H_/f");
        tree_->Branch("prefiringweightMuondownSyst2016H", &prefiringweightMuondownSyst2016H_, "prefiringweightMuondownSyst2016H_/f");
        tree_->Branch("prefiringweightMuonupStat2016H",   &prefiringweightMuonupStat2016H_,   "prefiringweightMuonupStat2016H_/f");
        tree_->Branch("prefiringweightMuondownStat2016H", &prefiringweightMuondownStat2016H_, "prefiringweightMuondownStat2016H_/f");
    }

    // gen level info
    tree_->Branch("decayMode", &decayMode_, "decayMode_/i");
    tree_->Branch("z_genMass", &z_genMass_,"z_genMass_/f");
    tree_->Branch("z_recoMass", &z_recoMass_,"z_recoMass_/f");

    tree_->Branch("muon_genVtxToPV", &muon_genVtxToPV_,"muon_genVtxToPV_/f");
    tree_->Branch("muon_genRecoMatches", &muon_genRecoMatches_,"muon_genRecoMatches_/i");
    tree_->Branch("muon_genRecoObj", &muon_genRecoObj_,"muon_genRecoObj_/I");
    tree_->Branch("muon_genPt", &muon_genPt_,"muon_genPt_/f");
    tree_->Branch("muon_genEta", &muon_genEta_,"muon_genEta_/f");
    tree_->Branch("muon_genPhi", &muon_genPhi_,"muon_genPhi_/f");

    tree_->Branch("antiMuon_genVtxToPV", &antiMuon_genVtxToPV_,"antiMuon_genVtxToPV_/f");
    tree_->Branch("antiMuon_genRecoMatches", &antiMuon_genRecoMatches_,"antiMuon_genRecoMatches_/i");
    tree_->Branch("antiMuon_genRecoObj", &antiMuon_genRecoObj_,"antiMuon_genRecoObj_/I");
    tree_->Branch("antiMuon_genPt", &antiMuon_genPt_,"antiMuon_genPt_/f");
    tree_->Branch("antiMuon_genEta", &antiMuon_genEta_,"antiMuon_genEta_/f");
    tree_->Branch("antiMuon_genPhi", &antiMuon_genPhi_,"antiMuon_genPhi_/f");

    // reco level info
    tree_->Branch("nMuon", &nMuon_,"nMuon_/s");
    tree_->Branch("Muon_charge", &muon_charge_);
    tree_->Branch("Muon_ID", &muon_ID_);
    tree_->Branch("Muon_genPartFlav", &muon_genPartFlav_);
    tree_->Branch("Muon_tkRelIso", &muon_tkRelIso_);
    tree_->Branch("Muon_pfRelIso04_all", &muon_pfRelIso04_all_);
    tree_->Branch("Muon_triggerBits", &muon_triggerBits_);
    tree_->Branch("Muon_dxy", &muon_dxy_);
    tree_->Branch("Muon_dz", &muon_dz_);
    tree_->Branch("Muon_pt", &muon_pt_);
    tree_->Branch("Muon_eta", &muon_eta_);
    tree_->Branch("Muon_phi", &muon_phi_);
    tree_->Branch("Muon_isFromPV", &muon_isFromPV_);

    tree_->Branch("Muon_isMedium", &muon_isMedium_);
    tree_->Branch("Muon_isStandalone", &muon_isStandalone_);
    tree_->Branch("Muon_isTracker", &muon_isTracker_);
    tree_->Branch("Muon_isGlobal", &muon_isGlobal_);
    tree_->Branch("Muon_isPFCand", &muon_isPFCand_);
    tree_->Branch("Muon_normChi2", &muon_normChi2_);
    tree_->Branch("Muon_nTrackHits", &muon_nTrackHits_);
    tree_->Branch("Muon_nStations", &muon_nStations_);
    tree_->Branch("Muon_SegmentCompatibility", &muon_SegmentCompatibility_);
    tree_->Branch("Muon_chi2LocalPosition", &muon_chi2LocalPosition_);
    tree_->Branch("Muon_validFraction", &muon_validFraction_);
    tree_->Branch("Muon_trkKink", &muon_trkKink_);
    tree_->Branch("Muon_nPixelHits", &muon_nPixelHits_);
    tree_->Branch("Muon_nTrackerLayers", &muon_nTrackerLayers_);

}

// ------------ method called once each job just after ending the event loop  ------------
void
ZCounting::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZCounting::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//--------------------------------------------------------------------------------------------------
void ZCounting::clearVariables(){
    LogDebug("ZCounting::clearVariables()");

    // Event info
    runNumber_ = 0;
    lumiBlock_ = 0;
    eventNumber_ = 0;

    eventweight_ = 0.;
    nPV_ = 0;
    nPU_ = 0;
    decayMode_ = 0;
    z_genMass_ = 0.;
    z_recoMass_ = 0;

    muon_genPt_ = 0.;
    muon_genEta_ = 0.;
    muon_genPhi_ = 0.;
    muon_genVtxToPV_ = -1;
    muon_genRecoMatches_ = 0;
    muon_genRecoObj_ = -1;

    antiMuon_genPt_ = 0.;
    antiMuon_genEta_ = 0.;
    antiMuon_genPhi_ = 0.;
    antiMuon_genVtxToPV_ = -1;
    antiMuon_genRecoMatches_ = 0;
    antiMuon_genRecoObj_ = -1;

    nMuon_ = 0;
    muon_charge_.clear();
    muon_genPartFlav_.clear();
    muon_ID_.clear();
    muon_tkRelIso_.clear();
    muon_pfRelIso04_all_.clear();
    muon_triggerBits_.clear();
    muon_dxy_.clear();
    muon_dz_.clear();
    muon_pt_.clear();
    muon_eta_.clear();
    muon_phi_.clear();
    muon_isFromPV_.clear();
    muon_isMedium_.clear();
    muon_isStandalone_.clear();
    muon_isTracker_.clear();
    muon_isGlobal_.clear();
    muon_isPFCand_.clear();
    muon_normChi2_.clear();
    muon_nTrackHits_.clear();
    muon_nStations_.clear();
    muon_SegmentCompatibility_.clear();
    muon_chi2LocalPosition_.clear();
    muon_trkKink_.clear();
    muon_nPixelHits_.clear();
    muon_nTrackerLayers_.clear();
    muon_validFraction_.clear();


}

//--------------------------------------------------------------------------------------------------
std::string ZCounting::get_triggerPath(std::string pattern, const edm::TriggerNames& triggerNames){
    std::string path = "";
    if (edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
        std::vector<std::vector<std::string>::const_iterator> matches =
            edm::regexMatch(triggerNames.triggerNames(), pattern);
        if (matches.empty()) {
            edm::LogWarning("ZCounting") << "requested pattern [" << pattern << "] does not match any HLT paths";
        }
        else {
            for (auto const& match : matches) {
                path = *match;
            }
        }
    }
    else {  // take full HLT path name given
        path = pattern;
    }
    return path;
}

//--------------------------------------------------------------------------------------------------
std::vector<pat::TriggerObjectStandAlone>
ZCounting::get_muonTriggerObjects(const std::vector<pat::TriggerObjectStandAlone> & tObjCol,
                                  const edm::TriggerNames& names,
                                  const std::string& path
                                  ){
    std::vector<pat::TriggerObjectStandAlone> muonTriggerObjects;
    for(pat::TriggerObjectStandAlone obj : tObjCol) {
        obj.unpackPathNames(names);
        const std::vector<std::string> pathNamesAll = obj.pathNames(false);

        // bool passedTrigger = false;
        for(const std::string iPath: pathNamesAll) {
            if(iPath == path){
                muonTriggerObjects.push_back(obj);
                continue;
            }
        }
    }

    return muonTriggerObjects;
}


//--------------------------------------------------------------------------------------------------
bool ZCounting::isTriggerObject(const std::vector<pat::TriggerObjectStandAlone> &tObjCol, const pat::Muon &muon){
    for(const pat::TriggerObjectStandAlone tObj : tObjCol) {
        if (reco::deltaR(muon.eta(), muon.phi(), tObj.eta(), tObj.phi()) < DRMAX) {
            return true;
        }
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
bool ZCounting::isGoodPV(const reco::Vertex &vtx){
    if (vtx.isFake())
        return false;
    if (vtx.tracksSize() < VtxNTracksFitCut_)
        return false;
    if (vtx.ndof() < VtxNdofCut_)
        return false;
    if (fabs(vtx.z()) > VtxAbsZCut_)
        return false;
    if (vtx.position().Rho() > VtxRhoCut_)
        return false;

    return true;
}

//--------------------------------------------------------------------------------------------------
bool ZCounting::isValidTrack(const reco::Track &trk){
    if(trk.hitPattern().trackerLayersWithMeasurement() >= 6 && trk.hitPattern().numberOfValidPixelHits() >= 1)
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool ZCounting::customIsTightMuon(const pat::Muon &mu){
    // From https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon releasing the PV criterias
    return (mu.isGlobalMuon()
        && mu.isPFMuon()
        && (mu.globalTrack()->normalizedChi2() < 10.)
        && (mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
        && (mu.numberOfMatchedStations() > 1)
        && (mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0)
        && (mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5));
}

//--------------------------------------------------------------------------------------------------
int ZCounting::getMuonID(const pat::Muon &mu, const reco::Vertex &vtx){
    if(mu.isTightMuon(vtx)) return 5;
    if(customIsTightMuon(mu)) return 4;
    if(mu.isGlobalMuon()) return 3;
    if(mu.isStandAloneMuon()) return 2;
    if(isValidTrack(*(mu.innerTrack()))) return 1;
    return 0;
}

//--------------------------------------------------------------------------------------------------
double ZCounting::tkIso(const pat::Muon &mu){
    return mu.isolationR03().sumPt / mu.pt();
}

//--------------------------------------------------------------------------------------------------
double ZCounting::pfIso(const pat::Muon &mu){
    return (mu.pfIsolationR04().sumChargedHadronPt + std::max(0., mu.pfIsolationR04().sumNeutralHadronEt
                                                                         + mu.pfIsolationR04().sumPhotonEt
                                                                         - 0.5 * mu.pfIsolationR04().sumPUPt)
                   ) / mu.pt();
}

//--------------------------------------------------------------------------------------------------
double ZCounting::dxy(const pat::Muon &mu, const reco::Vertex &vtx){
    return fabs(mu.muonBestTrack()->dxy(vtx.position()));
}

//--------------------------------------------------------------------------------------------------
double ZCounting::dz(const pat::Muon &mu, const reco::Vertex &vtx){
    return fabs(mu.muonBestTrack()->dz(vtx.position()));
}

//--------------------------------------------------------------------------------------------------
double ZCounting::pointsDistance(const reco::Candidate::Point &p1, const reco::Candidate::Point &p2){
    // computes the euclidean distance of two points
    return std::sqrt(std::pow(p1.x() - p2.x(),2) + std::pow(p1.y() - p2.y(),2) + std::pow(p1.z() - p2.z(),2));
}

//--------------------------------------------------------------------------------------------------
bool ZCounting::isPVClosestVertex(const std::vector<reco::Vertex> &vtxCol, const pat::Muon &mu){
    reco::Vertex vtx = *vtxCol.begin();
    double minDist = 999.;
    int nPV = 0;
    for (auto const& itVtx : vtxCol) {
        if(!isGoodPV(itVtx))
            continue;
        const float dist = pointsDistance(itVtx.position(), mu.vertex());
        if (dist < minDist){
            vtx = itVtx;
            minDist = dist;
            if(nPV > 0)
                return false;
        }
        nPV++;
    }
    return true;
}


//--------------------------------------------------------------------------------------------------
bool ZCounting::isTau(const reco::GenParticle* lepton)const
{
    return std::abs(lepton->pdgId()) == 15;
}

//--------------------------------------------------------------------------------------------------
const reco::GenParticle* ZCounting::tauDaughter(const reco::GenParticle* tau)const
{
    for(size_t iDaughter = 0; iDaughter < tau->numberOfDaughters(); ++iDaughter){
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(tau->daughter(iDaughter));
        if(std::abs(daughter->pdgId())==11 || std::abs(daughter->pdgId())==13) return daughter;
        else if(this->isTau(daughter)) return this->tauDaughter(daughter);
    }
    return tau;
}



//define this as a plug-in
DEFINE_FWK_MODULE(ZCounting);
