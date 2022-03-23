// -*- C++ -*-
//
// Package:    ZCountingAOD/ZCountAnalyze
// Class:      ZCountingAOD
//
/**\class ZCountingAOD ZCountingAOD.cc ZCountingAOD/ZCountAnalyze/plugins/ZCountingAOD.cc

 Description:
     Purpose is to derive MC corrections of Z counting efficiency in respect of true (MC) Z efficiency
     ZCountingAODAOD.cc is derived from ZCountingAOD.cc but adapted for AOD

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
#include <algorithm>
#include <iterator>

// CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// local framework includes
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"
#include "ZCounting/ZUtils/interface/triggertool.h"

#include "ZCounting/ZUtils/interface/Helper.h"
#include "ZCounting/ZUtils/interface/RoccoR.h"
#include "ZCounting/ZUtils/interface/getFilename.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"
#include <TLorentzVector.h>
//
// class declaration
//


class ZCountingAOD : //public edm::one::EDAnalyzer<edm::one::SharedResources>
    public edm::EDAnalyzer
{
public:
    explicit ZCountingAOD(const edm::ParameterSet&);
    ~ZCountingAOD();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void clearVariables();

    std::string get_triggerPath(std::string, const edm::TriggerNames&);

    int getMuonID(const reco::Muon&, const reco::Vertex&);
    
    bool isMuonTriggerObjEmulated(const double eta, const double phi, const long long unsigned eventNumber);


    // ----------member data ---------------------------
    const double MUON_MASS = 0.105658369;

    edm::Service<TFileService> fs;
    TTree *tree_;

    HLTConfigProvider hltConfigProvider_;
    triggertool *triggers;

    bool emulateTrigger_;
    TFile *_fileHLTEmulation = 0;
        
    // rocchester corrections
    RoccoR rc; 
    std::string roccorFile;   
    TRandom3 rand_;

    TLorentzVector vMuon;
    TLorentzVector vAntiMuon;

    // --- input
    const edm::InputTag triggerResultsInputTag_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> standaloneCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> pvCollection_;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZCollection_;
    edm::EDGetTokenT<TtGenEvent> genTtEvent_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoCollection_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;
    edm::EDGetTokenT<LHEEventProduct> lheEventInfo_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleCollection_;

    std::string era_;

    // Triggers
    bool hltChanged_;

    // MET trigger
    std::vector<std::string> metTriggerPatterns_;

    // Muon trigger
    std::vector<std::string> muonTriggerPatterns_;
    // max dR matching between muon and hlt object
    double DRMAX;

    // flags
    bool isData_;
    bool hasGenZ_;
    bool hasGenTt_;


    // --- output

    // ... for event info
    unsigned int runNumber_;
    unsigned int lumiBlock_;
    unsigned int eventNumber_;

    std::vector<int> v_genWeightIDs_;
    std::vector<int> v_pdfWeightIDs_;

    float eventweight_;
    std::vector<float> v_psWeight_;
    std::vector<float> v_meWeight_;
    std::vector<float> v_pdfWeight_;

    int nPU_;
    int nPV_;
    int decayMode_;
    float z_genMass_;
    float z_recoMass_;

    int met_triggerBits_;

    int muon_genRecoMatches_;
    int muon_genRecoObj_;
    int muon_genRecoTrackMatches_;
    int muon_genRecoTrackObj_;
    float muon_genPt_;
    float muon_genEta_;
    float muon_genPhi_;

    float muon_genVx_;
    float muon_genVy_;
    float muon_genVz_;
    // float muon_genDxyError_;
    // float muon_genDzError_;
    
    int antiMuon_genRecoMatches_;
    int antiMuon_genRecoObj_;
    int antiMuon_genRecoTrackMatches_;
    int antiMuon_genRecoTrackObj_;
    float antiMuon_genPt_;
    float antiMuon_genEta_;
    float antiMuon_genPhi_;

    float antiMuon_genVx_;
    float antiMuon_genVy_;
    float antiMuon_genVz_;
    // float antiMuon_genDxyError_;
    // float antiMuon_genDzError_;
    
    unsigned int nMuon_;
    std::vector<float> muon_pt_;
    std::vector<float> muon_eta_;
    std::vector<float> muon_phi_;
    std::vector<float> muon_ptTrk_;
    std::vector<float> muon_etaTrk_;
    std::vector<float> muon_phiTrk_;
    std::vector<float> muon_ptSta_;
    std::vector<float> muon_etaSta_;
    std::vector<float> muon_phiSta_;
    std::vector<int> muon_charge_;
    // reference point position
    std::vector<float> muon_dx_;
    std::vector<float> muon_dy_;
    std::vector<float> muon_dz_;
    // vertex position
    std::vector<float> muon_vx_;
    std::vector<float> muon_vy_;
    std::vector<float> muon_vz_;
    // uncertainties on position
    std::vector<float> muon_DxyError_;
    std::vector<float> muon_DzError_;
    
    std::vector<int> muon_matchValue_;
    std::vector<float> muon_tkRelIso_;
    std::vector<float> muon_pfRelIso04_all_;
    std::vector<int> muon_ID_;
    std::vector<int> muon_triggerBits_;

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
    std::vector<float> muon_trackAlgo_;

    std::vector<float> muon_ScaleCorr_;
    //std::vector<float> muon_ScaleCorr_stat_RMS_;
    //std::vector<float> muon_ScaleCorr_Zpt_;
    //std::vector<float> muon_ScaleCorr_Ewk_;
    //std::vector<float> muon_ScaleCorr_deltaM_;
    //std::vector<float> muon_ScaleCorr_Ewk2_;
    //std::vector<float> muon_ScaleCorr_Total_;

    unsigned int nStandalone_;
    std::vector<float> sta_pt_;
    std::vector<float> sta_eta_;
    std::vector<float> sta_phi_;
    std::vector<int> sta_charge_;
    
    std::vector<float> sta_dx_;
    std::vector<float> sta_dy_;
    std::vector<float> sta_dz_;
    std::vector<float> sta_DxyError_;
    std::vector<float> sta_DzError_;
    std::vector<float> sta_vx_;
    std::vector<float> sta_vy_;
    std::vector<float> sta_vz_;

    std::vector<float> sta_trackAlgo_;

    std::vector<float> sta_ScaleCorr_;

    unsigned int nTrack_;
    std::vector<float> track_pt_;
    std::vector<float> track_eta_;
    std::vector<float> track_phi_;
    std::vector<int> track_charge_;
    
    std::vector<float> track_dx_;
    std::vector<float> track_dy_;
    std::vector<float> track_dz_;
    std::vector<float> track_DxyError_;
    std::vector<float> track_DzError_;
    std::vector<float> track_vx_;
    std::vector<float> track_vy_;
    std::vector<float> track_vz_;

    std::vector<int> track_nPixelHits_;
    std::vector<int> track_nTrackerLayers_;
    std::vector<float> track_validFraction_;
    std::vector<float> track_trackAlgo_;

    std::vector<float> track_ScaleCorr_;
};

//
// constructors and destructor
//
ZCountingAOD::ZCountingAOD(const edm::ParameterSet& iConfig):
    triggerResultsInputTag_(iConfig.getParameter<edm::InputTag>("TriggerResults")),
    muonCollection_  (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("reco_muons"))),
    standaloneCollection_  (consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("reco_standalones"))),
    trackCollection_  (consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("reco_tracks"))),
    pvCollection_    (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("edmPVName"))),
    genZCollection_  (consumes<std::vector<GenZDecayProperties> > (iConfig.getParameter<edm::InputTag>("genZLeptonCollection"))),
    genTtEvent_  (consumes<TtGenEvent> (iConfig.getParameter<edm::InputTag>("genTtCollection"))),
    pileupInfoCollection_  (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfoCollection"))),
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfo_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    genParticleCollection_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles")))
{
    LogDebug("ZCountingAOD")<<"ZCountingAOD(...)";

    era_ = iConfig.getParameter<std::string>("era");

    isData_ = iConfig.getUntrackedParameter<bool>("isData");
    hasGenZ_ = iConfig.getUntrackedParameter<bool>("hasGenZ");
    hasGenTt_ = iConfig.getUntrackedParameter<bool>("hasGenTt");

    v_genWeightIDs_ = iConfig.getParameter<std::vector<int>>("genWeights");
    v_pdfWeightIDs_ = iConfig.getParameter<std::vector<int>>("pdfWeights");

    roccorFile = iConfig.getParameter<std::string>("roccorFile");

    hltChanged_ = true;
    emulateTrigger_ = iConfig.getUntrackedParameter<bool>("emulateTrigger");

    muonTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("muon_trigger_patterns");
    metTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("met_trigger_patterns");
    
    DRMAX = iConfig.getUntrackedParameter<double>("muon_trigger_DRMAX");


    triggers = new triggertool();
    triggers->setTriggerResultsToken(consumes<edm::TriggerResults>(triggerResultsInputTag_));
    triggers->setTriggerEventToken(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent")));
    triggers->setDRMAX(DRMAX);

    edm::LogVerbatim("ZCountingAOD") << "getInput: set trigger names";
    for(unsigned int i = 0; i < muonTriggerPatterns_.size(); ++i) {
        triggers->addTriggerRecord(muonTriggerPatterns_.at(i));
    }
    for(unsigned int i = 0; i < metTriggerPatterns_.size(); ++i) {
        triggers->addTriggerRecord(metTriggerPatterns_.at(i));
    }

}


ZCountingAOD::~ZCountingAOD()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZCountingAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    LogDebug("ZCountingAOD::analyze");

    this->clearVariables();

    triggers->readEvent(iEvent);
    
    // take data events only if one of the required triggers has fired
    if(iEvent.isRealData() && !(emulateTrigger_ || triggers->pass(metTriggerPatterns_) || triggers->pass(muonTriggerPatterns_))){
        return;
    }

    // Event info
    eventNumber_ = iEvent.id().event();

    const reco::GenParticle* genLepton = 0;
    const reco::GenParticle* genAntiLepton = 0;
    if(!iEvent.isRealData()){

        edm::Handle<GenEventInfoProduct> evt_info;
        iEvent.getByToken(genEventInfo_, evt_info);

        // nominal ME weights
        eventweight_ = evt_info->weight();

        edm::Handle<LHEEventProduct> lhe_info;
        try {iEvent.getByToken(lheEventInfo_, lhe_info);}
        catch (...) {;}
    	if(lhe_info.isValid()){
            // ME variations
            int pdf_begin=-2;
            int pdf_end=-2;
            for(size_t iwgt=0; iwgt<lhe_info->weights().size(); ++iwgt){
        		const LHEEventProduct::WGT& wgt = lhe_info->weights().at(iwgt);
                // check if weight is in the list of meWeightIDs
                bool exists = std::find(std::begin(v_genWeightIDs_), std::end(v_genWeightIDs_), std::stoi(wgt.id)) != std::end(v_genWeightIDs_);
                if(exists){
        		    v_meWeight_.push_back((float) (wgt.wgt/lhe_info->originalXWGTUP()));
        		}

                // store indices of pdf weights
                if(wgt.id == v_pdfWeightIDs_.at(0)){
                    pdf_begin = iwgt;
                }
                if(wgt.id == v_pdfWeightIDs_.at(1)){
                    pdf_end = iwgt;
                }
            }
            // store PDF weights
            if (pdf_begin > -2 && pdf_end > -2 && pdf_end > pdf_begin) {
        	    unsigned nPDFvars = pdf_end - pdf_begin;
                for (unsigned i=0; i <= nPDFvars; ++i) {
            		const LHEEventProduct::WGT& wgt = lhe_info->weights().at(pdf_begin+i);
            		v_pdfWeight_.push_back((float) wgt.wgt/lhe_info->originalXWGTUP());
        	    }
        	}
        }

        // --- PS weights
        std::vector<double> psWeights;
        psWeights = evt_info->weights();
        if(psWeights.size()>0) {
            const double nominal = psWeights.at(0);
            if(nominal == 0.0){
                throw cms::Exception("LogicError") << "@@@ ZCountingAOD::analyze -- "
                    << "nominal PS-Weight equal to zero (failed to normalize other PS-Weights)";
            }
            // store only the default set of variations (2 and 0.5)
            for(unsigned int i=6; i<10; ++i){
                v_psWeight_.push_back((float) (psWeights.at(i) / nominal));
            }
        }

        if(hasGenZ_){
            edm::Handle<std::vector<GenZDecayProperties> > genZCollection;
            iEvent.getByToken(genZCollection_, genZCollection);
            decayMode_ = (genZCollection->size() != 0 ? genZCollection->at(0).decayMode() : 0);
            if(genZCollection->size() > 1)
                decayMode_+= genZCollection->at(1).decayMode()*100000;
            genLepton = genZCollection->at(0).stableLepton();
            genAntiLepton = genZCollection->at(0).stableAntiLepton();

        }
        if(hasGenTt_){
            LogDebug("ZCountingAOD::analyze")<<"hasGenZ";

            edm::Handle<TtGenEvent> genTtEvent;
            iEvent.getByToken(genTtEvent_, genTtEvent);
            if(genTtEvent->lepton()){
                LogDebug("ZCountingAOD::analyze")<<"hasLepton";
                genLepton = genTtEvent->lepton();
                if(isTau(genLepton)){
                    decayMode_ += 150000;
                    genLepton = tauDaughter(genTtEvent->lepton());
                }

            }
            if(genTtEvent->leptonBar()){
                LogDebug("ZCountingAOD::analyze")<<"hasAntiLepton";
                genAntiLepton = genTtEvent->leptonBar();
                if(isTau(genAntiLepton)){
                    decayMode_ += 150000;
                    genAntiLepton = tauDaughter(genTtEvent->leptonBar());
                }
            }
            LogDebug("ZCountingAOD::analyze")<<"setDecayMode";
            if(genLepton)
                decayMode_ += 100*std::abs(genLepton->pdgId());
            if(genAntiLepton)
                decayMode_ += std::abs(genAntiLepton->pdgId());
        }
    }

    // --- PV selection
    edm::Handle<std::vector<reco::Vertex> > pvCollection;
    iEvent.getByToken(pvCollection_, pvCollection);
    const std::vector<reco::Vertex> *pvCol = pvCollection.product();
    reco::Vertex pv = *pvCol->begin();
    nPV_ = 0;
    for (const reco::Vertex &itVtx : *pvCol) {
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
        
        muon_genVx_ = genLepton->vx();
        muon_genVy_ = genLepton->vy();
        muon_genVz_ = genLepton->vz();
        // muon_genDxyError_ = genLepton->dxyError();
        // muon_genDzError_ = genLepton->dzError();

    }
    if(genAntiLepton){
        antiMuon_genPt_ = genAntiLepton->pt();
        antiMuon_genEta_ = genAntiLepton->eta();
        antiMuon_genPhi_ = genAntiLepton->phi();
        
        antiMuon_genVx_ = genAntiLepton->vx();
        antiMuon_genVy_ = genAntiLepton->vy();
        antiMuon_genVz_ = genAntiLepton->vz();
        // antiMuon_genDxyError_ = genAntiLepton->dxyError();
        // antiMuon_genDzError_ = genAntiLepton->dzError();        
        
    }
    if(genLepton && genAntiLepton){
        z_genMass_ = (genLepton->p4() + genAntiLepton->p4()).M();
    }

    for(unsigned j = 0; j < metTriggerPatterns_.size(); ++j){
        met_triggerBits_ += std::pow(2,j) * triggers->pass(metTriggerPatterns_.at(j));
    }

    // --- store all reco muons
    LogDebug("ZCountingAOD") << "find reco muons";
    edm::Handle<std::vector<reco::Muon>> muonCollection;
    iEvent.getByToken(muonCollection_, muonCollection);
    for (const reco::Muon &mu : *muonCollection){
        // reject muons if best track, outer track, and inner track fail acceptance cuts
        if (std::abs(mu.eta()) > 2.4
            && (mu.outerTrack().isNull() || mu.outerTrack()->eta() > 2.4)
            && (mu.innerTrack().isNull() || mu.innerTrack()->eta() > 2.4)        
        ) {
            continue;
        }
        if (mu.pt() < 20 
            && (mu.outerTrack().isNull() || mu.outerTrack()->pt() < 20)
            && (mu.innerTrack().isNull() || mu.innerTrack()->pt() < 20)
        ) {
            continue;
        }

        if (mu.outerTrack().isNonnull() && mu.outerTrack()->extra().isNonnull() && !(mu.isGlobalMuon())){
            // Check if muon is twice in this collection, once as standAloneMuon and once as globalMuon
            //    We want to keep the globalMuon doublicate, hence we only consider muons that are not global
            bool isDoublicate = false;
            for (const reco::Muon &mu2 : *muonCollection) {
                if ((&mu2 == &mu) || mu2.outerTrack().isNull() || mu2.outerTrack()->extra().isNull())
                    continue;
                
                if (mu2.outerTrack()->extra().get() == mu.outerTrack()->extra().get()) {
                    isDoublicate = true;
                    break;
                }
            }            
            if(isDoublicate)
                continue;
        }


        int match_value_ = -1;
        if(!iEvent.isRealData()){
            // Lepton truth matching
            edm::Handle<reco::GenParticleCollection> genParticleCollection;
            iEvent.getByToken(genParticleCollection_, genParticleCollection);

            float pt_gen = -999.;
            float eta_gen = -999.;
            float phi_gen = -999.;
            int id_gen = -999;
            int isPromptFinalState_gen = -999;
            int isDirectPromptTauDecayProductFinalState_gen = -999;
            match_value_ = doMatch(genParticleCollection, mu.pt(), mu.eta(), mu.phi(), -13*mu.charge(),
                pt_gen, eta_gen, phi_gen, id_gen,
                isPromptFinalState_gen, isDirectPromptTauDecayProductFinalState_gen);
        }
        muon_matchValue_.push_back(match_value_);
                
        // Look for the distance to (0,0,0)
        muon_dx_.push_back(mu.muonBestTrack()->referencePoint().x());
        muon_dy_.push_back(mu.muonBestTrack()->referencePoint().y());
        muon_dz_.push_back(mu.muonBestTrack()->referencePoint().z());
        // vertex position 
        muon_vx_.push_back(mu.vx());
        muon_vy_.push_back(mu.vy());
        muon_vz_.push_back(mu.vz());
        // uncertainty
        muon_DxyError_.push_back(mu.dxyError());
        muon_DzError_.push_back(mu.dzError());

        muon_pt_.push_back(mu.pt());
        muon_eta_.push_back(mu.eta());
        muon_phi_.push_back(mu.phi());        
        muon_charge_.push_back(mu.charge());
        muon_tkRelIso_.push_back(getTkIso(mu));
        muon_pfRelIso04_all_.push_back(getPFIso(mu));
        muon_ID_.push_back(getMuonID(mu, pv));

        muon_isMedium_.push_back(muon::isMediumMuon(mu));
        muon_isStandalone_.push_back(mu.isStandAloneMuon());
        muon_isTracker_.push_back(mu.isTrackerMuon());
        muon_isGlobal_.push_back(mu.isGlobalMuon());
        muon_isPFCand_.push_back(mu.isPFMuon());
        muon_nStations_.push_back(mu.numberOfMatchedStations());
        muon_SegmentCompatibility_.push_back(muon::segmentCompatibility(mu));
        muon_chi2LocalPosition_.push_back(mu.combinedQuality().chi2LocalPosition);
        muon_trkKink_.push_back(mu.combinedQuality().trkKink);

        if(mu.innerTrack().isNonnull()){
            muon_ptTrk_.push_back(mu.innerTrack()->pt());
            muon_etaTrk_.push_back(mu.innerTrack()->eta());
            muon_phiTrk_.push_back(mu.innerTrack()->phi());      
            
            muon_nPixelHits_.push_back(mu.innerTrack()->hitPattern().numberOfValidPixelHits());
            muon_nTrackerLayers_.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
            muon_validFraction_.push_back(mu.innerTrack()->validFraction());
            muon_trackAlgo_.push_back(mu.innerTrack()->originalAlgo());
        }
        else{
            muon_ptTrk_.push_back(-1);
            muon_etaTrk_.push_back(-1);
            muon_phiTrk_.push_back(-1);    
            
            muon_nPixelHits_.push_back(-1);
            muon_nTrackerLayers_.push_back(-1);
            muon_validFraction_.push_back(-1);
            muon_trackAlgo_.push_back(-1);
        }
        if(mu.outerTrack().isNonnull()){
            muon_ptSta_.push_back(mu.outerTrack()->pt());
            muon_etaSta_.push_back(mu.outerTrack()->eta());
            muon_phiSta_.push_back(mu.outerTrack()->phi());            
        }
        else{            
            muon_ptSta_.push_back(-1);
            muon_etaSta_.push_back(-1);
            muon_phiSta_.push_back(-1);                 
        }


        if(mu.isGlobalMuon()){
            muon_normChi2_.push_back(mu.globalTrack()->normalizedChi2());
            muon_nTrackHits_.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
        }
        else{
            muon_normChi2_.push_back(-1);
            muon_nTrackHits_.push_back(-1);
        }
        
        int bits_ = 0;    
        if(emulateTrigger_){
            bits_ = isMuonTriggerObjEmulated(mu.eta(), mu.phi(), eventNumber_);
        }
        else{
            for(unsigned j = 0; j < muonTriggerPatterns_.size(); ++j){
                bits_ += std::pow(2,j) * triggers->passObj(muonTriggerPatterns_.at(j), mu.eta(), mu.phi());
            }
        }
        muon_triggerBits_.push_back(bits_);            
    


        double roccorSF = 1.; // Rochester correction
        if(!iEvent.isRealData()){

            if(genLepton && mu.pdgId() == 13 && reco::deltaR(mu.eta(), mu.phi(), genLepton->eta(), genLepton->phi()) < 0.03){
                muon_genRecoMatches_++;
                if(muon_genRecoObj_ == -1
                    || std::abs(mu.pt() - muon_genPt_) < std::abs(muon_pt_[muon_genRecoObj_] - muon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    muon_genRecoObj_ = nMuon_;
                    roccorSF = rc.kSpreadMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), muon_genPt_);
                }
            }

            else if(genAntiLepton && mu.pdgId() == -13 && reco::deltaR(mu.eta(), mu.phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                antiMuon_genRecoMatches_++;
                if(antiMuon_genRecoObj_ == -1
                    || std::abs(mu.pt() - antiMuon_genPt_) < std::abs(muon_pt_[antiMuon_genRecoObj_] - antiMuon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    antiMuon_genRecoObj_ = nMuon_;
                    roccorSF = rc.kSpreadMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), antiMuon_genPt_);
                }
            }
            else{
                const int nl = mu.innerTrack().isNonnull() ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0;
                roccorSF = rc.kSmearMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), nl, rand_.Rndm());
            }

            // store muon Rochester corrections and uncertainties
            
            
            // muon_ScaleCorr_.push_back(mu.hasUserFloat("MuonEnergyCorr")                   ? mu.userFloat("MuonEnergyCorr")          : 1.0);
            // muon_ScaleCorr_stat_RMS_.push_back(mu.hasUserFloat("MuonEnergyCorr_stat_RMS") ? mu.userFloat("MuonEnergyCorr_stat_RMS") : 0.0);
            // muon_ScaleCorr_Zpt_.push_back(mu.hasUserFloat("MuonEnergyCorr_Zpt")           ? mu.userFloat("MuonEnergyCorr_Zpt")      : 0.0);
            // muon_ScaleCorr_Ewk_.push_back(mu.hasUserFloat("MuonEnergyCorr_Ewk")           ? mu.userFloat("MuonEnergyCorr_Ewk")      : 0.0);
            // muon_ScaleCorr_deltaM_.push_back(mu.hasUserFloat("MuonEnergyCorr_deltaM")     ? mu.userFloat("MuonEnergyCorr_deltaM")   : 0.0);
            // muon_ScaleCorr_Ewk2_.push_back(mu.hasUserFloat("MuonEnergyCorr_Ewk2")         ? mu.userFloat("MuonEnergyCorr_Ewk2")     : 0.0);
            // muon_ScaleCorr_Total_.push_back(mu.hasUserFloat("MuonEnergyCorr_Total")       ? mu.userFloat("MuonEnergyCorr_Total")    : 0.0);
        }
        else{
            roccorSF = rc.kScaleDT(mu.charge(), mu.pt(), mu.eta(), mu.phi());
        }
        
        // Rochester corrections 
        muon_ScaleCorr_.push_back(roccorSF);

        nMuon_++;
    }
    
    // --- store all reco tracks
    edm::Handle<std::vector<reco::Track>> trackCollection;
    iEvent.getByToken(trackCollection_, trackCollection);
    for (const reco::Track &trk : *trackCollection){        
        if(std::abs(trk.eta()) > 2.4) continue;
        if(trk.pt() < 20) continue;

        // Check if track is not a muon
        bool isMuon = false;
        for (const reco::Muon &mu : *muonCollection) {
            if (mu.innerTrack().isNonnull() && mu.innerTrack().get() == &trk) {
                isMuon = true;
                break;
            }
        }
        if(isMuon)
            continue;
        
        track_pt_.push_back(trk.pt());
        track_eta_.push_back(trk.eta());
        track_phi_.push_back(trk.phi());

        track_charge_.push_back(trk.charge());
        track_dx_.push_back(trk.referencePoint().x());
        track_dy_.push_back(trk.referencePoint().y());
        track_dz_.push_back(trk.referencePoint().z());

        // vertex position 
        track_vx_.push_back(trk.vx());
        track_vy_.push_back(trk.vy());
        track_vz_.push_back(trk.vz());
        // uncertainty
        track_DxyError_.push_back(trk.dxyError());
        track_DzError_.push_back(trk.dzError());
        
        track_nPixelHits_.push_back(trk.hitPattern().numberOfValidPixelHits());
        track_nTrackerLayers_.push_back(trk.hitPattern().trackerLayersWithMeasurement());
        track_validFraction_.push_back(trk.validFraction());
        track_trackAlgo_.push_back(trk.originalAlgo());

        double roccorSF = 1.; // Rochester correction
        if(!iEvent.isRealData()){

            if(genLepton && trk.charge() > 0 && reco::deltaR(trk.eta(), trk.phi(), genLepton->eta(), genLepton->phi()) < 0.03){
                muon_genRecoTrackMatches_++;
                if(muon_genRecoTrackObj_ == -1
                    || std::abs(trk.pt() - muon_genPt_) < std::abs(muon_pt_[muon_genRecoTrackObj_] - muon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    muon_genRecoTrackObj_ = nTrack_;
                    roccorSF = rc.kSpreadMC(trk.charge(), trk.pt(), trk.eta(), trk.phi(), muon_genPt_);
                }
            }

            else if(genAntiLepton && trk.charge() < 0 && reco::deltaR(trk.eta(), trk.phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                antiMuon_genRecoTrackMatches_++;
                if(antiMuon_genRecoTrackObj_ == -1
                    || std::abs(trk.pt() - antiMuon_genPt_) < std::abs(muon_pt_[antiMuon_genRecoTrackObj_] - antiMuon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    antiMuon_genRecoTrackObj_ = nTrack_;
                    roccorSF = rc.kSpreadMC(trk.charge(), trk.pt(), trk.eta(), trk.phi(), antiMuon_genPt_);
                }
            }
            else{
                roccorSF = rc.kSmearMC(trk.charge(), trk.pt(), trk.eta(), trk.phi(), trk.hitPattern().trackerLayersWithMeasurement(), rand_.Rndm());
            }
        }
        else{
            roccorSF = rc.kScaleDT(trk.charge(), trk.pt(), trk.eta(), trk.phi());
        }

        // Rochester corrections 
        track_ScaleCorr_.push_back(roccorSF);            
        nTrack_++;
    }    

    // --- store all standalone tracks
    edm::Handle<std::vector<reco::Track>> standaloneCollection;
    iEvent.getByToken(standaloneCollection_, standaloneCollection);
    for (const reco::Track &sta : *standaloneCollection){        
        if(std::abs(sta.eta()) > 2.4) continue;
        if(sta.pt() < 20) continue;

        // Check if standalone muon is not a muon
        bool isMuon = false;
        for (const reco::Muon &mu : *muonCollection) {
            if (mu.outerTrack().isNull())
                continue;
            if (mu.outerTrack().get() == &sta) {
                isMuon = true;
                break;
            }
            if (mu.outerTrack()->extra().isNonnull() && mu.outerTrack()->extra().get() == sta.extra().get()) {
                isMuon = true;
                break;                
            }
        }
        if(isMuon)
            continue;
        
        sta_pt_.push_back(sta.pt());
        sta_eta_.push_back(sta.eta());
        sta_phi_.push_back(sta.phi());

        sta_charge_.push_back(sta.charge());
        sta_dx_.push_back(sta.referencePoint().x());
        sta_dy_.push_back(sta.referencePoint().y());
        sta_dz_.push_back(sta.referencePoint().z());

        // vertex position 
        sta_vx_.push_back(sta.vx());
        sta_vy_.push_back(sta.vy());
        sta_vz_.push_back(sta.vz());
        // uncertainty
        sta_DxyError_.push_back(sta.dxyError());
        sta_DzError_.push_back(sta.dzError());
        
        sta_trackAlgo_.push_back(sta.originalAlgo());

        double roccorSF = 1.; // Rochester correction
        if(!iEvent.isRealData()){

            if(genLepton && sta.charge() > 0 && reco::deltaR(sta.eta(), sta.phi(), genLepton->eta(), genLepton->phi()) < 0.03){
                muon_genRecoTrackMatches_++;
                if(muon_genRecoTrackObj_ == -1
                    || std::abs(sta.pt() - muon_genPt_) < std::abs(muon_pt_[muon_genRecoTrackObj_] - muon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    muon_genRecoTrackObj_ = nStandalone_;
                    roccorSF = rc.kSpreadMC(sta.charge(), sta.pt(), sta.eta(), sta.phi(), muon_genPt_);
                }
            }

            else if(genAntiLepton && sta.charge() < 0 && reco::deltaR(sta.eta(), sta.phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                antiMuon_genRecoTrackMatches_++;
                if(antiMuon_genRecoTrackObj_ == -1
                    || std::abs(sta.pt() - antiMuon_genPt_) < std::abs(muon_pt_[antiMuon_genRecoTrackObj_] - antiMuon_genPt_)
                ){
                    // store index of reco match for the one with the closest pt in case of ambiguity
                    antiMuon_genRecoTrackObj_ = nStandalone_;
                    roccorSF = rc.kSpreadMC(sta.charge(), sta.pt(), sta.eta(), sta.phi(), antiMuon_genPt_);
                }
            }
            else{
                roccorSF = rc.kSmearMC(sta.charge(), sta.pt(), sta.eta(), sta.phi(), sta.hitPattern().trackerLayersWithMeasurement(), rand_.Rndm());
            }
        }
        else{
            roccorSF = rc.kScaleDT(sta.charge(), sta.pt(), sta.eta(), sta.phi());
        }

        // Rochester corrections 
        sta_ScaleCorr_.push_back(roccorSF);            
        nStandalone_++;
    }

    if(!iEvent.isRealData()){
        edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoCollection;
        iEvent.getByToken(pileupInfoCollection_, pileupInfoCollection);

        for(std::vector<PileupSummaryInfo>::const_iterator puI = pileupInfoCollection->begin(); puI != pileupInfoCollection->end(); ++puI)
        {
            if(puI->getBunchCrossing() == 0){
                nPU_ = puI->getTrueNumInteractions();
            }
        }

        if(antiMuon_genRecoObj_ != -1 && muon_genRecoObj_ != -1){
            vMuon.SetPtEtaPhiM(muon_pt_[muon_genRecoObj_], muon_eta_[muon_genRecoObj_], muon_phi_[muon_genRecoObj_], MUON_MASS);
            vAntiMuon.SetPtEtaPhiM(muon_pt_[antiMuon_genRecoObj_], muon_eta_[antiMuon_genRecoObj_], muon_phi_[antiMuon_genRecoObj_], MUON_MASS);

            z_recoMass_ = (vMuon + vAntiMuon).M();
        }
    }

    tree_->Fill();
}

// ------------ method called once each new LS  ------------
void
ZCountingAOD::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
    lumiBlock_ = iLumi.id().luminosityBlock();

    std::cout<<"ZCountingAOD::beginLuminosityBlock --- now at LS "<< lumiBlock_ <<std::endl;

    if(emulateTrigger_){
        // find and open file with emulated HLT information
        const std::string fNameHLT = getFilename(runNumber_, lumiBlock_);

        if(_fileHLTEmulation != 0)
            _fileHLTEmulation->Close();

        _fileHLTEmulation = TFile::Open(("root://xrootd-cms.infn.it//"+fNameHLT).c_str());        
        fs->cd();
    }
}

// ------------ method called once each run just before starting event loop  ------------
void ZCountingAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    // edm::LogVerbatim("ZCountingAOD") << "now at "<<iRun.id();
    std::cout<< "new run"<<std::endl;
    runNumber_ = iRun.id().run();
    hltChanged_ = true;

    edm::LogVerbatim("ZCountingAOD") << "now at "<<iRun.id();
    
    if (hltConfigProvider_.init(iRun, iSetup, triggerResultsInputTag_.process(), hltChanged_)) {
        edm::LogVerbatim("ZCountingAOD")<<" [TriggerObjMatchValueMapsProducer::beginRun] HLTConfigProvider initialized [processName() = \""
            << hltConfigProvider_.processName() << "\", tableName() = \"" << hltConfigProvider_.tableName()
            << "\", size() = " << hltConfigProvider_.size() << "]";
    } else {
        edm::LogError("ZCountingAOD") << "Initialization of HLTConfigProvider failed for Run=" << iRun.id() << " (process=\""
        << triggerResultsInputTag_.process() << "\") -> plugin will not produce outputs for this Run";
        return;
    }
    edm::LogVerbatim("ZCountingAOD") << "hlt: "<<hltChanged_;

    triggers->initHLTObjects(hltConfigProvider_);

}

// ------------ method called once each job just before starting event loop  ------------
void
ZCountingAOD::beginJob()
{
    LogDebug("ZCountingAOD")<<"beginJob()";

    rc.init(edm::FileInPath(roccorFile).fullPath());
    rand_ = TRandom3();
    rand_.SetSeed(1);
    
    if( !fs ){
        edm::LogError("ZCountingAOD") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    // Event info
    tree_->Branch("runNumber", &runNumber_, "runNumber/i");
    tree_->Branch("lumiBlock", &lumiBlock_,"lumiBlock/i");
    tree_->Branch("eventNumber", &eventNumber_, "eventNumber/i");

    tree_->Branch("nPV", &nPV_,"nPV_/i");
    if(!isData_){
        tree_->Branch("nPU", &nPU_,"nPU_/i");
    }

    // weights
    // ME weight
    if(!isData_){
        tree_->Branch("eventweight", &eventweight_, "eventweight_/f");
        tree_->Branch("MEWeight",  &v_meWeight_);
        tree_->Branch("PDFWeight", &v_pdfWeight_);

        // PS weights
        tree_->Branch("PSWeight", &v_psWeight_);
    }


    if(!isData_){
        // gen level info
        tree_->Branch("decayMode", &decayMode_, "decayMode_/i");
        tree_->Branch("z_genMass", &z_genMass_,"z_genMass_/f");
        tree_->Branch("z_recoMass", &z_recoMass_,"z_recoMass_/f");

        tree_->Branch("muon_genPt", &muon_genPt_,"muon_genPt_/f");
        tree_->Branch("muon_genEta", &muon_genEta_,"muon_genEta_/f");
        tree_->Branch("muon_genPhi", &muon_genPhi_,"muon_genPhi_/f");
        
        tree_->Branch("muon_vx", &muon_genVx_,"muon_genVx_/f");
        tree_->Branch("muon_vy", &muon_genVy_,"muon_genVy_/f");
        tree_->Branch("muon_vz", &muon_genVz_,"muon_genVz_/f");
        // tree_->Branch("muon_DxyError", &muon_genDxyError_,"muon_genDxyError_/f");
        // tree_->Branch("muon_DzError", &muon_genDzError_,"muon_genDzError_/f"); 
        
        tree_->Branch("muon_genRecoMatches", &muon_genRecoMatches_,"muon_genRecoMatches_/i");
        tree_->Branch("muon_genRecoTrackMatches", &muon_genRecoTrackMatches_,"muon_genRecoTrackMatches_/i");
        tree_->Branch("muon_genRecoObj", &muon_genRecoObj_,"muon_genRecoObj_/I");
        tree_->Branch("muon_genRecoTrackObj", &muon_genRecoTrackObj_,"muon_genRecoTrackObj_/I");

        tree_->Branch("antiMuon_genPt", &antiMuon_genPt_,"antiMuon_genPt_/f");
        tree_->Branch("antiMuon_genEta", &antiMuon_genEta_,"antiMuon_genEta_/f");
        tree_->Branch("antiMuon_genPhi", &antiMuon_genPhi_,"antiMuon_genPhi_/f");
    
        tree_->Branch("antiMuon_vx", &antiMuon_genVx_,"antiMuon_genVx_/f");
        tree_->Branch("antiMuon_vy", &antiMuon_genVy_,"antiMuon_genVy_/f");
        tree_->Branch("antiMuon_vz", &antiMuon_genVz_,"antiMuon_genVz_/f");        
        // tree_->Branch("antiMuon_DxyError", &antiMuon_genDxyError_,"antiMuon_genDxyError_/f");
        // tree_->Branch("antiMuon_DzError", &antiMuon_genDzError_,"antiMuon_genDzError_/f"); 
        
        tree_->Branch("antiMuon_genRecoMatches", &antiMuon_genRecoMatches_,"antiMuon_genRecoMatches_/i");
        tree_->Branch("antiMuon_genRecoTrackMatches", &antiMuon_genRecoTrackMatches_,"antiMuon_genRecoTrackMatches_/i");
        tree_->Branch("antiMuon_genRecoObj", &antiMuon_genRecoObj_,"antiMuon_genRecoObj_/I");
        tree_->Branch("antiMuon_genRecoTrackObj", &antiMuon_genRecoTrackObj_,"antiMuon_genRecoTrackObj_/I");
    }

    // reco level info
    tree_->Branch("MET_triggerBits", &met_triggerBits_);

    // muons
    tree_->Branch("nMuon", &nMuon_,"nMuon_/s");
    tree_->Branch("Muon_pt", &muon_pt_);
    tree_->Branch("Muon_eta", &muon_eta_);
    tree_->Branch("Muon_phi", &muon_phi_);
    tree_->Branch("Muon_ptTrk", &muon_ptTrk_);
    tree_->Branch("Muon_etaTrk", &muon_etaTrk_);
    tree_->Branch("Muon_phiTrk", &muon_phiTrk_);
    tree_->Branch("Muon_ptSta", &muon_ptSta_);
    tree_->Branch("Muon_etaSta", &muon_etaSta_);
    tree_->Branch("Muon_phiSta", &muon_phiSta_);
    tree_->Branch("Muon_charge", &muon_charge_);
    
    tree_->Branch("Muon_dx", &muon_dx_);
    tree_->Branch("Muon_dy", &muon_dy_);
    tree_->Branch("Muon_dz", &muon_dz_);
    tree_->Branch("Muon_vx", &muon_vx_);
    tree_->Branch("Muon_vy", &muon_vy_);
    tree_->Branch("Muon_vz", &muon_vz_);
    tree_->Branch("Muon_DxyError", &muon_DxyError_);
    tree_->Branch("Muon_DzError", &muon_DzError_);
    
    tree_->Branch("Muon_matchValue", &muon_matchValue_);
    tree_->Branch("Muon_ID", &muon_ID_);
    tree_->Branch("Muon_tkRelIso", &muon_tkRelIso_);
    tree_->Branch("Muon_pfRelIso04_all", &muon_pfRelIso04_all_);
    tree_->Branch("Muon_triggerBits", &muon_triggerBits_);

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
    tree_->Branch("Muon_trackAlgo", &muon_trackAlgo_);

    
    //Roccester corrections
    tree_->Branch("Muon_ScaleCorr",          &muon_ScaleCorr_);
    // tree_->Branch("Muon_ScaleCorr_stat_RMS", &muon_ScaleCorr_stat_RMS_);
    // tree_->Branch("Muon_ScaleCorr_Zpt",      &muon_ScaleCorr_Zpt_);
    // tree_->Branch("Muon_ScaleCorr_Ewk",      &muon_ScaleCorr_Ewk_);
    // tree_->Branch("Muon_ScaleCorr_deltaM",   &muon_ScaleCorr_deltaM_);
    // tree_->Branch("Muon_ScaleCorr_Ewk2",     &muon_ScaleCorr_Ewk2_);
    // tree_->Branch("Muon_ScaleCorr_Total",    &muon_ScaleCorr_Total_);

    // standalone muons
    tree_->Branch("nStandalone", &nStandalone_,"nStandalone_/s");
    tree_->Branch("Sta_pt", &sta_pt_);
    tree_->Branch("Sta_eta", &sta_eta_);
    tree_->Branch("Sta_phi", &sta_phi_);
    tree_->Branch("Sta_charge", &sta_charge_);
    
    tree_->Branch("Sta_dx", &sta_dx_);
    tree_->Branch("Sta_dy", &sta_dx_);
    tree_->Branch("Sta_dz", &sta_dz_);
    tree_->Branch("Sta_vx", &sta_vx_);
    tree_->Branch("Sta_vy", &sta_vy_);
    tree_->Branch("Sta_vz", &sta_vz_);
    tree_->Branch("Sta_DxyError", &sta_DxyError_);
    tree_->Branch("Sta_DzError",  &sta_DzError_);

    tree_->Branch("Sta_trackAlgo", &sta_trackAlgo_);


    //Roccester corrections
    tree_->Branch("Track_ScaleCorr", &track_ScaleCorr_);
    
    // tracks
    tree_->Branch("nTrack", &nTrack_,"nTrack_/s");
    tree_->Branch("Track_pt", &track_pt_);
    tree_->Branch("Track_eta", &track_eta_);
    tree_->Branch("Track_phi", &track_phi_);
    tree_->Branch("Track_charge", &track_charge_);
    
    tree_->Branch("Track_dx", &track_dx_);
    tree_->Branch("Track_dy", &track_dx_);
    tree_->Branch("Track_dz", &track_dz_);
    tree_->Branch("Track_vx", &track_vx_);
    tree_->Branch("Track_vy", &track_vy_);
    tree_->Branch("Track_vz", &track_vz_);
    tree_->Branch("Track_DxyError", &track_DxyError_);
    tree_->Branch("Track_DzError", &track_DzError_);

    tree_->Branch("Track_nPixelHits", &track_nPixelHits_);
    tree_->Branch("Track_nTrackerLayers", &track_nTrackerLayers_);
    tree_->Branch("Track_validFraction", &track_validFraction_);
    tree_->Branch("Track_trackAlgo", &track_trackAlgo_);

    //Roccester corrections
    tree_->Branch("Track_ScaleCorr", &track_ScaleCorr_);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZCountingAOD::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZCountingAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//--------------------------------------------------------------------------------------------------
void ZCountingAOD::clearVariables(){
    LogDebug("ZCountingAOD::clearVariables()");

    // Event info
    eventNumber_ = 0;

    eventweight_ = 0.;
    v_psWeight_.clear();
    v_meWeight_.clear();
    v_pdfWeight_.clear();

    nPV_ = 0;
    nPU_ = 0;
    decayMode_ = 0;
    z_genMass_ = 0.;
    z_recoMass_ = 0;
    
    met_triggerBits_ = 0;

    muon_genPt_ = 0.;
    muon_genEta_ = 0.;
    muon_genPhi_ = 0.;

    muon_genVx_ = 99.;
    muon_genVy_ = 99.;
    muon_genVz_ = 99.;
    // muon_genDxyError_ = 99.;
    // muon_genDzError_ = 99.;
        
    muon_genRecoMatches_ = 0;
    muon_genRecoObj_ = -1;
    muon_genRecoTrackMatches_ = 0;
    muon_genRecoTrackObj_ = -1;

    antiMuon_genPt_ = 0.;
    antiMuon_genEta_ = 0.;
    antiMuon_genPhi_ = 0.;

    antiMuon_genVx_ = 99.;
    antiMuon_genVy_ = 99.;
    antiMuon_genVz_ = 99.;
    // antiMuon_genDxyError_ = 99.;
    // antiMuon_genDzError_ = 99.;
        
    antiMuon_genRecoMatches_ = 0;
    antiMuon_genRecoObj_ = -1;
    antiMuon_genRecoTrackMatches_ = 0;
    antiMuon_genRecoTrackObj_ = -1;

    nMuon_ = 0;
    muon_pt_.clear();
    muon_eta_.clear();
    muon_phi_.clear();    
    muon_ptTrk_.clear();
    muon_etaTrk_.clear();
    muon_phiTrk_.clear(); 
    muon_ptSta_.clear();
    muon_etaSta_.clear();
    muon_phiSta_.clear();    
    muon_charge_.clear();
    muon_dx_.clear();
    muon_dy_.clear();
    muon_dz_.clear();

    muon_vx_.clear();
    muon_vy_.clear();
    muon_vz_.clear();
    muon_DxyError_.clear();
    muon_DzError_.clear();    
    
    muon_matchValue_.clear();
    muon_tkRelIso_.clear();
    muon_pfRelIso04_all_.clear();
    muon_ID_.clear();
    muon_triggerBits_.clear();
    
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
    muon_trackAlgo_.clear();

    muon_ScaleCorr_.clear();
    // muon_ScaleCorr_stat_RMS_.clear();
    // muon_ScaleCorr_Zpt_.clear();
    // muon_ScaleCorr_Ewk_.clear();
    // muon_ScaleCorr_deltaM_.clear();
    // muon_ScaleCorr_Ewk2_.clear();
    // muon_ScaleCorr_Total_.clear();

    nStandalone_ = 0;
    sta_pt_.clear();
    sta_eta_.clear();
    sta_phi_.clear();
    sta_charge_.clear();
    sta_dx_.clear();
    sta_dy_.clear();
    sta_dz_.clear();
    sta_vx_.clear();
    sta_vy_.clear();
    sta_vz_.clear();
    sta_DxyError_.clear();
    sta_DzError_.clear();  

    sta_trackAlgo_.clear();

    sta_ScaleCorr_.clear();

    nTrack_ = 0;
    track_pt_.clear();
    track_eta_.clear();
    track_phi_.clear();
    track_charge_.clear();
    track_dx_.clear();
    track_dy_.clear();
    track_dz_.clear();
    track_vx_.clear();
    track_vy_.clear();
    track_vz_.clear();
    track_DxyError_.clear();
    track_DzError_.clear();  

    track_nPixelHits_.clear();
    track_nTrackerLayers_.clear();
    track_validFraction_.clear();
    track_trackAlgo_.clear();

    track_ScaleCorr_.clear();

}

//--------------------------------------------------------------------------------------------------
std::string ZCountingAOD::get_triggerPath(std::string pattern, const edm::TriggerNames& triggerNames){
    std::string path = "";
    if (edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
        std::vector<std::vector<std::string>::const_iterator> matches =
            edm::regexMatch(triggerNames.triggerNames(), pattern);
        if (matches.empty()) {
            edm::LogWarning("ZCountingAOD") << "requested pattern [" << pattern << "] does not match any HLT paths";
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
int ZCountingAOD::getMuonID(const reco::Muon &mu, const reco::Vertex &vtx){
    if(muon::isTightMuon(mu, vtx) ) return 5;
    if(isCustomTightMuon(mu)) return 4;
    if(mu.isGlobalMuon()) return 3;
    if(mu.isStandAloneMuon()) return 2;
    if(isValidTrack(*(mu.innerTrack()))) return 1;
    return 0;
}

//--------------------------------------------------------------------------------------------------
// For trigger emulation in the 2017H (low PU) dataset
// We emulated the HLT_IsoMu24_v11 in separated samples and need to get the trigger objects from these separate samples

bool ZCountingAOD::isMuonTriggerObjEmulated(const double eta, const double phi, const long long unsigned eventNumber) {

    // filter tag for HLT_IsoMu24_v11:
    const edm::InputTag filterTag("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07", "", "HLTX");
    
    TTreeReader myReader("Events", _fileHLTEmulation);
    TTreeReaderValue<edm::EventAuxiliary> eventAuxiliary_(myReader, "EventAuxiliary");
    TTreeReaderValue<trigger::TriggerEvent> triggerEvent_(myReader, "triggerTriggerEvent_hltTriggerSummaryAOD__HLTX.obj");

    // std::cout<<"Look for event "<< eventNumber <<std::endl;

    while(myReader.Next()){
        // find event
        if(eventNumber != eventAuxiliary_->event())
            continue;

        // std::cout<<"Found event!"<<std::endl;
        // look for trigger objects
        if(triggerEvent_->filterIndex(filterTag) < triggerEvent_->sizeFilters()){
            const trigger::Keys& trigKeys = triggerEvent_->filterKeys(triggerEvent_->filterIndex(filterTag));
            const trigger::TriggerObjectCollection & trigObjColl(triggerEvent_->getObjects());
            //now loop of the trigger objects passing filter
            for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
                const trigger::TriggerObject& obj = trigObjColl[*keyIt];
                // std::cout<<"Trigger object(pt | eta) = "<<obj.pt()<< " | "<<obj.eta()<<std::endl;
                if (reco::deltaR(eta, phi, obj.eta(), obj.phi()) < DRMAX){
                    return true;
                }
            }
        }
        return false;
    }
    edm::LogWarning("isMuonTriggerObjEmulated")<<"Event was not found!"<<std::endl;
    return false;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ZCountingAOD);
