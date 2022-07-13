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
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// For CMSSW_12_X: 
// #include "CommonTools/Egamma/interface/ConversionTools.h"
// #include "CommonTools/Egamma/interface/EffectiveAreas.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

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
    
    float getMinDxy(const reco::Track &trk, const std::vector<reco::Vertex> &pvCol);
    float getMinDz(const reco::Track &trk, const std::vector<reco::Vertex> &pvCol);

    int getMuonID(const reco::Muon&, const reco::Vertex&);
        
    bool isMuonTriggerObjEmulated(const double eta, const double phi, const long long unsigned eventNumber);

    float getElectronIso(const reco::GsfElectron &el, const double rho);

    // ----------member data ---------------------------
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

    // --- input
    const edm::InputTag triggerResultsInputTag_;
    edm::EDGetTokenT<std::vector<reco::Muon>> muonCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> standaloneCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> standaloneUpdatedCollection_;
    edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronCollection_;
    edm::EDGetTokenT<std::vector<reco::SuperCluster>> superclusterCollection_;
    edm::EDGetTokenT<double> rho_;
    edm::EDGetTokenT<reco::BeamSpot> beamspot_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> conversionCollection_;
    edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> pvCollection_;
    edm::EDGetTokenT<std::vector<GenZDecayProperties>> genZCollection_;
    edm::EDGetTokenT<TtGenEvent> genTtEvent_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoCollection_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;
    edm::EDGetTokenT<LHEEventProduct> lheEventInfo_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleCollection_;

    // Effective area constants
    EffectiveAreas effectiveAreas_;

    std::string era_;

    // Triggers
    bool hltChanged_;

    std::vector<std::string> muonTriggerPatterns_;      // Muon trigger    
    std::vector<std::string> electronTriggerPatterns_;  // Electron trigger    
    std::vector<std::string> metTriggerPatterns_;       // MET trigger

    // max dR matching between muon and hlt object
    double DRMAX;

    // flags
    bool isData_;
    bool hasGenZ_;
    bool hasGenTt_;
    
    bool store_muons_;
    bool store_electrons_;

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

    float lepton_genPt_;
    float lepton_genEta_;
    float lepton_genPhi_;

    float lepton_genVx_;
    float lepton_genVy_;
    float lepton_genVz_;
    
    float antiLepton_genPt_;
    float antiLepton_genEta_;
    float antiLepton_genPhi_;

    float antiLepton_genVx_;
    float antiLepton_genVy_;
    float antiLepton_genVz_;
    
    unsigned int nMuon_;
    std::vector<float> muon_pt_;
    std::vector<float> muon_eta_;
    std::vector<float> muon_phi_;
    std::vector<int> muon_charge_;
    
    std::vector<float> muon_ptTrk_;
    std::vector<float> muon_etaTrk_;
    std::vector<float> muon_phiTrk_;
    std::vector<int> muon_chargeTrk_;
    
    std::vector<float> muon_ptStaReg_;
    std::vector<float> muon_etaStaReg_;
    std::vector<float> muon_phiStaReg_;
    std::vector<int> muon_chargeStaReg_;

    std::vector<int> muon_nStationsStaReg_;

    std::vector<float> muon_ptStaUpd_;
    std::vector<float> muon_etaStaUpd_;
    std::vector<float> muon_phiStaUpd_;
    std::vector<int> muon_chargeStaUpd_;
    
    std::vector<bool> muon_useUpdated_;

    // impact parameter from (first good) primary vertex
    std::vector<float> muon_dxyPV_;
    std::vector<float> muon_dzPV_;
    
    std::vector<float> muon_dxyPVmin_;
    std::vector<float> muon_dzPVmin_;
    
    // reference point position
    std::vector<float> muon_dx_;
    std::vector<float> muon_dy_;
    std::vector<float> muon_dz_;
    // uncertainties on position
    std::vector<float> muon_DxyError_;
    std::vector<float> muon_DzError_;
    
    // std::vector<bool> muon_isDoublicate_;
    std::vector<float> muon_tkRelIso_;
    std::vector<float> muon_pfRelIso04_all_;
    std::vector<int> muon_ID_;
    std::vector<int> muon_triggerBits_;

    std::vector<bool> muon_isMedium_;
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

    unsigned int nTrack_;
    std::vector<float> track_pt_;
    std::vector<float> track_eta_;
    std::vector<float> track_phi_;
    std::vector<int> track_charge_;

    std::vector<float> track_dxyPV_;
    std::vector<float> track_dzPV_;    
    std::vector<float> track_dxyPVmin_;
    std::vector<float> track_dzPVmin_;
    std::vector<float> track_dx_;
    std::vector<float> track_dy_;
    std::vector<float> track_dz_;
    std::vector<float> track_DxyError_;
    std::vector<float> track_DzError_;

    std::vector<int> track_nPixelHits_;
    std::vector<int> track_nTrackerLayers_;
    std::vector<float> track_validFraction_;
    std::vector<float> track_trackAlgo_;

    std::vector<float> track_ScaleCorr_;

    unsigned int nElectron_;
    std::vector<float> electron_pt_;
    std::vector<float> electron_eta_;
    std::vector<float> electron_phi_;
    std::vector<int> electron_charge_;

    std::vector<float> electron_dxyPV_;
    std::vector<float> electron_dzPV_;    
    std::vector<float> electron_dxyPVmin_;
    std::vector<float> electron_dzPVmin_;
        
    std::vector<float> electron_superclusterEnergy_;
    std::vector<float> electron_superclusterEta_;
    std::vector<float> electron_superclusterPhi_;
    std::vector<float> electron_dEtaInSeed_;
    std::vector<float> electron_pfRelIso_;
    std::vector<bool> electron_hasMatchedConversion_;
    std::vector<float> electron_full5x5_sigmaIetaIeta_;
    std::vector<float> electron_deltaPhiSuperClusterTrackAtVtx_;
    std::vector<float> electron_hadronicOverEm_;
    std::vector<float> electron_eSuperClusterOverP_;
    std::vector<float> electron_ecalEnergy_;
    std::vector<int> electron_gsfTrackMissingInnerHits_;

    std::vector<int> electron_triggerBits_;

    unsigned int nSupercluster_;
    std::vector<float> supercluster_energy_;
    std::vector<float> supercluster_eta_;
    std::vector<float> supercluster_phi_;
    
};

//
// constructors and destructor
//
ZCountingAOD::ZCountingAOD(const edm::ParameterSet& iConfig):
    triggerResultsInputTag_(iConfig.getParameter<edm::InputTag>("TriggerResults")),
    muonCollection_  (consumes<std::vector<reco::Muon>> (iConfig.getParameter<edm::InputTag>("reco_muons"))),
    standaloneCollection_  (consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("reco_standalones"))),
    standaloneUpdatedCollection_  (consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("reco_standalonesUpdated"))),
    electronCollection_  (consumes<std::vector<reco::GsfElectron>> (iConfig.getParameter<edm::InputTag>("reco_electrons"))),
    superclusterCollection_  (consumes<std::vector<reco::SuperCluster>> (iConfig.getParameter<edm::InputTag>("reco_superclusters"))),
    rho_ (consumes<double> (iConfig.getParameter<edm::InputTag>("rhoName"))),
    beamspot_ (consumes<reco::BeamSpot> (iConfig.getParameter<edm::InputTag>("beamspotName"))),
    conversionCollection_ (consumes<std::vector<reco::Conversion>> (iConfig.getParameter<edm::InputTag>("conversionsName"))),
    trackCollection_  (consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("reco_tracks"))),
    pvCollection_    (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("edmPVName"))),
    genZCollection_  (consumes<std::vector<GenZDecayProperties> > (iConfig.getParameter<edm::InputTag>("genZLeptonCollection"))),
    genTtEvent_  (consumes<TtGenEvent> (iConfig.getParameter<edm::InputTag>("genTtCollection"))),
    pileupInfoCollection_  (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfoCollection"))),
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfo_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    genParticleCollection_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath())
{
    LogDebug("ZCountingAOD")<<"ZCountingAOD(...)";

    era_ = iConfig.getParameter<std::string>("era");

    isData_ = iConfig.getUntrackedParameter<bool>("isData");
    hasGenZ_ = iConfig.getUntrackedParameter<bool>("hasGenZ");
    hasGenTt_ = iConfig.getUntrackedParameter<bool>("hasGenTt");
    
    store_muons_ = iConfig.getUntrackedParameter<bool>("store_muons");
    store_electrons_ = iConfig.getUntrackedParameter<bool>("store_electrons");

    v_genWeightIDs_ = iConfig.getParameter<std::vector<int>>("genWeights");
    v_pdfWeightIDs_ = iConfig.getParameter<std::vector<int>>("pdfWeights");

    roccorFile = iConfig.getParameter<std::string>("roccorFile");

    hltChanged_ = true;
    emulateTrigger_ = iConfig.getUntrackedParameter<bool>("emulateTrigger");

    muonTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("muon_trigger_patterns");
    electronTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("electron_trigger_patterns");
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
    for(unsigned int i = 0; i < electronTriggerPatterns_.size(); ++i) {
        triggers->addTriggerRecord(electronTriggerPatterns_.at(i));
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
    
    // -------------------------------------------------------------------------
    // Generator information
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

        // --- Pileup 
        edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoCollection;
        iEvent.getByToken(pileupInfoCollection_, pileupInfoCollection);
        
        for(std::vector<PileupSummaryInfo>::const_iterator puI = pileupInfoCollection->begin(); puI != pileupInfoCollection->end(); ++puI)
        {
            if(puI->getBunchCrossing() == 0){
                nPU_ = puI->getTrueNumInteractions();
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

        if(genLepton){
            lepton_genPt_  = genLepton->pt();
            lepton_genEta_ = genLepton->eta();
            lepton_genPhi_ = genLepton->phi();
            
            lepton_genVx_ = genLepton->vx();
            lepton_genVy_ = genLepton->vy();
            lepton_genVz_ = genLepton->vz();
        }
        if(genAntiLepton){
            antiLepton_genPt_  = genAntiLepton->pt();
            antiLepton_genEta_ = genAntiLepton->eta();
            antiLepton_genPhi_ = genAntiLepton->phi();
            
            antiLepton_genVx_ = genAntiLepton->vx();
            antiLepton_genVy_ = genAntiLepton->vy();
            antiLepton_genVz_ = genAntiLepton->vz();        
        }
        if(genLepton && genAntiLepton){
            z_genMass_ = (genLepton->p4() + genAntiLepton->p4()).M();
        }
    }

    // --- PV selection
    edm::Handle<std::vector<reco::Vertex> > pvCollection;
    iEvent.getByToken(pvCollection_, pvCollection);
    const std::vector<reco::Vertex> *pvCol = pvCollection.product();
    const reco::Vertex *pv = 0;
    nPV_ = 0;
    for (const reco::Vertex &itVtx : *pvCol) {
        if(!isGoodPV(itVtx))
            continue;
        if (nPV_ == 0)
            pv = &itVtx;
        nPV_++;
    }

    for(unsigned j = 0; j < metTriggerPatterns_.size(); ++j){
        met_triggerBits_ += std::pow(2,j) * triggers->pass(metTriggerPatterns_.at(j));
    }
    
    // -------------------------------------------------------------------------
    // Muons
    if(store_muons_){
        // --- store all standAloneMuons
        LogDebug("ZCountingAOD") << "find reco muons";
        edm::Handle<std::vector<reco::Muon>> muonCollection;
        edm::Handle<std::vector<reco::Track>> standaloneCollection;
        edm::Handle<std::vector<reco::Track>> standaloneUpdatedCollection;
        iEvent.getByToken(muonCollection_, muonCollection);
        iEvent.getByToken(standaloneCollection_, standaloneCollection);
        iEvent.getByToken(standaloneUpdatedCollection_, standaloneUpdatedCollection);
        for (const reco::Track &StaReg : *standaloneCollection){
            
            if (StaReg.extra().isNull())
                edm::LogError("ZCountingAOD") << "Track from standAloneMuons has no TrackExtra ref.";

            const reco::Track *StaUpd = 0;
            // standAloneMuons:UpdatedAtVtx is a subset of standAloneMuons
            for (const reco::Track &trk : *standaloneUpdatedCollection){
                // check if a candidate in standAloneMuons:UpdatedAtVtx has been reconstructed
                if (trk.extra().isNull())
                    edm::LogError("ZCountingAOD") << "Track from standAloneMuons:UpdatedAtVtx has no TrackExtra ref.";
            
                if (StaReg.extra().get() == trk.extra().get()) {
                    // we found the corresponding candidate
                    StaUpd = &trk;
                    break;                
                }
            }

            bool useUpdated = false;
            if ((StaUpd != 0) && (StaUpd->eta() * StaReg.eta() >= 0))
                useUpdated = true;
            
            // globalMuons is a subset of standAloneMuons
            //    it uses the track from standAloneMuons:UpdatedAtVtx if it is available and
            //    the sign of eta is the same as of the corresponding track from standAloneMuons
            //    otherwise the track from standAloneMuons is used
            const reco::Muon *mu = 0;
            for (const reco::Muon &muon : *muonCollection) {
                if (muon.outerTrack().isNull())
                    continue;

                if (useUpdated && muon.outerTrack().get() == StaUpd) {
                    // we found the corresponding candidate through common TrackExtra from standAloneMuons:UpdatedAtVtx
                    mu = &muon;
                    break;                
                }      
                else if (!useUpdated && muon.outerTrack().get() == &StaReg) {
                    // we found the corresponding candidate through direct link to standAloneMuons
                    mu = &muon;
                    break;    
                }           
            }    
            
            // inner track of muon
            const reco::Track *trk = 0;
            if (mu!=0 && mu->innerTrack().isNonnull())
                trk = mu->innerTrack().get();
            
            // reject muons only if best track, both outer tracks, and inner track fail acceptance cuts
            if (   (             std::abs(StaReg.eta()) > 2.4  || StaReg.pt() < 20)
                && (mu==0     || std::abs(mu->eta()) > 2.4     || mu->pt() < 20)
                && (trk==0    || std::abs(trk->eta()) > 2.4    || trk->pt() < 20)
                && (StaUpd==0 || std::abs(StaUpd->eta()) > 2.4 || StaUpd->pt() < 20)
            ) {
                continue;
            }

            muon_ptStaReg_.push_back(StaReg.pt());
            muon_etaStaReg_.push_back(StaReg.eta());
            muon_phiStaReg_.push_back(StaReg.phi());     
            muon_chargeStaReg_.push_back(StaReg.charge());

            muon_nStationsStaReg_.push_back(StaReg.hitPattern().muonStationsWithValidHits());
            
            if(StaUpd != 0 && StaReg.hitPattern().muonStationsWithValidHits() != StaUpd->hitPattern().muonStationsWithValidHits()){
                std::cout<< StaReg.hitPattern().muonStationsWithValidHits() <<" | "<< StaUpd->hitPattern().muonStationsWithValidHits() <<std::endl;
            }
            
            muon_useUpdated_.push_back(useUpdated);

            if(StaUpd != 0){
                muon_ptStaUpd_.push_back(StaUpd->pt());
                muon_etaStaUpd_.push_back(StaUpd->eta());
                muon_phiStaUpd_.push_back(StaUpd->phi());            
                muon_chargeStaUpd_.push_back(StaUpd->charge());

            }
            else{            
                muon_ptStaUpd_.push_back(-1);
                muon_etaStaUpd_.push_back(-1);
                muon_phiStaUpd_.push_back(-1);                 
                muon_chargeStaUpd_.push_back(0);
            }

            if(trk != 0){
                muon_ptTrk_.push_back(trk->pt());
                muon_etaTrk_.push_back(trk->eta());
                muon_phiTrk_.push_back(trk->phi());               
                muon_chargeTrk_.push_back(trk->charge());
            
                muon_nPixelHits_.push_back(trk->hitPattern().numberOfValidPixelHits());
                muon_nTrackerLayers_.push_back(trk->hitPattern().trackerLayersWithMeasurement());
                muon_validFraction_.push_back(trk->validFraction());
                muon_trackAlgo_.push_back(trk->originalAlgo());
            }
            else{
                muon_ptTrk_.push_back(-1);
                muon_etaTrk_.push_back(-1);
                muon_phiTrk_.push_back(-1);    
                muon_chargeTrk_.push_back(0);
        
                muon_nPixelHits_.push_back(-1);
                muon_nTrackerLayers_.push_back(-1);
                muon_validFraction_.push_back(-1);
                muon_trackAlgo_.push_back(-1);

            }

            if(mu!=0){
                muon_pt_.push_back(mu->pt());
                muon_eta_.push_back(mu->eta());
                muon_phi_.push_back(mu->phi());          
                muon_charge_.push_back(mu->charge());
                
                // Look for the distance to the first good primary vertex
                muon_dxyPV_.push_back(pv != nullptr ? mu->muonBestTrack()->dxy(pv->position()) : 999);
                muon_dzPV_.push_back(pv != nullptr ? mu->muonBestTrack()->dz(pv->position()) : 999);

                // Look for the distance to the first good primary vertex
                muon_dxyPVmin_.push_back(getMinDxy(*(mu->muonBestTrack()), *pvCol));
                muon_dzPVmin_.push_back(getMinDz(*(mu->muonBestTrack()), *pvCol));
                
                // Look for the distance to (0,0,0)
                muon_dx_.push_back(mu->muonBestTrack()->referencePoint().x());
                muon_dy_.push_back(mu->muonBestTrack()->referencePoint().y());
                muon_dz_.push_back(mu->muonBestTrack()->referencePoint().z());
                // uncertainty
                muon_DxyError_.push_back(mu->dxyError());
                muon_DzError_.push_back(mu->dzError());
                
                muon_tkRelIso_.push_back(getTkIso(*mu));
                muon_pfRelIso04_all_.push_back(getPFIso(*mu));
                muon_ID_.push_back(getMuonID(*mu, *pv));
                
                muon_isMedium_.push_back(muon::isMediumMuon(*mu));
                muon_isTracker_.push_back(mu->isTrackerMuon());
                muon_isGlobal_.push_back(mu->isGlobalMuon());
                muon_isPFCand_.push_back(mu->isPFMuon());
                
                muon_nStations_.push_back(mu->numberOfMatchedStations());
                muon_SegmentCompatibility_.push_back(muon::segmentCompatibility(*mu));
                muon_chi2LocalPosition_.push_back(mu->combinedQuality().chi2LocalPosition);
                muon_trkKink_.push_back(mu->combinedQuality().trkKink);    

                if(mu->isGlobalMuon()){
                    muon_normChi2_.push_back(mu->globalTrack()->normalizedChi2());
                    muon_nTrackHits_.push_back(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
                }
                else{
                    muon_normChi2_.push_back(-1);
                    muon_nTrackHits_.push_back(-1);
                }

                int bits_ = 0;    
                if(emulateTrigger_){
                    bits_ = isMuonTriggerObjEmulated(mu->eta(), mu->phi(), eventNumber_);
                }
                else{
                    for(unsigned j = 0; j < muonTriggerPatterns_.size(); ++j){
                        bits_ += std::pow(2,j) * triggers->passObj(muonTriggerPatterns_.at(j), mu->eta(), mu->phi());
                    }
                }
                muon_triggerBits_.push_back(bits_);   

                double roccorSF = 1.; // Rochester correction
                if(!iEvent.isRealData()){
                
                    if(genLepton && mu->pdgId() == 13 && reco::deltaR(mu->eta(), mu->phi(), genLepton->eta(), genLepton->phi()) < 0.03){
                        roccorSF = rc.kSpreadMC(mu->charge(), mu->pt(), mu->eta(), mu->phi(), lepton_genPt_);
                    }
                    else if(genAntiLepton && mu->pdgId() == -13 && reco::deltaR(mu->eta(), mu->phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                        roccorSF = rc.kSpreadMC(mu->charge(), mu->pt(), mu->eta(), mu->phi(), antiLepton_genPt_);
                    }
                    else{
                        const int nl = mu->innerTrack().isNonnull() ? mu->innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0;
                        roccorSF = rc.kSmearMC(mu->charge(), mu->pt(), mu->eta(), mu->phi(), nl, rand_.Rndm());
                    }
                
                    // store muon Rochester corrections and uncertainties
                
                
                    // muon_ScaleCorr_.push_back(mu->hasUserFloat("MuonEnergyCorr")                   ? mu->userFloat("MuonEnergyCorr")          : 1.0);
                    // muon_ScaleCorr_stat_RMS_.push_back(mu->hasUserFloat("MuonEnergyCorr_stat_RMS") ? mu->userFloat("MuonEnergyCorr_stat_RMS") : 0.0);
                    // muon_ScaleCorr_Zpt_.push_back(mu->hasUserFloat("MuonEnergyCorr_Zpt")           ? mu->userFloat("MuonEnergyCorr_Zpt")      : 0.0);
                    // muon_ScaleCorr_Ewk_.push_back(mu->hasUserFloat("MuonEnergyCorr_Ewk")           ? mu->userFloat("MuonEnergyCorr_Ewk")      : 0.0);
                    // muon_ScaleCorr_deltaM_.push_back(mu->hasUserFloat("MuonEnergyCorr_deltaM")     ? mu->userFloat("MuonEnergyCorr_deltaM")   : 0.0);
                    // muon_ScaleCorr_Ewk2_.push_back(mu->hasUserFloat("MuonEnergyCorr_Ewk2")         ? mu->userFloat("MuonEnergyCorr_Ewk2")     : 0.0);
                    // muon_ScaleCorr_Total_.push_back(mu->hasUserFloat("MuonEnergyCorr_Total")       ? mu->userFloat("MuonEnergyCorr_Total")    : 0.0);
                }
                else{
                    roccorSF = rc.kScaleDT(mu->charge(), mu->pt(), mu->eta(), mu->phi());
                }
                
                // Rochester corrections 
                muon_ScaleCorr_.push_back(roccorSF);
            }
            else{
                // if no muon candidate was reconstructed, use the parameters of the standalone muon that was considered
                if(useUpdated){
                    muon_pt_.push_back(StaUpd->pt());
                    muon_eta_.push_back(StaUpd->eta());
                    muon_phi_.push_back(StaUpd->phi());     
                    muon_charge_.push_back(StaUpd->charge());  
                    
                    // Look for distance to first good primary vertex
                    muon_dxyPV_.push_back(pv != nullptr ? StaUpd->dxy(pv->position()) : 999);          
                    muon_dzPV_.push_back(pv != nullptr ? StaUpd->dz(pv->position()) : 999);          

                    // Look for the distance to the first good primary vertex
                    muon_dxyPVmin_.push_back(getMinDxy(*StaUpd, *pvCol));
                    muon_dzPVmin_.push_back(getMinDz(*StaUpd, *pvCol));

                    // Look for the distance to (0,0,0)
                    muon_dx_.push_back(StaUpd->referencePoint().x());
                    muon_dy_.push_back(StaUpd->referencePoint().y());
                    muon_dz_.push_back(StaUpd->referencePoint().z());
                    // uncertainty
                    muon_DxyError_.push_back(StaUpd->dxyError());
                    muon_DzError_.push_back(StaUpd->dzError());
                }
                else{
                    muon_pt_.push_back(StaReg.pt());
                    muon_eta_.push_back(StaReg.eta());
                    muon_phi_.push_back(StaReg.phi());     
                    muon_charge_.push_back(StaReg.charge());

                    // Look for distance to first good primary vertex
                    muon_dxyPV_.push_back(pv != nullptr ? StaReg.dxy(pv->position()): 999);          
                    muon_dzPV_.push_back(pv != nullptr ? StaReg.dz(pv->position()): 999);    

                    // Look for the distance to the first good primary vertex
                    muon_dxyPVmin_.push_back(getMinDxy(StaReg, *pvCol));
                    muon_dzPVmin_.push_back(getMinDz(StaReg, *pvCol));

                    // Look for the distance to (0,0,0)
                    muon_dx_.push_back(StaReg.referencePoint().x());
                    muon_dy_.push_back(StaReg.referencePoint().y());
                    muon_dz_.push_back(StaReg.referencePoint().z());
                    // uncertainty
                    muon_DxyError_.push_back(StaReg.dxyError());
                    muon_DzError_.push_back(StaReg.dzError());
                }
            
                muon_tkRelIso_.push_back(-1);
                muon_pfRelIso04_all_.push_back(-1);
                muon_ID_.push_back(-1);
                
                muon_isMedium_.push_back(0);
                muon_isTracker_.push_back(0);
                muon_isGlobal_.push_back(0);
                muon_isPFCand_.push_back(0);
                
                muon_nStations_.push_back(-1);
                muon_SegmentCompatibility_.push_back(-1);
                muon_chi2LocalPosition_.push_back(-1);
                muon_trkKink_.push_back(-1);       
                
                muon_normChi2_.push_back(-1);
                muon_nTrackHits_.push_back(-1);     

                muon_triggerBits_.push_back(0);       

                muon_ScaleCorr_.push_back(1.);
            }
            
            nMuon_++;
        }

        // --- store all reco tracks
        edm::Handle<std::vector<reco::Track>> trackCollection;
        iEvent.getByToken(trackCollection_, trackCollection);
        for (const reco::Track &trk : *trackCollection){        
            if(std::abs(trk.eta()) > 2.4) continue;
            if(trk.pt() < 20) continue;
        
            // Check if track is already in muons
            bool isMuon = false;
            for (const reco::Muon &mu : *muonCollection) {
                // only standalone muons are in muon collection
                if(!mu.isStandAloneMuon())
                    continue;
                    
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
            
            // Look for distance to first good primary vertex
            track_dxyPV_.push_back(pv != nullptr ? trk.dxy(pv->position()) : 999);
            track_dzPV_.push_back(pv != nullptr ? trk.dz(pv->position()) : 999);

            // Look for the distance to the first good primary vertex
            track_dxyPVmin_.push_back(getMinDxy(trk, *pvCol));
            track_dzPVmin_.push_back(getMinDz(trk, *pvCol));
            
            // Look for the distance to (0,0,0)
            track_dx_.push_back(trk.referencePoint().x());
            track_dy_.push_back(trk.referencePoint().y());
            track_dz_.push_back(trk.referencePoint().z());
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
                    roccorSF = rc.kSpreadMC(trk.charge(), trk.pt(), trk.eta(), trk.phi(), lepton_genPt_);
                }
        
                else if(genAntiLepton && trk.charge() < 0 && reco::deltaR(trk.eta(), trk.phi(), genAntiLepton->eta(), genAntiLepton->phi()) < 0.03){
                    roccorSF = rc.kSpreadMC(trk.charge(), trk.pt(), trk.eta(), trk.phi(), antiLepton_genPt_);
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

    }
    
    // -------------------------------------------------------------------------
    // Electrons
    if(store_electrons_){
        // --- store all gsf electrons
        LogDebug("ZCountingAOD") << "find reco gsf electrons";
        // Get Rho
        edm::Handle<double> rho;
        edm::Handle<reco::BeamSpot> beamspot;
        edm::Handle<std::vector<reco::Conversion>> conversionCollection;
        edm::Handle<std::vector<reco::GsfElectron>> electronCollection;
        
        iEvent.getByToken(rho_, rho);
        iEvent.getByToken(beamspot_, beamspot);
        iEvent.getByToken(conversionCollection_, conversionCollection);
        iEvent.getByToken(electronCollection_, electronCollection);
        for (const reco::GsfElectron &el : *electronCollection){
            
            const reco::SuperCluster *sc = 0;
            if(el.superCluster().isNonnull())
                sc = el.superCluster().get();
                
            const float sc_pt = sc->energy() * sqrt(1 - std::pow(std::tanh(sc->eta()), 2));

            if (                (std::abs(el.eta()) > 2.5 || el.pt() < 20)
                && (sc==0    || std::abs(sc->eta()) > 2.5 || sc_pt < 20))          
            {
                continue;
            }
            

            electron_pt_.push_back(el.pt());
            electron_eta_.push_back(el.eta());
            electron_phi_.push_back(el.phi());     
            electron_charge_.push_back(el.charge());
            
            electron_full5x5_sigmaIetaIeta_.push_back(el.full5x5_sigmaIetaIeta());
            electron_deltaPhiSuperClusterTrackAtVtx_.push_back(el.deltaPhiSuperClusterTrackAtVtx());
            electron_hadronicOverEm_.push_back(el.hadronicOverEm());
            electron_eSuperClusterOverP_.push_back(el.eSuperClusterOverP());
            electron_ecalEnergy_.push_back(el.ecalEnergy());
            
            electron_pfRelIso_.push_back(getElectronIso(el, *rho) / el.pt());
            electron_hasMatchedConversion_.push_back(ConversionTools::hasMatchedConversion(el, *conversionCollection, beamspot->position()));

            if(sc != 0){
                electron_superclusterEnergy_.push_back(sc->energy());
                electron_superclusterEta_.push_back(sc->eta());
                electron_superclusterPhi_.push_back(sc->phi());
            }
            else{
                electron_superclusterEnergy_.push_back(0);
                electron_superclusterEta_.push_back(999);
                electron_superclusterPhi_.push_back(999);
            }
            if(sc != 0 && sc->seed().isNonnull()){
                 electron_dEtaInSeed_.push_back(el.deltaEtaSuperClusterTrackAtVtx() - sc->eta() + sc->seed()->eta());
            }
            else{
                 electron_dEtaInSeed_.push_back(999);
            }

            const reco::Track *gsf = 0;
            if(el.gsfTrack().isNonnull())
                gsf = el.gsfTrack().get();
            
            if(gsf != 0){
                // Look for distance to first good primary vertex
                electron_dxyPV_.push_back(pv != nullptr ? gsf->dxy(pv->position()) : 999);
                electron_dzPV_.push_back(pv != nullptr ? gsf->dz(pv->position()) : 999);

                // Look for the distance to the first good primary vertex
                electron_dxyPVmin_.push_back(getMinDxy(*gsf, *pvCol));
                electron_dzPVmin_.push_back(getMinDz(*gsf, *pvCol));

                electron_gsfTrackMissingInnerHits_.push_back(gsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
            }
            else{
                electron_dxyPV_.push_back(999);
                electron_dzPV_.push_back(999);
                electron_dxyPVmin_.push_back(999);
                electron_dzPVmin_.push_back(999);
                
                electron_gsfTrackMissingInnerHits_.push_back(-1);
            }
            
            int bits_ = 0;
            for(unsigned j = 0; j < electronTriggerPatterns_.size(); ++j){
                bits_ += std::pow(2,j) * triggers->passObj(electronTriggerPatterns_.at(j), el.eta(), el.phi());
            }
            electron_triggerBits_.push_back(bits_);   

            nElectron_++;
        }
        
        // store all superclusters that are not already stored as electrons
        edm::Handle<std::vector<reco::SuperCluster>> superclusterCollection;
        iEvent.getByToken(superclusterCollection_, superclusterCollection);

        for (const reco::SuperCluster &sc : *superclusterCollection){
            
            const float sc_pt = sc.energy() * sqrt(1 - std::pow(std::tanh(sc.eta()), 2));
            
            if (std::abs(sc.eta()) > 2.5 || sc_pt < 20){
                continue;
            }

            // Check if track is already in electrons
            bool isElectron = false;
            for (const reco::GsfElectron &el : *electronCollection){
                
                if(el.superCluster().isNonnull() && el.superCluster().get() == &sc){
                    isElectron = true;
                    break;
                }
            }
            if(isElectron)
                continue;

            supercluster_energy_.push_back(sc.energy());
            supercluster_eta_.push_back(sc.eta());
            supercluster_phi_.push_back(sc.phi());
            
            nSupercluster_++;
            
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
        tree_->Branch("decayMode",  &decayMode_,  "decayMode_/i");
        tree_->Branch("z_genMass",  &z_genMass_,  "z_genMass_/f");
        tree_->Branch("z_recoMass", &z_recoMass_, "z_recoMass_/f");
        
        // muon
        tree_->Branch("lepton_genPt",  &lepton_genPt_,  "lepton_genPt_/f");
        tree_->Branch("lepton_genEta", &lepton_genEta_, "lepton_genEta_/f");
        tree_->Branch("lepton_genPhi", &lepton_genPhi_, "lepton_genPhi_/f");
        
        tree_->Branch("lepton_vx", &lepton_genVx_, "lepton_genVx_/f");
        tree_->Branch("lepton_vy", &lepton_genVy_, "lepton_genVy_/f");
        tree_->Branch("lepton_vz", &lepton_genVz_, "lepton_genVz_/f");

        // anti muon
        tree_->Branch("antiLepton_genPt",  &antiLepton_genPt_,  "antiLepton_genPt_/f");
        tree_->Branch("antiLepton_genEta", &antiLepton_genEta_, "antiLepton_genEta_/f");
        tree_->Branch("antiLepton_genPhi", &antiLepton_genPhi_, "antiLepton_genPhi_/f");
    
        tree_->Branch("antiLepton_vx", &antiLepton_genVx_, "antiLepton_genVx_/f");
        tree_->Branch("antiLepton_vy", &antiLepton_genVy_, "antiLepton_genVy_/f");
        tree_->Branch("antiLepton_vz", &antiLepton_genVz_, "antiLepton_genVz_/f");        
    }

    // reco level info
    tree_->Branch("MET_triggerBits", &met_triggerBits_);
    
    if(store_muons_){
        // muons
        tree_->Branch("nMuon", &nMuon_,"nMuon_/s");
        tree_->Branch("Muon_pt", &muon_pt_);
        tree_->Branch("Muon_eta", &muon_eta_);
        tree_->Branch("Muon_phi", &muon_phi_);
        tree_->Branch("Muon_charge", &muon_charge_);

        tree_->Branch("Muon_ptTrk", &muon_ptTrk_);
        tree_->Branch("Muon_etaTrk", &muon_etaTrk_);
        tree_->Branch("Muon_phiTrk", &muon_phiTrk_);
        tree_->Branch("Muon_chargeTrk", &muon_chargeTrk_);
        
        tree_->Branch("Muon_ptStaReg", &muon_ptStaReg_);
        tree_->Branch("Muon_etaStaReg", &muon_etaStaReg_);
        tree_->Branch("Muon_phiStaReg", &muon_phiStaReg_);
        tree_->Branch("Muon_chargeStaReg", &muon_chargeStaReg_);

        tree_->Branch("Muon_nStationsStaReg", &muon_nStationsStaReg_);

        tree_->Branch("Muon_ptStaUpd", &muon_ptStaUpd_);
        tree_->Branch("Muon_etaStaUpd", &muon_etaStaUpd_);
        tree_->Branch("Muon_phiStaUpd", &muon_phiStaUpd_);
        tree_->Branch("Muon_chargeStaUpd", &muon_chargeStaUpd_);
        
        tree_->Branch("Muon_useUpdated", &muon_useUpdated_);

        tree_->Branch("Muon_dxyPV", &muon_dxyPV_);
        tree_->Branch("Muon_dzPV", &muon_dzPV_);  
        
        tree_->Branch("Muon_dxyPVmin", &muon_dxyPVmin_);  
        tree_->Branch("Muon_dzPVmin", &muon_dzPVmin_);  

        tree_->Branch("Muon_dx", &muon_dx_);
        tree_->Branch("Muon_dy", &muon_dy_);
        tree_->Branch("Muon_dz", &muon_dz_);
        tree_->Branch("Muon_DxyError", &muon_DxyError_);
        tree_->Branch("Muon_DzError", &muon_DzError_);
        
        tree_->Branch("Muon_ID", &muon_ID_);
        tree_->Branch("Muon_tkRelIso", &muon_tkRelIso_);
        tree_->Branch("Muon_pfRelIso04_all", &muon_pfRelIso04_all_);
        tree_->Branch("Muon_triggerBits", &muon_triggerBits_);

        tree_->Branch("Muon_isMedium", &muon_isMedium_);
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
    }

    if(store_electrons_){
        // electrons
        tree_->Branch("nElectron", &nElectron_,"nElectron_/s");
        tree_->Branch("Electron_pt", &electron_pt_);
        tree_->Branch("Electron_eta", &electron_eta_);
        tree_->Branch("Electron_phi", &electron_phi_);
        tree_->Branch("Electron_charge", &electron_charge_);

        tree_->Branch("Electron_dxyPV", &electron_dxyPV_);
        tree_->Branch("Electron_dzPV", &electron_dzPV_);
        tree_->Branch("Electron_dxyPVmin", &electron_dxyPVmin_);  
        tree_->Branch("Electron_dzPVmin", &electron_dzPVmin_);  
        
        tree_->Branch("Electron_superclusterEnergy", &electron_superclusterEnergy_);
        tree_->Branch("Electron_superclusterEta", &electron_superclusterEta_);
        tree_->Branch("Electron_superclusterPhi", &electron_superclusterPhi_);
        
        tree_->Branch("Electron_dEtaInSeed", &electron_dEtaInSeed_);
        tree_->Branch("Electron_full5x5_sigmaIetaIeta", &electron_full5x5_sigmaIetaIeta_);
        tree_->Branch("Electron_deltaPhiSuperClusterTrackAtVtx", &electron_deltaPhiSuperClusterTrackAtVtx_);
        tree_->Branch("Electron_hadronicOverEm", &electron_hadronicOverEm_);
        tree_->Branch("Electron_eSuperClusterOverP", &electron_eSuperClusterOverP_);
        tree_->Branch("Electron_ecalEnergy", &electron_ecalEnergy_);
        tree_->Branch("Electron_gsfTrackMissingInnerHits", &electron_gsfTrackMissingInnerHits_);
        tree_->Branch("Electron_pfRelIso", &electron_pfRelIso_);
        tree_->Branch("Electron_hasMatchedConversion", &electron_hasMatchedConversion_);
        
        tree_->Branch("Electron_triggerBits", &electron_triggerBits_);

        tree_->Branch("nSupercluster", &nSupercluster_,"nSupercluster_/s");
        tree_->Branch("Supercluster_energy", &supercluster_energy_);
        tree_->Branch("Supercluster_eta", &supercluster_eta_);
        tree_->Branch("Supercluster_phi", &supercluster_phi_);
    }
    
    // tracks
    tree_->Branch("nTrack", &nTrack_,"nTrack_/s");
    tree_->Branch("Track_pt", &track_pt_);
    tree_->Branch("Track_eta", &track_eta_);
    tree_->Branch("Track_phi", &track_phi_);
    tree_->Branch("Track_charge", &track_charge_);
    
    tree_->Branch("Track_dxyPV", &track_dxyPV_);
    tree_->Branch("Track_dzPV", &track_dzPV_);       
    tree_->Branch("Track_dxyPVmin", &track_dxyPVmin_);  
    tree_->Branch("Track_dzPVmin", &track_dzPVmin_);  
    
    tree_->Branch("Track_dx", &track_dx_);
    tree_->Branch("Track_dy", &track_dy_);
    tree_->Branch("Track_dz", &track_dz_);
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

    lepton_genPt_  = 0.;
    lepton_genEta_ = 0.;
    lepton_genPhi_ = 0.;

    lepton_genVx_ = 99.;
    lepton_genVy_ = 99.;
    lepton_genVz_ = 99.;

    antiLepton_genPt_  = 0.;
    antiLepton_genEta_ = 0.;
    antiLepton_genPhi_ = 0.;

    antiLepton_genVx_ = 99.;
    antiLepton_genVy_ = 99.;
    antiLepton_genVz_ = 99.;

    nMuon_ = 0;
    muon_pt_.clear();
    muon_eta_.clear();
    muon_phi_.clear();    
    muon_charge_.clear();
    
    muon_ptTrk_.clear();
    muon_etaTrk_.clear();
    muon_phiTrk_.clear(); 
    muon_chargeTrk_.clear();
    
    muon_ptStaReg_.clear();
    muon_etaStaReg_.clear();
    muon_phiStaReg_.clear();    
    muon_chargeStaReg_.clear();

    muon_nStationsStaReg_.clear();

    muon_ptStaUpd_.clear();
    muon_etaStaUpd_.clear();
    muon_phiStaUpd_.clear();    
    muon_chargeStaUpd_.clear();
    
    muon_useUpdated_.clear();
    
    muon_dxyPV_.clear();
    muon_dzPV_.clear();
    muon_dxyPVmin_.clear();    
    muon_dzPVmin_.clear();
    
    muon_dx_.clear();
    muon_dy_.clear();
    muon_dz_.clear();
    muon_DxyError_.clear();
    muon_DzError_.clear();    
    
    muon_tkRelIso_.clear();
    muon_pfRelIso04_all_.clear();
    muon_ID_.clear();
    muon_triggerBits_.clear();
    
    muon_isMedium_.clear();
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

    nTrack_ = 0;
    track_pt_.clear();
    track_eta_.clear();
    track_phi_.clear();
    track_charge_.clear();
    
    track_dxyPV_.clear();
    track_dzPV_.clear();
    track_dxyPVmin_.clear();    
    track_dzPVmin_.clear();
    
    track_dx_.clear();
    track_dy_.clear();
    track_dz_.clear();
    track_DxyError_.clear();
    track_DzError_.clear();  

    track_nPixelHits_.clear();
    track_nTrackerLayers_.clear();
    track_validFraction_.clear();
    track_trackAlgo_.clear();

    track_ScaleCorr_.clear();

    nElectron_ = 0;
    electron_pt_.clear();
    electron_eta_.clear();
    electron_phi_.clear();    
    electron_charge_.clear();

    electron_dxyPV_.clear();
    electron_dzPV_.clear();
    electron_dxyPVmin_.clear();    
    electron_dzPVmin_.clear();
        
    electron_superclusterEnergy_.clear();
    electron_superclusterEta_.clear();
    electron_superclusterPhi_.clear();
    electron_dEtaInSeed_.clear();
    electron_full5x5_sigmaIetaIeta_.clear();
    electron_deltaPhiSuperClusterTrackAtVtx_.clear();
    electron_hadronicOverEm_.clear();
    electron_eSuperClusterOverP_.clear();
    electron_ecalEnergy_.clear();
    electron_gsfTrackMissingInnerHits_.clear();
    electron_pfRelIso_.clear();
    electron_hasMatchedConversion_.clear();
    
    electron_triggerBits_.clear();

    nSupercluster_ = 0;
    supercluster_energy_.clear();
    supercluster_eta_.clear();
    supercluster_phi_.clear(); 

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
// Look for closest distance in xy to any good primary vertex
float ZCountingAOD::getMinDxy(const reco::Track &trk, const std::vector<reco::Vertex> &pvCol){
    float min_dxy = 999;
    for (const reco::Vertex &itVtx : pvCol) {
        if(!isGoodPV(itVtx))
            continue;
        if(std::abs(trk.dxy(itVtx.position())) < min_dxy)
            min_dxy = std::abs(trk.dxy(itVtx.position()));
    }
    return min_dxy;
}
//--------------------------------------------------------------------------------------------------
// Look for closest distance in z to any good primary vertex
float ZCountingAOD::getMinDz(const reco::Track &trk, const std::vector<reco::Vertex> &pvCol){
    float min_dz = 999;
    for (const reco::Vertex &itVtx : pvCol) {
        if(!isGoodPV(itVtx))
            continue;
        if(std::abs(trk.dz(itVtx.position())) < min_dz)
            min_dz = std::abs(trk.dz(itVtx.position()));
    }
    return min_dz;
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

//--------------------------------------------------------------------------------------------------
float ZCountingAOD::getElectronIso(const reco::GsfElectron &el, const double rho){
    
    const reco::GsfElectron::PflowIsolationVariables& pfIso = el.pfIsolationVariables();
    const double chad = pfIso.sumChargedHadronPt;
    const double nhad = pfIso.sumNeutralHadronEt;
    const double pho = pfIso.sumPhotonEt;    
    const double eA = effectiveAreas_.getEffectiveArea(fabs(el.superCluster()->eta()));
        
    return chad + std::max(nhad + pho - rho * eA, 0.0);;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZCountingAOD);
