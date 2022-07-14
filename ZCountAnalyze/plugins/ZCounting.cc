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
#include <algorithm>
#include <iterator>

// CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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
#include "ZCounting/ZUtils/interface/Helper.h"

// ROOT includes
#include "TTree.h"
#include "Math/Vector4D.h"
#include <TLorentzVector.h>

//
// class declaration
//


class ZCounting : public edm::one::EDAnalyzer<>
{
public:
    explicit ZCounting(const edm::ParameterSet&);
    ~ZCounting() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    void beginRun(const edm::Run&, const edm::EventSetup&);
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

    int getMuonID(const pat::Muon&, const reco::Vertex&);
    double dxy(const pat::Muon&, const reco::Vertex&);
    double dz(const pat::Muon&, const reco::Vertex&);


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
    edm::EDGetTokenT<std::vector<pat::MET> >  metPFCollection_;
    edm::EDGetTokenT<std::vector<pat::MET> >  metPuppiCollection_;
    edm::EDGetTokenT<std::vector<GenZDecayProperties> > genZCollection_;
    edm::EDGetTokenT<TtGenEvent> genTtEvent_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoCollection_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;
    edm::EDGetTokenT<LHEEventProduct> lheEventInfo_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleCollection_;

    edm::EDGetTokenT< double > prefweightECAL_token;
    edm::EDGetTokenT< double > prefweightupECAL_token;
    edm::EDGetTokenT< double > prefweightdownECAL_token;

    edm::EDGetTokenT< double > prefweightECAL2017H_token;
    edm::EDGetTokenT< double > prefweightupECAL2017H_token;
    edm::EDGetTokenT< double > prefweightdownECAL2017H_token;

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

    // Muon trigger
    bool hltChanged_;
    std::vector<std::string> muonTriggerPatterns_;
    std::vector<std::string> muonTriggerPaths_;
    // MET trigger
    std::vector<std::string> metTriggerPatterns_;
    std::vector<std::string> metTriggerPatternsExt_;
    std::vector<std::string> metTriggerPaths_;
    std::vector<std::string> metTriggerPathsExt_;

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

    int met_triggerBits_;
    int met_triggerBitsExt_;

    float met_PF_pt_;
    float met_PF_phi_;
    float met_Puppi_pt_;
    float met_Puppi_phi_;

    std::vector<int> v_genWeightIDs_;
    std::vector<int> v_pdfWeightIDs_;

    float eventweight_;
    std::vector<float> v_psWeight_;
    std::vector<float> v_meWeight_;
    std::vector<float> v_pdfWeight_;

    float prefiringweightECAL_;
    float prefiringweightECALup_;
    float prefiringweightECALdown_;

    float prefiringweightECAL2017H_;
    float prefiringweightECALup2017H_;
    float prefiringweightECALdown2017H_;

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
    std::vector<int> muon_matchValue_;
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

    std::vector<float> v_muon_ScaleCorr_;
    std::vector<float> v_muon_ScaleCorr_stat_RMS_;
    std::vector<float> v_muon_ScaleCorr_Zpt_;
    std::vector<float> v_muon_ScaleCorr_Ewk_;
    std::vector<float> v_muon_ScaleCorr_deltaM_;
    std::vector<float> v_muon_ScaleCorr_Ewk2_;
    std::vector<float> v_muon_ScaleCorr_Total_;
};

//
// constructors and destructor
//
ZCounting::ZCounting(const edm::ParameterSet& iConfig):
    triggerBits_     (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigger_bits"))),
    triggerObjects_  (consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("trigger_objects"))),
    muonCollection_  (consumes<pat::MuonCollection> (iConfig.getParameter<edm::InputTag>("pat_muons"))),
    pvCollection_    (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("edmPVName"))),
    metPFCollection_   (consumes<std::vector<pat::MET> >  (iConfig.getParameter<edm::InputTag>("pat_met_PF"))),
    metPuppiCollection_   (consumes<std::vector<pat::MET> >  (iConfig.getParameter<edm::InputTag>("pat_met_puppi"))),
    genZCollection_  (consumes<std::vector<GenZDecayProperties> > (iConfig.getParameter<edm::InputTag>("genZLeptonCollection"))),
    genTtEvent_  (consumes<TtGenEvent> (iConfig.getParameter<edm::InputTag>("genTtCollection"))),
    pileupInfoCollection_  (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfoCollection"))),
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfo_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    genParticleCollection_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles")))
{
    LogDebug("ZCounting")<<"ZCounting(...)";

    era_ = iConfig.getParameter<std::string>("era");

    isData_ = iConfig.getUntrackedParameter<bool>("isData");
    hasGenZ_ = iConfig.getUntrackedParameter<bool>("hasGenZ");
    hasGenTt_ = iConfig.getUntrackedParameter<bool>("hasGenTt");

    v_genWeightIDs_ = iConfig.getParameter<std::vector<int>>("genWeights");
    v_pdfWeightIDs_ = iConfig.getParameter<std::vector<int>>("pdfWeights");

    metTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("met_trigger_patterns");
    metTriggerPatternsExt_ = iConfig.getParameter<std::vector<std::string>>("met_trigger_patterns_ext");

    muonTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("muon_trigger_patterns");
    DRMAX = iConfig.getUntrackedParameter<double>("muon_trigger_DRMAX");

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
    else if(era_ == "2017"){
        std::cout<<"2017 low PU"<<std::endl;

        prefweightECAL2017H_token = consumes< double >(edm::InputTag("prefiringweight2017H:nonPrefiringProbECAL"));
        prefweightupECAL2017H_token = consumes< double >(edm::InputTag("prefiringweight2017H:nonPrefiringProbECALUp"));
        prefweightdownECAL2017H_token = consumes< double >(edm::InputTag("prefiringweight2017H:nonPrefiringProbECALDown"));
    }


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
    prefiringweightMuonupSyst_ =(*theprefweightupsystMuon);
    prefiringweightMuondownSyst_ =(*theprefweightdownsystMuon);
    prefiringweightMuonupStat_ =(*theprefweightupstatMuon);
    prefiringweightMuondownStat_ =(*theprefweightdownstatMuon);

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
    else if(era_ == "2017"){
        edm::Handle< double > theprefweightECAL2017H;
        edm::Handle< double > theprefweightupECAL2017H;
        edm::Handle< double > theprefweightdownECAL2017H;

        iEvent.getByToken(prefweightECAL2017H_token, theprefweightECAL2017H);
        iEvent.getByToken(prefweightupECAL2017H_token, theprefweightupECAL2017H);
        iEvent.getByToken(prefweightdownECAL2017H_token, theprefweightdownECAL2017H);

        prefiringweightECAL2017H_ =(*theprefweightECAL2017H);
        prefiringweightECALup2017H_ =(*theprefweightupECAL2017H);
        prefiringweightECALdown2017H_ =(*theprefweightdownECAL2017H);
    }
    // <<< ---

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
                throw cms::Exception("LogicError") << "@@@ ZCounting::analyze -- "
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
            LogDebug("ZCounting::analyze")<<"hasGenZ";

            edm::Handle<TtGenEvent> genTtEvent;
            iEvent.getByToken(genTtEvent_, genTtEvent);
            if(genTtEvent->lepton()){
                LogDebug("ZCounting::analyze")<<"hasLepton";
                genLepton = genTtEvent->lepton();
                if(isTau(genLepton)){
                    decayMode_ += 150000;
                    genLepton = tauDaughter(genTtEvent->lepton());
                }

            }
            if(genTtEvent->leptonBar()){
                LogDebug("ZCounting::analyze")<<"hasAntiLepton";
                genAntiLepton = genTtEvent->leptonBar();
                if(isTau(genAntiLepton)){
                    decayMode_ += 150000;
                    genAntiLepton = tauDaughter(genTtEvent->leptonBar());
                }
            }
            LogDebug("ZCounting::analyze")<<"setDecayMode";
            if(genLepton)
                decayMode_ += 100*std::abs(genLepton->pdgId());
            if(genAntiLepton)
                decayMode_ += std::abs(genAntiLepton->pdgId());
        }
    }

    // --- get trigger names
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    LogDebug("ZCounting") << "get trigger names";
    const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
    if(hltChanged_ == true){
        muonTriggerPaths_.clear();
        metTriggerPaths_.clear();
        metTriggerPathsExt_.clear();
        for(const std::string pattern: muonTriggerPatterns_){
            muonTriggerPaths_.push_back(get_triggerPath(pattern, triggerNames));
        }
        for(const std::string pattern: metTriggerPatterns_){
            metTriggerPaths_.push_back(get_triggerPath(pattern, triggerNames));
        }
        for(const std::string pattern: metTriggerPatternsExt_){
            metTriggerPathsExt_.push_back(get_triggerPath(pattern, triggerNames));
        }
        hltChanged_=false;
    }

    // --- check if met trigger has fired
    const unsigned int nTrigger = triggerBits->size();
    // loop over all trigger of the event
    for(unsigned int iTrigger = 0; iTrigger < nTrigger; ++iTrigger){
        // check if this trigger has fired
        if(!(triggerBits.product()->accept(iTrigger))) continue;

        const std::string& triggerName = triggerNames.triggerName(iTrigger);

        auto it = std::find(metTriggerPaths_.begin(), metTriggerPaths_.end(), triggerName);
        if(it != metTriggerPaths_.end()){
            // found trigger
            met_triggerBits_ += std::pow(2, std::distance(metTriggerPaths_.begin(), it));
        }

        it = std::find(metTriggerPathsExt_.begin(), metTriggerPathsExt_.end(), triggerName);
        if(it != metTriggerPathsExt_.end()){
            // found trigger
            met_triggerBitsExt_ += std::pow(2, std::distance(metTriggerPathsExt_.begin(), it));
        }

    }

    // --- PF MET
    edm::Handle<std::vector<pat::MET> > metPFCollection;
    iEvent.getByToken(metPFCollection_, metPFCollection);

    met_PF_pt_ = metPFCollection->at(0).shiftedP4(pat::MET::NoShift, pat::MET::Type1).Pt();
    met_PF_phi_= metPFCollection->at(0).shiftedP4(pat::MET::NoShift, pat::MET::Type1).Phi();

    // --- Puppi MET
    edm::Handle<std::vector<pat::MET> > metPuppiCollection;
    iEvent.getByToken(metPuppiCollection_, metPuppiCollection);

    met_Puppi_pt_ = metPuppiCollection->at(0).shiftedP4(pat::MET::NoShift, pat::MET::Type1).Pt();
    met_Puppi_phi_= metPuppiCollection->at(0).shiftedP4(pat::MET::NoShift, pat::MET::Type1).Phi();

    // --- PV selection
    edm::Handle<std::vector<reco::Vertex> > pvCollection;
    iEvent.getByToken(pvCollection_, pvCollection);
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
    edm::Handle<pat::MuonCollection> muonCollection;
    iEvent.getByToken(muonCollection_, muonCollection);
    for (pat::Muon mu : *muonCollection){
        if(std::abs(mu.eta()) > 2.4) continue;
        if(mu.pt() < 15) continue;

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

        if(mu.innerTrack().isNonnull()){
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

        muon_tkRelIso_.push_back(getTkIso(mu));
        muon_pfRelIso04_all_.push_back(getPFIso(mu));
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

            // store muon Rochester corrections and uncertainties
            v_muon_ScaleCorr_.push_back(mu.hasUserFloat("MuonEnergyCorr")                   ? mu.userFloat("MuonEnergyCorr")          : 1.0);
            v_muon_ScaleCorr_stat_RMS_.push_back(mu.hasUserFloat("MuonEnergyCorr_stat_RMS") ? mu.userFloat("MuonEnergyCorr_stat_RMS") : 0.0);
            v_muon_ScaleCorr_Zpt_.push_back(mu.hasUserFloat("MuonEnergyCorr_Zpt")           ? mu.userFloat("MuonEnergyCorr_Zpt")      : 0.0);
            v_muon_ScaleCorr_Ewk_.push_back(mu.hasUserFloat("MuonEnergyCorr_Ewk")           ? mu.userFloat("MuonEnergyCorr_Ewk")      : 0.0);
            v_muon_ScaleCorr_deltaM_.push_back(mu.hasUserFloat("MuonEnergyCorr_deltaM")     ? mu.userFloat("MuonEnergyCorr_deltaM")   : 0.0);
            v_muon_ScaleCorr_Ewk2_.push_back(mu.hasUserFloat("MuonEnergyCorr_Ewk2")         ? mu.userFloat("MuonEnergyCorr_Ewk2")     : 0.0);
            v_muon_ScaleCorr_Total_.push_back(mu.hasUserFloat("MuonEnergyCorr_Total")       ? mu.userFloat("MuonEnergyCorr_Total")    : 0.0);
        }

        nMuon_++;

    }
    
    // we don't need data events with less than two muons
    if(iEvent.isRealData() && nMuon_ < 2)
        return;

    if(antiMuon_genRecoObj_ != -1 && muon_genRecoObj_ != -1){
        vMuon.SetPtEtaPhiM(muon_pt_[muon_genRecoObj_], muon_eta_[muon_genRecoObj_], muon_phi_[muon_genRecoObj_], MUON_MASS);
        vAntiMuon.SetPtEtaPhiM(muon_pt_[antiMuon_genRecoObj_], muon_eta_[antiMuon_genRecoObj_], muon_phi_[antiMuon_genRecoObj_], MUON_MASS);

        z_recoMass_ = (vMuon + vAntiMuon).M();
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

    tree_->Branch("MET_triggerBits",    &met_triggerBits_,      "met_triggerBits_/I");
    tree_->Branch("MET_triggerBitsExt", &met_triggerBitsExt_,   "met_triggerBitsExt_/I");
    tree_->Branch("MET_PF_pt",    &met_PF_pt_,  "met_PF_pt_/f");
    tree_->Branch("MET_PF_phi",   &met_PF_phi_, "met_PF_phi_/f");
    tree_->Branch("MET_Puppi_pt",    &met_Puppi_pt_,  "met_Puppi_pt_/f");
    tree_->Branch("MET_Puppi_phi",   &met_Puppi_phi_, "met_Puppi_phi_/f");

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
    else if(era_ == "2017"){
        tree_->Branch("prefiringweightECAL2017H",         &prefiringweightECAL2017H_,         "prefiringweightECAL2017H_/f");
        tree_->Branch("prefiringweightECALup2017H",       &prefiringweightECALup2017H_,       "prefiringweightECALup2017H_/f");
        tree_->Branch("prefiringweightECALdown2017H",     &prefiringweightECALdown2017H_,     "prefiringweightECALdown2017H_/f");
    }

    if(!isData_){
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
    }

    // reco level info
    tree_->Branch("nMuon", &nMuon_,"nMuon_/s");
    tree_->Branch("Muon_matchValue", &muon_matchValue_);
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

    if(!isData_){
        //Roccester corrections
        tree_->Branch("Muon_ScaleCorr",          &v_muon_ScaleCorr_);
        tree_->Branch("Muon_ScaleCorr_stat_RMS", &v_muon_ScaleCorr_stat_RMS_);
        tree_->Branch("Muon_ScaleCorr_Zpt",      &v_muon_ScaleCorr_Zpt_);
        tree_->Branch("Muon_ScaleCorr_Ewk",      &v_muon_ScaleCorr_Ewk_);
        tree_->Branch("Muon_ScaleCorr_deltaM",   &v_muon_ScaleCorr_deltaM_);
        tree_->Branch("Muon_ScaleCorr_Ewk2",     &v_muon_ScaleCorr_Ewk2_);
        tree_->Branch("Muon_ScaleCorr_Total",    &v_muon_ScaleCorr_Total_);
    }
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

    met_triggerBits_ = 0;
    met_triggerBitsExt_ = 0;

    met_PF_pt_ = 0;
    met_PF_phi_ = 0;
    met_Puppi_pt_ = 0;
    met_Puppi_phi_ = 0;

    eventweight_ = 0.;
    v_psWeight_.clear();
    v_meWeight_.clear();
    v_pdfWeight_.clear();

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
    muon_matchValue_.clear();
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

    v_muon_ScaleCorr_.clear();
    v_muon_ScaleCorr_stat_RMS_.clear();
    v_muon_ScaleCorr_Zpt_.clear();
    v_muon_ScaleCorr_Ewk_.clear();
    v_muon_ScaleCorr_deltaM_.clear();
    v_muon_ScaleCorr_Ewk2_.clear();
    v_muon_ScaleCorr_Total_.clear();
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
        // default: pathNames(pathLastFilterAccept=false, pathL3FilterAccept=true)
        // this means also trigger objects that do not pass the last filter are obtained
        // but we want only last filters
        const std::vector<std::string> pathNamesAll = obj.pathNames(true, true);

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
int ZCounting::getMuonID(const pat::Muon &mu, const reco::Vertex &vtx){
    if(mu.isTightMuon(vtx)) return 5;
    if(isCustomTightMuon(mu)) return 4;
    if(mu.isGlobalMuon()) return 3;
    if(mu.isStandAloneMuon()) return 2;
    if(isValidTrack(*(mu.innerTrack()))) return 1;
    return 0;
}

//--------------------------------------------------------------------------------------------------
double ZCounting::dxy(const pat::Muon &mu, const reco::Vertex &vtx){
    return fabs(mu.muonBestTrack()->dxy(vtx.position()));
}

//--------------------------------------------------------------------------------------------------
double ZCounting::dz(const pat::Muon &mu, const reco::Vertex &vtx){
    return fabs(mu.muonBestTrack()->dz(vtx.position()));
}



//define this as a plug-in
DEFINE_FWK_MODULE(ZCounting);
