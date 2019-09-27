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
        - at least one good primary vertex
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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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

#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

// ROOT includes
#include "TTree.h"
#include "Math/Vector4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

//
// class declaration
//


class ZCounting : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit ZCounting(const edm::ParameterSet&);
    ~ZCounting();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void clearVariables();

    std::string get_triggerPath(std::string, const edm::TriggerNames&);
    std::vector<pat::TriggerObjectStandAlone> get_muonTriggerObjects(const std::vector<pat::TriggerObjectStandAlone> &,
                                                                     const edm::TriggerNames&,
                                                                     std::vector<std::string>
                                                                     );
    std::vector<pat::TriggerObjectStandAlone> get_muonTriggerObjects(const std::vector<pat::TriggerObjectStandAlone> &,
                                                                     const edm::TriggerNames&,
                                                                     std::string
                                                                     );
    bool isTriggerObject(const std::vector<pat::TriggerObjectStandAlone>&, const pat::Muon&);
    bool isGoodPV(const reco::Vertex&);
    bool isValidTrack(const reco::Track&);
    bool customIsTightMuon(const pat::Muon&);
    double pfIso(const pat::Muon&);
    double tkIso(const pat::Muon&);
    int getMuonID(const pat::Muon&, const reco::Vertex&);

    // ----------member data ---------------------------

    edm::Service<TFileService> fs;
    TTree *tree_;

    const double DRMAX = 0.2;

    // --- input
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
    edm::EDGetTokenT<pat::MuonCollection> muonCollection_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > pvCollection_;
    edm::EDGetTokenT<std::vector<GenZDecayProperties> > genZCollection_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoCollection_;

    std::vector<std::string> muonTriggerPatterns_;

    // primary vertex cuts
    double VtxNTracksFitCut_;
    double VtxNdofCut_;
    double VtxAbsZCut_;
    double VtxRhoCut_;

    // --- output
    int nPU_;
    int nPV_;

    int muon_recoMatches_;
    bool muon_hasRecoObj_;
    float muon_tkIso_;
    float muon_pfIso_;
    int muon_ID_;
    int muon_triggerBits_;
    float muon_recoPt_;
    float muon_recoEta_;
    float muon_recoPhi_;
    float muon_genPt_;
    float muon_genEta_;
    float muon_genPhi_;

    int antiMuon_recoMatches_;
    bool antiMuon_hasRecoObj_;
    float antiMuon_tkIso_;
    float antiMuon_pfIso_;
    int antiMuon_ID_;
    int antiMuon_triggerBits_;
    float antiMuon_recoPt_;
    float antiMuon_recoEta_;
    float antiMuon_recoPhi_;
    float antiMuon_genPt_;
    float antiMuon_genEta_;
    float antiMuon_genPhi_;

    float z_genMass_;
    float z_recoMass_;
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
    pileupInfoCollection_  (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryInfoCollection")))
{
    edm::LogInfo("ZCounting")<<"ZCounting(...)";

    muonTriggerPatterns_ = iConfig.getParameter<std::vector<std::string>>("muon_trigger_patterns");

    VtxNTracksFitCut_ = iConfig.getUntrackedParameter<double>("VtxNTracksFitMin");
    VtxNdofCut_       = iConfig.getUntrackedParameter<double>("VtxNdofMin");
    VtxAbsZCut_       = iConfig.getUntrackedParameter<double>("VtxAbsZMax");
    VtxRhoCut_        = iConfig.getUntrackedParameter<double>("VtxRhoMax");

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
    edm::LogInfo("ZCounting")<<"analyze(...)";

    this->clearVariables();

    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
    edm::Handle<pat::MuonCollection> muonCollection;
    edm::Handle<std::vector<reco::Vertex> > pvCollection;
    edm::Handle<std::vector<GenZDecayProperties> > genZCollection;
    edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoCollection;

    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(muonCollection_, muonCollection);
    iEvent.getByToken(pvCollection_, pvCollection);
    iEvent.getByToken(genZCollection_, genZCollection);
    iEvent.getByToken(pileupInfoCollection_, pileupInfoCollection);

    // only keep Z -> mu mu
    if(!genZCollection->size() || genZCollection->at(0).decayMode() != 13) return;

    // PV selection
    const std::vector<reco::Vertex> *pvCol = pvCollection.product();
    reco::Vertex pv = *pvCol->begin();
    nPV_ = 0;
    for (auto const& itVtx : *pvCollection) {
        if(!isGoodPV(itVtx))
            continue;
        if (nPV_ == 0)
            pv = itVtx;
        nPV_++;
    }

    edm::LogVerbatim("ZCounting") << "get trigger objects";
    const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);

    // --- gen muons from Z
    edm::LogVerbatim("ZCounting") << "get gen muons";
    const reco::GenParticle *ZLepton = genZCollection->at(0).stableLepton();
    const reco::GenParticle *ZAntiLepton = genZCollection->at(0).stableAntiLepton();

    if(!ZLepton || ! ZAntiLepton) return;

    muon_genPt_ = ZLepton->pt();
    muon_genEta_ = ZLepton->eta();
    muon_genPhi_ = ZLepton->phi();

    antiMuon_genPt_ = ZAntiLepton->pt();
    antiMuon_genEta_ = ZAntiLepton->eta();
    antiMuon_genPhi_ = ZAntiLepton->phi();

    z_genMass_ = (ZLepton->p4() + ZAntiLepton->p4()).M();

    // --- find reco muons corresponding to the gen muons from Z
    edm::LogVerbatim("ZCounting") << "find reco muons";
    LorentzVector muon_reco;
    LorentzVector antiMuon_reco;
    for (pat::Muon mu : *muonCollection){

        if(reco::deltaR(mu.eta(), mu.phi(), ZLepton->eta(), ZLepton->phi()) < 0.03 && mu.pdgId() == 13){
            muon_recoMatches_++;
            muon_reco = mu.p4();
            muon_hasRecoObj_ = true;
            muon_recoPt_ = mu.pt();
            muon_recoEta_ = mu.eta();
            muon_recoPhi_ = mu.phi();
            muon_ID_ = getMuonID(mu, pv);
            muon_tkIso_ = tkIso(mu);
            muon_pfIso_ = pfIso(mu);

            for(unsigned j = 0, m = muonTriggerPatterns_.size(); j < m; ++j){
                std::vector<pat::TriggerObjectStandAlone> tObjCol = get_muonTriggerObjects(*triggerObjects, triggerNames, muonTriggerPatterns_.at(j));
                muon_triggerBits_ += std::pow(2,j) * isTriggerObject(tObjCol, mu);
            }

        }
        if(reco::deltaR(mu.eta(), mu.phi(), ZAntiLepton->eta(), ZAntiLepton->phi()) < 0.03 && mu.pdgId() == -13){
            antiMuon_recoMatches_++;
            antiMuon_reco = mu.p4();
            antiMuon_hasRecoObj_ = true;
            antiMuon_recoPt_ = mu.pt();
            antiMuon_recoEta_ = mu.eta();
            antiMuon_recoPhi_ = mu.phi();
            antiMuon_ID_ = getMuonID(mu, pv);
            antiMuon_tkIso_ = tkIso(mu);
            antiMuon_pfIso_ = pfIso(mu);

            for(unsigned j = 0, m = muonTriggerPatterns_.size(); j < m; ++j){
                std::vector<pat::TriggerObjectStandAlone> tObjCol = get_muonTriggerObjects(*triggerObjects, triggerNames, muonTriggerPatterns_.at(j));
                antiMuon_triggerBits_ += std::pow(2,j) * isTriggerObject(tObjCol, mu);
            }
        }
    }

    if(muon_hasRecoObj_ && antiMuon_hasRecoObj_)
        z_recoMass_ = (muon_reco + antiMuon_reco).M();

    for(std::vector<PileupSummaryInfo>::const_iterator puI = pileupInfoCollection->begin(); puI != pileupInfoCollection->end(); ++puI)
    {
        if(puI->getBunchCrossing() == 0) nPU_ = puI->getTrueNumInteractions();
    }

    tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZCounting::beginJob()
{
    edm::LogInfo("ZCounting")<<"beginJob()";

    if( !fs ){
        edm::LogError("ZCounting") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    tree_->Branch("nPV", &nPV_,"nPV_/i");
    tree_->Branch("nPU", &nPU_,"nPU_/i");

    tree_->Branch("muon_recoMatches", &muon_recoMatches_,"muon_recoMatches_/i");
    tree_->Branch("muon_hasRecoObj", &muon_hasRecoObj_,"muon_hasRecoObj_/b");
    tree_->Branch("muon_ID", &muon_ID_,"muon_ID_/i");
    tree_->Branch("muon_tkIso", &muon_tkIso_,"muon_tkIso_/f");
    tree_->Branch("muon_pfIso", &muon_pfIso_,"muon_pfIso_/f");
    tree_->Branch("muon_triggerBits", &muon_triggerBits_,"muon_triggerBits_/i");
    tree_->Branch("muon_recoPt", &muon_recoPt_,"muon_recoPt_/f");
    tree_->Branch("muon_recoEta", &muon_recoEta_,"muon_recoEta_/f");
    tree_->Branch("muon_recoPhi", &muon_recoPhi_,"muon_recoPhi_/f");
    tree_->Branch("muon_genPt", &muon_genPt_,"muon_genPt_/f");
    tree_->Branch("muon_genEta", &muon_genEta_,"muon_genEta_/f");
    tree_->Branch("muon_genPhi", &muon_genPhi_,"muon_genPhi_/f");

    tree_->Branch("antiMuon_recoMatches", &antiMuon_recoMatches_,"antiMuon_recoMatches_/i");
    tree_->Branch("antiMuon_hasRecoObj", &antiMuon_hasRecoObj_,"antiMuon_hasRecoObj_/b");
    tree_->Branch("antiMuon_ID", &antiMuon_ID_,"antiMuon_ID_/i");
    tree_->Branch("antiMuon_tkIso", &antiMuon_tkIso_,"antiMuon_tkIso_/f");
    tree_->Branch("antiMuon_pfIso", &antiMuon_pfIso_,"antiMuon_pfIso_/f");
    tree_->Branch("antiMuon_triggerBits", &antiMuon_triggerBits_,"antiMuon_triggerBits_/i");
    tree_->Branch("antiMuon_recoPt", &antiMuon_recoPt_,"antiMuon_recoPt_/f");
    tree_->Branch("antiMuon_recoEta", &antiMuon_recoEta_,"antiMuon_recoEta_/f");
    tree_->Branch("antiMuon_recoPhi", &antiMuon_recoPhi_,"antiMuon_recoPhi_/f");
    tree_->Branch("antiMuon_genPt", &antiMuon_genPt_,"antiMuon_genPt_/f");
    tree_->Branch("antiMuon_genEta", &antiMuon_genEta_,"antiMuon_genEta_/f");
    tree_->Branch("antiMuon_genPhi", &antiMuon_genPhi_,"antiMuon_genPhi_/f");

    tree_->Branch("z_genMass", &z_genMass_,"z_genMass_/f");
    tree_->Branch("z_recoMass", &z_recoMass_,"z_recoMass_/f");

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
    edm::LogInfo("ZCounting")<<"clearVariables()";

    nPV_ = 0;
    nPU_ = 0;

    muon_recoMatches_ = 0;
    muon_hasRecoObj_ = false;
    muon_ID_ = 0;
    muon_tkIso_ = 0.;
    muon_pfIso_ = 0.;
    muon_triggerBits_ = 0;
    muon_recoPt_ = 0.;
    muon_recoEta_ = 0.;
    muon_recoPhi_ = 0.;
    muon_genPt_ = 0.;
    muon_genEta_ = 0.;
    muon_genPhi_ = 0.;

    antiMuon_recoMatches_ = 0;
    antiMuon_hasRecoObj_ = false;
    antiMuon_ID_ = 0;
    antiMuon_tkIso_ = 0.;
    antiMuon_pfIso_ = 0.;
    antiMuon_triggerBits_ = 0;
    antiMuon_recoPt_ = 0.;
    antiMuon_recoEta_ = 0.;
    antiMuon_recoPhi_ = 0.;
    antiMuon_genPt_ = 0.;
    antiMuon_genEta_ = 0.;
    antiMuon_genPhi_ = 0.;

    z_genMass_ = 0.;
    z_recoMass_ = 0.;
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
                                  std::string pattern
                                  ){

    std::vector<pat::TriggerObjectStandAlone> muonTriggerObjects;
    for (pat::TriggerObjectStandAlone obj : tObjCol) {
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesAll = obj.pathNames(false);
        std::string path = get_triggerPath(pattern, names);
        bool passedTrigger = false;
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            if(pathNamesAll[h] == path)
                passedTrigger = true;
        }
        if(!passedTrigger) continue;
        muonTriggerObjects.push_back(obj);
    }

    return muonTriggerObjects;
}

//--------------------------------------------------------------------------------------------------
std::vector<pat::TriggerObjectStandAlone>
ZCounting::get_muonTriggerObjects(const std::vector<pat::TriggerObjectStandAlone> & tObjCol,
                                  const edm::TriggerNames& names,
                                  std::vector<std::string> patterns
                                  ){

    std::vector<pat::TriggerObjectStandAlone> muonTriggerObjects;
    for (pat::TriggerObjectStandAlone obj : tObjCol) {
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesAll = obj.pathNames(false);
        for(unsigned j = 0, m = patterns.size(); j < m; ++j){
            std::string path = get_triggerPath(patterns.at(j), names);
            bool passedTrigger = false;
            for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
                if(pathNamesAll[h] == path)
                    passedTrigger = true;
            }
            if(!passedTrigger) continue;
            muonTriggerObjects.push_back(obj);
        }
    }

    return muonTriggerObjects;
}

//--------------------------------------------------------------------------------------------------
bool ZCounting::isTriggerObject(const std::vector<pat::TriggerObjectStandAlone> &tObjCol, const pat::Muon &muon){
    for(pat::TriggerObjectStandAlone tObj : tObjCol) {
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



//define this as a plug-in
DEFINE_FWK_MODULE(ZCounting);
