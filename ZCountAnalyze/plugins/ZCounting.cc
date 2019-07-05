//
// Original Author:  David Walter
//         Created:  Sat, 29 Jun 2019 13:26:06 GMT
//
//

// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// includes from ZCounting project
#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

#include "ZCounting/ZCountAnalyze/plugins/ZCounting.h"

//
// constructors and destructor
//
ZCounting::ZCounting(const edm::ParameterSet& iConfig)
{
    isData_ = iConfig.getUntrackedParameter<bool>("isData");

    zcount_trigger *triggerModule = new zcount_trigger();
    triggerModule->setTriggerResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults")));
    triggerModule->setTriggerEventToken(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent")));

    addModule(triggerModule);

    zcount_PV* pvModule = new zcount_PV();
    pvModule->setPVToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("edmPVName")));

    addModule(pvModule);

    zcount_muons* muonModule = new zcount_muons();
    muonModule->setMuonToken(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("edmMuonName")));
    muonModule->setTrackToken(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("edmTrackName")));
    muonModule->setPVModule(*pvModule);
    muonModule->setTriggerModule(*triggerModule);

    addModule(muonModule);

    addModule(new zcount_eventInfo());

    if(!isData_){
        zcount_genInfo* genModule = new zcount_genInfo();
        genModule->setGenInfoToken(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")));
        genModule->setGenZInfoToken(consumes<std::vector<GenZDecayProperties>>(iConfig.getParameter<edm::InputTag>("genZCollection")));

        addModule(genModule);
    }

    for(auto& m: modules_){
        m->getInput(iConfig);
    }
}


ZCounting::~ZCounting()
{
    for(auto& m:modules_)
        delete m;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZCounting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::LogInfo("ZCounting")<<"analyze";


    for(auto& m:modules_){
        m->readSetup(iSetup);
        if(!m->readEvent(iEvent)) return;
        m->fillBranches();
    }


    tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZCounting::beginJob()
{
    if( !fs ){
        edm::LogError("ZCounting") << "TFile Service is not registered in cfg file";
        return;
    }

    tree_=(fs->make<TTree>("tree" ,"tree" ));

    for(auto& m:modules_)
        m->initBranches(tree_);

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

