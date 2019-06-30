//
// Original Author:  David Walter
//         Created:  Sat, 29 Jun 2019 13:26:06 GMT
//
//

// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "ZCounting/ZCountAnalyze/plugins/ZCounting.h"

//
// constructors and destructor
//
ZCounting::ZCounting(const edm::ParameterSet& iConfig)
{

    zcount_PV* pvModule = new zcount_PV();
    pvModule->setPVToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("edmPVName")));

    addModule(pvModule);

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


    for(auto& m:modules_){
        m->readSetup(iSetup);
        if(!m->readEvent(iEvent)) return;
        m->fillBranches();
    }

    edm::LogInfo("ZCounting") << "test ";

    tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ZCounting::beginJob()
{
    if( !fs ){
        edm::LogError("ZCounting") << "TFile Service is not registered in cfg file" << std::endl;
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

