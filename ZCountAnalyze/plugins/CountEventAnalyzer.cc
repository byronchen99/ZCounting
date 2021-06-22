#include <iostream>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>
#include <TH1I.h>
#include <TH1D.h>

class CountEventAnalyzer : public edm::EDAnalyzer
{

public:

    CountEventAnalyzer(const edm::ParameterSet&);
    ~CountEventAnalyzer();

protected:
    void beginJob();
    void analyze(const edm::Event&, const edm::EventSetup&);
    void endJob();

private:
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;

    TH1* unweightedEvents_;
    TH1* weightedEvents_;
};

CountEventAnalyzer::CountEventAnalyzer(const edm::ParameterSet& iConfig):
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")))
{}

CountEventAnalyzer::~CountEventAnalyzer()
{}

void CountEventAnalyzer::beginJob()
{
    edm::Service<TFileService> fs;
    if(!fs){
        throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
    }
    unweightedEvents_ = fs->make<TH1I>("unweightedEvents", "N of Events", 1, 0.5, 1.5);
    weightedEvents_ = fs->make<TH1D>("weightedEvents", "N of weighted events", 1, 0.5, 1.5);

}

void CountEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    unweightedEvents_->Fill(1);
    if(iEvent.isRealData()) return;

    edm::Handle<GenEventInfoProduct> evt_info;
    iEvent.getByToken(genEventInfo_, evt_info);

    weightedEvents_->Fill(1, evt_info->weight());

}

void CountEventAnalyzer::endJob()
{}

DEFINE_FWK_MODULE(CountEventAnalyzer);
