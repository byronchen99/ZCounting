#include <iostream>
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <TH1.h>
#include <TH1I.h>
#include <TH1D.h>

class CountEventAnalyzer : public edm::one::EDAnalyzer<>
{

public:

    explicit CountEventAnalyzer(const edm::ParameterSet&);
    virtual ~CountEventAnalyzer() {};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

private:
    edm::EDGetTokenT<GenEventInfoProduct> genEventInfo_;

    TH1* unweightedEvents_;
    TH1* weightedEvents_;
};

CountEventAnalyzer::CountEventAnalyzer(const edm::ParameterSet& iConfig):
    genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo")))
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

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CountEventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CountEventAnalyzer);


