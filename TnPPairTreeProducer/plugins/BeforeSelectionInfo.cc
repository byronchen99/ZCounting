// CMSSW includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "ZCounting/TnPPairTreeProducer/plugins/BeforeSelectionInfo.h"


BeforeSelectionInfo::BeforeSelectionInfo(const edm::ParameterSet& iConfig):
    // genEventInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    src_pileupInfo_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("src")))
{
  // must state that we are using the TFileService
  // usesResource("TFileService");
}

void BeforeSelectionInfo::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
    edm::Handle<std::vector<PileupSummaryInfo> > puInfos;
    iEvent.getByToken(src_pileupInfo_, puInfos);

    h_LS_->Fill(iEvent.id().luminosityBlock());

    return;
}

void BeforeSelectionInfo::beginJob()
{
    edm::Service<TFileService> fs;
    if(!fs){
        throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
    }


}

void BeforeSelectionInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    const unsigned int run_ = iRun.run();

    std::cout<<"BeforeSelectionInfo::beginRun now at run "<<run_<< std::endl;
    edm::Service<TFileService> fs;
    if(!fs)
    {
    throw edm::Exception(edm::errors::Configuration, "no TFileService registered in cfg file");
    }

    h_LS_ = fs->make<TH1I>(
        ("LS_run"+std::to_string(run_)).c_str(),
        ("LS_run"+std::to_string(run_)).c_str(), 10000, 0, 10000);

    // h_TrueNInt_outOfTimeEarly_ = fs->make<TH1D>("MC_TrueNInt_outOfTimeEarly_run"+std::to_string(run_), "MC_TrueNInt_outOfTimeEarly_run"+std::to_string(iRun.id()), 120, 0., 120.);
    // h_LS_vs_PU_ = fs->make<TH2D>(
    //     ("LS_vs_PU_run"+std::to_string(run_)).c_str(),
    //     ("LS_vs_PU_run"+std::to_string(run_)).c_str(), 10000, 0, 10000, 120, 0., 120.);

    // h_TrueNInt_outOfTimeLate_  = fs->make<TH1D>("MC_TrueNInt_outOfTimeLate_run"+std::to_string(run_),  "MC_TrueNInt_outOfTimeLate_run"+std::to_string(iRun.id()),  120, 0., 120.);
}

void BeforeSelectionInfo::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    std::cout<<"BeforeSelectionInfo::endRun now at run"<<iRun.run()<< std::endl;

    // h_TrueNInt_outOfTimeEarly_->Delete();
    h_LS_->Write();
    h_LS_->Delete();
    // h_LS_vs_PU_->Write();
    // h_LS_vs_PU_->Delete();
    // h_TrueNInt_outOfTimeLate_->Delete();
}
