#include <iostream>
#include <FWCore/Framework/interface/one/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include <TH1D.h>

//
// Multi-Thread EDAnalyzer using TFileService
// REF: https://twiki.cern.ch/twiki/bin/view/CMSPublic/FWMultithreadedAnalysisEDAnalyzer
//
class PileupMCTemplateMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
 public:
  explicit PileupMCTemplateMaker(const edm::ParameterSet&);
  virtual ~PileupMCTemplateMaker() override = default;

 protected:
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > src_pileupInfo_;

  TH1D* h_TrueNInt_outOfTimeEarly_;
  TH1D* h_TrueNInt_inTime_;
  TH1D* h_TrueNInt_outOfTimeLate_;

 private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
};

PileupMCTemplateMaker::PileupMCTemplateMaker(const edm::ParameterSet& iConfig)
{
  // must state that we are using the TFileService
  usesResource("TFileService");

  src_pileupInfo_ = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("src"));
}

void PileupMCTemplateMaker::beginJob()
{
  edm::Service<TFileService> fs;
  if(!fs)
  {
    throw edm::Exception(edm::errors::Configuration, "no TFileService registered in cfg file");
  }

  h_TrueNInt_outOfTimeEarly_ = fs->make<TH1D>("MC_TrueNInt_outOfTimeEarly", "MC_TrueNInt_outOfTimeEarly", 120, 0., 120.);
  h_TrueNInt_inTime_         = fs->make<TH1D>("MC_TrueNInt_inTime",         "MC_TrueNInt_inTime",         120, 0., 120.);
  h_TrueNInt_outOfTimeLate_  = fs->make<TH1D>("MC_TrueNInt_outOfTimeLate",  "MC_TrueNInt_outOfTimeLate",  120, 0., 120.);
}

void PileupMCTemplateMaker::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
  // MC True Num. of Interactions
  if(iEvent.isRealData() == false)
  {
    edm::Handle<std::vector<PileupSummaryInfo> > puInfos;
    iEvent.getByToken(src_pileupInfo_, puInfos);

    for(std::vector<PileupSummaryInfo>::const_iterator puI = puInfos->begin(); puI != puInfos->end(); ++puI)
    {
      if     (puI->getBunchCrossing() <  0){ h_TrueNInt_outOfTimeEarly_->Fill(puI->getTrueNumInteractions()); }
      else if(puI->getBunchCrossing() == 0){ h_TrueNInt_inTime_        ->Fill(puI->getTrueNumInteractions()); }
      else if(puI->getBunchCrossing() >  0){ h_TrueNInt_outOfTimeLate_ ->Fill(puI->getTrueNumInteractions()); }
    }
  }
  // ----------------------------

  return;
}

DEFINE_FWK_MODULE(PileupMCTemplateMaker);
