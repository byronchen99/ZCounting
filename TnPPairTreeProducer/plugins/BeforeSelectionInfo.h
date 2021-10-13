#include <iostream>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include <TH1I.h>
// #include <TH2D.h>


class BeforeSelectionInfo :
    public edm::EDAnalyzer
{
    public:
        explicit BeforeSelectionInfo(const edm::ParameterSet&);
        virtual ~BeforeSelectionInfo() {}

    protected:
        edm::EDGetTokenT<std::vector<PileupSummaryInfo> > src_pileupInfo_;

        // TH1D* h_TrueNInt_outOfTimeEarly_;
        // TH2D* h_TrueNInt_inTime_;
        // TH1D* h_TrueNInt_outOfTimeLate_;

        TH1I* h_LS_;
        // TH2D* h_LS_vs_PU_;

    private:
        void beginJob();
        void beginRun(const edm::Run&, const edm::EventSetup&) override;
        void endRun(const edm::Run&, const edm::EventSetup&) override;
        void analyze(const edm::Event&, const edm::EventSetup&) override;
};


DEFINE_FWK_MODULE(BeforeSelectionInfo);
