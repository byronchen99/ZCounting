// -*- C++ -*-
//
// Package:    ZCounting/ZCounting
// Class:      ZCounting
// 
/**\class ZCounting ZCounting.h ZCounting/ZCountAnalyze/plugins/ZCounting.h

 Description: 
    Base Class for Z Counting Studies

 Implementation:
    The ZCounting Analyzer steers various modules from the base class zcount_module
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// ZCounting Project includes
#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_PV.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_muons.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_eventInfo.h"
#include "ZCounting/ZCountAnalyze/interface/zcount_genInfo.h"



// ROOT includes
#include "TTree.h"


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

    // ----------member data ---------------------------

    edm::Service<TFileService> fs;
    TTree *tree_;

    bool isData_;
    bool selectEvents_;

    std::vector<zcount_module* > modules_;
    zcount_module * addModule(zcount_module *m){
        modules_.push_back(m);
        return m;
    }

};


//define this as a plug-in
DEFINE_FWK_MODULE(ZCounting);
