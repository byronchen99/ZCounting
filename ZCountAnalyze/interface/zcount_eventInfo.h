/**
 * *
 *    
 * Created on: 29 June 2019
 *   Author: David Walter
 *  
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_EVENTINFO_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_EVENTINFO_H_

#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class zcount_eventInfo: public zcount_module{
public:
    zcount_eventInfo();
    ~zcount_eventInfo();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

private:

    // input
    
    // output
    int runNumber_;
    int lumiSection_;
    int eventNumber_;
};

#endif

