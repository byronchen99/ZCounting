/**
 *  * Base class for modules to inherit from.
 *   
 * Created on: 29 June 2019
 *    Author: David Walter
 *
 *  *  *
 * Functions:
 *    getInput
 *        - executed right after construction 
 *        - set constants like pt cut, ...
 *    initBranches
 *        - executed after constructing the TTree
 *        - set branches for variables that should be stored in the output root file
 *        - only use 'addBranch'
 *    readSetup
 *        - executed for each event
 *        - set event specific variables
 *    readEvent
 *        - executed right after readSetup
 *        - determine, if the event should be stored(return true) or dicarded(return false)
 *    fillBranches
 *        - executed after readEvent
 *        - do final calculations for variables to be stored
 *  *  *
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_MODULE_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_MODULE_H_

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TTree.h"

#include "TString.h"


class zcount_module{
public:
    zcount_module():read_(false){}
    virtual ~zcount_module();

    virtual void getInput(const edm::ParameterSet& iConfig){}
    virtual void initBranches(TTree* )=0;
    virtual void readSetup(const edm::EventSetup& iSetup){}
    virtual bool readEvent(const edm::Event& iEvent)=0;
    virtual void fillBranches()=0;

protected:
    template <class T>
    void addBranch(TTree* t, const char* name, T*, const char* leaflist=0);

    // definitions
    enum MuonIDTypes { NoneID, LooseID, MediumID, TightID };
    enum MuonIsoTypes { NoneIso, TrackerIso, PFIso };
    enum MuonCategory { cNone, cTrk, cSta, cGlo, cSel, cHLT};

    const double MUON_MASS = 0.105658369;
    const double MUON_MAX_DELTAR_MATCH = 0.5;   // maximum delta R between gen and reco muon by matching
    const double MUON_MAX_DRELPT_MATCH = 0.5;   // maximum relative pt difference between gen and reco muon by matching

private:
    std::vector<TString> allbranches_;
    bool read_;
};

template <class T>
void zcount_module::addBranch(TTree* t, const char* name,  T* address, const char* leaflist)
{

    if(read_ ){
        t->SetBranchAddress(name,address);
    }
    else{
        if(leaflist)
            t->Branch(name  ,address  ,leaflist );
        else
            t->Branch(name  ,address);
    }
    allbranches_.push_back((TString)name);
}


#endif
