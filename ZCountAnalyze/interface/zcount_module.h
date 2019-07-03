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
