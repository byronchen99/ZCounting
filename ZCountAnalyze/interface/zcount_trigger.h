/**
 * * Trigger handling
 *    
 * Created on: 01 July 2019
 *   Author: David Walter
 *  
 */

#ifndef ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_TRIGGER_H_
#define ZCOUNTING_ZCOUNTANALYZER_INTERFACE_ZCOUNT_TRIGGER_H_

#include "ZCounting/ZCountAnalyze/interface/zcount_module.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <bitset>

const unsigned int kNTrigBit = 128;
typedef std::bitset<kNTrigBit> TriggerBits;
const unsigned int kNTrigObjectBit = 256;
typedef std::bitset<kNTrigObjectBit> TriggerObjectBits;


class zcount_trigger: public zcount_module{
public:
    zcount_trigger();
    ~zcount_trigger();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

    void setTriggerResultsToken(edm::EDGetTokenT<edm::TriggerResults> token) {
        fHLTTag_token = token;
    }

    void setTriggerEventToken(edm::EDGetTokenT<trigger::TriggerEvent> token) {
        fHLTObjTag_token = token;
    }

    void addTriggerRecord(const std::string &name = "",
                          const std::string &objName = ""){
        Record rec;
        rec.hltPattern = name;
        rec.hltObjName = objName;
        records.push_back(rec);
    }

    bool pass(const std::vector<std::string> &iNames) const;
    bool passObj(const std::vector<std::string> &iNames, const std::vector<std::string> &iObjNames, const double eta, const double phi) const;

private:

    struct Record {
        std::string hltPattern;                         // HLT path name/pattern (wildcards allowed: *,?)
        std::string hltPathName = "";                   // HLT path name in trigger menu
        unsigned int hltPathIndex = (unsigned int)-1;   // HLT path index in trigger menu
        std::string hltObjName = "";                    // trigger object name in trigger menu
    };
    std::vector<Record> records;

    edm::EDGetTokenT<edm::TriggerResults>   fHLTTag_token;
    edm::EDGetTokenT<trigger::TriggerEvent> fHLTObjTag_token;

    edm::Handle<edm::TriggerResults>   hTrgRes;
    edm::Handle<trigger::TriggerEvent> hTrgEvt;

    edm::ParameterSetID fTriggerNamesID;

    // initialization from HLT menu; needs to be called on every change in HLT menu
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);

    TriggerObjectBits matchHLT(const double eta, const double phi) const;

    int getTriggerBit(const std::string &iName) const;
    int getTriggerObjectBit(const std::string &iName, const std::string &iObjName) const;

    TriggerBits triggerBits;

    // input

    // output

};

#endif

