#ifndef ZCOUNTING_TNPPAIRTREEPRODUCER_INTERFACE_TRIGGERTOOL_H
#define ZCOUNTING_TNPPAIRTREEPRODUCER_INTERFACE_TRIGGERTOOL_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <bitset>

const unsigned int kNTrigBit = 128;
typedef std::bitset<kNTrigBit> TriggerBits;
const unsigned int kNTrigObjectBit = 256;
typedef std::bitset<kNTrigObjectBit> TriggerObjectBits;


class triggertool{
public:
    triggertool(){};
    ~triggertool(){};

    bool readEvent(const edm::Event& iEvent);

    void fillBranches();

    void setTriggerResultsToken(edm::EDGetTokenT<edm::TriggerResults> token) {
        fHLTTag_token = token;
    }

    void setTriggerEventToken(edm::EDGetTokenT<trigger::TriggerEvent> token) {
        fHLTObjTag_token = token;
    }

    void addTriggerRecord(const std::string &name){
        Record rec;
        rec.hltPattern = name;
        records.push_back(rec);
    }
    void addTriggerRecord(const std::string &name,
                          const std::string &objName){
        Record rec;
        rec.hltPattern = name;
        rec.hltObjName = objName;
        records.push_back(rec);
    }

    void initHLTObjects(const HLTConfigProvider& hltConfigProvider_);

    bool pass(const std::vector<std::string> &iNames) const;
    bool pass(const std::string &iName) const;
    bool passObj(const std::string &iName, const double eta, const double phi) const;
    bool passObj(const std::vector<std::string> &iNames, const double eta, const double phi) const;

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
    void initPathNames(const edm::TriggerResults&, const std::vector<std::string>&);
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);

    TriggerObjectBits matchHLT(const double eta, const double phi) const;

    int getTriggerBit(const std::string &iName) const;

    TriggerBits triggerBits;

};

#endif
