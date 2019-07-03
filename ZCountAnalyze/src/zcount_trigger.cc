#include "ZCounting/ZCountAnalyze/interface/zcount_trigger.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <algorithm>

zcount_trigger::zcount_trigger(){
}

zcount_trigger::~zcount_trigger(){
}

void zcount_trigger::getInput(const edm::ParameterSet& iConfig){
    edm::LogVerbatim("zcount_trigger") << " getInput";
}

void zcount_trigger::initBranches(TTree* tree){
}

bool zcount_trigger::readEvent(const edm::Event& iEvent){
    edm::LogVerbatim("zcount_triggers") << "zcount_triggers: readEvent";

    iEvent.getByToken(fHLTTag_token, hTrgRes);
    if (!hTrgRes.isValid())
        return false;

    iEvent.getByToken(fHLTObjTag_token, hTrgEvt);

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*hTrgRes);
    if (fTriggerNamesID != triggerNames.parameterSetID()) {
        fTriggerNamesID = triggerNames.parameterSetID();
        initHLT(*hTrgRes, triggerNames);
    }

    triggerBits.reset();
    for (unsigned int i = 0; i < records.size(); i++) {
        if (records.at(i).hltPathIndex == (unsigned int)-1)
            continue;
        if (hTrgRes->accept(records.at(i).hltPathIndex)) {
            triggerBits[i] = true;
        }
    }
    edm::LogVerbatim("zcount_triggers") << "zcount_triggers: bitset = "<<triggerBits[1]<<triggerBits[0];

    return true;
}

void zcount_trigger::fillBranches(){
    return;
}

//--------------------------------------------------------------------------------------------------
int zcount_trigger::getTriggerBit(const std::string &iName) const {
    int lId = -1;
    for (unsigned int i = 0; i < records.size(); i++) {
        if (iName == records.at(i).hltPattern)
            lId = i;
    }
    if (lId == -1)
        edm::LogWarning("ZCounting") << "=== Missing Trigger ==" << iName ;
    return lId;
}

//--------------------------------------------------------------------------------------------------
int zcount_trigger::getTriggerObjectBit(const std::string &iName, const std::string &iObjName) const {
    int lId = getTriggerBit(iName);
    if (lId == -1)
        return -1;

    for (unsigned int i = 0; i < records.size(); i++) {
        if (iObjName != records.at(i).hltObjName)
            continue;
        return i;
    }
    return -1;
}

//--------------------------------------------------------------------------------------------------
bool zcount_trigger::pass(const std::vector<std::string> &iNames) const {
    // return true, if one of the triggers, given by iNames has fired
    int lId = -1;
    for(unsigned int i = 0; i < iNames.size(); i++){
        lId = getTriggerBit(iNames.at(i));
        if(lId != -1 && triggerBits[lId])
            return true;
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
bool zcount_trigger::passObj(const std::vector<std::string> &iNames,
                             const std::vector<std::string> &iObjNames,
                             const double eta,
                             const double phi) const {
    // return true, if one of the triggers, given by iNames and iObjNames has fired an Object, which is matched with
    //   the given eta, phi direction
    const TriggerObjectBits iTrigObj = matchHLT(eta, phi);

    int lId = -1;
    for(unsigned int i = 0; i < iNames.size(); i++){
        lId = getTriggerObjectBit(iNames.at(i), iObjNames.at(i));
        if(lId != -1 && iTrigObj[lId])
            return true;
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
void zcount_trigger::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames) {
    edm::LogVerbatim("zcount_triggers") << " initHLT" ;
    for (auto &iRec: records) {
        iRec.hltPathName = "";
        iRec.hltPathIndex = (unsigned int)-1;
        const std::string pattern = iRec.hltPattern;
        if (edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
            std::vector<std::vector<std::string>::const_iterator> matches =
                edm::regexMatch(triggerNames.triggerNames(), pattern);
            if (matches.empty()) {
                edm::LogWarning("ZCounting") << "requested pattern [" << pattern << "] does not match any HLT paths";
            }
            else {
                for (auto const& match : matches) {
                    iRec.hltPathName = *match;
                }
            }
        }
        else {  // take full HLT path name given
             iRec.hltPathName = pattern;
        }
        // Retrieve index in trigger menu corresponding to HLT path
        unsigned int index = triggerNames.triggerIndex(iRec.hltPathName);
        if (index < result.size()) {  // check for valid index
            iRec.hltPathIndex = index;
        }
    }
}

//--------------------------------------------------------------------------------------------------
TriggerObjectBits zcount_trigger::matchHLT(const double eta, const double phi) const {
    const double dRMax = 0.2;

    TriggerObjectBits matchBits;
    for (unsigned int i = 0; i < records.size(); i++) {
        const std::string& filterName = records.at(i).hltObjName;
        const unsigned int filterBit = i;

        edm::InputTag filterTag(filterName, "", "HLT");
        // filterIndex must be less than the size of trgEvent or you get a CMSException: _M_range_check
        if (hTrgEvt->filterIndex(filterTag) < hTrgEvt->sizeFilters()) {
            const trigger::TriggerObjectCollection& toc(hTrgEvt->getObjects());
            const trigger::Keys& keys(hTrgEvt->filterKeys(hTrgEvt->filterIndex(filterTag)));

            for (unsigned int hlto = 0; hlto < keys.size(); hlto++) {
                trigger::size_type hltf = keys[hlto];
                const trigger::TriggerObject& tobj(toc[hltf]);
                if (reco::deltaR(eta, phi, tobj.eta(), tobj.phi()) < dRMax) {
                    matchBits[filterBit] = true;
                }
            }
        }
    }
    return matchBits;
}