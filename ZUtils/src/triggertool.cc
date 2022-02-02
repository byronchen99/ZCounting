#include "ZCounting/ZUtils/interface/triggertool.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <algorithm>



//--------------------------------------------------------------------------------------------------
bool triggertool::readEvent(const edm::Event& iEvent){
    edm::LogVerbatim("triggertools") << "triggertools: readEvent";

    iEvent.getByToken(fHLTTag_token, hTrgRes);
    if (!hTrgRes.isValid()){
        edm::LogWarning("triggertool") << "No valid trigger result product found" ;
    }

    iEvent.getByToken(fHLTObjTag_token, hTrgEvt);
    if (!hTrgEvt.isValid()){
        edm::LogWarning("triggertool") << "No valid trigger event product found" ;
    }

    triggerBits.reset();
    for (unsigned int i = 0; i < records.size(); i++) {
        if (records.at(i).hltPathIndex == (unsigned int)-1){
            edm::LogWarning("triggertools")<<"hltPathIndex has not been set"<<std::endl;
            continue;
        }
        if (hTrgRes->accept(records.at(i).hltPathIndex)) {
            triggerBits[i] = true;
        }
    }
    edm::LogVerbatim("triggertools") << "triggertools: bitset = "<<triggerBits[1]<<triggerBits[0];

    return true;
}

//--------------------------------------------------------------------------------------------------
int triggertool::getTriggerBit(const std::string &iName) const {
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
bool triggertool::pass(const std::string &iName) const {
    const int lId = getTriggerBit(iName);
    if(lId != -1 && triggerBits[lId])
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool triggertool::pass(const std::vector<std::string> &iNames) const {
    // return true, if one of the triggers, given by iNames has fired
    for(unsigned int i = 0; i < iNames.size(); i++){
        if(pass(iNames.at(i)))
            return true;
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
bool triggertool::passObj(const std::string &iName,
                          const double eta,
                          const double phi) const {
    // return true, if the triggers, given by iName and iObjName has fired an Object, which is matched with
    //   the given eta, phi direction
    const TriggerObjectBits iTrigObj = matchHLT(eta, phi);

    const int lId = getTriggerBit(iName);
    if(lId != -1 && iTrigObj[lId])
        return true;
    return false;
}

//--------------------------------------------------------------------------------------------------
bool triggertool::passObj(const std::vector<std::string> &iNames,
                             const double eta,
                             const double phi) const {
    // return true, if one of the triggers, given by iNames and iObjNames has fired an Object, which is matched with
    //   the given eta, phi direction
    const TriggerObjectBits iTrigObj = matchHLT(eta, phi);

    int lId = -1;
    for(unsigned int i = 0; i < iNames.size(); i++){
        lId = getTriggerBit(iNames.at(i));
        if(lId != -1 && iTrigObj[lId])
            return true;
    }
    return false;
}

//--------------------------------------------------------------------------------------------------
void triggertool::initHLTObjects(const HLTConfigProvider& hltConfigProvider_){
    /*
        execture each run to initialize the last filter of each trigger corresponding to the corresponding object that has fired the trigger
    */
    edm::LogVerbatim("triggertools")<<"initHLTObjects";
    const std::vector<std::string>& triggerNames(hltConfigProvider_.triggerNames());

    initPathNames(triggerNames);

    for (auto &iRec: records) {

        std::vector<std::string> hltFiltersWithTags_;

        for (auto const& iPathName : triggerNames) {

            edm::LogVerbatim("triggertools")<<"trigger name"<<iPathName;

            if(iPathName != iRec.hltPathName){
                continue;
            }

            iRec.hltPathIndex = hltConfigProvider_.triggerIndex(iPathName);

            auto const& moduleLabels(hltConfigProvider_.moduleLabels(iRec.hltPathIndex));

            for(int idx=moduleLabels.size()-1; idx >= 0; --idx){
                auto const& moduleLabel(moduleLabels.at(idx));

                auto const& moduleEDMType(hltConfigProvider_.moduleEDMType(moduleLabel));
                if(moduleEDMType != "EDFilter"){
                    continue;
                }

                auto const& moduleType(hltConfigProvider_.moduleType(moduleLabel));
                if((moduleType == "HLTTriggerTypeFilter") or (moduleType == "HLTBool") or (moduleType == "HLTPrescaler")){
                    continue;
                }

                if(!hltConfigProvider_.saveTags(moduleLabel)){
                    continue;
                }
                edm::LogVerbatim("triggertools")<<"new hlt object name: "<< moduleLabel;

                iRec.hltObjName = moduleLabel;
                break;
            }

            break;

        }
    }
}

//--------------------------------------------------------------------------------------------------
void triggertool::initPathNames(const std::vector<std::string>& triggerNames) {
    /*
        init HLT path every run (e.g. versions can change)
    */
    edm::LogVerbatim("triggertools") << " initHLT" ;
    for (auto &iRec: records) {
        iRec.hltPathName = "";
        iRec.hltPathIndex = (unsigned int)-1;
        const std::string pattern = iRec.hltPattern;
        if (edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
            std::vector<std::vector<std::string>::const_iterator> matches =
                edm::regexMatch(triggerNames, pattern);
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
    }
}

//--------------------------------------------------------------------------------------------------
TriggerObjectBits triggertool::matchHLT(const double eta, const double phi) const {

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
                if (reco::deltaR(eta, phi, tobj.eta(), tobj.phi()) < DRMAX) {
                    matchBits[filterBit] = true;
                }
            }
        }
    }
    return matchBits;
}
