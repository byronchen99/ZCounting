#include "ZCounting/ZUtils/interface/GenZDecayProperties.h"

#include "DataFormats/Common/interface/Wrapper.h"




namespace {
  struct dictionary {
    GenZDecayProperties genZDecayProperties;
    std::vector<GenZDecayProperties> v_genZDecayProperties;
    edm::Wrapper<GenZDecayProperties> w_genZDecayProperties;
    edm::RefProd<GenZDecayProperties> rp_genZDecayProperties;
    edm::Wrapper<std::vector<GenZDecayProperties> > w_v_genZDecayProperties;
  };
}
