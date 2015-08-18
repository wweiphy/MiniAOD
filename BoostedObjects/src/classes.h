#include "DataFormats/Common/interface/Wrapper.h"

//Add includes for your classes here
#include "MiniAOD/BoostedObjects/interface/BoostedJet.h"
#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"

#include <vector>

namespace {

  struct Boosted_Collections {

    //add 'dummy' Wrapper variable for each class type you put into the Event
    boosted::BoostedJet boostedjetdummy0;
    edm::Wrapper<boosted::BoostedJet> boostedjetdummy1;
    std::vector<boosted::BoostedJet> boostedjetdummy2;
    edm::Wrapper<std::vector<boosted::BoostedJet> > boostedjetdummy3;
    
    boosted::HTTTopJet htttopjetdummy0;
    edm::Wrapper<boosted::HTTTopJet> htttopjetdummy1;
    std::vector<boosted::HTTTopJet> htttopjetdummy2;
    edm::Wrapper<std::vector<boosted::HTTTopJet> > htttopjetdummy3;
    
    boosted::SubFilterJet subfilterjetdummy0;
    edm::Wrapper<boosted::SubFilterJet> subfilterjetdummy1;
    std::vector<boosted::SubFilterJet> subfilterjetdummy2;
    edm::Wrapper<std::vector<boosted::SubFilterJet> > subfilterjetdummy3;
  };
}
