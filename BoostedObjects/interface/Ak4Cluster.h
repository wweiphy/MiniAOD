#ifndef MINIAOD_BOOSTEDOBJECTS_AK4CLUSTER_H
#define MINIAOD_BOOSTEDOBJECTS_AK4CLUSTER_H

#include <vector>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TMath.h"

namespace boosted {

  class Ak4Cluster {

    public:

      Ak4Cluster() :
        fatjet(math::XYZTLorentzVector()){};

      // Fatjet or cluster
      math::XYZTLorentzVector fatjet;

      // Subjets of the cluster
      std::vector<pat::Jet> ak4jets;

  };

  typedef std::vector<Ak4Cluster> Ak4ClusterCollection;

}

#endif
