#include <string>

#include "MiniAOD/MiniAODHelper/interface/PUJetID.h"

PUJetID::WP PUJetID::get(const std::string& wp_name) {
    if(wp_name=="loose") return PUJetID::loose;
    else if(wp_name=="medium") return PUJetID::medium;
    else if(wp_name=="tight") return PUJetID::tight;
    else return PUJetID::none;
}

std::string PUJetID::toString(const PUJetID::WP wp) {
    if(wp==PUJetID::loose) return "loose";
    else if(wp==PUJetID::medium) return "medium";
    else if(wp==PUJetID::tight) return "tight";
    else return "none";
}

int PUJetID::toInt(const PUJetID::WP wp) {
    if(wp==PUJetID::loose) return 4;
    else if(wp==PUJetID::medium) return 6;
    else if(wp==PUJetID::tight) return 7;
    else return 0;
}