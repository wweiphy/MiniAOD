#ifndef PUJETID_H
#define PUJETID_H

#include <string>

class PUJetID {
public:
  enum WP {
    none,
    loose,
    medium,
    tight
  };
  static PUJetID::WP get(const std::string& wp_name);
  static std::string toString(const PUJetID::WP wp);
  static int toInt(const PUJetID::WP wp);
  
};
#endif
    