#ifndef __RICO_REGISTRY_H__
#define __RICO_REGISTRY_H__

#include "Rico.h"
#include <map>

////////////////////////////////////////////////////////////////////////////////////
// Objects such as BasicRV's may have id codes from here, but not be registered.

typedef int RegistryID;

class RicoRegistry {
 private: 
  static RegistryID registry_id;
  static map<RegistryID, RicoPtr> * rico_registry;

 private: // do not instantiate
  RicoRegistry() {}
  ~RicoRegistry() {}

 public:
  static RicoPtr getRicoPtr (RegistryID id) {return (*rico_registry)[id];}

  static RegistryID registerRicoPtr(RicoPtr rico) {
    (*rico_registry)[registry_id] = rico; 
    return registry_id++;
  }

};

////////////////////////////////////////////////////////////////////////////////////
#endif // __RICO_REGISTRY_H__

