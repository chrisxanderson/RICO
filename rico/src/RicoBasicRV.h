#ifndef __RICO_BASICRV_H__
#define __RICO_BASICRV_H__

#include "Rico.h"

/////////////////////////////////////////////////////////////////////////

class RicoBasicRV : public Rico {
 private:
  static RicoID rico_id;
  static map<RicoID, const RicoBasicRV*> * basic_rv_map;

 private:
  const RicoID m_id;

 public:
  RicoBasicRV();
  ~RicoBasicRV() {}
  
 public:
  virtual RicoID id() const {return m_id;}

  virtual ppPair getPair() const {
    cout << "ERROR: RicoNumber::getPair() undefined." << endl;
    return ppPair(Rico::One, Rico::One);
  }

 public: // Interface
  static const RicoBasicRV& getBasicRV(RicoID id) {return *((*basic_rv_map)[id]);}
  virtual string  getDefinition() const = 0;
  virtual string  getExpression() const;

  virtual RicoPtr differentiate(RicoPtr This, RicoID Wrt) const;

 public: // Util
  virtual void collectBasicRVs(set<RicoID>& rv_set) const {rv_set.insert(id());}

 public: // Operators
  virtual bool operator==(RicoPtr That) const;

}; // RicoBasicRV

/////////////////////////////////////////////////////////////////////////////////////////
#endif  // __RICO_BASICRV_H__
