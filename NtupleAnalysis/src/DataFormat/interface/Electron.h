// -*- c++ -*-
#ifndef DataFormat_Electron_h
#define DataFormat_Electron_h

#include "DataFormat/interface/ElectronGenerated.h"
#include "DataFormat/interface/ParticleIterator.h"

class Electron;

class ElectronCollection: public ElectronGeneratedCollection, public ParticleIteratorAdaptor<ElectronCollection> {
public:
  using value_type = Electron;

  ElectronCollection() {}
  ElectronCollection(const std::string& prefix): ElectronGeneratedCollection(prefix) {}
  ~ElectronCollection() {}

  void setupBranches(BranchManager& mgr);

  Electron operator[](size_t i) const;
  std::vector<Electron> toVector() const;

  friend class Electron;
  friend class ElectronGenerated<ElectronCollection>;
  friend class Particle<ElectronCollection>;

protected:
};

class Electron: public ElectronGenerated<ElectronCollection> {
public:
  Electron() {}
  Electron(const ElectronCollection* coll, size_t index): ElectronGenerated(coll, index) {}
  ~Electron() {}
};

inline
Electron ElectronCollection::operator[](size_t i) const {
  return Electron(this, i);
}

inline
std::vector<Electron> ElectronCollection::toVector() const {
  return ParticleCollectionBase::toVector(*this);
}

#endif