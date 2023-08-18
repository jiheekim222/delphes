//////////////////////////////////////////////////////////////
// 08/15/2023 Jihee Kim (jkim11@bnl.gov)
// Generate HepMC file
//////////////////////////////////////////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include <TMath.h>
#include <cmath>
#include <iostream>
#include <math.h>
#include <random>

using namespace HepMC3;

void generate_hepmc(int n_events = 1e3, double pdgid = 11., double e_start = 5.0, double e_end = 5.0, const char* out_fname = "test.hepmc") {
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Particles [pdgid, mass in GeV]
  double particles[15][2] = {
    {11.,   0.51099895e-3}, // electron  
    {13.,   0.1056583745}, // muon-  
    {-13.,  0.1056583745}, // muon+  
    {2212., 0.938272}, // proton
    {22.,   0.0}, // photon
    {-11.,  0.51099895e-3}, // positron
    {2112., 0.939565}, // neutron
    {111.,  0.1349768}, // pion0
    {211.,  0.13957039}, // pion+
    {-211., 0.13957039}, // pion-
    {311.,  0.497611}, // kaon0
    {321.,  0.493677}, // kaon+
    {-321., 0.493677}, // kaon-
    {130.,  0.497611}, // kaon0 long
    {310.,  0.497611}  // kaon0 short 
  };

  // Look for mass of desired particle
  double mass;
  for (int i = 0; i < 11; ++i)
  {
    if (particles[i][0] == pdgid)
      mass = particles[i][1];
  } 
  
  // Random number generator
  TRandom* r1 = new TRandom();

  // Constraining the solid angle, but larger than that subtended by the detector
  // Backward ECAL: 155 < theta_angle < 177 (-3.6 < eta < -1.5) 
  // Barrel ECAL:
  // Forward ECAL: 
  double cos_theta_min = std::cos(M_PI * (45.0 / 180.0));
  double cos_theta_max = std::cos(M_PI * (135.0 / 180.0));

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector((px,py,pz,e),pdgid,status)
    // type 4 beam: 10 GeV electron (pdgid=11) with proton (pdgid=2212) at rest
    GenParticlePtr p1 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum
    Double_t p        = r1->Uniform(e_start, e_end);
    Double_t phi      = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t costheta = r1->Uniform(cos_theta_min, cos_theta_max);
    Double_t theta    = std::acos(costheta);
    Double_t px       = p * std::cos(phi) * std::sin(theta);
    Double_t py       = p * std::sin(phi) * std::sin(theta);
    Double_t pz       = p * std::cos(theta);

    // type 1 final state
    GenParticlePtr p3 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt((p * p) + (mass * mass))), (int)pdgid, 1);

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    v1->add_particle_out(p3);
    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
