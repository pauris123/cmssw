#ifndef RECOMTD_TIMINGIDTOOLS_MTDTRACKQUALITYMVA
#define RECOMTD_TIMINGIDTOOLS_MTDTRACKQUALITYMVA

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include <memory>
#include <vector>
#include <map>

#define MTDTRACKQUALITYMVA_VARS(MTDDNNVAR)       \
  MTDDNNVAR(Track_pt)                            \
  MTDDNNVAR(Track_eta)                           \
  MTDDNNVAR(Track_phi)                           \
  MTDDNNVAR(Track_dz)                            \
  MTDDNNVAR(Track_dxy)                           \
  MTDDNNVAR(Track_chi2)                          \
  MTDDNNVAR(Track_ndof)                          \
  MTDDNNVAR(Track_npixBarrelValidHits)           \
  MTDDNNVAR(Track_npixEndcapValidHits)           \
  MTDDNNVAR(Track_BTLchi2)                       \
  MTDDNNVAR(Track_BTLtime_chi2)                  \
  MTDDNNVAR(Track_ETLchi2)                       \
  MTDDNNVAR(Track_ETLtime_chi2)                  \
  MTDDNNVAR(Track_Tmtd)                          \
  MTDDNNVAR(Track_sigmaTmtd)                     \
  MTDDNNVAR(Track_length)                        \
  MTDDNNVAR(Track_lHitPos)                        

#define MTDDNNVAR_ENUM(ENUM) ENUM,
#define MTDDNNVAR_STRING(STRING) #STRING,

class MTDTrackQualityMVA {
public:
  //---ctors---
  MTDTrackQualityMVA(std::string model_file);

  enum class VarID { MTDTRACKQUALITYMVA_VARS(MTDDNNVAR_ENUM) };

  //---getters---
  // 4D
  float evaluate(const reco::TrackRef& trk,
                   const reco::BeamSpot& beamspot,
                   const edm::ValueMap<int>& npixBarrels,
                   const edm::ValueMap<int>& npixEndcaps,
                   const edm::ValueMap<float>& btl_chi2s,
                   const edm::ValueMap<float>& btl_time_chi2s,
                   const edm::ValueMap<float>& etl_chi2s,
                   const edm::ValueMap<float>& etl_time_chi2s,
                   const edm::ValueMap<float>& tmtds,
                   const edm::ValueMap<float>& sigmatmtds,
                   const edm::ValueMap<float>& trk_lengths,
                   const edm::ValueMap<float>& trk_lhitpos) const;

private:
  // Variables
  std::vector<std::string> vars_; // Names of variables

  // ONNX runtime components
  std::unique_ptr<Ort::Session> session_;
  Ort::Env env_;
  Ort::SessionOptions session_options_;
};

#endif
