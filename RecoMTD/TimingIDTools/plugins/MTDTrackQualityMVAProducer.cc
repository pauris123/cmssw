#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


#include <vector>
#include <memory>
#include <map>

using namespace std;
using namespace edm;

class MTDTrackQualityMVAProducer : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit MTDTrackQualityMVAProducer(const edm::ParameterSet& iConfig, const cms::Ort::ONNXRuntime* cache);
  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet& iConfig);
  static void globalEndJob(const cms::Ort::ONNXRuntime* cache);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  template <class H, class T>
  void fillValueMap(edm::Event& iEvent,
                    const edm::Handle<H>& handle,
                    const std::vector<T>& vec,
                    const std::string& name) const;

  void produce(edm::Event& iEvent, const edm::EventSetup&iSetup) final; /// was override instead of final

private:

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksMTDToken_;
  edm::EDGetTokenT<reco::BeamSpot> RecBeamSpotToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> mtdTimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmamtdTimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixBarrelToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixEndcapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> outermostHitPositionToken_;

  const std::vector<float> means_ = {1.11761003e+00, -4.92049884e-03, 7.26265140e-04, 3.43678598e+00,
                                     2.09002652e-02, 3.52629064e+01, 2.90593730e+01, 2.00602628e+00,
                                     4.62898690e+00, 1.25299188e+00, -3.33957107e-01, 4.21849305e+00,
                                     1.86598480e-01, 9.00257806e+00, 3.21656671e-02, 2.65071352e+02,
                                     3.34686673e+01};
  const std::vector<float> scales_ = {1.52643100e+00, 1.83874287e+00, 1.81416222e+00, 2.62584658e+00,
                                      6.25817336e-02, 2.96932933e+01, 1.08452474e+01, 1.51489321e+00,
                                      4.37572451e+00, 6.16078599e+00, 1.21758077e+00, 9.20447762e+00,
                                      1.35221298e+00, 2.42142681e+00, 7.63215421e-03, 7.23910626e+01,
                                      1.93965525e+02};
};

std::unique_ptr<cms::Ort::ONNXRuntime> MTDTrackQualityMVAProducer::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  
  std::string modelPath = iConfig.getParameter<std::string>("DNN_model");

  // Create session options (optional - customize based on your needs)
  Ort::SessionOptions sessionOptions;

  // Example: Set graph optimization level (this is optional, depending on your needs)
  sessionOptions.SetGraphOptimizationLevel(ORT_ENABLE_EXTENDED);

  // Optionally, set other session options if required
  // sessionOptions.SetLogVerbosityLevel(4); // Example of log verbosity setting
  // sessionOptions.SetIntraOpNumThreads(4); // Set number of threads for computation, if needed

  // Create the ONNXRuntime object with the model path and session options
  return std::make_unique<cms::Ort::ONNXRuntime>(modelPath, &sessionOptions);
  
  //return std::make_unique<cms::Ort::ONNXRuntime>(iConfig.getParameter<std::string>("RecoMTD/TimingIDTools/data/DNN_model_final__bugfix_test.onnx"));
}

void MTDTrackQualityMVAProducer::globalEndJob(const cms::Ort::ONNXRuntime* cache) {}

MTDTrackQualityMVAProducer::MTDTrackQualityMVAProducer(const edm::ParameterSet& iConfig, const cms::Ort::ONNXRuntime* cache)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      RecBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBS"))),
      btlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2Src"))),
      btlMatchTimeChi2Token_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2Src"))),
      etlMatchChi2Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2Src"))),
      etlMatchTimeChi2Token_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2Src"))),
      mtdTimeToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mtdTimeSrc"))),
      sigmamtdTimeToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmamtdTimeSrc"))),
      pathLengthToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"))),
      npixBarrelToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixBarrelSrc"))),
      npixEndcapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixEndcapSrc"))),
      outermostHitPositionToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("outermostHitPositionSrc"))) {
      //modelPath(iConfig.getParameter<edm::FileInPath>("DNN_file").fullPath()) {
  produces<edm::ValueMap<float>>("mtdQualMVA");
}

void MTDTrackQualityMVAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"))->setComment("Input tracks collection");
  desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchChi2"))
      ->setComment("BTL Chi2 Matching value Map");
  desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"))
      ->setComment("BTL Chi2 Matching value Map");
  desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchChi2"))
      ->setComment("ETL Chi2 Matching value Map");
  desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"))
      ->setComment("ETL Chi2 Matching value Map");
  desc.add<edm::InputTag>("mtdTimeSrc", edm::InputTag("trackExtenderWithMTD", "generalTracktmtd"))
      ->setComment("MTD TIme value Map");
  desc.add<edm::InputTag>("sigmamtdTimeSrc", edm::InputTag("trackExtenderWithMTD", "generalTracksigmatmtd"))
      ->setComment("sigma MTD TIme value Map");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD", "generalTrackPathLength"))
      ->setComment("MTD PathLength value Map");
  desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"))
      ->setComment("# of Barrel pixel associated to refitted tracks");
  desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"))
      ->setComment("# of Endcap pixel associated to refitted tracks");
  desc.add<edm::InputTag>("outermostHitPositionSrc", edm::InputTag("trackExtenderWithMTD", "generalTrackOutermostHitPosition"));
  desc.add<edm::InputTag>("offlineBS", edm::InputTag("offlineBeamSpot"));
  desc.add<std::string>("DNN_model", std::string("/afs/cern.ch/user/n/nstrautn/CMSSW_15_0_0_pre3/src/RecoMTD/TimingIDTools/data/DNN_model_final__bugfix_test.onnx"));
  descriptions.add("MTDTrackQualityMVAProducer", desc);
}

template <class H, class T>
void MTDTrackQualityMVAProducer::fillValueMap(edm::Event& iEvent,
                                              const edm::Handle<H>& handle,
                                              const std::vector<T>& vec,
                                              const std::string& name) const {
  auto out = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*out);
  filler.insert(handle, vec.begin(), vec.end());
  filler.fill();
  iEvent.put(std::move(out), name);
}

void MTDTrackQualityMVAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //auto qualityScores = std::make_unique<std::vector<float>>();
  std::vector<float> qualityScores;

  edm::Handle<reco::TrackCollection> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> BeamSpotH;
  iEvent.getByToken(RecBeamSpotToken_, BeamSpotH);
  beamSpot = *BeamSpotH;

  const auto& btlMatchChi2 = iEvent.get(btlMatchChi2Token_);
  const auto& btlMatchTimeChi2 = iEvent.get(btlMatchTimeChi2Token_);
  const auto& etlMatchChi2 = iEvent.get(etlMatchChi2Token_);
  const auto& etlMatchTimeChi2 = iEvent.get(etlMatchTimeChi2Token_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& npixBarrel = iEvent.get(npixBarrelToken_);
  const auto& npixEndcap = iEvent.get(npixEndcapToken_);
  const auto& mtdTime = iEvent.get(mtdTimeToken_);
  const auto& sigmamtdTime = iEvent.get(sigmamtdTimeToken_);
  const auto& lHitPos = iEvent.get(outermostHitPositionToken_);

  std::vector<std::string> inputNames = {"dense_input"};  // Input tensor name
  std::vector<std::string> outputNames = {"dense_2"};  // Output tensor name

  

  for(unsigned int itrack = 0; itrack < tracks.size(); ++itrack) {
    const reco::TrackRef trackref(tracksH, itrack);
    if (pathLength[trackref] == -1.)
      //qualityScores->push_back(-1.);
      qualityScores.push_back(-1.);
    else {
      std::vector<float> inputFeatures = {
        static_cast<float>(trackref->pt()),
        static_cast<float>(trackref->eta()),
        static_cast<float>(trackref->phi()),
        static_cast<float>(trackref->dz(beamSpot.position())),
        static_cast<float>(trackref->dxy(beamSpot.position())),
        static_cast<float>(trackref->chi2()),
        static_cast<float>(trackref->ndof()),
        static_cast<float>(npixBarrel[trackref]),
        static_cast<float>(npixEndcap[trackref]),
        static_cast<float>(btlMatchChi2[trackref]),
        static_cast<float>(btlMatchTimeChi2[trackref]),
        static_cast<float>(etlMatchChi2[trackref]),
        static_cast<float>(etlMatchTimeChi2[trackref]),
        static_cast<float>(mtdTime[trackref]),
        static_cast<float>(sigmamtdTime[trackref]),
        static_cast<float>(pathLength[trackref]),
        static_cast<float>(lHitPos[trackref])
      };

      //std::vector<long unsigned int> inputShapes = {1, inputFeatures.size()}; // long unsigned int
      //std::vector<std::vector<int64_t>> inputShapesWrapped = { inputShapes };

      std::vector<int64_t> inputShape = {1, static_cast<int64_t>(inputFeatures.size())};
      std::vector<std::vector<int64_t>> inputShapesWrapped = {inputShape};

      std::vector<std::vector<float>> wrappedInput = { inputFeatures };
      cms::Ort::FloatArrays inputData(wrappedInput);

      for (size_t j = 0; j < inputFeatures.size(); ++j) {
        inputFeatures[j] = (inputFeatures[j] - means_[j]) / scales_[j];
      }
      //std::vector<float> output = globalCache()->run(inputNames, inputData, inputShapesWrapped, outputNames, 1);
      std::vector<std::vector<float>> output = globalCache()->run(inputNames, inputData, inputShapesWrapped, outputNames, 1);
      //qualityScores->push_back(output[0]);
      for (float score : output[0]) {
        //qualityScores->push_back(score);
        qualityScores.push_back(score);
      }
    }  
  }
  fillValueMap(iEvent, tracksH, qualityScores, "mtdQualMVA");
  
  //iEvent.put(std::move(qualityScores), "mtdQualMVA");
}

DEFINE_FWK_MODULE(MTDTrackQualityMVAProducer);
