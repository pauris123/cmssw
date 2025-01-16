#include "RecoMTD/TimingIDTools/interface/MTDTrackQualityMVA.h"



MTDTrackQualityMVA::MTDTrackQualityMVA(std::string model_file)
    : env_(ORT_LOGGING_LEVEL_WARNING, "MTDTrackQualityMVA"), session_options_() {
  // Configure ONNX runtime session
  session_options_.SetIntraOpNumThreads(1);
  session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  // Load the ONNX model
  session_ = std::make_unique<Ort::Session>(env_, model_file.c_str(), session_options_);

  // Initialize variable names
  std::string vars_array[] = {MTDTRACKQUALITYMVA_VARS(MTDDNNVAR_STRING)};
  int nvars = sizeof(vars_array) / sizeof(vars_array[0]);
  vars_.assign(vars_array, vars_array + nvars);
}


float MTDTrackQualityMVA::evaluate(const reco::TrackRef& trk,
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
                                     const edm::ValueMap<float>& trk_lhitpos) const {
  constexpr float minPtForDNN = 0.5;
  if (trk->pt() < minPtForDNN) return -1;

  if (tmtds[trk] <= 0) return -1;

  std::vector<float> input(vars_.size());
  input[int(VarID::Track_pt)] = trk->pt();
  input[int(VarID::Track_eta)] = trk->eta();
  input[int(VarID::Track_phi)] = trk->phi();
  input[int(VarID::Track_dz)] = trk->dz(beamspot.position());
  input[int(VarID::Track_dxy)] = trk->dxy(beamspot.position());
  input[int(VarID::Track_chi2)] = trk->chi2();
  input[int(VarID::Track_ndof)] = trk->ndof();
  input[int(VarID::Track_npixBarrelValidHits)] = npixBarrels[trk];
  input[int(VarID::Track_npixEndcapValidHits)] = npixEndcaps[trk];
  input[int(VarID::Track_BTLchi2)] = btl_chi2s[trk];
  input[int(VarID::Track_BTLtime_chi2)] = btl_time_chi2s[trk];
  input[int(VarID::Track_ETLchi2)] = etl_chi2s[trk];
  input[int(VarID::Track_ETLtime_chi2)] = etl_time_chi2s[trk];
  input[int(VarID::Track_Tmtd)] = tmtds[trk];
  input[int(VarID::Track_sigmaTmtd)] = sigmatmtds[trk];
  input[int(VarID::Track_length)] = trk_lengths[trk];
  input[int(VarID::Track_lHitPos)] = trk_lhitpos[trk];

  // Scaling: Apply StandardScaler (mean, std) taken from python DNN training code (hardcoded values)
  const std::vector<float> means = {1.11761003e+00, -4.92049884e-03,  7.26265140e-04,  3.43678598e+00,
        2.09002652e-02,  3.52629064e+01,  2.90593730e+01,  2.00602628e+00,
        4.62898690e+00,  1.25299188e+00, -3.33957107e-01,  4.21849305e+00,
        1.86598480e-01,  9.00257806e+00,  3.21656671e-02,  2.65071352e+02,
        3.34686673e+01};
  const std::vector<float> scales = {1.52643100e+00, 1.83874287e+00, 1.81416222e+00, 2.62584658e+00,
       6.25817336e-02, 2.96932933e+01, 1.08452474e+01, 1.51489321e+00,
       4.37572451e+00, 6.16078599e+00, 1.21758077e+00, 9.20447762e+00,
       1.35221298e+00, 2.42142681e+00, 7.63215421e-03, 7.23910626e+01,
       1.93965525e+02};


  // Scale input features
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = (input[i] - means[i]) / scales[i];
  }

  // Run ONNX inference

  std::vector<int64_t> inputShape = {1, static_cast<int64_t>(input.size())};
  Ort::Value inputTensor = Ort::Value::CreateTensor<float>(
        Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault),
        input.data(), input.size(),
        inputShape.data(), inputShape.size());

  // Perform inference
  const char* inputNames[] = {"dense_input"};
  const char* outputNames[] = {"dense_2"};
  auto outputTensors = session_->Run(Ort::RunOptions{nullptr},
                                        inputNames, &inputTensor, 1,
                                        outputNames, 1);

  // Retrieve the result
  const float* outputData = outputTensors[0].GetTensorData<float>();
  return outputData[0];  // Assuming the model outputs a single value

}
