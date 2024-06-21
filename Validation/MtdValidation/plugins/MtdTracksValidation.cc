/**/
#define EDM_ML_DEBUG
#include <string>

#include "TTree.h"
#include "TFile.h"
#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h> // Important to add <CommonTools/UtilAlgos> in buildfile.xml file!!!
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDSectorForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetSector.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"

#include "RecoLocalFastTime/FTLClusterizer/interface/MTDClusterParameterEstimator.h"
#include "RecoLocalFastTime/Records/interface/MTDCPERecord.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "MTDHit.h"

class MtdTracksValidation : public DQMEDAnalyzer {
public:
  explicit MtdTracksValidation(const edm::ParameterSet&);
  ~MtdTracksValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  bool isSameCluster(const FTLCluster&, const FTLCluster&);

  const std::pair<bool, bool> checkAcceptance(
      const reco::Track&, const edm::Event&, const edm::EventSetup&, size_t&, float&, float&, float&, float&);

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaTPSelLV(const TrackingParticle&);
  const bool mvaTPSelAll(const TrackingParticle&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&, const bool&);
  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);

  const unsigned long int uniqueId(const uint32_t x, const EncodedEventId& y) {
    const uint64_t a = static_cast<uint64_t>(x);
    const uint64_t b = static_cast<uint64_t>(y.rawId());

    if (x < y.rawId())
      return (b << 32) | a;
    else
      return (a << 32) | b;
  }

  bool isETL(const double eta) const { return (std::abs(eta) > trackMinEtlEta_) && (std::abs(eta) < trackMaxEtlEta_); }

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;
  const bool OutputForBDTinput_;

  static constexpr double etacutGEN_ = 4.;               // |eta| < 4;
  static constexpr double etacutREC_ = 3.;               // |eta| < 3;
  static constexpr double pTcut_ = 0.7;                  // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;              // dz separation 1 mm
  static constexpr double deltaPTcut_ = 0.05;            // dPT < 5%
  static constexpr double deltaDRcut_ = 0.03;            // DeltaR separation
  static constexpr double depositBTLthreshold_ = 1;      // threshold for energy deposit in BTL cell [MeV]
  static constexpr double depositETLthreshold_ = 0.001;  // threshold for energy deposit in ETL cell [MeV]
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  static constexpr double etaMatchCut_ = 0.05;
  static constexpr double cluDRradius_ = 0.05;  // to cluster rechits around extrapolated track

  // cuts for BDT training input
  static constexpr double BDT_track_eta_cut = 3.0;
  static constexpr double BDT_track_pt_cut = 0.5;
  static constexpr double BDT_track_dz_cut = 1.0; // 1cm cut for BDT training sample

  // Parameters to be saved for BDT output in TTree
  TTree* track_BDT_tree;
  std::vector<double> Ttrack_pt, Ttrack_eta, Ttrack_phi, Ttrack_dz, Ttrack_dxy, Ttrack_chi2, Ttrack_BTLchi2, Ttrack_BTLtime_chi2, Ttrack_ETLchi2, Ttrack_ETLtime_chi2, Ttrack_t0, Ttrack_sigmat0, Ttrack_Tmtd, Ttrack_sigmaTmtd, Ttrack_lenght, Ttrack_MtdMVA, TtrackTP_pt, TtrackTP_eta, TtrackTP_phi;
  std::vector<int> Ttrack_ndof, Ttrack_nValidHits, Ttrack_npixBarrelValidHits, Ttrack_npixEndcapValidHits, TtrackTP_nValidHits, Ttrack_det, Ttrack_subdet, Ttrack_subdetmtd,Ttrack_nSimClust;
  std::vector<bool> Ttrack_Signal, Ttrack_Associated, Ttrack_HasTP, Ttrack_TPHasSimClu, Ttrack_TPSimCluDirect, Ttrack_hitMTD, Ttrack_hit, Ttrack_hitMatchMTD, Ttrack_hasRecoClu, Ttrack_RecoSimLink, Ttrack_FullMatch;
  std::vector<double> Ttrack_CluX1, Ttrack_CluX2, Ttrack_CluY1, Ttrack_CluY2, Ttrack_CluTime1, Ttrack_CluTime2;


  const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> mtdTrackingHitToken_;

  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> r2sAssociationMapToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> tp2SimAssociationMapToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit>> btlSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit>> etlSimHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixBarrelToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> npixEndcapToken_;

  edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> mtdlayerToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  MonitorElement* meBTLTrackRPTime_;
  MonitorElement* meBTLTrackEffEtaTot_;
  MonitorElement* meBTLTrackEffPhiTot_;
  MonitorElement* meBTLTrackEffPtTot_;
  MonitorElement* meBTLTrackEffEtaMtd_;
  MonitorElement* meBTLTrackEffPhiMtd_;
  MonitorElement* meBTLTrackEffPtMtd_;
  MonitorElement* meBTLTrackPtRes_;

  MonitorElement* meETLTrackRPTime_;
  MonitorElement* meETLTrackEffEtaTot_[2];
  MonitorElement* meETLTrackEffPhiTot_[2];
  MonitorElement* meETLTrackEffPtTot_[2];
  MonitorElement* meETLTrackEffEtaMtd_[2];
  MonitorElement* meETLTrackEffPhiMtd_[2];
  MonitorElement* meETLTrackEffPtMtd_[2];
  MonitorElement* meETLTrackEffEta2Mtd_[2];
  MonitorElement* meETLTrackEffPhi2Mtd_[2];
  MonitorElement* meETLTrackEffPt2Mtd_[2];
  MonitorElement* meETLTrackPtRes_;

  MonitorElement* meTracktmtd_;
  MonitorElement* meTrackt0Src_;
  MonitorElement* meTrackSigmat0Src_;
  MonitorElement* meTrackt0Pid_;
  MonitorElement* meTrackSigmat0Pid_;
  MonitorElement* meTrackt0SafePid_;
  MonitorElement* meTrackSigmat0SafePid_;
  MonitorElement* meTrackNumHits_;
  MonitorElement* meTrackNumHitsNT_;
  MonitorElement* meTrackMVAQual_;
  MonitorElement* meTrackPathLenghtvsEta_;

  MonitorElement* meTrackSigmaTof_[3];
  MonitorElement* meTrackSigmaTofvsP_[3];

  MonitorElement* meTrackPtTot_;
  MonitorElement* meMVATrackEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtMtd_;
  MonitorElement* meExtraPtMtd_;
  MonitorElement* meExtraPtEtl2Mtd_;

  MonitorElement* meBTLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResMtd_;
  MonitorElement* meBTLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTP2PtRatioGen_;
  MonitorElement* meBTLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTP2PtRatioMtd_;
  MonitorElement* meBTLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResvsPtMtd_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtGen_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtMtd_;

  MonitorElement* meTrackMatchedTPEffPtTot_;
  MonitorElement* meTrackMatchedTPEffPtMtd_;
  MonitorElement* meTrackMatchedTPEffPtEtl2Mtd_;
  MonitorElement* meTrackMatchedTPmtdEffPtTot_;
  MonitorElement* meTrackMatchedTPmtdEffPtMtd_;
  MonitorElement* meTrackEtaTot_;
  MonitorElement* meMVATrackEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaMtd_;
  MonitorElement* meExtraEtaMtd_;
  MonitorElement* meExtraEtaEtl2Mtd_;
  MonitorElement* meTrackMatchedTPEffEtaTot_;
  MonitorElement* meTrackMatchedTPEffEtaMtd_;
  MonitorElement* meTrackMatchedTPEffEtaEtl2Mtd_;
  MonitorElement* meTrackMatchedTPmtdEffEtaTot_;
  MonitorElement* meTrackMatchedTPmtdEffEtaMtd_;
  MonitorElement* meMVATrackResTot_;
  MonitorElement* meMVATrackPullTot_;
  MonitorElement* meMVATrackZposResTot_;

  MonitorElement* meExtraPhiAtBTL_;
  MonitorElement* meExtraPhiAtBTLmatched_;
  MonitorElement* meExtraBTLeneInCone_;
  MonitorElement* meExtraMTDfailExtenderEta_;
  MonitorElement* meExtraMTDfailExtenderPt_;
};

bool MtdTracksValidation::isSameCluster(const FTLCluster& clu1, const FTLCluster& clu2) {
  return clu1.id() == clu2.id() && clu1.size() == clu2.size() && clu1.x() == clu2.x() && clu1.y() == clu2.y() &&
         clu1.time() == clu2.time();
}

// ------------ constructor and destructor --------------
MtdTracksValidation::MtdTracksValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")),
      OutputForBDTinput_(iConfig.getParameter<bool>("OutputForBDTinput")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  tp2SimAssociationMapToken_ =
      consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMapTag"));
  btlSimHitsToken_ = consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("btlSimHits"));
  etlSimHitsToken_ = consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("etlSimHits"));
  btlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("btlRecHits"));
  etlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("etlRecHits"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  SigmaTofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofPi"));
  SigmaTofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofK"));
  SigmaTofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofP"));
  btlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchChi2Src"));
  btlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2Src"));
  etlMatchChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchChi2Src"));
  etlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2Src"));
  npixBarrelToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixBarrelSrc"));
  npixEndcapToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("npixEndcapSrc"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  btlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagBTL"));
  etlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagETL"));
  mtdTrackingHitToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("trkHitTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  mtdgeoToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  mtdlayerToken_ = esConsumes<MTDDetLayerGeometry, MTDRecoGeometryRecord>();
  magfieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();

  edm::Service<TFileService> fs;
  track_BDT_tree = fs->make<TTree>("Events","Events");
  // // Track kinematic/tracker variables
  track_BDT_tree->Branch("Track_pt",&Ttrack_pt);
  track_BDT_tree->Branch("Track_eta",&Ttrack_eta);
  track_BDT_tree->Branch("Track_phi",&Ttrack_phi);
  track_BDT_tree->Branch("Track_dz",&Ttrack_dz);
  track_BDT_tree->Branch("Track_dxy",&Ttrack_dxy);
  track_BDT_tree->Branch("Track_chi2",&Ttrack_chi2);
  track_BDT_tree->Branch("Track_ndof",&Ttrack_ndof);
  track_BDT_tree->Branch("Track_nValidHits",&Ttrack_nValidHits);

  track_BDT_tree->Branch("TrackTP_pt",&TtrackTP_pt);
  track_BDT_tree->Branch("TrackTP_eta",&TtrackTP_eta);
  track_BDT_tree->Branch("TrackTP_phi",&TtrackTP_phi);
  track_BDT_tree->Branch("TrackTP_nValidHits",&TtrackTP_nValidHits);

  track_BDT_tree->Branch("Track_npixBarrelValidHits",&Ttrack_npixBarrelValidHits);
  track_BDT_tree->Branch("Track_npixEndcapValidHits",&Ttrack_npixEndcapValidHits);
  track_BDT_tree->Branch("Track_Signal",&Ttrack_Signal);
  track_BDT_tree->Branch("Track_Associated",&Ttrack_Associated);
  // Track MTD variables
  track_BDT_tree->Branch("Track_BTLchi2",&Ttrack_BTLchi2);
  track_BDT_tree->Branch("Track_BTLtime_chi2",&Ttrack_BTLtime_chi2);
  track_BDT_tree->Branch("Track_ETLchi2",&Ttrack_ETLchi2);
  track_BDT_tree->Branch("Track_ETLtime_chi2",&Ttrack_ETLtime_chi2);
  track_BDT_tree->Branch("Track_t0",&Ttrack_t0);
  track_BDT_tree->Branch("Track_sigmat0",&Ttrack_sigmat0);
  track_BDT_tree->Branch("Track_Tmtd",&Ttrack_Tmtd);
  track_BDT_tree->Branch("Track_MtdMVA",&Ttrack_MtdMVA);
  track_BDT_tree->Branch("Track_sigmaTmtd",&Ttrack_sigmaTmtd);
  track_BDT_tree->Branch("Track_lenght",&Ttrack_lenght);

  // Track boolean check parameters ->  Ttrack_HasTP, Ttrack_TPHasSimClu, Ttrack_TPSimCluDirect, Ttrack_hitMTD, Ttrack_hitMatchMTD, Ttrack_hasRecoClu, Ttrack_RecoSimLink;

  track_BDT_tree->Branch("Track_HasTP",&Ttrack_HasTP);
  track_BDT_tree->Branch("Track_TPHasSimClu",&Ttrack_TPHasSimClu);
  track_BDT_tree->Branch("Track_TPSimCluDirect",&Ttrack_TPSimCluDirect);
  track_BDT_tree->Branch("Track_hitMTD",&Ttrack_hitMTD);
  track_BDT_tree->Branch("Track_hit",&Ttrack_hit);
  track_BDT_tree->Branch("Track_hitMatchMTD",&Ttrack_hitMatchMTD);
  track_BDT_tree->Branch("Track_hasRecoClu",&Ttrack_hasRecoClu);
  track_BDT_tree->Branch("Track_RecoSimLink",&Ttrack_RecoSimLink);
  track_BDT_tree->Branch("Track_FullMatch",&Ttrack_FullMatch); 

  track_BDT_tree->Branch("Track_det",&Ttrack_det);
  track_BDT_tree->Branch("Track_subdet",&Ttrack_subdet);
  track_BDT_tree->Branch("Track_subdetmtd",&Ttrack_subdetmtd);
  track_BDT_tree->Branch("Track_nSimClust",&Ttrack_nSimClust);

  track_BDT_tree->Branch("Track_CluX1",&Ttrack_CluX1);
  track_BDT_tree->Branch("Track_CluX2",&Ttrack_CluX2);
  track_BDT_tree->Branch("Track_CluY1",&Ttrack_CluY1);
  track_BDT_tree->Branch("Track_CluY2",&Ttrack_CluY2);
  track_BDT_tree->Branch("Track_CluTime1",&Ttrack_CluTime1);
  track_BDT_tree->Branch("Track_CluTime2",&Ttrack_CluTime2);
}

MtdTracksValidation::~MtdTracksValidation() {}

// ------------ method called for each event  ------------
void MtdTracksValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));

  auto btlRecCluHandle = makeValid(iEvent.getHandle(btlRecCluToken_));
  auto etlRecCluHandle = makeValid(iEvent.getHandle(etlRecCluToken_));
  auto mtdTrkHitHandle = makeValid(iEvent.getHandle(mtdTrackingHitToken_));

  std::unordered_map<uint32_t, MTDHit> m_btlHits;
  std::unordered_map<uint32_t, MTDHit> m_etlHits;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_btlTrkPerCell;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_etlTrkPerCell;
  const auto& tp2SimAssociationMap = iEvent.get(tp2SimAssociationMapToken_);

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& SigmaTofPi = iEvent.get(SigmaTofPiToken_);
  const auto& SigmaTofK = iEvent.get(SigmaTofKToken_);
  const auto& SigmaTofP = iEvent.get(SigmaTofPToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);

  const auto& primRecoVtx = *(RecVertexHandle.product()->begin());
  std::vector<reco::Vertex> vertices = *RecVertexHandle;

  const auto& btlMatchChi2 = iEvent.get(btlMatchChi2Token_);
  const auto& btlMatchTimeChi2 = iEvent.get(btlMatchTimeChi2Token_);
  const auto& etlMatchChi2 = iEvent.get(etlMatchChi2Token_);
  const auto& etlMatchTimeChi2 = iEvent.get(etlMatchTimeChi2Token_);
  const auto& npixBarrel = iEvent.get(npixBarrelToken_);
  const auto& npixEndcap = iEvent.get(npixEndcapToken_);
  const auto& r2sAssociationMap = iEvent.get(r2sAssociationMapToken_);

  // clearing vector for BDT output
  if(OutputForBDTinput_){
   
    Ttrack_pt.clear();
    Ttrack_eta.clear();
    Ttrack_phi.clear();
    Ttrack_dz.clear();
    Ttrack_dxy.clear();
    Ttrack_chi2.clear();
    Ttrack_ndof.clear();
    Ttrack_nValidHits.clear();

    TtrackTP_pt.clear();
    TtrackTP_eta.clear();
    TtrackTP_phi.clear();
    TtrackTP_nValidHits.clear();

    Ttrack_npixBarrelValidHits.clear();
    Ttrack_npixEndcapValidHits.clear();
    Ttrack_Signal.clear();
    Ttrack_Associated.clear();

    Ttrack_BTLchi2.clear();
    Ttrack_BTLtime_chi2.clear();
    Ttrack_ETLchi2.clear();
    Ttrack_ETLtime_chi2.clear();
    Ttrack_t0.clear();
    Ttrack_sigmat0.clear();
    Ttrack_Tmtd.clear();
    Ttrack_MtdMVA.clear();
    Ttrack_sigmaTmtd.clear();
    Ttrack_lenght.clear();

    Ttrack_HasTP.clear();
    Ttrack_TPHasSimClu.clear();
    Ttrack_TPSimCluDirect.clear();
    Ttrack_hitMTD.clear();
    Ttrack_hit.clear();
    Ttrack_hitMatchMTD.clear();
    Ttrack_hasRecoClu.clear();
    Ttrack_RecoSimLink.clear();
    Ttrack_FullMatch.clear();

    Ttrack_det.clear();
    Ttrack_subdet.clear();
    Ttrack_subdetmtd.clear();
    Ttrack_nSimClust.clear();

    Ttrack_CluX1.clear();
    Ttrack_CluX2.clear();
    Ttrack_CluY1.clear();
    Ttrack_CluY2.clear();
    Ttrack_CluTime1.clear();
    Ttrack_CluTime2.clear();

  } 

  // generator level information (HepMC format)
  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

  auto pdt = iSetup.getHandle(particleTableToken_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  auto simToRecoH = makeValid(iEvent.getHandle(simToRecoAssociationToken_));
  s2r_ = simToRecoH.product();

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  unsigned int index = 0;

  // flag to select events with reco vertex close to true simulated primary vertex, or PV fake (particle guns)
  const bool isGoodVtx = std::abs(primRecoVtx.z() - zsim) < deltaZcut_ || primRecoVtx.isFake();

  // additional check for PV (For BDT)
  reco::Vertex Vtx_chosen;
  // ~1% of cases, the 1st vertex is the "good" one, instead of 0th.
  for (int iVtx = 0; iVtx < (int)vertices.size(); iVtx++) {
    const reco::Vertex& vertex = vertices.at(iVtx);
    if (!vertex.isFake() && vertex.ndof() >= 4) {
      Vtx_chosen = vertex;
      break;
    }
  }

  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogInfo("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track& track = *mtdTrackref;

    bool isBTL = false;
    bool isETL = false;
    bool twoETLdiscs = false;
    bool noCrack = std::abs(trackGen.eta()) < trackMaxBtlEta_ || std::abs(trackGen.eta()) > trackMinEtlEta_;

    if(OutputForBDTinput_){

      double track_vtx_dz_check = std::abs(trackGen.dz(Vtx_chosen.position())); // let's see if this works

      if(std::abs(trackGen.eta()) < BDT_track_eta_cut && trackGen.pt() > BDT_track_pt_cut && track_vtx_dz_check < BDT_track_dz_cut){
      
        // Inputs for BDT
        

        bool good_association = false;

        const reco::TrackBaseRef trkrefb(trackref);
        auto found = r2s_->find(trkrefb); // Find TP!
        if (found != r2s_->end()){

          Ttrack_HasTP.push_back(true);


          Ttrack_pt.push_back(trackGen.pt());
          Ttrack_phi.push_back(trackGen.phi());
          Ttrack_eta.push_back(trackGen.eta());
          Ttrack_dz.push_back(std::abs(trackGen.dz(Vtx_chosen.position())));
          Ttrack_dxy.push_back(std::abs(trackGen.dxy(Vtx_chosen.position())));
          Ttrack_chi2.push_back(trackGen.chi2());
          Ttrack_ndof.push_back(trackGen.ndof());
          Ttrack_nValidHits.push_back(trackGen.numberOfValidHits());

          Ttrack_npixBarrelValidHits.push_back(npixBarrel[trackref]);
          Ttrack_npixEndcapValidHits.push_back(npixEndcap[trackref]);
          Ttrack_BTLchi2.push_back(btlMatchChi2[trackref]);
          Ttrack_BTLtime_chi2.push_back(btlMatchTimeChi2[trackref]);
          Ttrack_ETLchi2.push_back(etlMatchChi2[trackref]);
          Ttrack_ETLtime_chi2.push_back(etlMatchTimeChi2[trackref]);

          Ttrack_t0.push_back(t0Src[trackref]);
          Ttrack_sigmat0.push_back(Sigmat0Src[trackref]);
          Ttrack_Tmtd.push_back(tMtd[trackref]);
          Ttrack_sigmaTmtd.push_back(SigmatMtd[trackref]);
          Ttrack_lenght.push_back(pathLength[trackref]);
          Ttrack_MtdMVA.push_back(mtdQualMVA[trackref]);

          const auto& tp = (found->val)[0]; // almost all tracks have just one TP, a few have 2.  (can scan through with "for(const auto& tp : found->val)")
          
          TtrackTP_pt.push_back(tp.first->pt());
          TtrackTP_eta.push_back(tp.first->eta());
          TtrackTP_phi.push_back(tp.first->phi());
          TtrackTP_nValidHits.push_back(tp.first->numberOfHits());

          auto simClustersRefs = tp2SimAssociationMap.find(tp.first); // finds a simClusterReference!!
          const bool withMTD = (simClustersRefs != tp2SimAssociationMap.end()); 

          //bool directHit = false; // We find if Sim cluster linked to TP is from direct hits. For good association selection.
          //bool MTD_hit = false;
          //bool recH_to_Mtdhit_match = false;
          //bool hit_clust_match_btl = false;
          //bool hit_clust_match_etl = false;
          //int nGoodCluster = -1;

          // 1) Link track RecHit to MTdTrackingRecHit (I know which RecHits, hit MTD)
          // 2) Get the MTD Reco Cluster from MTDTrackingRecHit info
          // 3) Find the MTD sim cluster that is linked to MTD reco cluster in the previous step
          // 4) Check if the MTD sim cluster found in previous step is the same as MTD Sim cluster that is linked to TP.

          if (withMTD){ // TP link to MTDsimCluster

            // In test file, all TPs had only 1 simCluster linked to them
            Ttrack_TPHasSimClu.push_back(true);
            const auto& SimCluRefs = (simClustersRefs->val)[0];
            if((*SimCluRefs).trackIdOffset() == 0){ // SimCluster linked to TP is from DirectHit!!!
              Ttrack_TPSimCluDirect.push_back(true);
              //Ttrack_nSimClust.push_back(nGoodCluster);
              //directHit = true;

              //for (const auto& hit : trackGen.recHits()) { // scan through track recHits until we find MTD hit // General Tracks
              for (const auto& hit : track.recHits()) { // Extended track with MTD
                if (good_association)
                  continue; // if goodd assoc found, do not go through all the following checks.
                if (hit->isValid() == false)
                  continue;

                Ttrack_hit.push_back(true);

                MTDDetId Hit = hit->geographicalId();
                //DetId hHit = hit->geographicalId();
                Ttrack_det.push_back(Hit.det());
                Ttrack_subdet.push_back(Hit.subdetId());
                Ttrack_subdetmtd.push_back(Hit.mtdSubDetector());
                
                if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1 || Hit.mtdSubDetector() == 2)) { // trackingRecHit is a hit in MTD
                  //MTD_hit = true;
                  Ttrack_hitMTD.push_back(true);

                  const MTDTrackingRecHit* mtdhit1 = static_cast<const MTDTrackingRecHit*>(hit); // Why I can't I access the mtdcluster info directly from TrackingRecHit?
                  const FTLCluster& hit_cluster_check = mtdhit1->mtdCluster();

                  if (abs(track.eta()) < 1.5){ // Should be a BTL cluster
                    for (const auto& DetSetCluBTL : *btlRecCluHandle) { // BTL check
                      if (good_association)
                        break;
                      for (const auto& clusterBTL : DetSetCluBTL) { // Scan throguh btl reco clusters to find a match
                        if (good_association)
                          break;
                        if(isSameCluster(hit_cluster_check , clusterBTL)){ // find the reco Cluster inside the recoCluster collections
                          //hit_clust_match_btl = true;
                          Ttrack_hasRecoClu.push_back(true);

                          edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> clusterRefBTL = edmNew::makeRefTo(btlRecCluHandle, &clusterBTL); // get the reference to reco cluster inside the collections
                          auto itp = r2sAssociationMap.equal_range(clusterRefBTL); // find the linked simCluster
                          if (itp.first != itp.second) { // find the linked simCluster
                            std::vector<MtdSimLayerClusterRef> simClustersRefs_RecoMatchBTL = (*itp.first).second;  // the range of itp.first, itp.second should be always 1
                            Ttrack_RecoSimLink.push_back(true);
                            
                            for (unsigned int i = 0; i < simClustersRefs_RecoMatchBTL.size(); i++) {
                              auto simClusterRef_RecoMatchBTL = simClustersRefs_RecoMatchBTL[i];

                              //Ttrack_CluX1.push_back((*SimCluRefs).x());
                              //Ttrack_CluX2.push_back((*simClusterRef_RecoMatchBTL).x());
                              Ttrack_CluTime1.push_back((*SimCluRefs).simLCTime());
                              Ttrack_CluTime2.push_back((*simClusterRef_RecoMatchBTL).simLCTime());
                              Ttrack_nSimClust.push_back(i);

                              //if(simClusterRef_RecoMatchBTL == SimCluRefs){
                              if((*simClusterRef_RecoMatchBTL).simLCTime() == (*SimCluRefs).simLCTime()){ // check if the sim cluster linked to reco cluster is the same as the one linked to TP.
                                good_association = true;
                                Ttrack_FullMatch.push_back(true);
                                break;
                              }
                            }
                          }
                        }else{
                          continue;
                        } // mtd hit matched to btl reco cluster
                      } // loop through BTL reco clusters
                    } // loop thorugh set of BTL reco clusters  
                  }else{ // Should be an ETL cluster
                    for (const auto& DetSetCluETL : *etlRecCluHandle) { // ETL check
                      if (good_association)
                        break;
                      for (const auto& clusterETL : DetSetCluETL) { // Scan throguh etl reco clusters to find a match
                        if (good_association)
                          break; 
                        if(isSameCluster(hit_cluster_check , clusterETL)){
                          //hit_clust_match_etl = true;
                          Ttrack_hasRecoClu.push_back(true);

                          edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> clusterRefETL = edmNew::makeRefTo(etlRecCluHandle, &clusterETL);
                          auto itp = r2sAssociationMap.equal_range(clusterRefETL);
                          if (itp.first != itp.second) {
                            std::vector<MtdSimLayerClusterRef> simClustersRefs_RecoMatchETL = (*itp.first).second;  // the range of itp.first, itp.second should be always 1
                            Ttrack_RecoSimLink.push_back(true);
                            
                            for (unsigned int i = 0; i < simClustersRefs_RecoMatchETL.size(); i++) {
                              auto simClusterRef_RecoMatchETL = simClustersRefs_RecoMatchETL[i];

                              //Ttrack_CluX1.push_back((*SimCluRefs).x());
                              //Ttrack_CluX2.push_back((*simClusterRef_RecoMatchETL).x());
                              //Ttrack_CluY1.push_back((*SimCluRefs).y());
                              //Ttrack_CluY2.push_back((*simClusterRef_RecoMatchETL).y());
                              Ttrack_CluTime1.push_back((*SimCluRefs).simLCTime());
                              Ttrack_CluTime2.push_back((*simClusterRef_RecoMatchETL).simLCTime());
                              Ttrack_nSimClust.push_back(i);

                              //if(simClusterRef_RecoMatchETL == SimCluRefs){
                              if((*simClusterRef_RecoMatchETL).simLCTime() == (*SimCluRefs).simLCTime()){
                                good_association = true;
                                Ttrack_FullMatch.push_back(true);
                                break;
                              }
                            }
                          }
                        }else{
                          continue;
                        } // mtd hit matched to etl reco cluster
                      } // loop through ETL reco clusters
                    } // loop thorugh set of ETL reco clusters 
                  } // BTL/ETL cluster search split 
                }else{ // trackingRecHit is a hit in MTD
                  continue;
                } // Hits in MTD
              } // Loop through trackHits 
            }
          } // TP linked to a MTD sim Cluster

          if (tp.first->eventId().bunchCrossing() == 0 && tp.first->eventId().event() == 0){ // Signal vs PU seperation
            Ttrack_Signal.push_back(true); // Signal track
          }else{
            Ttrack_Signal.push_back(false); // PU track?
          }

          Ttrack_Associated.push_back(good_association);

        } // Found TP that is matched to the GTrack

        

      } // basic track eta/pT/dz cuts
    } // BDT output section


    if (track.pt() >= trackMinPt_ && std::abs(track.eta()) <= trackMaxEtlEta_) {
      meTracktmtd_->Fill(tMtd[trackref]);
      if (std::round(SigmatMtd[trackref] - Sigmat0Pid[trackref]) != 0) {
        LogWarning("mtdTracks")
            << "TimeError associated to refitted track is different from TimeError stored in tofPID "
               "sigmat0 ValueMap: this should not happen";
      }

      meTrackt0Src_->Fill(t0Src[trackref]);
      meTrackSigmat0Src_->Fill(Sigmat0Src[trackref]);

      meTrackt0Pid_->Fill(t0Pid[trackref]);
      meTrackSigmat0Pid_->Fill(Sigmat0Pid[trackref]);
      meTrackt0SafePid_->Fill(t0Safe[trackref]);
      meTrackSigmat0SafePid_->Fill(Sigmat0Safe[trackref]);
      meTrackMVAQual_->Fill(mtdQualMVA[trackref]);

      meTrackSigmaTof_[0]->Fill(SigmaTofPi[trackref] * 1e3);  //save as ps
      meTrackSigmaTof_[1]->Fill(SigmaTofK[trackref] * 1e3);
      meTrackSigmaTof_[2]->Fill(SigmaTofP[trackref] * 1e3);
      meTrackSigmaTofvsP_[0]->Fill(track.p(), SigmaTofPi[trackref] * 1e3);
      meTrackSigmaTofvsP_[1]->Fill(track.p(), SigmaTofK[trackref] * 1e3);
      meTrackSigmaTofvsP_[2]->Fill(track.p(), SigmaTofP[trackref] * 1e3);

      meTrackPathLenghtvsEta_->Fill(std::abs(track.eta()), pathLength[trackref]);

      if (std::abs(track.eta()) < trackMaxBtlEta_) {
        // --- all BTL tracks (with and without hit in MTD) ---
        meBTLTrackEffEtaTot_->Fill(track.eta());
        meBTLTrackEffPhiTot_->Fill(track.phi());
        meBTLTrackEffPtTot_->Fill(track.pt());

        bool MTDBtl = false;
        int numMTDBtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
            MTDBtl = true;
            numMTDBtlvalidhits++;
          }
        }
        meTrackNumHits_->Fill(numMTDBtlvalidhits);

        // --- keeping only tracks with last hit in MTD ---
        if (MTDBtl == true) {
          isBTL = true;
          meBTLTrackEffEtaMtd_->Fill(track.eta());
          meBTLTrackEffPhiMtd_->Fill(track.phi());
          meBTLTrackEffPtMtd_->Fill(track.pt());
          meBTLTrackRPTime_->Fill(track.t0());
          meBTLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
        }
        if (isBTL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(numMTDBtlvalidhits);
        }
      }  //loop over (geometrical) BTL tracks

      else {
        // --- all ETL tracks (with and without hit in MTD) ---
        if ((track.eta() < -trackMinEtlEta_) && (track.eta() > -trackMaxEtlEta_)) {
          meETLTrackEffEtaTot_[0]->Fill(track.eta());
          meETLTrackEffPhiTot_[0]->Fill(track.phi());
          meETLTrackEffPtTot_[0]->Fill(track.pt());
        }

        if ((track.eta() > trackMinEtlEta_) && (track.eta() < trackMaxEtlEta_)) {
          meETLTrackEffEtaTot_[1]->Fill(track.eta());
          meETLTrackEffPhiTot_[1]->Fill(track.phi());
          meETLTrackEffPtTot_[1]->Fill(track.pt());
        }

        bool MTDEtlZnegD1 = false;
        bool MTDEtlZnegD2 = false;
        bool MTDEtlZposD1 = false;
        bool MTDEtlZposD2 = false;
        int numMTDEtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
            isETL = true;
            ETLDetId ETLHit = hit->geographicalId();

            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZnegD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZposD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
          }
        }
        meTrackNumHits_->Fill(-numMTDEtlvalidhits);
        if (isETL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(-numMTDEtlvalidhits);
        }

        // --- keeping only tracks with last hit in MTD ---
        if ((track.eta() < -trackMinEtlEta_) && (track.eta() > -trackMaxEtlEta_)) {
          twoETLdiscs = (MTDEtlZnegD1 == true) && (MTDEtlZnegD2 == true);
          if ((MTDEtlZnegD1 == true) || (MTDEtlZnegD2 == true)) {
            meETLTrackEffEtaMtd_[0]->Fill(track.eta());
            meETLTrackEffPhiMtd_[0]->Fill(track.phi());
            meETLTrackEffPtMtd_[0]->Fill(track.pt());
            if (twoETLdiscs) {
              meETLTrackEffEta2Mtd_[0]->Fill(track.eta());
              meETLTrackEffPhi2Mtd_[0]->Fill(track.phi());
              meETLTrackEffPt2Mtd_[0]->Fill(track.pt());
            }
          }
        }
        if ((track.eta() > trackMinEtlEta_) && (track.eta() < trackMaxEtlEta_)) {
          twoETLdiscs = (MTDEtlZposD1 == true) && (MTDEtlZposD2 == true);
          if ((MTDEtlZposD1 == true) || (MTDEtlZposD2 == true)) {
            meETLTrackEffEtaMtd_[1]->Fill(track.eta());
            meETLTrackEffPhiMtd_[1]->Fill(track.phi());
            meETLTrackEffPtMtd_[1]->Fill(track.pt());
            if (twoETLdiscs) {
              meETLTrackEffEta2Mtd_[1]->Fill(track.eta());
              meETLTrackEffPhi2Mtd_[1]->Fill(track.phi());
              meETLTrackEffPt2Mtd_[1]->Fill(track.pt());
            }
          }
        }
      }

      LogDebug("MtdTracksValidation") << "Track p/pt = " << track.p() << " " << track.pt() << " eta " << track.eta()
                                      << " BTL " << isBTL << " ETL " << isETL << " 2disks " << twoETLdiscs;

      // TrackingParticle based matching

      const reco::TrackBaseRef trkrefb(trackref);
      auto tp_info = getMatchedTP(trkrefb);

      meTrackPtTot_->Fill(trackGen.pt());
      meTrackEtaTot_->Fill(std::abs(trackGen.eta()));
      if (tp_info != nullptr && mvaTPSelAll(**tp_info)) {
        if (track.pt() < 12.) {
          if (isBTL) {
            meBTLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                               std::abs(trackGen.pt() - (*tp_info)->pt()));
            meBTLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meBTLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meBTLTrackMatchedTPPtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meBTLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                 (trackGen.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
            meBTLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(), (track.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
          }
          if (isETL && !twoETLdiscs) {
            meETLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                               std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meETLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meETLTrackMatchedTPPtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                 (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            meETLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                 (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
          }
          if (isETL && twoETLdiscs) {
            meETLTrackMatchedTP2PtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                                std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTP2PtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meETLTrackMatchedTP2PtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meETLTrackMatchedTP2PtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTP2DPtvsPtGen_->Fill((*tp_info)->pt(),
                                                  (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            meETLTrackMatchedTP2DPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                  (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
          }
        }
        auto simClustersRefs = tp2SimAssociationMap.find(*tp_info);
        const bool withMTD = (simClustersRefs != tp2SimAssociationMap.end());
        if (noCrack) {
          meTrackMatchedTPEffPtTot_->Fill(trackGen.pt());
          if (withMTD) {
            meTrackMatchedTPmtdEffPtTot_->Fill(trackGen.pt());
          }
        }
        meTrackMatchedTPEffEtaTot_->Fill(std::abs(trackGen.eta()));
        if (withMTD) {
          meTrackMatchedTPmtdEffEtaTot_->Fill(std::abs(trackGen.eta()));
        }
        if (isBTL || isETL) {
          if (noCrack) {
            meTrackMatchedTPEffPtMtd_->Fill(trackGen.pt());
            if (isBTL || twoETLdiscs) {
              meTrackMatchedTPEffPtEtl2Mtd_->Fill(trackGen.pt());
            }
            if (withMTD) {
              meTrackMatchedTPmtdEffPtMtd_->Fill(trackGen.pt());
            }
          }
          meTrackMatchedTPEffEtaMtd_->Fill(std::abs(trackGen.eta()));
          if (isBTL || twoETLdiscs) {
            meTrackMatchedTPEffEtaEtl2Mtd_->Fill(std::abs(trackGen.eta()));
          }
          if (withMTD) {
            meTrackMatchedTPmtdEffEtaMtd_->Fill(std::abs(trackGen.eta()));
          }
        }

        // detailed extrapolation check only on tracks associated to TP from signal event
        if (!mvaTPSelLV(**tp_info)) {
          continue;
        }
        size_t nlayers(0);
        float extrho(0.);
        float exteta(0.);
        float extphi(0.);
        float selvar(0.);
        auto accept = checkAcceptance(trackGen, iEvent, iSetup, nlayers, extrho, exteta, extphi, selvar);
        if (accept.first && std::abs(exteta) < trackMaxBtlEta_) {
          meExtraPhiAtBTL_->Fill(angle_units::operators::convertRadToDeg(extphi));
          meExtraBTLeneInCone_->Fill(selvar);
        }
        if (accept.second) {
          if (std::abs(exteta) < trackMaxBtlEta_) {
            meExtraPhiAtBTLmatched_->Fill(angle_units::operators::convertRadToDeg(extphi));
          }
          if (noCrack) {
            meExtraPtMtd_->Fill(trackGen.pt());
            if (nlayers == 2) {
              meExtraPtEtl2Mtd_->Fill(trackGen.pt());
            }
          }
          meExtraEtaMtd_->Fill(std::abs(trackGen.eta()));
          if (nlayers == 2) {
            meExtraEtaEtl2Mtd_->Fill(std::abs(trackGen.eta()));
          }
          if (accept.first && accept.second && !(isBTL || isETL)) {
            edm::LogInfo("MtdTracksValidation")
                << "MtdTracksValidation: extender fail in " << iEvent.id().run() << " " << iEvent.id().event()
                << " pt= " << trackGen.pt() << " eta= " << trackGen.eta();
            meExtraMTDfailExtenderEta_->Fill(std::abs(trackGen.eta()));
            if (noCrack) {
              meExtraMTDfailExtenderPt_->Fill(trackGen.pt());
            }
          }
        }

      }  // TP matching
    }

    if (isGoodVtx) {
      const bool vtxFake = primRecoVtx.isFake();

      if (mvaRecSel(trackGen, primRecoVtx, t0Safe[trackref], Sigmat0Safe[trackref])) {
        // reco-gen matching used for MVA quality flag

        if (noCrack) {
          meMVATrackEffPtTot_->Fill(trackGen.pt());
        }
        meMVATrackEffEtaTot_->Fill(std::abs(trackGen.eta()));

        double dZ = trackGen.vz() - zsim;
        double dT(-9999.);
        double pullT(-9999.);
        if (Sigmat0Safe[trackref] != -1.) {
          dT = t0Safe[trackref] - tsim;
          pullT = dT / Sigmat0Safe[trackref];
        }
        for (const auto& genP : mc->particle_range()) {
          // select status 1 genParticles and match them to the reconstructed track

          float charge = pdTable->particle(HepPDT::ParticleID(genP->pdg_id())) != nullptr
                             ? pdTable->particle(HepPDT::ParticleID(genP->pdg_id()))->charge()
                             : 0.f;
          if (mvaGenSel(*genP, charge)) {
            if (mvaGenRecMatch(*genP, zsim, trackGen, vtxFake)) {
              meMVATrackZposResTot_->Fill(dZ);
              if (noCrack) {
                meMVATrackMatchedEffPtTot_->Fill(trackGen.pt());
              }
              meMVATrackMatchedEffEtaTot_->Fill(std::abs(trackGen.eta()));
              if (isBTL || isETL) {
                meMVATrackResTot_->Fill(dT);
                meMVATrackPullTot_->Fill(pullT);
                if (noCrack) {
                  meMVATrackMatchedEffPtMtd_->Fill(trackGen.pt());
                }
                meMVATrackMatchedEffEtaMtd_->Fill(std::abs(trackGen.eta()));
              }
              break;
            }
          }
        }
      }
    }  // MC truth matich analysis for good PV
  }    //RECO tracks loop
  if(OutputForBDTinput_){
    track_BDT_tree->Fill();
  }
}

const std::pair<bool, bool> MtdTracksValidation::checkAcceptance(const reco::Track& track,
                                                                 const edm::Event& iEvent,
                                                                 edm::EventSetup const& iSetup,
                                                                 size_t& nlayers,
                                                                 float& extrho,
                                                                 float& exteta,
                                                                 float& extphi,
                                                                 float& selvar) {
  bool isMatched(false);
  nlayers = 0;
  extrho = 0.;
  exteta = -999.;
  extphi = -999.;
  selvar = 0.;

  auto geometryHandle = iSetup.getTransientHandle(mtdgeoToken_);
  const MTDGeometry* geom = geometryHandle.product();
  auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();

  auto layerHandle = iSetup.getTransientHandle(mtdlayerToken_);
  const MTDDetLayerGeometry* layerGeo = layerHandle.product();

  auto magfieldHandle = iSetup.getTransientHandle(magfieldToken_);
  const MagneticField* mfield = magfieldHandle.product();

  auto ttrackBuilder = iSetup.getTransientHandle(builderToken_);

  auto tTrack = ttrackBuilder->build(track);
  TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
  float theMaxChi2 = 500.;
  float theNSigma = 10.;
  std::unique_ptr<MeasurementEstimator> theEstimator =
      std::make_unique<Chi2MeasurementEstimator>(theMaxChi2, theNSigma);
  SteppingHelixPropagator prop(mfield, anyDirection);

  auto btlRecHitsHandle = makeValid(iEvent.getHandle(btlRecHitsToken_));
  auto etlRecHitsHandle = makeValid(iEvent.getHandle(etlRecHitsToken_));

  edm::LogVerbatim("MtdTracksValidation")
      << "MtdTracksValidation: extrapolating track, pt= " << track.pt() << " eta= " << track.eta();

  //try BTL
  bool inBTL = false;
  float eneSum(0.);
  const std::vector<const DetLayer*>& layersBTL = layerGeo->allBTLLayers();
  for (const DetLayer* ilay : layersBTL) {
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inBTL) {
      inBTL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
      edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at BTL surface, rho= " << extrho
                                              << " eta= " << exteta << " phi= " << extphi;
    }
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check energy in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *btlRecHitsHandle) {
        BTLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(0., 0., 0.);
          local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " ene= " << recHit.energy()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            eneSum += recHit.energy();
            //extrho = detWithState.second.globalPosition().perp();
            //exteta = detWithState.second.globalPosition().eta();
            //extphi = detWithState.second.globalPosition().phi();
          }
        }
      }
    }
    if (eneSum > depositBTLthreshold_) {
      nlayers++;
      selvar = eneSum;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: BTL matched, energy= " << eneSum << " #layers= " << nlayers;
    }
  }
  if (inBTL) {
    return std::make_pair(inBTL, isMatched);
  }

  //try ETL
  bool inETL = false;
  const std::vector<const DetLayer*>& layersETL = layerGeo->allETLLayers();
  for (const DetLayer* ilay : layersETL) {
    size_t hcount(0);
    const BoundDisk& disk = static_cast<const MTDSectorForwardDoubleLayer*>(ilay)->specificSurface();
    const double diskZ = disk.position().z();
    if (tsos.globalPosition().z() * diskZ < 0)
      continue;  // only propagate to the disk that's on the same side
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inETL) {
      inETL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
    }
    edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at ETL surface, rho= " << extrho
                                            << " eta= " << exteta << " phi= " << extphi;
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check hits in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *etlRecHitsHandle) {
        ETLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId();
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(topo.localX(recHit.row()), topo.localY(recHit.column()), 0.);
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " time= " << recHit.time()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            hcount++;
            if (hcount == 1) {
              //extrho = detWithState.second.globalPosition().perp();
              //exteta = detWithState.second.globalPosition().eta();
              //extphi = detWithState.second.globalPosition().phi();
            }
          }
        }
      }
    }
    if (hcount > 0) {
      nlayers++;
      selvar = (float)hcount;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: ETL matched, counts= " << hcount << " #layers= " << nlayers;
    }
  }

  if (!inBTL && !inETL) {
    edm::LogVerbatim("MtdTracksValidation")
        << "MtdTracksValidation: track not extrapolating to MTD: pt= " << track.pt() << " eta= " << track.eta()
        << " phi= " << track.phi() << " vz= " << track.vz()
        << " vxy= " << std::sqrt(track.vx() * track.vx() + track.vy() * track.vy());
  }
  return std::make_pair(inETL, isMatched);
}

// ------------ method for histogram booking ------------
void MtdTracksValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meBTLTrackRPTime_ = ibook.book1D("TrackBTLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meBTLTrackEffEtaTot_ = ibook.book1D("TrackBTLEffEtaTot", "Track efficiency vs eta (Tot);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiTot_ =
      ibook.book1D("TrackBTLEffPhiTot", "Track efficiency vs phi (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtTot_ = ibook.book1D("TrackBTLEffPtTot", "Track efficiency vs pt (Tot);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackEffEtaMtd_ = ibook.book1D("TrackBTLEffEtaMtd", "Track efficiency vs eta (Mtd);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiMtd_ =
      ibook.book1D("TrackBTLEffPhiMtd", "Track efficiency vs phi (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtMtd_ = ibook.book1D("TrackBTLEffPtMtd", "Track efficiency vs pt (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackPtRes_ =
      ibook.book1D("TrackBTLPtRes", "Track pT resolution  ;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);
  meETLTrackRPTime_ = ibook.book1D("TrackETLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meETLTrackEffEtaTot_[0] =
      ibook.book1D("TrackETLEffEtaTotZneg", "Track efficiency vs eta (Tot) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaTot_[1] =
      ibook.book1D("TrackETLEffEtaTotZpos", "Track efficiency vs eta (Tot) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiTot_[0] =
      ibook.book1D("TrackETLEffPhiTotZneg", "Track efficiency vs phi (Tot) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiTot_[1] =
      ibook.book1D("TrackETLEffPhiTotZpos", "Track efficiency vs phi (Tot) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtTot_[0] =
      ibook.book1D("TrackETLEffPtTotZneg", "Track efficiency vs pt (Tot) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtTot_[1] =
      ibook.book1D("TrackETLEffPtTotZpos", "Track efficiency vs pt (Tot) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEtaMtd_[0] =
      ibook.book1D("TrackETLEffEtaMtdZneg", "Track efficiency vs eta (Mtd) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaMtd_[1] =
      ibook.book1D("TrackETLEffEtaMtdZpos", "Track efficiency vs eta (Mtd) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiMtd_[0] =
      ibook.book1D("TrackETLEffPhiMtdZneg", "Track efficiency vs phi (Mtd) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiMtd_[1] =
      ibook.book1D("TrackETLEffPhiMtdZpos", "Track efficiency vs phi (Mtd) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtMtd_[0] =
      ibook.book1D("TrackETLEffPtMtdZneg", "Track efficiency vs pt (Mtd) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtMtd_[1] =
      ibook.book1D("TrackETLEffPtMtdZpos", "Track efficiency vs pt (Mtd) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEta2Mtd_[0] =
      ibook.book1D("TrackETLEffEta2MtdZneg", "Track efficiency vs eta (Mtd 2 hit) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEta2Mtd_[1] =
      ibook.book1D("TrackETLEffEta2MtdZpos", "Track efficiency vs eta (Mtd 2 hit) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhi2Mtd_[0] = ibook.book1D(
      "TrackETLEffPhi2MtdZneg", "Track efficiency vs phi (Mtd 2 hit) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhi2Mtd_[1] = ibook.book1D(
      "TrackETLEffPhi2MtdZpos", "Track efficiency vs phi (Mtd 2 hit) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPt2Mtd_[0] =
      ibook.book1D("TrackETLEffPt2MtdZneg", "Track efficiency vs pt (Mtd 2 hit) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPt2Mtd_[1] =
      ibook.book1D("TrackETLEffPt2MtdZpos", "Track efficiency vs pt (Mtd 2 hit) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackPtRes_ =
      ibook.book1D("TrackETLPtRes", "Track pT resolution;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);

  meTracktmtd_ = ibook.book1D("Tracktmtd", "Track time from TrackExtenderWithMTD;tmtd [ns]", 150, 1, 16);
  meTrackt0Src_ = ibook.book1D("Trackt0Src", "Track time from TrackExtenderWithMTD;t0Src [ns]", 100, -1.5, 1.5);
  meTrackSigmat0Src_ =
      ibook.book1D("TrackSigmat0Src", "Time Error from TrackExtenderWithMTD; #sigma_{t0Src} [ns]", 100, 0, 0.1);

  meTrackt0Pid_ = ibook.book1D("Trackt0Pid", "Track t0 as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0Pid_ = ibook.book1D("TrackSigmat0Pid", "Sigmat0 as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackt0SafePid_ = ibook.book1D("Trackt0SafePID", "Track t0 Safe as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0SafePid_ =
      ibook.book1D("TrackSigmat0SafePID", "Sigmat0 Safe as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackNumHits_ = ibook.book1D("TrackNumHits", "Number of valid MTD hits per track ; Number of hits", 10, -5, 5);
  meTrackNumHitsNT_ = ibook.book1D(
      "TrackNumHitsNT", "Number of valid MTD hits per track no time associated; Number of hits", 10, -5, 5);
  meTrackMVAQual_ = ibook.book1D("TrackMVAQual", "Track MVA Quality as stored in Value Map ; MVAQual", 100, 0, 1);
  meTrackPathLenghtvsEta_ = ibook.bookProfile(
      "TrackPathLenghtvsEta", "MTD Track pathlength vs MTD track Eta;|#eta|;Pathlength", 100, 0, 3.2, 100.0, 400.0, "S");

  meTrackSigmaTof_[0] =
      ibook.book1D("TrackSigmaTof_Pion", "Sigma(TOF) for pion hypothesis; #sigma_{t0} [ps]", 10, 0, 5);
  meTrackSigmaTof_[1] =
      ibook.book1D("TrackSigmaTof_Kaon", "Sigma(TOF) for kaon hypothesis; #sigma_{t0} [ps]", 25, 0, 25);
  meTrackSigmaTof_[2] =
      ibook.book1D("TrackSigmaTof_Proton", "Sigma(TOF) for proton hypothesis; #sigma_{t0} [ps]", 50, 0, 50);

  meTrackSigmaTofvsP_[0] = ibook.bookProfile("TrackSigmaTofvsP_Pion",
                                             "Sigma(TOF) for pion hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");
  meTrackSigmaTofvsP_[1] = ibook.bookProfile("TrackSigmaTofvsP_Kaon",
                                             "Sigma(TOF) for kaon hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");
  meTrackSigmaTofvsP_[2] = ibook.bookProfile("TrackSigmaTofvsP_Proton",
                                             "Sigma(TOF) for proton hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");

  meMVATrackEffPtTot_ = ibook.book1D("MVAEffPtTot", "Pt of tracks associated to LV; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtTot_ =
      ibook.book1D("MVAMatchedEffPtTot", "Pt of tracks associated to LV matched to GEN; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtMtd_ = ibook.book1D(
      "MVAMatchedEffPtMtd", "Pt of tracks associated to LV matched to GEN with time; track pt [GeV] ", 110, 0., 11.);

  meExtraPtMtd_ = ibook.book1D("ExtraPtMtd", "Pt of tracks extrapolated to hits; track pt [GeV] ", 110, 0., 11.);
  meExtraPtEtl2Mtd_ =
      ibook.book1D("ExtraPtEtl2Mtd", "Pt of tracks extrapolated to hits, 2 ETL layers; track pt [GeV] ", 110, 0., 11.);

  meTrackPtTot_ = ibook.book1D("TrackPtTot", "Pt of tracks ; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtTot_ =
      ibook.book1D("MatchedTPEffPtTot", "Pt of tracks  matched to TP; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtMtd_ =
      ibook.book1D("MatchedTPEffPtMtd", "Pt of tracks  matched to TP with time; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtEtl2Mtd_ = ibook.book1D(
      "MatchedTPEffPtEtl2Mtd", "Pt of tracks  matched to TP with time, 2 ETL hits; track pt [GeV] ", 110, 0., 11.);

  meBTLTrackMatchedTPPtResMtd_ = ibook.book1D(
      "TrackMatchedTPBTLPtResMtd",
      "Pt resolution of tracks matched to TP-BTL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meETLTrackMatchedTPPtResMtd_ = ibook.book1D(
      "TrackMatchedTPETLPtResMtd",
      "Pt resolution of tracks matched to TP-ETL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meETLTrackMatchedTP2PtResMtd_ = ibook.book1D(
      "TrackMatchedTPETL2PtResMtd",
      "Pt resolution of tracks matched to TP-ETL 2hits  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meBTLTrackMatchedTPPtRatioGen_ = ibook.book1D(
      "TrackMatchedTPBTLPtRatioGen", "Pt ratio of Gentracks (BTL)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meETLTrackMatchedTPPtRatioGen_ = ibook.book1D(
      "TrackMatchedTPETLPtRatioGen", "Pt ratio of Gentracks (ETL 1hit)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meETLTrackMatchedTP2PtRatioGen_ = ibook.book1D(
      "TrackMatchedTPETL2PtRatioGen", "Pt ratio of Gentracks (ETL 2hits)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meBTLTrackMatchedTPPtRatioMtd_ = ibook.book1D("TrackMatchedTPBTLPtRatioMtd",
                                                "Pt ratio of tracks matched to TP-BTL hits  ;pT_{MTDtrack}/pT_{truth} ",
                                                100,
                                                0.9,
                                                1.1);
  meETLTrackMatchedTPPtRatioMtd_ = ibook.book1D("TrackMatchedTPETLPtRatioMtd",
                                                "Pt ratio of tracks matched to TP-ETL hits  ;pT_{MTDtrack}/pT_{truth} ",
                                                100,
                                                0.9,
                                                1.1);
  meETLTrackMatchedTP2PtRatioMtd_ =
      ibook.book1D("TrackMatchedTPETL2PtRatioMtd",
                   "Pt ratio of tracks matched to TP-ETL 2hits  ;pT_{MTDtrack}/pT_{truth} ",
                   100,
                   0.9,
                   1.1);
  meBTLTrackMatchedTPPtResvsPtMtd_ = ibook.bookProfile("TrackMatchedTPBTLPtResvsPtMtd",
                                                       "Pt resolution of tracks matched to TP-BTL hit vs Pt;pT_{truth} "
                                                       "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       0.,
                                                       4.,
                                                       "s");
  meETLTrackMatchedTPPtResvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETLPtResvsPtMtd",
                                                       "Pt resolution of tracks matched to TP-ETL hit vs Pt;pT_{truth} "
                                                       "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       0.,
                                                       4.,
                                                       "s");
  meETLTrackMatchedTP2PtResvsPtMtd_ =
      ibook.bookProfile("TrackMatchedTPETL2PtResvsPtMtd",
                        "Pt resolution of tracks matched to TP-ETL 2hits Pt pT;pT_{truth} "
                        "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                        20,
                        0.7,
                        10.,
                        0.,
                        4.,
                        "s");
  meBTLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPBTLDPtvsPtGen",
      "Pt relative difference of Gentracks (BTL) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meETLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPETLDPtvsPtGen",
      "Pt relative difference of Gentracks (ETL 1hit) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meETLTrackMatchedTP2DPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPETL2DPtvsPtGen",
      "Pt relative difference  of Gentracks (ETL 2hits) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meBTLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPBTLDPtvsPtMtd",
                                                     "Pt relative difference of tracks matched to TP-BTL hits vs "
                                                     "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                     20,
                                                     0.7,
                                                     10.,
                                                     -0.1,
                                                     0.1,
                                                     "s");
  meETLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETLDPtvsPtMtd",
                                                     "Pt relative difference of tracks matched to TP-ETL hits vs "
                                                     "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                     20,
                                                     0.7,
                                                     10.,
                                                     -0.1,
                                                     0.1,
                                                     "s");
  meETLTrackMatchedTP2DPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETL2DPtvsPtMtd",
                                                      "Pt relative difference of tracks matched to TP-ETL 2hits vs "
                                                      "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                      20,
                                                      0.7,
                                                      10.,
                                                      -0.1,
                                                      0.1,
                                                      "s");

  meTrackMatchedTPmtdEffPtTot_ =
      ibook.book1D("MatchedTPmtdEffPtTot", "Pt of tracks  matched to TP-mtd hit; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPmtdEffPtMtd_ = ibook.book1D(
      "MatchedTPmtdEffPtMtd", "Pt of tracks  matched to TP-mtd hit with time; track pt [GeV] ", 110, 0., 11.);

  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Eta of tracks associated to LV; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaTot_ =
      ibook.book1D("MVAMatchedEffEtaTot", "Eta of tracks associated to LV matched to GEN; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaMtd_ = ibook.book1D(
      "MVAMatchedEffEtaMtd", "Eta of tracks associated to LV matched to GEN with time; track eta ", 66, 0., 3.3);

  meExtraEtaMtd_ = ibook.book1D("ExtraEtaMtd", "Eta of tracks extrapolated to hits; track eta ", 66, 0., 3.3);
  meExtraEtaEtl2Mtd_ =
      ibook.book1D("ExtraEtaEtl2Mtd", "Eta of tracks extrapolated to hits, 2 ETL layers; track eta ", 66, 0., 3.3);

  meTrackEtaTot_ = ibook.book1D("TrackEtaTot", "Eta of tracks ; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaTot_ =
      ibook.book1D("MatchedTPEffEtaTot", "Eta of tracks  matched to TP; track eta ", 66, 0., 3.3);
  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Eta of tracks ; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaMtd_ =
      ibook.book1D("MatchedTPEffEtaMtd", "Eta of tracks  matched to TP with time; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaEtl2Mtd_ = ibook.book1D(
      "MatchedTPEffEtaEtl2Mtd", "Eta of tracks  matched to TP with time, 2 ETL hits; track eta ", 66, 0., 3.3);

  meTrackMatchedTPmtdEffEtaTot_ =
      ibook.book1D("MatchedTPmtdEffEtaTot", "Eta of tracks  matched to TP-mtd hit; track eta ", 66, 0., 3.3);
  meTrackMatchedTPmtdEffEtaMtd_ =
      ibook.book1D("MatchedTPmtdEffEtaMtd", "Eta of tracks  matched to TP-mtd hit with time; track eta ", 66, 0., 3.3);

  meMVATrackResTot_ = ibook.book1D(
      "MVATrackRes", "t_{rec} - t_{sim} for LV associated tracks; t_{rec} - t_{sim} [ns] ", 120, -0.15, 0.15);
  meMVATrackPullTot_ =
      ibook.book1D("MVATrackPull", "Pull for associated tracks; (t_{rec}-t_{sim})/#sigma_{t}", 50, -5., 5.);
  meMVATrackZposResTot_ = ibook.book1D(
      "MVATrackZposResTot", "Z_{PCA} - Z_{sim} for associated tracks;Z_{PCA} - Z_{sim} [cm] ", 100, -0.1, 0.1);

  meExtraPhiAtBTL_ =
      ibook.book1D("ExtraPhiAtBTL", "Phi at BTL surface of extrapolated tracks; phi [deg]", 720, -180., 180.);
  meExtraPhiAtBTLmatched_ = ibook.book1D("ExtraPhiAtBTLmatched",
                                         "Phi at BTL surface of extrapolated tracksi matched with BTL hits; phi [deg]",
                                         720,
                                         -180.,
                                         180.);
  meExtraBTLeneInCone_ = ibook.book1D(
      "ExtraBTLeneInCone", "BTL reconstructed energy in cone arounnd extrapolated track; E [MeV]", 100, 0., 50.);
  meExtraMTDfailExtenderEta_ =
      ibook.book1D("ExtraMTDfailExtenderEta",
                   "Eta of tracks extrapolated to MTD with no track extender match to hits; track eta",
                   66,
                   0.,
                   3.3);
  ;
  meExtraMTDfailExtenderPt_ =
      ibook.book1D("ExtraMTDfailExtenderPt",
                   "Pt of tracks extrapolated to MTD with no track extender match to hits; track pt [GeV] ",
                   110,
                   0.,
                   11.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdTracksValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tp2SimAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("btlSimHits", edm::InputTag("mix", "g4SimHitsFastTimerHitsBarrel"));
  desc.add<edm::InputTag>("etlSimHits", edm::InputTag("mix", "g4SimHitsFastTimerHitsEndcap"));
  desc.add<edm::InputTag>("btlRecHits", edm::InputTag("mtdRecHits", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecHits", edm::InputTag("mtdRecHits", "FTLEndcap"));
  desc.add<edm::InputTag>("recCluTagBTL", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("recCluTagETL", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("trkHitTag", edm::InputTag("mtdTrackingRecHits"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD:btlMatchChi2"));
  desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD:btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD:etlMatchChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD:etlMatchTimeChi2"));
  desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD:npixBarrel"));
  desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD:npixEndcap"));
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("sigmaTofPi", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmaTofK", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmaTofP", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<double>("trackMinimumPt", 0.7);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);
  desc.addUntracked<bool>("optionalPlots", true);
  desc.add<bool>("OutputForBDTinput",true);

  descriptions.add("mtdTracksValid", desc);
}

const bool MtdTracksValidation::mvaGenSel(const HepMC::GenParticle& gp, const float& charge) {
  bool match = false;
  if (gp.status() != 1) {
    return match;
  }
  match = charge != 0.f && gp.momentum().perp() > pTcut_ && std::abs(gp.momentum().eta()) < etacutGEN_;
  return match;
}

const bool MtdTracksValidation::mvaTPSelLV(const TrackingParticle& tp) {
  bool match = (tp.status() != 1) ? false : true;
  return match;
}

const bool MtdTracksValidation::mvaTPSelAll(const TrackingParticle& tp) {
  bool match = false;

  auto x_pv = tp.parentVertex()->position().x();
  auto y_pv = tp.parentVertex()->position().y();
  auto z_pv = tp.parentVertex()->position().z();

  auto r_pv = std::sqrt(x_pv * x_pv + y_pv * y_pv);

  match =
      tp.charge() != 0 && tp.pt() > pTcut_ && std::abs(tp.eta()) < etacutGEN_ && r_pv < rBTL_ && std::abs(z_pv) < zETL_;
  return match;
}

const bool MtdTracksValidation::mvaRecSel(const reco::TrackBase& trk,
                                          const reco::Vertex& vtx,
                                          const double& t0,
                                          const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ &&
          (std::abs(trk.vz() - vtx.z()) <= deltaZcut_ || vtx.isFake());
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}

const bool MtdTracksValidation::mvaGenRecMatch(const HepMC::GenParticle& genP,
                                               const double& zsim,
                                               const reco::TrackBase& trk,
                                               const bool& vtxFake) {
  bool match = false;
  double dR = reco::deltaR(genP.momentum(), trk.momentum());
  double genPT = genP.momentum().perp();
  match = std::abs(genPT - trk.pt()) < trk.pt() * deltaPTcut_ && dR < deltaDRcut_ &&
          (std::abs(trk.vz() - zsim) < deltaZcut_ || vtxFake);
  return match;
}

const edm::Ref<std::vector<TrackingParticle>>* MtdTracksValidation::getMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP associated to in time events
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0)
      return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

DEFINE_FWK_MODULE(MtdTracksValidation);

/**/
