//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TLorentzVector.h>
//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

/////////////
//// Fast Jet
#include "RecoJets/JetProducers/plugins/VirtualJetProducer.h"
//#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
//#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
//#include <fastjet/ClusterSequenceAreaBase.hh>
//#include "fastjet/ClusterSequenceArea.hh"
//#include "RecoJets/JetProducers/plugins/FixedGridRhoProducerFastjet.h"
// #include <fastjet/ClusterSequenceAreaBase.hh>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackNtupleMaker : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackNtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~L1TrackNtupleMaker();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void FillJets(std::vector<fastjet::PseudoJet>  JetInputs_, bool Prompt, bool TrueTP, edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle);  
  virtual double FillRecoPrimaryVtx();
  virtual double FillPrimaryVtx(edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle,edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle, const TrackerTopology* const tTopo, const TrackerGeometry* const theTrackerGeom);
  virtual void FillClusteringMaps();  
  virtual void FindEtaClusters();  
  //virtual void FindL1L2Clusters();  
  virtual void FindClusters();  
  virtual float FindMaxHTSlice();  
  virtual void MergePromptCandidates();  
protected:
  
private:
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)
  int L1Tk_nPar;        // use 4 or 5 parameter track fit? 
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt 
  double TP_maxEta;     // save TPs with |eta| < maxEta 
  double TP_maxZ0;      // save TPs with |z0| < maxZ0 
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  
  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag TrackingVertexInputTag;
  edm::InputTag TrueVertexInputTag;
  
  edm::InputTag GenJetAK4;

  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > ttClusterToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > ttClusterMCTruthToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT< vector<reco::GenParticle> > HEPMCVertexToken_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > TrackingVertexToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> >GenJetCollectionToken_;
  std::vector<fastjet::PseudoJet>  RecoJetInputs_;
  std::vector<fastjet::PseudoJet>  TrueJetInputs_;
  std::vector<fastjet::PseudoJet>  JetInputs_;
  std::vector<fastjet::PseudoJet>  JetInputsPU1_;
  std::vector<fastjet::PseudoJet>  JetInputsPU2_;
  std::vector<fastjet::PseudoJet>  JetInputsPU3_;
  std::vector<fastjet::PseudoJet>  JetInputsPU4_;
  std::vector<fastjet::PseudoJet>  JetInputsPU5_;
  std::vector<fastjet::PseudoJet>  JetInputsPU6_;
  std::vector<fastjet::PseudoJet>  JetInputsPU7_;
  std::vector<fastjet::PseudoJet>  JetInputsPU8_;
  std::vector<fastjet::PseudoJet>  JetInputsPU9_;
  std::vector<fastjet::PseudoJet>  JetInputsPU10_;
  std::vector<fastjet::PseudoJet>  JetInputsPU11_;
  std::vector<fastjet::PseudoJet>  JetInputsPU12_;
  std::vector<fastjet::PseudoJet>  JetInputsPU13_;
  std::vector<fastjet::PseudoJet>  JetInputsPU14_;
  std::vector<fastjet::PseudoJet>  JetInputsPU15_;
  std::vector<fastjet::PseudoJet>  JetInputsPU16_;
  std::vector<fastjet::PseudoJet>  JetInputsPU17_;
  std::vector<fastjet::PseudoJet>  JetInputsPU18_;
  std::vector<fastjet::PseudoJet>  JetInputsPU19_;
  std::vector<fastjet::PseudoJet>  JetInputsPU20_;


  std::vector<fastjet::PseudoJet> JetOutputs_;
  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;
  // primary vertex
  std::vector<float>* m_pv_HTz;
  std::vector<float>* m_pv_L1reco;
  std::vector<float>* m_pv_L1;
  std::vector<float>* m_pv_MC;
  std::vector<float>* m_pv_MC_vr;
  std::vector<int>* m_MC_lep;
  // all L1 tracks

  std::vector<float>* m_trk_p;
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_d0;   // (filled if L1Tk_nPar==5, else 999)
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2; 
  std::vector<int>*   m_trk_nstub;
  std::vector<int>*   m_trk_genuine;
  std::vector<int>*   m_trk_loose;
  std::vector<int>*   m_trk_unknown;
  std::vector<int>*   m_trk_combinatoric;
  std::vector<int>*   m_trk_fake; //0 fake, 1 track from primary interaction, 2 secondary track
  std::vector<int>*   m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_dxy;

  // all tracking particles
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_dxy;
  std::vector<float>* m_tp_d0;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0_prod;
  std::vector<float>* m_tp_z0_prod;
  std::vector<int>*   m_tp_pdgid;
  std::vector<int>*   m_tp_nmatch;
  std::vector<int>*   m_tp_nstub;
  std::vector<int>*   m_tp_nstublayers;
  std::vector<int>*   m_tp_eventid;
  std::vector<int>*   m_tp_charge;
  // *L1 track* properties if m_tp_nmatch > 0
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2; 
  std::vector<int>*   m_matchtrk_nstub;

  // ALL stubs
  std::vector<float>* m_allstub_x;
  std::vector<float>* m_allstub_y;
  std::vector<float>* m_allstub_z;

  std::vector<int>*   m_allstub_isBarrel; // stub is in barrel (1) or in disk (0)
  std::vector<int>*   m_allstub_layer;
  std::vector<int>*   m_allstub_isPSmodule;

  std::vector<float>* m_allstub_trigDisplace;
  std::vector<float>* m_allstub_trigOffset;
  std::vector<float>* m_allstub_trigPos;
  std::vector<float>* m_allstub_trigBend;

  // stub associated with tracking particle ?
  std::vector<int>*   m_allstub_matchTP_pdgid; // -999 if not matched
  std::vector<float>* m_allstub_matchTP_pt;    // -999 if not matched
  std::vector<float>* m_allstub_matchTP_eta;   // -999 if not matched
  std::vector<float>* m_allstub_matchTP_phi;   // -999 if not matched

  std::vector<int>*   m_allstub_genuine;
  std::vector<float>* m_genak4jet_phi;
  std::vector<float>* m_genak4jet_neufrac;
  std::vector<float>* m_genak4jet_chgfrac;
  std::vector<float>* m_genak4jet_metfrac;
  std::vector<float>* m_genak4jet_eta;
  std::vector<float>* m_genak4jet_pt;
  std::vector<float>* m_genak4jet_p;
  std::vector<float>* m_genak4chgjet_phi;
  std::vector<float>* m_genak4chgjet_eta;
  std::vector<float>* m_genak4chgjet_pt;
  std::vector<float>* m_genak4chgjet_p;
  std::vector<float>* m_jet_vz;
  std::vector<float>* m_jet_p;
  std::vector<float>* m_jet_phi;
  std::vector<float>* m_jet_eta;
  std::vector<int>* m_jet_ntracks;
  std::vector<float>* m_jet_tp_sumpt;
  std::vector<float>* m_jet_truetp_sumpt;
  std::vector<float>* m_jet_pt;

  std::vector<float>* m_recojet_vz;
  std::vector<float>* m_recojet_p;
  std::vector<float>* m_recojet_phi;
  std::vector<float>* m_recojet_eta;
  std::vector<int>* m_recojet_ntracks;
  std::vector<float>* m_recojet_tp_sumpt;
  std::vector<float>* m_recojet_truetp_sumpt;
  std::vector<float>* m_recojet_pt;
  std::vector<float>* m_recoPromptjet_vz;
  std::vector<float>* m_recoPromptjet_p;
  std::vector<float>* m_recoPromptjet_phi;
  std::vector<float>* m_recoPromptjet_eta;
  std::vector<int>* m_recoPromptjet_ntracks;
  std::vector<float>* m_recoPromptjet_pt;

  std::vector<float>* m_recoClusjet_vz;
  std::vector<float>* m_recoClusjet_p;
  std::vector<float>* m_recoClusjet_phi;
  std::vector<float>* m_recoClusjet_eta;
  std::vector<int>* m_recoClusjet_ntracks;
  std::vector<float>* m_recoClusjet_avgZNum;
  std::vector<float>* m_recoClusjet_truetp_sumpt;
  std::vector<float>* m_recoClusjet_pt;
  std::vector<float>* m_recoClusjet_seedpt;
  std::vector<float>* m_recoClusjet_seedeta;
  std::vector<float>* m_recoClusjet_seedphi;
  std::vector<int>* m_recoClusjet_seedzbin;
  std::vector<float>* m_truejet_vz;
  std::vector<float>* m_truejet_p;
  std::vector<float>* m_truejet_phi;
  std::vector<float>* m_truejet_eta;
  std::vector<int>* m_truejet_ntracks;
  std::vector<float>* m_truejet_tp_sumpt;
  std::vector<float>* m_truejet_truetp_sumpt;
  std::vector<float>* m_truejet_pt;

  std::vector<float>* m_jetPU_vz;
  std::vector<float>* m_jetPU_p;
  std::vector<float>* m_jetPU_phi;
  std::vector<float>* m_jetPU_eta;
  std::vector<int>* m_jetPU_ntracks;
  std::vector<float>* m_jetPU_pt;
  std::vector<float>* m_jetPU_tp_sumpt;
  std::vector<float>* m_jetPU_truetp_sumpt;
  std::vector<float>* m_jet_matchtrk_sumpt;
  std::vector<float>* m_jet_loosematchtrk_sumpt;
  std::vector<float>* m_jet_trk_sumpt;
        TH1F* htmp;
        TH1F* htmp_weight;
	TH1F*ZbinShiftIndicies;
	TH1F*ZbinIndicies;
	TH2F*Cluster3x3_shiftzbins[30];
	TH2F*Cluster3x3_zbins[30];
	
};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackNtupleMaker::L1TrackNtupleMaker(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{

  MyProcess        = iConfig.getParameter< int >("MyProcess");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  SaveAllTracks    = iConfig.getParameter< bool >("SaveAllTracks");
  SaveStubs        = iConfig.getParameter< bool >("SaveStubs");
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");
  TP_minNStub      = iConfig.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = iConfig.getParameter< double >("TP_minPt");
  TP_maxEta        = iConfig.getParameter< double >("TP_maxEta");
  TP_maxZ0         = iConfig.getParameter< double >("TP_maxZ0");
  L1TrackInputTag      = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub    = iConfig.getParameter< int >("L1Tk_minNStub");

  TrueVertexInputTag = iConfig.getParameter<edm::InputTag>("HepMCContent");
  L1StubInputTag      = iConfig.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingVertexInputTag = iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag");
  GenJetAK4=iConfig.getParameter<edm::InputTag>("GenJetAK4");
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);
  HEPMCVertexToken_=consumes< std::vector<reco::GenParticle> >(TrueVertexInputTag);
  ttStubToken_ = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);
   
  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  TrackingVertexToken_ = consumes< std::vector< TrackingVertex > >(TrackingVertexInputTag);
  GenJetCollectionToken_=consumes< std::vector<reco::GenJet > >(GenJetAK4);
 int nbins=600;
 float xmin = -30 ;
 float xmax = +30 ;

  htmp = new TH1F("htmp",";z (cm); Tracks",nbins,xmin,xmax);
  htmp_weight = new TH1F("htmp_weight",";z (cm); Tracks",nbins,xmin,xmax);
 for (unsigned int i=0; i<30 ; ++i)Cluster3x3_zbins[i]=new TH2F(TString::Format("Cluster3x3_zbins%d", i), "", 24, -2.4, 2.4, 28, -3.14, 3.14);
 for (unsigned int i=0; i<30 ; ++i)Cluster3x3_shiftzbins[i]=new TH2F(TString::Format("Cluster3x3_zbins%d", i), "", 24, -2.4, 2.4, 28, -3.14, 3.14);
  ZbinIndicies=new TH1F("ZbinIndicies", "", 30, -15, 15);
  ZbinShiftIndicies=new TH1F("ZbinShiftIndicies", "", 30, -14.5, 15.5);
}

/////////////
// DESTRUCTOR
L1TrackNtupleMaker::~L1TrackNtupleMaker()
{
}  

//////////
// END JOB
void L1TrackNtupleMaker::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackNtupleMaker::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackNtupleMaker::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackNtupleMaker::beginJob" << endl;

  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;


  // initilize
  m_trk_p    = new std::vector<float>;
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_genuine      = new std::vector<int>;
  m_trk_loose        = new std::vector<int>;
  m_trk_unknown      = new std::vector<int>;
  m_trk_combinatoric = new std::vector<int>;
  m_trk_fake = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt = new std::vector<float>;
  m_trk_matchtp_eta = new std::vector<float>;
  m_trk_matchtp_phi = new std::vector<float>;
  m_trk_matchtp_z0 = new std::vector<float>;
  m_trk_matchtp_dxy = new std::vector<float>;

  m_tp_pt     = new std::vector<float>;
  m_tp_eta    = new std::vector<float>;
  m_tp_phi    = new std::vector<float>;
  m_tp_dxy    = new std::vector<float>;
  m_tp_d0     = new std::vector<float>;
  m_tp_z0     = new std::vector<float>;
  m_tp_d0_prod = new std::vector<float>;
  m_tp_z0_prod = new std::vector<float>;
  m_tp_pdgid  = new std::vector<int>;
  m_tp_nmatch = new std::vector<int>;
  m_tp_nstub  = new std::vector<int>;
  m_tp_nstublayers  = new std::vector<int>;
  m_tp_eventid = new std::vector<int>;
  m_tp_charge = new std::vector<int>;

  m_matchtrk_pt    = new std::vector<float>;
  m_matchtrk_eta   = new std::vector<float>;
  m_matchtrk_phi   = new std::vector<float>;
  m_matchtrk_z0    = new std::vector<float>;
  m_matchtrk_d0    = new std::vector<float>;
  m_matchtrk_chi2  = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  
  m_allstub_x = new std::vector<float>;
  m_allstub_y = new std::vector<float>;
  m_allstub_z = new std::vector<float>;

  m_allstub_isBarrel = new std::vector<int>;
  m_allstub_layer    = new std::vector<int>;
  m_allstub_isPSmodule = new std::vector<int>;
  m_allstub_trigDisplace = new std::vector<float>;
  m_allstub_trigOffset   = new std::vector<float>;
  m_allstub_trigPos      = new std::vector<float>;
  m_allstub_trigBend     = new std::vector<float>;

  m_allstub_matchTP_pdgid = new std::vector<int>;
  m_allstub_matchTP_pt    = new std::vector<float>;
  m_allstub_matchTP_eta   = new std::vector<float>;
  m_allstub_matchTP_phi   = new std::vector<float>;

  m_allstub_genuine = new std::vector<int>;

  m_pv_L1reco = new std::vector<float>;
  m_pv_L1 = new std::vector<float>;
  m_pv_HTz = new std::vector<float>;
  m_pv_MC = new std::vector<float>;
  m_pv_MC_vr = new std::vector<float>;
  m_MC_lep=new std::vector<int>;  
m_genak4jet_phi = new std::vector<float>;
  m_genak4jet_neufrac = new std::vector<float>;
  m_genak4jet_chgfrac = new std::vector<float>;
  m_genak4jet_metfrac = new std::vector<float>;
  m_genak4jet_eta = new std::vector<float>;
  m_genak4jet_pt = new std::vector<float>;
  m_genak4jet_p = new std::vector<float>;

  m_genak4chgjet_phi = new std::vector<float>;
  m_genak4chgjet_eta = new std::vector<float>;
  m_genak4chgjet_pt = new std::vector<float>;
  m_genak4chgjet_p = new std::vector<float>;


 m_jetPU_eta = new std::vector<float>;
  m_jetPU_vz = new std::vector<float>;
  m_jetPU_phi = new std::vector<float>;
  m_jetPU_p = new std::vector<float>;
  m_jetPU_pt = new std::vector<float>;
  m_jetPU_ntracks = new std::vector<int>;
  m_jetPU_tp_sumpt = new std::vector<float>;
  m_jetPU_truetp_sumpt = new std::vector<float>;

  m_recojet_eta = new std::vector<float>;
  m_recojet_vz = new std::vector<float>;
  m_recojet_phi = new std::vector<float>;
  m_recojet_p = new std::vector<float>;
  m_recojet_pt = new std::vector<float>;
  m_recojet_ntracks = new std::vector<int>;
  m_recojet_truetp_sumpt = new std::vector<float>;
  m_recoClusjet_eta = new std::vector<float>;
  m_recoClusjet_vz = new std::vector<float>;
  m_recoClusjet_phi = new std::vector<float>;
  m_recoClusjet_p = new std::vector<float>;
  m_recoClusjet_pt = new std::vector<float>;
  m_recoClusjet_seedpt = new std::vector<float>;
  m_recoClusjet_seedzbin = new std::vector<int>;
  m_recoClusjet_seedphi = new std::vector<float>;
  m_recoClusjet_seedeta = new std::vector<float>;
  m_recoClusjet_ntracks = new std::vector<int>;
  m_recoClusjet_truetp_sumpt = new std::vector<float>;
  m_recoClusjet_avgZNum= new std::vector<float>;

  m_recoPromptjet_vz=new std::vector<float>;;
   m_recoPromptjet_p=new std::vector<float>;
  m_recoPromptjet_phi=new std::vector<float>;
   m_recoPromptjet_eta=new std::vector<float>;
  m_recoPromptjet_ntracks=new std::vector<int>;
  m_recoPromptjet_pt=new std::vector<float>;

  m_jet_eta = new std::vector<float>;
  m_jet_vz = new std::vector<float>;
  m_jet_phi = new std::vector<float>;
  m_jet_p = new std::vector<float>;
  m_jet_pt = new std::vector<float>;
  m_jet_ntracks = new std::vector<int>;
  m_jet_tp_sumpt = new std::vector<float>;
  m_jet_truetp_sumpt = new std::vector<float>;
  m_jet_matchtrk_sumpt = new std::vector<float>;
  m_jet_loosematchtrk_sumpt = new std::vector<float>;
  m_jet_trk_sumpt = new std::vector<float>;
  m_truejet_eta = new std::vector<float>;
  m_truejet_vz = new std::vector<float>;
  m_truejet_phi = new std::vector<float>;
  m_truejet_p = new std::vector<float>;
  m_truejet_pt = new std::vector<float>;
  m_truejet_ntracks = new std::vector<int>;
  m_truejet_tp_sumpt = new std::vector<float>;
  m_truejet_truetp_sumpt = new std::vector<float>;
  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");
  
  if (SaveAllTracks) {
    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_d0",    &m_trk_d0);
    eventTree->Branch("trk_z0",    &m_trk_z0);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);

    eventTree->Branch("trk_genuine",      &m_trk_genuine);
    eventTree->Branch("trk_loose",        &m_trk_loose);
    eventTree->Branch("trk_unknown",      &m_trk_unknown);
    eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
    eventTree->Branch("trk_fake", &m_trk_fake);
    eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
    eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
    eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
    eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
    eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
    eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);
  }

  eventTree->Branch("tp_pt",     &m_tp_pt);
  eventTree->Branch("tp_eta",    &m_tp_eta);
  eventTree->Branch("tp_phi",    &m_tp_phi);
  eventTree->Branch("tp_dxy",    &m_tp_dxy);
  eventTree->Branch("tp_d0",     &m_tp_d0);
  eventTree->Branch("tp_z0",     &m_tp_z0);
  eventTree->Branch("tp_d0_prod",&m_tp_d0_prod);
  eventTree->Branch("tp_z0_prod",&m_tp_z0_prod);
  eventTree->Branch("tp_pdgid",  &m_tp_pdgid);
  eventTree->Branch("tp_nmatch", &m_tp_nmatch);
  eventTree->Branch("tp_nstub", &m_tp_nstub);
  eventTree->Branch("tp_nstublayers", &m_tp_nstublayers);
  eventTree->Branch("tp_eventid",&m_tp_eventid);
  eventTree->Branch("tp_charge",&m_tp_charge);

  eventTree->Branch("matchtrk_pt",      &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta",     &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi",     &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0",      &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0",      &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2",    &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_nstub",   &m_matchtrk_nstub);
    eventTree->Branch("pv_L1reco", &m_pv_L1reco);
    eventTree->Branch("pv_L1", &m_pv_L1);
    eventTree->Branch("MC_lep", &m_MC_lep);
    eventTree->Branch("pv_HTz", &m_pv_HTz);
    eventTree->Branch("pv_MC", &m_pv_MC);
    eventTree->Branch("pv_MC_vr", &m_pv_MC_vr);
    eventTree->Branch("truejet_eta", &m_truejet_eta);
    eventTree->Branch("truejet_vz", &m_truejet_vz);
    eventTree->Branch("truejet_p", &m_truejet_p);
    eventTree->Branch("truejet_pt", &m_truejet_pt);
    eventTree->Branch("truejet_phi", &m_truejet_phi);
    eventTree->Branch("truejet_ntracks", &m_truejet_ntracks);
    eventTree->Branch("truejet_tp_sumpt", &m_truejet_tp_sumpt);
    eventTree->Branch("truejet_truetp_sumpt", &m_truejet_truetp_sumpt);
    eventTree->Branch("recojet_eta", &m_recojet_eta);
    eventTree->Branch("recojet_vz", &m_recojet_vz);
    eventTree->Branch("recojet_p", &m_recojet_p);
    eventTree->Branch("recojet_pt", &m_recojet_pt);
    eventTree->Branch("recojet_phi", &m_recojet_phi);
    eventTree->Branch("recojet_ntracks", &m_recojet_ntracks);
    eventTree->Branch("recojet_truetp_sumpt", m_recojet_truetp_sumpt);
    eventTree->Branch("recoClusjet_eta", &m_recoPromptjet_eta);
    eventTree->Branch("recoClusjet_vz", &m_recoPromptjet_vz);
    eventTree->Branch("recoClusjet_p", &m_recoPromptjet_p);
    eventTree->Branch("recoClusjet_pt", &m_recoPromptjet_pt);
    eventTree->Branch("recoClusjet_seedpt", &m_recoClusjet_seedpt);
    eventTree->Branch("recoClusjet_seedzbin", &m_recoClusjet_seedzbin);
    eventTree->Branch("recoClusjet_seedphi", &m_recoClusjet_seedphi);
    eventTree->Branch("recoClusjet_seedeta", &m_recoClusjet_seedeta);
    eventTree->Branch("recoClusjet_phi", &m_recoPromptjet_phi);
    eventTree->Branch("recoClusjet_ntracks", &m_recoPromptjet_ntracks);
    eventTree->Branch("jet_eta", &m_jet_eta);
    eventTree->Branch("jet_vz", &m_jet_vz);
    eventTree->Branch("jet_p", &m_jet_p);
    eventTree->Branch("jet_pt", &m_jet_pt);
    eventTree->Branch("jet_phi", &m_jet_phi);
    eventTree->Branch("jet_ntracks", &m_jet_ntracks);
    eventTree->Branch("jet_tp_sumpt", &m_jet_tp_sumpt);
    eventTree->Branch("jet_truetp_sumpt", &m_jet_truetp_sumpt);
    eventTree->Branch("jetPU_eta", &m_jetPU_eta);
    eventTree->Branch("jetPU_vz", &m_jetPU_vz);
    eventTree->Branch("jetPU_p", &m_jetPU_p);
    eventTree->Branch("jetPU_pt", &m_jetPU_pt);
    eventTree->Branch("jetPU_phi", &m_jetPU_phi);
    eventTree->Branch("jetPU_ntracks", &m_jetPU_ntracks);
    eventTree->Branch("jetPU_tp_sumpt", &m_jetPU_tp_sumpt);
    eventTree->Branch("jetPU_truetp_sumpt", &m_jetPU_truetp_sumpt);
    eventTree->Branch("genjetak4_neufrac", &m_genak4jet_neufrac);
    eventTree->Branch("genjetak4_chgfrac", &m_genak4jet_chgfrac);
    eventTree->Branch("genjetak4_metfrac", &m_genak4jet_metfrac);
    eventTree->Branch("genjetak4_eta", &m_genak4jet_eta);
    eventTree->Branch("genjetak4_phi", &m_genak4jet_phi);
    eventTree->Branch("genjetak4_p", &m_genak4jet_p);
    eventTree->Branch("genjetak4_pt", &m_genak4jet_pt);
    eventTree->Branch("genjetchgak4_eta", &m_genak4chgjet_eta);
    eventTree->Branch("genjetchgak4_phi", &m_genak4chgjet_phi);
    eventTree->Branch("genjetchgak4_p", &m_genak4chgjet_p);
    eventTree->Branch("genjetchgak4_pt", &m_genak4chgjet_pt);
  if (SaveStubs) {
    eventTree->Branch("allstub_x", &m_allstub_x);
    eventTree->Branch("allstub_y", &m_allstub_y);
    eventTree->Branch("allstub_z", &m_allstub_z);

    eventTree->Branch("allstub_isBarrel", &m_allstub_isBarrel);
    eventTree->Branch("allstub_layer", &m_allstub_layer);
    eventTree->Branch("allstub_isPSmodule", &m_allstub_isPSmodule);
    
    eventTree->Branch("allstub_trigDisplace", &m_allstub_trigDisplace);
    eventTree->Branch("allstub_trigOffset", &m_allstub_trigOffset);
    eventTree->Branch("allstub_trigPos", &m_allstub_trigPos);
    eventTree->Branch("allstub_trigBend", &m_allstub_trigBend);

    eventTree->Branch("allstub_matchTP_pdgid", &m_allstub_matchTP_pdgid);
    eventTree->Branch("allstub_matchTP_pt", &m_allstub_matchTP_pt);
    eventTree->Branch("allstub_matchTP_eta", &m_allstub_matchTP_eta);
    eventTree->Branch("allstub_matchTP_phi", &m_allstub_matchTP_phi);

    eventTree->Branch("allstub_genuine", &m_allstub_genuine);
    eventTree->Branch("pv_L1", &m_pv_L1);
    eventTree->Branch("pv_MC", &m_pv_MC);
  }

}


//////////
// ANALYZE
void L1TrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    return;
  }

  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    return;
  }
 htmp -> Reset();
 htmp_weight -> Reset();
  // clear variables
  if (SaveAllTracks) {
    m_trk_p->clear();
    m_trk_pt->clear();
    m_trk_eta->clear();
    m_trk_phi->clear();
    m_trk_d0->clear();
    m_trk_z0->clear();
    m_trk_chi2->clear();
    m_trk_nstub->clear();
    m_trk_genuine->clear();
    m_trk_loose->clear();
    m_trk_unknown->clear();
    m_trk_combinatoric->clear();
    m_trk_fake->clear();
    m_trk_matchtp_pdgid->clear();
    m_trk_matchtp_pt->clear();
    m_trk_matchtp_eta->clear();
    m_trk_matchtp_phi->clear();
    m_trk_matchtp_z0->clear();
    m_trk_matchtp_dxy->clear();
  }
  
  m_tp_pt->clear();
  m_tp_eta->clear();
  m_tp_phi->clear();
  m_tp_dxy->clear();
  m_tp_d0->clear();
  m_tp_z0->clear();
  m_tp_d0_prod->clear();
  m_tp_z0_prod->clear();
  m_tp_pdgid->clear();
  m_tp_nmatch->clear();
  m_tp_nstub->clear();
  m_tp_nstublayers->clear();
  m_tp_eventid->clear();
  m_tp_charge->clear();

  m_matchtrk_pt->clear();
  m_matchtrk_eta->clear();
  m_matchtrk_phi->clear();
  m_matchtrk_z0->clear();
  m_matchtrk_d0->clear();
  m_matchtrk_chi2->clear();
  m_matchtrk_nstub->clear();

  if (SaveStubs) {
    m_allstub_x->clear();
    m_allstub_y->clear();
    m_allstub_z->clear();

    m_allstub_isBarrel->clear();
    m_allstub_layer->clear();
    m_allstub_isPSmodule->clear();

    m_allstub_trigDisplace->clear();
    m_allstub_trigOffset->clear();
    m_allstub_trigPos->clear();
    m_allstub_trigBend->clear();

    m_allstub_matchTP_pdgid->clear();
    m_allstub_matchTP_pt->clear();
    m_allstub_matchTP_eta->clear();
    m_allstub_matchTP_phi->clear();

    m_allstub_genuine->clear();
  }
  m_genak4chgjet_phi->clear();
  m_genak4chgjet_eta->clear() ;
  m_genak4chgjet_pt->clear() ;
  m_genak4chgjet_p->clear() ;
  m_genak4jet_phi->clear();
  m_genak4jet_neufrac->clear() ;
  m_genak4jet_chgfrac->clear() ;
  m_genak4jet_metfrac->clear() ;
  m_genak4jet_eta->clear() ;
  m_genak4jet_pt->clear() ;
  m_genak4jet_p->clear() ;
  m_recojet_eta->clear();
  m_recojet_pt->clear();
  m_recojet_vz->clear();
  m_recojet_phi->clear();
  m_recojet_p->clear();
  m_recojet_ntracks->clear();
  m_recojet_truetp_sumpt->clear();

  m_recoClusjet_eta->clear();
  m_recoClusjet_seedeta->clear();
  m_recoClusjet_seedpt->clear();
  m_recoClusjet_seedzbin->clear();
  m_recoClusjet_seedphi->clear();
  m_recoClusjet_pt->clear();
  m_recoClusjet_vz->clear();
  m_recoClusjet_phi->clear();
  m_recoClusjet_p->clear();
  m_recoClusjet_ntracks->clear();
  m_recoClusjet_truetp_sumpt->clear();
m_recoClusjet_avgZNum->clear();
  m_recoPromptjet_vz->clear();
  m_recoPromptjet_p->clear();
  m_recoPromptjet_phi->clear();
  m_recoPromptjet_eta->clear();
  m_recoPromptjet_ntracks->clear();
  m_recoPromptjet_pt->clear();

m_jet_eta->clear();
  m_jet_pt->clear();
  m_jet_vz->clear();
  m_jet_phi->clear();
  m_jet_p->clear();
  m_jet_ntracks->clear();
  m_jet_tp_sumpt->clear();
  m_jet_truetp_sumpt->clear();
  m_jet_matchtrk_sumpt->clear();
  m_jet_loosematchtrk_sumpt->clear();
  m_jet_trk_sumpt->clear();
  m_truejet_eta->clear();
  m_truejet_pt->clear();
  m_truejet_vz->clear();
  m_truejet_phi->clear();
  m_truejet_p->clear();
  m_truejet_ntracks->clear();
  m_truejet_tp_sumpt->clear();
  m_truejet_truetp_sumpt->clear();

  m_jetPU_eta->clear();
  m_jetPU_pt->clear();
  m_jetPU_vz->clear();
  m_jetPU_phi->clear();
  m_jetPU_p->clear();
  m_jetPU_ntracks->clear();
  m_jetPU_tp_sumpt->clear();
  m_jetPU_truetp_sumpt->clear();
  m_pv_L1reco->clear();
  m_pv_L1->clear();
  m_pv_HTz->clear();
  m_pv_MC->clear();
  m_pv_MC_vr->clear();
  m_MC_lep->clear();

  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);
    edm::Handle<std::vector<reco::GenJet> >GenJetsAK4Handle;
  iEvent.getByToken(GenJetCollectionToken_,GenJetsAK4Handle); 
  // L1 stubs
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  if (SaveStubs) iEvent.getByToken(ttStubToken_, TTStubHandle);


   edm::Handle< std::vector< reco::GenParticle> > GenParticleHandle;
   iEvent.getByToken(HEPMCVertexToken_,GenParticleHandle);
  // MC truth association maps
  edm::Handle< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  iEvent.getByToken(TrackingVertexToken_, TrackingVertexHandle);


  // -----------------------------------------------------------------------------------------------
  // more for TTStubs
  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();



     float zvtx_gen = -999;
     float rvtx_gen=-999;

     // if (GenParticleHandle.isValid() ) {
        vector<reco::GenParticle>::const_iterator genpartIter ;
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
           int status = genpartIter -> status() ;
	   //if (status< 21 || status>29) continue;
	   //if (status< 51 || status>59) continue;
	   if (status!=1) continue;
	   if ( genpartIter -> numberOfMothers() == 0) continue;   // the incoming hadrons
	   //if (abs(genpartIter ->pdgId())!=5)continue;
	   //float part_zvertex = genpartIter ->daughter(0)-> vz() ;
	   float part_zvertex = genpartIter -> vz() ;
	   zvtx_gen = part_zvertex ;
	   //rvtx_gen=sqrt(genpartIter ->daughter(0)->daughter(0) ->vx() * genpartIter ->daughter(0)->daughter(0) ->vx() + genpartIter ->daughter(0) ->vy()*genpartIter ->daughter(0) ->vy());
	//   rvtx_gen=sqrt(genpartIter ->daughter(0)->daughter(0) ->vx() * genpartIter ->daughter(0)->daughter(0) ->vx() + genpartIter ->daughter(0)->daughter(0) ->vy()*genpartIter ->daughter(0)->daughter(0) ->vy());
	   //std::cout<<"PdgId "<<genpartIter ->daughter(0)->daughter(0)->pdgId()<<std::endl;
	   //std::cout<<"R vertex "<<rvtx_gen<<std::endl;
	   //std::cout<<"MC z "<<part_zvertex<<std::endl;
	   break;	// 
	}
    // }
m_pv_MC->push_back(zvtx_gen);

    int leptonicCount=0;
/*
    bool leptonic=false;
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
		if (abs(genpartIter ->pdgId())!=24)continue;
		//std::cout<<"W-boson mother "<<genpartIter ->mother(0)->pdgId()<<std::endl;
		if(abs(genpartIter ->mother(0)->pdgId())!=6 && abs(genpartIter ->mother(0)->pdgId())!=24 )continue;
		if( abs(genpartIter ->daughter(0)->pdgId())==11  ||  abs(genpartIter ->daughter(0)->pdgId())==13 ||  abs(genpartIter ->daughter(0)->pdgId())==15){leptonic=true;break;}
	}
m_pv_MC_vr->push_back(rvtx_gen);
if(leptonic)m_MC_lep->push_back(1);
else m_MC_lep->push_back(0);
*/
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
		if (abs(genpartIter ->pdgId())!=24)continue;
		//std::cout<<"W-boson mother "<<genpartIter ->mother(0)->pdgId()<<std::endl;
		if(abs(genpartIter ->mother(0)->pdgId())!=6 && abs(genpartIter ->mother(0)->pdgId())!=24 )continue;
		if( abs(genpartIter ->daughter(0)->pdgId())==11  ||  abs(genpartIter ->daughter(0)->pdgId())==13 ||  abs(genpartIter ->daughter(0)->pdgId())==15){++leptonicCount;}
	}
m_MC_lep->push_back(leptonicCount);
std::cout<<"Leptonic W count "<<leptonicCount<<std::endl;
  // ----------------------------------------------------------------------------------------------
  // loop over L1 stubs
  // ----------------------------------------------------------------------------------------------

  if (SaveStubs) {

    for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
    
      DetId detid = (*gd)->geographicalId();
      if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue;
      DetId stackDetid = tTopo->stack(detid); // Stub module detid

      if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;

      // Get the DetSets of the Clusters
      edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );
    
      // loop over stubs
      for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) {
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > >
	  tempStubPtr = edmNew::makeRefTo( TTStubHandle, stubIter );
	
	int isBarrel = 0;
	int layer=-999999;
	if ( detid.subdetId()==StripSubdetector::TOB ) {
	  isBarrel = 1;
	  layer  = static_cast<int>(tTopo->layer(detid));
	}
	else if ( detid.subdetId()==StripSubdetector::TID ) {
	  isBarrel = 0;
	  layer  = static_cast<int>(tTopo->layer(detid));
	}
	else {
	  cout << "WARNING -- neither TOB or TID stub, shouldn't happen..." << endl;
	  layer = -1;
	}

	int isPSmodule=0;
	if (topol->nrows() == 960) isPSmodule=1;

	MeasurementPoint coords = tempStubPtr->getClusterRef(0)->findAverageLocalCoordinatesCentered();      
	LocalPoint clustlp = topol->localPosition(coords);
	GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);

	double tmp_stub_x=posStub.x();
	double tmp_stub_y=posStub.y();
	double tmp_stub_z=posStub.z();

	float trigDisplace = tempStubPtr->getTriggerDisplacement();
	float trigOffset = tempStubPtr->getTriggerOffset();
	float trigPos = tempStubPtr->getTriggerPosition();
	float trigBend = tempStubPtr->getTriggerBend();

	m_allstub_x->push_back(tmp_stub_x);
	m_allstub_y->push_back(tmp_stub_y);
	m_allstub_z->push_back(tmp_stub_z);

	m_allstub_isBarrel->push_back(isBarrel);
	m_allstub_layer->push_back(layer);
	m_allstub_isPSmodule->push_back(isPSmodule);

	m_allstub_trigDisplace->push_back(trigDisplace);
	m_allstub_trigOffset->push_back(trigOffset);
	m_allstub_trigPos->push_back(trigPos);
	m_allstub_trigBend->push_back(trigBend);

	// matched to tracking particle? 
	edm::Ptr< TrackingParticle > my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);

	int myTP_pdgid = -999;
	float myTP_pt  = -999;
	float myTP_eta = -999;
	float myTP_phi = -999;

	if (my_tp.isNull() == false) {
	  //int tmp_eventid = my_tp->eventId().event();

	 // if (tmp_eventid > 0) continue; // this means stub from pileup track

	  
	  myTP_pdgid = my_tp->pdgId();
	  myTP_pt = my_tp->p4().pt();
	  myTP_eta = my_tp->p4().eta();
	  myTP_phi = my_tp->p4().phi();
	}

	m_allstub_matchTP_pdgid->push_back(myTP_pdgid);
	m_allstub_matchTP_pt->push_back(myTP_pt);
	m_allstub_matchTP_eta->push_back(myTP_eta);
	m_allstub_matchTP_phi->push_back(myTP_phi);
	
	int tmp_stub_genuine = 0;
	if (MCTruthTTStubHandle->isGenuine(tempStubPtr)) tmp_stub_genuine = 1;

	m_allstub_genuine->push_back(tmp_stub_genuine);

      }
      
    }

  }


  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------
    m_pv_L1->push_back(FillPrimaryVtx(TrackingParticleHandle,  MCTruthTTStubHandle,tTopo,  theTrackerGeom));

  if (SaveAllTracks) {
    
    if (DebugMode) {
      cout << endl << "Loop over L1 tracks!" << endl;
      cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
    }
    
    int this_l1track = 0;
    std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
	RecoJetInputs_.clear();

    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
      
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
      this_l1track++;
      
      float tmp_trk_p   = iterL1Track->getMomentum(L1Tk_nPar).mag();
      float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
      float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
      float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
      float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm

      float tmp_trk_d0 = -999;
      if (L1Tk_nPar == 5) {
	float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
	float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();	
	tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
      }

      float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
      int tmp_trk_nstub  = (int) iterL1Track->getStubRefs().size();


      // ----------------------------------------------------------------------------------------------
      // loop over stubs on tracks
      if (DebugMode && SaveStubs) {

	// loop over stubs
	std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
	for (int is=0; is<tmp_trk_nstub; is++) {
	  
	  //detID of stub
	  DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	  MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();      
	  const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
	  Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );
	  
	  double x=posStub.x();
	  double y=posStub.y();
	  double z=posStub.z();
	  
	  int layer=-999999;
	  if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    cout << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
	  }
	  else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    cout << "   stub in disk " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
	  }	  
	  
	}//end loop over stubs
      }
      // ----------------------------------------------------------------------------------------------


      int tmp_trk_genuine = 0;
      int tmp_trk_loose = 0;
      int tmp_trk_unknown = 0;
      int tmp_trk_combinatoric = 0;
      if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
      if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
      if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;
      
      if (DebugMode) {
	cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi 
	     << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
	if (tmp_trk_genuine) cout << " (is genuine)" << endl; 
	if (tmp_trk_unknown) cout << " (is unknown)" << endl; 
	if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl; 
      }
      m_trk_p ->push_back(tmp_trk_p); 
      m_trk_pt ->push_back(tmp_trk_pt);
      m_trk_eta->push_back(tmp_trk_eta);
      m_trk_phi->push_back(tmp_trk_phi);
      m_trk_z0 ->push_back(tmp_trk_z0);
      if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
      else m_trk_d0->push_back(999.);
      m_trk_chi2 ->push_back(tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar));
      m_trk_nstub->push_back(tmp_trk_nstub);
      m_trk_genuine->push_back(tmp_trk_genuine);
      m_trk_loose->push_back(tmp_trk_loose);
      m_trk_unknown->push_back(tmp_trk_unknown);
      m_trk_combinatoric->push_back(tmp_trk_combinatoric);
      

      // ----------------------------------------------------------------------------------------------
      // for studying the fake rate
      // ----------------------------------------------------------------------------------------------

      edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

      int myFake = 0;

      int myTP_pdgid = -999;
      float myTP_pt = -999;
      float myTP_eta = -999;
      float myTP_phi = -999;
      float myTP_z0 = -999;
      float myTP_dxy = -999;


      if (my_tp.isNull()) myFake = 0;
      else {
	int tmp_eventid = my_tp->eventId().event();

	if (tmp_eventid > 0) myFake = 2;
	else myFake = 1;

	myTP_pdgid = my_tp->pdgId();
	myTP_pt = my_tp->p4().pt();
	myTP_eta = my_tp->p4().eta();
	myTP_phi = my_tp->p4().phi();
	myTP_z0 = my_tp->vertex().z();
	
	float myTP_x0 = my_tp->vertex().x();
	float myTP_y0 = my_tp->vertex().y();
	myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

	if (DebugMode) {
	  cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
	       << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
	       << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
	}
      }

      m_trk_fake->push_back(myFake);

      m_trk_matchtp_pdgid->push_back(myTP_pdgid);
      m_trk_matchtp_pt->push_back(myTP_pt);
      m_trk_matchtp_eta->push_back(myTP_eta);
      m_trk_matchtp_phi->push_back(myTP_phi);
      m_trk_matchtp_z0->push_back(myTP_z0);
      m_trk_matchtp_dxy->push_back(myTP_dxy);

    }//end track loop

  }//end if SaveAllTracks
  m_pv_L1reco->push_back(FillRecoPrimaryVtx());
  RecoJetInputs_.clear();
   for(unsigned int t=0; t<m_trk_z0->size(); ++t){
    if(m_trk_pt->at(t)>2 && fabs(m_trk_z0->at(t)-m_pv_L1->at(0))<0.5 && m_trk_chi2->at(t)<5	&& m_trk_nstub->at(t)>4
	){
    //if(m_trk_pt->at(t)>2 && m_trk_chi2->at(t)<5){
      TLorentzVector temp;
      //if(m_trk_pt->at(t)<200)
	temp.SetPtEtaPhiE(m_trk_pt->at(t), m_trk_eta->at(t), m_trk_phi->at(t), m_trk_p->at(t));
      //else temp.SetPtEtaPhiE(200, m_trk_eta->at(t), m_trk_phi->at(t), cosh(m_trk_eta->at(t))*200);
      fastjet::PseudoJet psuedoJet(temp.Px(), temp.Py(), temp.Pz(), temp.P()); 
      RecoJetInputs_.push_back( psuedoJet);
      RecoJetInputs_.back().set_user_index(t);
    }
}
//std::cout<<" Total Tracks "<<m_trk_z0->size()<<"PV L1 "<<m_pv_L1->at(0)<<" Reco Tracks "<<m_pv_L1reco->at(0)<<std::endl;
TrueJetInputs_.clear();

if(GenJetsAK4Handle.isValid()){
        for (unsigned iGenJet = 0; iGenJet < GenJetsAK4Handle->size(); ++iGenJet) {

         const reco::GenJet& genJet = (*GenJetsAK4Handle) [iGenJet];
   if(fabs(genJet.eta())<2.4 ){
          m_genak4jet_p->push_back(genJet.energy());
          m_genak4jet_pt->push_back(genJet.pt());
	  m_genak4jet_metfrac->push_back(genJet.invisibleEnergy()/genJet.energy());
          m_genak4jet_eta->push_back(genJet.eta());
          m_genak4jet_phi->push_back(genJet.phi());
	float NeuEnergy=0;
	float ChgEnergy=0;

	for(unsigned g=0; g<genJet.getGenConstituents().size(); ++g){
	if(genJet.getGenConstituent(g)->charge()!=0 && genJet.getGenConstituent(g)->pt()>2){
	fastjet::PseudoJet psuedoJet(genJet.getGenConstituent(g)->px(), genJet.getGenConstituent(g)->py(), genJet.getGenConstituent(g)->pz(), genJet.getGenConstituent(g)->p());
	TrueJetInputs_.push_back(psuedoJet);
	//TrueJetInputs_.back().set_user_index(this_tp);
	///	std::cout<<"Gen Jet Const vertex "<<"Pdg Id "<<genJet.getGenConstituent(g)->pdgId()<<" "<<genJet.getGenConstituent(g)->vx()<<", "<<genJet.getGenConstituent(g)->vy()<<std::endl;
		ChgEnergy=ChgEnergy+genJet.getGenConstituent(g)->energy();
	}
	else NeuEnergy=NeuEnergy+genJet.getGenConstituent(g)->energy();
        }
	m_genak4jet_chgfrac->push_back(ChgEnergy);
        m_genak4jet_neufrac->push_back(NeuEnergy);
	//if(m_MC_lep->at(0)<1  && ChgEnergy/genJet.energy()>0.98)std::cout<<"getGenConstituent Print "<<genJet.print()<<std::endl;
	//std::cout<<"Charge Energy Fraction "<<ChgEnergy/genJet.energy()<<std::endl;
//	if(ChgEnergy/genJet.energy()>0.6 && genJet.pt()>30){
//		std::cout<<"Gen Jet "<<genJet.eta()<<", "<<genJet.phi()<<std::endl;
//	}
    }
  }
}

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------


  if (DebugMode) cout << endl << "Loop over tracking particles!" << endl;
JetInputs_.clear();
JetInputsPU1_.clear();
JetInputsPU2_.clear();
JetInputsPU3_.clear();
JetInputsPU4_.clear();
JetInputsPU5_.clear();
JetInputsPU6_.clear();
JetInputsPU7_.clear();
JetInputsPU8_.clear();
JetInputsPU9_.clear();
JetInputsPU10_.clear();
JetInputsPU11_.clear();
JetInputsPU12_.clear();
JetInputsPU13_.clear();
JetInputsPU14_.clear();
JetInputsPU15_.clear();
JetInputsPU16_.clear();
JetInputsPU17_.clear();
JetInputsPU18_.clear();
JetInputsPU19_.clear();
JetInputsPU20_.clear();
  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;

  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
 
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int) theStubRefs.size();

    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (unsigned int is=0; is<theStubRefs.size(); is++) {
	  
	  //detID of stub
	  DetId detIdStub = theTrackerGeom->idToDet( (theStubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	  
         int layer = -1;
	  if ( detIdStub.subdetId()==StripSubdetector::TOB ) {

	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer-1;
	   //if(fabs(tp_ptr->eta())>1.35 && fabs(tp_ptr->eta())<1.45)std::cout<<"I am in the TOB layer"<<layer <<std::endl;
	  }
	
	  else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer+5;
	   //if(fabs(tp_ptr->eta())>1.35 && fabs(tp_ptr->eta())<1.45)std::cout<<"I am in the TID layer"<<layer <<std::endl;
	  }	  
	  else if ( detIdStub.subdetId()==StripSubdetector::TIB ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	   // std::cout<<"I am in the TIB layer"<<layer<<std::endl;

	  }	  
	  else if ( detIdStub.subdetId()==StripSubdetector::TEC ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	   // std::cout<<"I am in the TEC layer"<<layer<<std::endl;
	  }
	else{
		std::cout<<"I am a stub I don't know where i am; Layer "<<static_cast<int>(tTopo->layer(detIdStub))<<std::endl;
	}	  
         if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)hasStubInLayer[layer] = 1;
         else hasStubInLayer[layer] = 2;
}

    //                                                          //treat genuine stubs separately (==2 is genuine, ==1 is not)
    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum=0; isum<11; isum++) {
        if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
        if ( hasStubInLayer[isum] == 2) nStubLayerTP_g += 1;
   }
        int tmp_eventid = iterTP->eventId().event();
    float  tmp_tp_vx=tp_ptr->vx();
    float  tmp_tp_vy=tp_ptr->vy();
    float  tmp_tp_vz=tp_ptr->vz();
    float tmp_tp_eta=tp_ptr->eta();
    float tmp_tp_phi=tp_ptr->phi();
    float tmp_tp_pt=tp_ptr->pt();
    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;
	
    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;
    
    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;
	
    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));
	
    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);
    //if(tmp_eventid>0)continue;
/*
    if(fabs(iterTP->eta())<2.4 && tmp_eventid==0 && tp_ptr->pt() > TP_minPt && fabs(tmp_tp_z0-m_pv_L1->at(0))<0.5 ){
	fastjet::PseudoJet psuedoJet(tp_ptr->px(), tp_ptr->py(), tp_ptr->pz(), tp_ptr->p());
	TrueJetInputs_.push_back(psuedoJet);
	TrueJetInputs_.back().set_user_index(this_tp);
    }
*/
    if(fabs(iterTP->eta())<2.4 && tp_ptr->pt() > TP_minPt && fabs(tmp_tp_z0-m_pv_L1->at(0))<0.5 && nStubLayerTP>=4){ 
    fastjet::PseudoJet psuedoJet(tp_ptr->px(), tp_ptr->py(), tp_ptr->pz(), tp_ptr->p()); 
    JetInputs_.push_back(psuedoJet);	
    //JetInputs_.back().set_user_index(tp_index);
    JetInputs_.back().set_user_index(this_tp);
    }
    this_tp++;

    //if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)



  /* 
    float tmp_tp_pt  = iterTP->pt();
    float tmp_tp_eta = iterTP->eta();
    float tmp_tp_phi = iterTP->phi(); 
    float tmp_tp_vz  = iterTP->vz();
    float tmp_tp_vx  = iterTP->vx();
    float tmp_tp_vy  = iterTP->vy();
*/
    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = -tmp_tp_vx*sin(tmp_tp_phi) + tmp_tp_vy*cos(tmp_tp_phi);

    if (MyProcess==13 && abs(tmp_tp_pdgid) != 13) continue;
    if (MyProcess==11 && abs(tmp_tp_pdgid) != 11) continue;
    if ((MyProcess==6 || MyProcess==15 || MyProcess==211) && abs(tmp_tp_pdgid) != 211) continue;


    if (tmp_tp_pt < TP_minPt) continue;
    if (fabs(tmp_tp_eta) > TP_maxEta) continue;


    // ----------------------------------------------------------------------------------------------
    // get d0/z0 propagated back to the IP
    // ----------------------------------------------------------------------------------------------
    
    if (fabs(tmp_tp_z0) > TP_maxZ0) continue;


    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
    float tmp_tp_dxy = dxy;
    if (MyProcess==6 && (dxy > 1.0)) continue;

    if (DebugMode) cout << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta << " phi: " << tmp_tp_phi 
			<< " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0 
			<< " z_prod: " << tmp_tp_z0_prod << " d_prod: " << tmp_tp_d0_prod 
			<< " pdgid: " << tmp_tp_pdgid << " eventID: " << iterTP->eventId().event()
			<< " ttclusters " << MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() 
			<< " ttstubs " << MCTruthTTStubHandle->findTTStubRefs(tp_ptr).size()
			<< " tttracks " << MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr).size() << endl;


    // ----------------------------------------------------------------------------------------------
    // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)
    
    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() < 1) {
      if (DebugMode) cout << "No matching TTClusters for TP, continuing..." << endl;
      continue;
    }


    //std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    //int nStubTP = (int) theStubRefs.size(); 

    if (TP_minNStub > 0) {
      if (DebugMode) cout << "Only consider TPs with >= " << TP_minNStub << " stubs" << endl;
      if (nStubTP < TP_minNStub) {
	if (DebugMode) cout << "TP fails minimum nbr stubs requirement! Continuing..." << endl;
	continue;
      }
    }


    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle

    std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);
    
    int nMatch = 0;
    int i_track = -1;
    float i_chi2dof = 99999;

    if (matchedTracks.size() > 0) { 
    
      if (DebugMode && (matchedTracks.size()>1)) cout << "TrackingParticle has more than one matched L1 track!" << endl;

      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters 
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack

      for (int it=0; it<(int)matchedTracks.size(); it++) {

	bool tmp_trk_genuine = false;
	if (MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it))) tmp_trk_genuine = true;
	if (!tmp_trk_genuine) continue;


	if (DebugMode) {
	  if (MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it)).isNull()) {
	    cout << "track matched to TP is NOT uniquely matched to a TP" << endl;
	  }
	  else {
	    edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	    cout << "TP matched to track matched to TP ... tp pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
		 << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z() << endl;
	  }
	  cout << "   ... matched L1 track has pt = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).perp() 
	       << " eta = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).eta()
	       << " phi = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).phi()
	       << " chi2 = " << matchedTracks.at(it)->getChi2(L1Tk_nPar) 
	       << " consistency = " << matchedTracks.at(it)->getStubPtConsistency(L1Tk_nPar) 
	       << " z0 = " << matchedTracks.at(it)->getPOCA(L1Tk_nPar).z() 
	       << " nstub = " << matchedTracks.at(it)->getStubRefs().size();
	  if (tmp_trk_genuine) cout << " (genuine!) " << endl;
	}


	// ----------------------------------------------------------------------------------------------
	// further require L1 track to be (loosely) genuine, that there is only one TP matched to the track
	// + have >= L1Tk_minNStub stubs for it to be a valid match (only relevant is your track collection
	// e.g. stores 3-stub tracks but at plot level you require >= 4 stubs (--> tracklet case)

	int tmp_trk_nstub = matchedTracks.at(it)->getStubRefs().size();

	if (tmp_trk_nstub < L1Tk_minNStub) continue;
	
	float dmatch_pt  = 999;
	float dmatch_eta = 999;
	float dmatch_phi = 999;
	int match_id = 999;

	edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	dmatch_pt  = fabs(my_tp->p4().pt() - tmp_tp_pt);
	dmatch_eta = fabs(my_tp->p4().eta() - tmp_tp_eta);
	dmatch_phi = fabs(my_tp->p4().phi() - tmp_tp_phi);
	match_id = my_tp->pdgId();

	float tmp_trk_chi2dof = (matchedTracks.at(it)->getChi2(L1Tk_nPar)) / (2*tmp_trk_nstub - L1Tk_nPar);
	
	// ensure that track is uniquely matched to the TP we are looking at!
	if (dmatch_pt<0.1 && dmatch_eta<0.1 && dmatch_phi<0.1 && tmp_tp_pdgid==match_id) { 
	  nMatch++;
	  if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
	    i_track = it;
	    i_chi2dof = tmp_trk_chi2dof;
	  }
	}

      }// end loop over matched L1 tracks

    }// end has at least 1 matched L1 track
    // ----------------------------------------------------------------------------------------------


    float tmp_matchtrk_pt   = -999;
    float tmp_matchtrk_eta  = -999;
    float tmp_matchtrk_phi  = -999;
    float tmp_matchtrk_z0   = -999;
    float tmp_matchtrk_d0   = -999;
    float tmp_matchtrk_chi2 = -999;
    //float tmp_matchtrk_p=-999;
    int tmp_matchtrk_nstub  = -999;


    if (nMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to genuine L1 tracks ***" << endl;

    if (nMatch > 0) {
      tmp_matchtrk_pt   = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).perp();
      tmp_matchtrk_eta  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).eta();
      tmp_matchtrk_phi  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).phi();
      tmp_matchtrk_z0   = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).z();
      //tmp_matchtrk_p = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).mag();
      if (L1Tk_nPar == 5) {
	float tmp_matchtrk_x0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).x();
	float tmp_matchtrk_y0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).y();
	tmp_matchtrk_d0 = -tmp_matchtrk_x0*sin(tmp_matchtrk_phi) + tmp_matchtrk_y0*cos(tmp_matchtrk_phi);
      }

      tmp_matchtrk_chi2 = matchedTracks.at(i_track)->getChi2(L1Tk_nPar);
      tmp_matchtrk_nstub  = (int) matchedTracks.at(i_track)->getStubRefs().size();
    }


    m_tp_pt->push_back(tmp_tp_pt);
    m_tp_eta->push_back(tmp_tp_eta);
    m_tp_phi->push_back(tmp_tp_phi);
    m_tp_dxy->push_back(tmp_tp_dxy);
    m_tp_z0->push_back(tmp_tp_z0);
    m_tp_d0->push_back(tmp_tp_d0);
    m_tp_z0_prod->push_back(tmp_tp_z0_prod);
    m_tp_d0_prod->push_back(tmp_tp_d0_prod);
    m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_nmatch->push_back(nMatch);
    m_tp_nstub->push_back(nStubTP);
    m_tp_nstublayers->push_back(nStubLayerTP);
    m_tp_eventid->push_back(tmp_eventid);
    m_tp_charge->push_back(tmp_tp_charge);

    m_matchtrk_pt ->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0 ->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0 ->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2 ->push_back(tmp_matchtrk_chi2);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);
  } //end loop tracking particles

//  FillJets(TrueJetInputs_, true, true, TrackingParticleHandle);
  FillJets(JetInputs_,true, false, TrackingParticleHandle);
fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
fastjet::ClusterSequence csGen(TrueJetInputs_,jet_def);
JetOutputs_.clear();
JetOutputs_=fastjet::sorted_by_pt(csGen.inclusive_jets(0));
for (unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {

  m_genak4chgjet_phi->push_back(JetOutputs_[ijet].phi_std());
  m_genak4chgjet_eta->push_back(JetOutputs_[ijet].eta());
  m_genak4chgjet_pt->push_back(JetOutputs_[ijet].pt());
  m_genak4chgjet_p->push_back(JetOutputs_[ijet].modp());
}
fastjet::ClusterSequence csReco(RecoJetInputs_,jet_def);
JetOutputs_.clear();
JetOutputs_=fastjet::sorted_by_pt(csReco.inclusive_jets(0));
for (unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {
float avgZ=0;
float sumTrackpT=0;
float genuinepT=0;
std::vector<fastjet::PseudoJet> fjConstituents =fastjet::sorted_by_pt(csReco.constituents(JetOutputs_[ijet]));
for(unsigned int i=0; i<fjConstituents.size(); ++i){
	auto index =fjConstituents[i].user_index();
	sumTrackpT=sumTrackpT+m_trk_pt->at(index);
	avgZ=avgZ+m_trk_z0->at(index)*m_trk_pt->at(index);
	if(m_trk_genuine->at(index)>0)genuinepT=genuinepT+m_trk_pt->at(index);	
/*
*/
	}
avgZ=avgZ/sumTrackpT;
m_recojet_truetp_sumpt->push_back(genuinepT/sumTrackpT);
m_recojet_vz->push_back(avgZ);
  m_recojet_ntracks->push_back(fjConstituents.size());
  m_recojet_phi->push_back(JetOutputs_[ijet].phi_std());
  m_recojet_eta->push_back(JetOutputs_[ijet].eta());
  m_recojet_pt->push_back(JetOutputs_[ijet].pt());
  m_recojet_p->push_back(JetOutputs_[ijet].modp());
//if(JetOutputs_[ijet].pt()>10)std::cout<<"Fast Jet Eta Phi "<<JetOutputs_[ijet].eta()<<", "<<JetOutputs_[ijet].phi_std()<<std::endl;
}
eventTree->Fill();

} // end of analyze()
void L1TrackNtupleMaker::FillJets(std::vector<fastjet::PseudoJet>  JetInputs_, bool Prompt,bool TrueTP,edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle){
fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
fastjet::ClusterSequence csPU(JetInputs_,jet_def);
JetOutputs_.clear();
JetOutputs_=fastjet::sorted_by_pt(csPU.inclusive_jets(0));
for (unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {
double avgZ=0;
double sumpt=0; 
double ptThresTracks=0;
double trueSumpt=0;
//fjConstituents.resize(0);
std::vector<fastjet::PseudoJet> fjConstituentsPU =fastjet::sorted_by_pt(csPU.constituents(JetOutputs_[ijet]));
for(unsigned int i=0; i<fjConstituentsPU.size(); ++i){
	auto index =fjConstituentsPU[i].user_index();
	edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, index);
//	if(tp_ptr->pt()>3){
	sumpt=sumpt+tp_ptr->pt();
    float  tmp_tp_vx=tp_ptr->vx();
    float  tmp_tp_vy=tp_ptr->vy();
    float  tmp_tp_vz=tp_ptr->vz();
    float tmp_tp_eta=tp_ptr->eta();
    float tmp_tp_phi=tp_ptr->phi();
    float tmp_tp_pt=tp_ptr->pt();
    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));
	
    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;

    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;

    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;

    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));


    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);

	avgZ=avgZ+tp_ptr->pt()*tmp_tp_z0;
//	}
	int eventId=tp_ptr->eventId().event();
	if(eventId<=0)trueSumpt=trueSumpt+tp_ptr->pt();
	++ptThresTracks;
//	std::cout<<"Track Z PU"<<m_tp_z0->at(index)<<std::endl;

}
//std::cout<<"True Sum pT "<<trueSumpt/sumpt<<std::endl;
/*
m_jetPU_pt->push_back(JetOutputsPU_[ijet].pt());
m_jetPU_p->push_back(JetOutputsPU_[ijet].modp());
m_jetPU_eta->push_back(JetOutputsPU_[ijet].eta());
m_jetPU_phi->push_back(JetOutputsPU_[ijet].phi_std());
*/
avgZ=avgZ/sumpt;
//iterate:
//
TLorentzVector temp;
temp.SetPtEtaPhiE(JetOutputs_[ijet].pt(),JetOutputs_[ijet].eta(),JetOutputs_[ijet].phi_std(),JetOutputs_[ijet].modp());
//temp.SetPz(temp.Pz()-avgZ);
if(TrueTP){
	m_truejet_pt->push_back(temp.Pt());
	m_truejet_p->push_back(temp.P());
	m_truejet_eta->push_back(temp.Eta());
	m_truejet_phi->push_back(temp.Phi());
	m_truejet_ntracks->push_back(ptThresTracks);
	m_truejet_tp_sumpt->push_back(sumpt);
	m_truejet_truetp_sumpt->push_back(trueSumpt/sumpt);
	m_truejet_vz->push_back(avgZ);
	for(unsigned int g=0; g<m_genak4jet_eta->size(); ++g){
		if(m_genak4jet_chgfrac->at(g)/m_genak4jet_p->at(g)>0.98){
		float deta=m_genak4jet_eta->at(g)-temp.Eta();
		float dphi=m_genak4jet_phi->at(g)-temp.Phi();
	        float dR=sqrt((deta*deta)+(dphi*dphi));
		//std::cout<<"Jet eta, phi "<<temp.Eta()<<", "<<temp.Phi()<<"Jet deta, dphi, dR "<<deta<<", "<<dphi<<", "<<dR <<" Jet pT "<<temp.Pt()<<std::endl;	
		}
	}
}
if(Prompt && !TrueTP){
	m_jet_pt->push_back(temp.Pt());
	m_jet_p->push_back(temp.P());
	m_jet_eta->push_back(temp.Eta());
	m_jet_phi->push_back(temp.Phi());
	m_jet_ntracks->push_back(ptThresTracks);
	m_jet_tp_sumpt->push_back(sumpt);
	m_jet_truetp_sumpt->push_back(trueSumpt/sumpt);
	m_jet_vz->push_back(avgZ);
}
/*
if(!Prompt){
	m_jetPU_pt->push_back(temp.Pt());
	m_jetPU_p->push_back(temp.P());
	m_jetPU_eta->push_back(temp.Eta());
	m_jetPU_phi->push_back(temp.Phi());
	m_jetPU_ntracks->push_back(ptThresTracks);
	m_jetPU_tp_sumpt->push_back(sumpt);
	m_jetPU_truetp_sumpt->push_back(trueSumpt/sumpt);
	m_jetPU_vz->push_back(avgZ); 
    }
*/
  }
}
double L1TrackNtupleMaker::FillRecoPrimaryVtx(){
double vtxZ=0;
htmp->Reset();
htmp_weight->Reset();
for(unsigned int z=0; z<m_trk_z0->size(); ++z){
     //if(m_trk_nstub->at(z)<5)continue;
     if(m_trk_pt->at(z)<2)continue;
     if(m_trk_pt->at(z)!=m_trk_pt->at(z))continue;
     if(fabs(m_trk_z0->at(z))>15.)continue;
     if(m_trk_z0->at(z)!=m_trk_z0->at(z))continue;
     //if(m_trk_chi2->at(z)>5)continue;

     float pt=m_trk_pt->at(z);
     if(m_trk_pt->at(z)>200)pt=200;
     htmp -> Fill( m_trk_z0->at(z) );
     htmp_weight -> Fill( m_trk_z0->at(z), pt );
}
  //std::cout<<"Total TP "<<TrackingParticleHandle->size()<<std::endl;
  float zvtx_sliding = -999;
  float sigma_max = -999;
  int nb = htmp -> GetNbinsX();
  for (int i=2; i <= nb-1; i++) {
     float a0 = htmp -> GetBinContent(i-1);
     float a1 = htmp -> GetBinContent(i);
     float a2 = htmp -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp -> GetBinCenter(i-1);
        float z1 = htmp -> GetBinCenter(i);
        float z2 = htmp -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
 zvtx_sliding = -999;
 sigma_max = -999;
 for (int i=2; i <= nb-1; i++) {
     float a0 = htmp_weight -> GetBinContent(i-1);
     float a1 = htmp_weight -> GetBinContent(i);
     float a2 = htmp_weight -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp_weight -> GetBinCenter(i-1);
        float z1 = htmp_weight -> GetBinCenter(i);
        float z2 = htmp_weight -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
vtxZ=zvtx_sliding;
//std::cout<<"htmp Integral "<<htmp_weight->Integral()<<" Mean "<<htmp_weight->GetMean()<<std::endl;
//std::cout<<"Vtx Z "<<vtxZ<<std::endl;
return vtxZ;
}

double L1TrackNtupleMaker::FillPrimaryVtx(edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle,edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle,
const TrackerTopology* const tTopo, const TrackerGeometry* const theTrackerGeom
){
double vtxZ=0;
htmp->Reset();
htmp_weight->Reset();
int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
     ++this_tp;

    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    //int nStubTP = (int) theStubRefs.size();

    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (unsigned int is=0; is<theStubRefs.size(); is++) {
	  
	  //detID of stub
	  DetId detIdStub = theTrackerGeom->idToDet( (theStubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	  
         int layer = -1;
	  if ( detIdStub.subdetId()==StripSubdetector::TOB ) {

	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer-1;
	  }

	  else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer+5;
	  }	  
         if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)hasStubInLayer[layer] = 1;
         else hasStubInLayer[layer] = 2;
}

    //                                                          //treat genuine stubs separately (==2 is genuine, ==1 is not)
    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum=0; isum<11; isum++) {
        if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
        if ( hasStubInLayer[isum] == 2) nStubLayerTP_g += 1;
   }

   if(tp_ptr->pt()>200)continue;
   if(tp_ptr->pt()<2)continue;
   if(fabs(tp_ptr->vz())>25.)continue;
   if(nStubLayerTP<4)continue;
    float  tmp_tp_vx=tp_ptr->vx();
    float  tmp_tp_vy=tp_ptr->vy();
    float  tmp_tp_vz=tp_ptr->vz();
    float tmp_tp_eta=tp_ptr->eta();
    float tmp_tp_phi=tp_ptr->phi();
    float tmp_tp_pt=tp_ptr->pt();
    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));
	
    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;

    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;

    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;

    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));


    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);

     float z=tmp_tp_z0;
     float pt=tp_ptr->pt();
     htmp -> Fill( z );
     htmp_weight -> Fill( z, pt );
}
  //std::cout<<"Total TP "<<TrackingParticleHandle->size()<<std::endl;
  float zvtx_sliding = -999;
  float sigma_max = -999;
  int nb = htmp -> GetNbinsX();
  for (int i=2; i <= nb-1; i++) {
     float a0 = htmp -> GetBinContent(i-1);
     float a1 = htmp -> GetBinContent(i);
     float a2 = htmp -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp -> GetBinCenter(i-1);
        float z1 = htmp -> GetBinCenter(i);
        float z2 = htmp -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
 zvtx_sliding = -999;
 sigma_max = -999;
 for (int i=2; i <= nb-1; i++) {
     float a0 = htmp_weight -> GetBinContent(i-1);
     float a1 = htmp_weight -> GetBinContent(i);
     float a2 = htmp_weight -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp_weight -> GetBinCenter(i-1);
        float z1 = htmp_weight -> GetBinCenter(i);
        float z2 = htmp_weight -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
vtxZ=zvtx_sliding;
//std::cout<<"htmp Integral "<<htmp_weight->Integral()<<" Mean "<<htmp_weight->GetMean()<<std::endl;
//std::cout<<"Vtx Z "<<vtxZ<<std::endl;
return vtxZ;
}
void L1TrackNtupleMaker::FillClusteringMaps(){
for(unsigned int i=0; i<30; ++i)Cluster3x3_zbins[i]->Reset();
for(unsigned int i=0; i<30; ++i)Cluster3x3_shiftzbins[i]->Reset();
for(unsigned int z=0; z<m_trk_z0->size(); ++z){
	//if(m_trk_loose->at(z)<1)continue;
	if(fabs(m_trk_z0->at(z))>15 || m_trk_chi2->at(z)>5 || m_trk_pt->at(z)<2)continue;
	unsigned int zbin=ZbinIndicies->GetXaxis()->FindBin(m_trk_z0->at(z))-1;
       	unsigned int zbinOverlap=ZbinShiftIndicies->GetXaxis()->FindBin(m_trk_z0->at(z))-1;
	//if(zbin==zbinOverlap)std::cout<<"Overlap bin "<<m_trk_z0->at(z)<<" "<<ZbinIndicies->GetXaxis()->GetBinLowEdge(zbin+1)<<" "<<ZbinShiftIndicies->GetXaxis()->GetBinLowEdge(zbinOverlap+1)<<std::endl;	
	Cluster3x3_zbins[zbin]->Fill(m_trk_eta->at(z), m_trk_phi->at(z),m_trk_pt->at(z));
	if(zbinOverlap>1 && zbinOverlap<31)Cluster3x3_shiftzbins[zbinOverlap]->Fill(m_trk_eta->at(z), m_trk_phi->at(z),m_trk_pt->at(z));	
	}
//Handle the overlaps
}
void L1TrackNtupleMaker::FindEtaClusters(){
TH1D*EtaClusterSeeds[60];
float seedThresh=5;
//LAYER 1
for(unsigned int i=0; i<30; ++i){
	TH1D*EtaClusters=(TH1D*)Cluster3x3_zbins[i]->ProjectionX();
	EtaClusterSeeds[i]=(TH1D*)EtaClusters->Clone("EtaClusterSeeds");
	EtaClusterSeeds[i]->Reset();
	bool aboveThresh=true;
	while(aboveThresh){	
	unsigned int Maxetabin=EtaClusters->GetMaximumBin();
	if(EtaClusters->GetBinContent(Maxetabin)<seedThresh){aboveThresh=false;break;}
	EtaClusterSeeds[i]->SetBinContent(Maxetabin, EtaClusters->GetBinContent(Maxetabin));
	//std::cout<<"z-bin "<<i<<" Max Eta Bin "<< Maxetabin<<" Sum pT "<<EtaClusters->GetXaxis()->GetBinCenter(Maxetabin)<<", "<<EtaClusters->GetBinContent(Maxetabin)<<std::endl;
	//if(Maxetabin>1)EtaClusterSeeds[i]->SetBinContent(Maxetabin-1, EtaClusters->GetBinContent(Maxetabin-1));
	//if(Maxetabin<24)EtaClusterSeeds[i]->SetBinContent(Maxetabin+1, EtaClusters->GetBinContent(Maxetabin+1));
	EtaClusters->SetBinContent(Maxetabin, 0);
	if(Maxetabin>1)EtaClusters->SetBinContent(Maxetabin-1, 0);
	if(Maxetabin<24)EtaClusters->SetBinContent(Maxetabin+1, 0);
	}
	EtaClusters=(TH1D*)Cluster3x3_shiftzbins[i]->ProjectionX();
        EtaClusterSeeds[i+30]=(TH1D*)EtaClusters->Clone("EtaClusterSeeds");
        EtaClusterSeeds[i+30]->Reset();
	aboveThresh=true;

	while(aboveThresh){	
	unsigned int Maxetabin=EtaClusters->GetMaximumBin();
	if(EtaClusters->GetBinContent(Maxetabin)<seedThresh){aboveThresh=false;break;}
	EtaClusterSeeds[i+30]->SetBinContent(Maxetabin, EtaClusters->GetBinContent(Maxetabin));
	//if(Maxetabin>1)EtaClusterSeeds[i+30]->SetBinContent(Maxetabin-1, EtaClusters->GetBinContent(Maxetabin-1));
	//if(Maxetabin<24)EtaClusterSeeds[i+30]->SetBinContent(Maxetabin+1, EtaClusters->GetBinContent(Maxetabin+1));
//	std::cout<<"z-bin "<<i+30<<" Max Eta Bin "<< Maxetabin<<" Sum pT "<<EtaClusters->GetXaxis()->GetBinCenter(Maxetabin)<<", "<<EtaClusters->GetBinContent(Maxetabin)<<std::endl;
	
	EtaClusters->SetBinContent(Maxetabin, 0);
	
	if(Maxetabin>1)EtaClusters->SetBinContent(Maxetabin-1, 0);
	if(Maxetabin<24)EtaClusters->SetBinContent(Maxetabin+1, 0);
	}
    }
for(unsigned int j=0; j<60; ++j){
	unsigned int Maxetabin=EtaClusterSeeds[j]->GetMaximumBin();
	
	if(EtaClusterSeeds[j]->GetBinContent(Maxetabin)<seedThresh)continue;
	float sumPtEta=EtaClusterSeeds[j]->GetBinContent(Maxetabin)*EtaClusterSeeds[j]->GetXaxis()->GetBinCenter(Maxetabin);
	float sumPt=EtaClusterSeeds[j]->GetBinContent(Maxetabin);

	if(Maxetabin>1){sumPt=sumPt+EtaClusterSeeds[j]->GetBinContent(Maxetabin-1);sumPtEta=sumPtEta+EtaClusterSeeds[j]->GetBinContent(Maxetabin-1)*EtaClusterSeeds[j]->GetXaxis()->GetBinCenter(Maxetabin-1);}
	if(Maxetabin<24){sumPt=sumPt+EtaClusterSeeds[j]->GetBinContent(Maxetabin+1);sumPtEta=sumPtEta+EtaClusterSeeds[j]->GetBinContent(Maxetabin+1)*EtaClusterSeeds[j]->GetXaxis()->GetBinCenter(Maxetabin+1); }
	//m_recoClusjet_seedz->push_back(j+1);
	int etaMin=Maxetabin;
	int etaMax=Maxetabin;	
	if(Maxetabin>1)etaMin=Maxetabin-1;
	if(Maxetabin<24)etaMax=Maxetabin+1;
	if(j<30){
		TH1D*PhiClusters=(TH1D*)Cluster3x3_zbins[j]->ProjectionY("PhiClusters", etaMin,etaMax );
		 bool aboveThresh=true;
		 while(aboveThresh){
		
			unsigned int Maxphibin=PhiClusters->GetMaximumBin();
			if(PhiClusters->GetBinContent(Maxphibin)<seedThresh || Maxphibin<1){aboveThresh=false;break;}
			int minPhi=Maxphibin;
			int maxPhi=Maxphibin;
			if(Maxphibin==1)minPhi=28;//wraparound	
			if(Maxphibin==28)maxPhi=1;//wraparound	
			float phiSumPt=PhiClusters->GetBinContent(Maxphibin);
			phiSumPt=PhiClusters->GetBinContent(minPhi)+phiSumPt;
			phiSumPt=PhiClusters->GetBinContent(maxPhi)+phiSumPt;
			float phiPtWeight=PhiClusters->GetBinContent(Maxphibin)*PhiClusters->GetXaxis()->GetBinCenter(Maxphibin);
			phiPtWeight=phiPtWeight+PhiClusters->GetBinContent(minPhi)*PhiClusters->GetXaxis()->GetBinCenter(minPhi);
			phiPtWeight=phiPtWeight+PhiClusters->GetBinContent(maxPhi)*PhiClusters->GetXaxis()->GetBinCenter(maxPhi);
			phiPtWeight=phiPtWeight/phiSumPt;
			PhiClusters->SetBinContent(minPhi, 0);	
			PhiClusters->SetBinContent(maxPhi, 0);	
			PhiClusters->SetBinContent(Maxphibin, 0);	
			std::cout<<"z-bin "<<j<<"Eta, Phi bin "<<sumPtEta/sumPt<<" "<< phiPtWeight<<" sumPT "<<Cluster3x3_zbins[j]->GetBinContent(Maxetabin,Maxphibin)<<std::endl;
			m_recoClusjet_seedzbin->push_back(j);
			float sumpT3x3=0;
			//sumpT3x3=Cluster3x3_zbins[j]->GetBinContent(Maxetabin,Maxphibin)+Cluster3x3_zbins[j]
			m_recoClusjet_seedpt->push_back(phiSumPt);
			m_recoClusjet_seedeta->push_back(sumPtEta/sumPt);
			m_recoClusjet_seedphi->push_back(phiPtWeight);
		}
	    }
	  else{
		 bool aboveThresh=true;
		TH1D*PhiClusters=(TH1D*)Cluster3x3_shiftzbins[j-30]->ProjectionY("PhiClusters", etaMin,etaMax );
		
		 while(aboveThresh){
		
			unsigned int Maxphibin=PhiClusters->GetMaximumBin();
			if(PhiClusters->GetBinContent(Maxphibin)<seedThresh || Maxphibin<1){aboveThresh=false;break;}
			//std::cout<<"z-bin "<<j<<"Eta, Phi bin "<<sumPtEta/sumPt<<", phi bin "<<Maxphibin<<" "<< PhiClusters->GetXaxis()->GetBinCenter(Maxphibin) <<std::endl;
			
			int minPhi=Maxphibin;
			int maxPhi=Maxphibin;
			if(Maxphibin==1)minPhi=28;//wraparound	
			if(Maxphibin==28)maxPhi=1;//wraparound	
			float phiSumPt=PhiClusters->GetBinContent(Maxphibin);
			phiSumPt=PhiClusters->GetBinContent(minPhi)+phiSumPt;
			phiSumPt=PhiClusters->GetBinContent(maxPhi)+phiSumPt;
			float phiPtWeight=PhiClusters->GetBinContent(Maxphibin)*PhiClusters->GetXaxis()->GetBinCenter(Maxphibin);
			phiPtWeight=phiPtWeight+PhiClusters->GetBinContent(minPhi)*PhiClusters->GetXaxis()->GetBinCenter(minPhi);
			phiPtWeight=phiPtWeight+PhiClusters->GetBinContent(maxPhi)*PhiClusters->GetXaxis()->GetBinCenter(maxPhi);
			phiPtWeight=phiPtWeight/phiSumPt;
			std::cout<<"z-bin "<<j<<"Eta, Phi bin "<<sumPtEta/sumPt<<" "<< phiPtWeight<<" sumPT "<<Cluster3x3_shiftzbins[j-30]->GetBinContent(Maxetabin,Maxphibin)<<std::endl;
			PhiClusters->SetBinContent(minPhi, 0);	
			PhiClusters->SetBinContent(maxPhi, 0);	
			PhiClusters->SetBinContent(Maxphibin, 0);	
			//std::cout<<"z-bin "<<j<<"Eta, Phi bin "<<sumPtEta/sumPt<<" "<< phiPtWeight<<std::endl;	
			m_recoClusjet_seedzbin->push_back(j);
			m_recoClusjet_seedpt->push_back(phiSumPt);
			m_recoClusjet_seedeta->push_back(sumPtEta/sumPt);
			m_recoClusjet_seedphi->push_back(phiPtWeight);
				
		}
	     }
	}
// there will be two copies of seeds for the overlap z-bins: choose based on largest sumPT
	float maxSumpt=0;
	int primaryseed=-1;
	for(unsigned int s=0; s<m_recoClusjet_seedzbin->size(); ++s){
		int zindex=m_recoClusjet_seedzbin->at(s);
		if(maxSumpt<m_recoClusjet_seedpt->at(s)){;
		int seedindex=s;
		maxSumpt=m_recoClusjet_seedpt->at(s);
		primaryseed=s;
		}
	}
	int zindex=-1;
	if(primaryseed>-1){
	std::cout<<"Z bin "<<m_recoClusjet_seedzbin->at(primaryseed)<<std::endl;
	zindex=m_recoClusjet_seedzbin->at(primaryseed);	

	}
	//grab all the tracks
	for(unsigned int s=0; s<m_recoClusjet_seedzbin->size(); ++s){
		if(m_recoClusjet_seedzbin->at(s)!=zindex)continue;
	        float seedEta=m_recoClusjet_seedeta->at(s);
		float seedPhi=m_recoClusjet_seedphi->at(s);
		TLorentzVector JetCandidate;	
		 for(unsigned int z=0; z<m_trk_z0->size(); ++z){
			if(m_trk_pt->at(z)<2)continue;
			if(m_trk_chi2->at(z)>5)continue;
			if(fabs(m_trk_z0->at(z))>15)continue;
			//if(primaryseed<30)        
			
			float zbinCenter=ZbinIndicies->GetXaxis()->GetBinCenter(ZbinIndicies->GetXaxis()->FindBin(m_trk_z0->at(z)));

			if(primaryseed>30)zbinCenter=ZbinShiftIndicies->GetXaxis()->GetBinCenter(ZbinShiftIndicies->GetXaxis()->FindBin(m_trk_z0->at(z)));
			float dz=zbinCenter-m_trk_z0->at(z);
			float deta=seedEta-m_trk_eta->at(z);
			float dphi=seedPhi-m_trk_phi->at(z);
			float dR=sqrt((deta*deta)+(dphi*dphi));
			if(dR>0.4 || fabs(dz)>0.5)continue;
			TLorentzVector temp;
			temp.SetPtEtaPhiE(m_trk_pt->at(z), m_trk_eta->at(z), m_trk_phi->at(z), m_trk_p->at(z));
			JetCandidate=JetCandidate+temp;	
		}	
		std::cout<<"Jet Candidtate "<<JetCandidate.Pt()<<", "<<JetCandidate.Eta()<<", "<<JetCandidate.Phi()<<std::endl;
	}
	

	//	for(unsigned int z=0; z<SeedBinZ.size(); ++z)std::cout<<"Max Seed z -bins "<<SeedBinZ[z]<<std::endl;
}
void L1TrackNtupleMaker::FindClusters(){
for(unsigned int i=0; i<30; ++i){
	bool Seeds=true;
	float seedThres=10;
	std::vector<float>tmpEta;
	std::vector<float>tmpPhi;
	std::vector<float>tmpPt;
	std::vector<float>tmpP;
	std::vector<float>tmpZ;
	  for(unsigned int z=0; z<m_trk_z0->size(); ++z){
	        //if(m_trk_loose->at(z)<1)continue;

		if(m_trk_chi2->at(z)<5 && m_trk_pt->at(z)>2 && m_trk_z0->at(z)>=ZbinIndicies->GetXaxis()->GetBinLowEdge(i+1) && (m_trk_z0->at(z)<ZbinIndicies->GetXaxis()->GetBinUpEdge(i+1))){
		tmpEta.push_back(m_trk_eta->at(z));
		tmpP.push_back(m_trk_p->at(z));
		tmpPt.push_back(m_trk_pt->at(z));
		tmpPhi.push_back(m_trk_phi->at(z));
		tmpZ.push_back(m_trk_z0->at(z));
		//std::cout<<"All Tracks Pt, eta, phi"<<m_trk_pt->at(z)<<", "<<	m_trk_eta->at(z)<<", "<<m_trk_phi->at(z)<<std::endl;
		}
	}
	while(Seeds){
	unsigned int ix=Cluster3x3_zbins[i]->ProjectionX()->GetMaximumBin();
	unsigned int iy=Cluster3x3_zbins[i]->ProjectionY()->GetMaximumBin();
	if(Cluster3x3_zbins[i]->GetBinContent(ix,iy)<seedThres){Seeds=false;break;}
	  TLorentzVector JetCandidate;
	 float avgZ=0;
	float sumpT=0;
	int ntracks=0;
	std::vector<int>UnUsedtracks;
	  for(unsigned int z=0; z<tmpP.size(); ++z){
			float deta=tmpEta[z]-Cluster3x3_zbins[i]->ProjectionX()->GetBinCenter(ix);				
			float dphi=tmpPhi[z]-Cluster3x3_zbins[i]->ProjectionY()->GetBinCenter(iy);				
			float dR=sqrt((deta*deta)+(dphi*dphi));
			if(dR>0.3){UnUsedtracks.push_back(z);continue;}
			++ntracks;
			sumpT=sumpT+tmpPt[z];
			avgZ=avgZ+tmpPt[z]*tmpZ[z];
			TLorentzVector temp;
			temp.SetPtEtaPhiE(tmpPt[z], tmpEta[z],tmpPhi[z], tmpP[z]);
			JetCandidate=JetCandidate+temp;
			//std::cout<<"Seed "<<Cluster3x3_zbins[i]->GetBinContent(ix,iy)<<" eta, phi "<<Cluster3x3_zbins[i]->ProjectionX()->GetBinCenter(ix)<<", "<<Cluster3x3_zbins[i]->ProjectionY()->GetBinCenter(iy)<<" Matched Track "<<temp.Eta()<<", "<<temp.Phi()<<std::endl;	
	  }
        std::vector<float>tmp2Eta;
        std::vector<float>tmp2Phi;
        std::vector<float>tmp2Pt;
        std::vector<float>tmp2P;
        std::vector<float>tmp2Z;
 if(tmpP.size()>0){ 
m_recoClusjet_avgZNum->push_back(avgZ);
  m_recoClusjet_eta->push_back(JetCandidate.Eta());
  m_recoClusjet_seedpt->push_back(Cluster3x3_zbins[i]->GetBinContent(ix,iy));
	avgZ=avgZ/sumpT;
  m_recoClusjet_vz->push_back(avgZ);
  m_recoClusjet_phi->push_back(JetCandidate.Phi());
  m_recoClusjet_pt->push_back(JetCandidate.Pt());
  m_recoClusjet_p->push_back(JetCandidate.P());
  m_recoClusjet_ntracks->push_back(ntracks);
//if(ntracks==1)std::cout<<"Jet Candidate for single Track for the seed "<<JetCandidate.Pt()<<", "<<JetCandidate.Eta()<<", "<<JetCandidate.Phi()<<" z position "<<avgZ<<" True Z "<<m_pv_L1->at(0)<<std::endl;
//else std::cout<<"Jet Candidate Track for the seed "<<JetCandidate.Pt()<<", "<<JetCandidate.Eta()<<", "<<JetCandidate.Phi()<<" z position "<<avgZ<<" True Z "<<m_pv_L1->at(0)<<std::endl;
}
	 for(unsigned int k=0; k<UnUsedtracks.size(); ++k){
		int index=UnUsedtracks[k];
		tmp2Eta.push_back(tmpEta[index]);
		tmp2Phi.push_back(tmpPhi[index]);
		tmp2Pt.push_back(tmpPt[index]);
		tmp2P.push_back(tmpP[index]);
		tmp2Z.push_back(tmpZ[index]);
	}
	tmpEta=tmp2Eta;	
	tmpPhi=tmp2Phi;
	tmpPt=tmp2Pt;
        tmpP=tmp2P;
        tmpZ=tmp2Z;

	//std::cout<<"Total Tracks "<<m_trk_z0->size()<<" Used Tracks "<<usedTrack.size()<<std::endl;
		//if(sumpT<0.00001)continue;

		//mask the seeds in the 3x3
	for(int ibin=-1; ibin<2; ++ibin){
			for(int j=-1; j<2; ++j){
				
				Cluster3x3_zbins[i]->SetBinContent(ix+ibin, iy+j, 0);
			}
		}
	}
   }

}
float L1TrackNtupleMaker::FindMaxHTSlice(){
TH1F*MaxHTSlice=(TH1F*)ZbinIndicies->Clone("MaxHTSlice");
for(unsigned int z=0; z<m_recoClusjet_vz->size(); ++z){
	if(m_recoClusjet_ntracks->at(z)>1)MaxHTSlice->Fill(m_recoClusjet_vz->at(z),m_recoClusjet_pt->at(z));
}
float vtxZ=0;
int highestHTBin=MaxHTSlice->GetMaximumBin();
float zLowEdge=MaxHTSlice->GetXaxis()->GetBinLowEdge(highestHTBin);
float maxJetpT=0;
int maxJetIndex=-1;
for(unsigned int z=0; z<m_recoClusjet_vz->size(); ++z){
	//std::cout<<"Z edge for High HT bin "<<zLowEdge<<", "<<zLowEdge+1.0<<std::endl;
	if(m_recoClusjet_vz->at(z)>zLowEdge && m_recoClusjet_vz->at(z)<(zLowEdge+1.0 ) && maxJetpT<m_recoClusjet_pt->at(z)){maxJetpT=m_recoClusjet_pt->at(z); maxJetIndex=z;}
	
}
if(maxJetIndex<0)return -999.;
return m_recoClusjet_vz->at(maxJetIndex);
}
void L1TrackNtupleMaker::MergePromptCandidates(){

for(unsigned int z=0; z<m_recoClusjet_vz->size(); ++z){
	bool Merged=false;
	if(fabs(m_pv_HTz->at(0)-m_recoClusjet_vz->at(z))>0.5)continue;      
//	int zindex=ZbinIndicies->GetXaxis()->FindBin(m_pv_HTz->at(0));
	//if( m_recoClusjet_vz->at(z)<ZbinIndicies->GetXaxis()->GetBinLowEdge(zindex)-1.0 || m_recoClusjet_vz->at(z)>ZbinIndicies->GetXaxis()->GetBinUpEdge(zindex)+1.0)continue;	
//	if( m_recoClusjet_vz->at(z)<ZbinIndicies->GetXaxis()->GetBinLowEdge(zindex)-1.0 ||  m_recoClusjet_vz->at(z)>ZbinIndicies->GetXaxis()->GetBinUpEdge(zindex)+1.0)continue;
//	std::cout<<"Prompt Candidate Jets "<<m_recoClusjet_eta->at(z)<<", "<<m_recoClusjet_phi->at(z)<<std::endl;
		for(unsigned int z2=1; z2<m_recoClusjet_vz->size(); ++z2){
		if(z==z2 || z>z2)continue;
		if(fabs(m_pv_HTz->at(0)-m_recoClusjet_vz->at(z2))>0.5)continue;      
		float deta=m_recoClusjet_eta->at(z)-m_recoClusjet_eta->at(z2);		
		float dphi=m_recoClusjet_phi->at(z)-m_recoClusjet_phi->at(z2);
		float dz=fabs(m_recoClusjet_vz->at(z)-m_recoClusjet_vz->at(z2));	
		if(sqrt(deta*deta+dphi*dphi)<0.3){
		
		//std::cout<<"Merge these "<<dz<<std::endl;
		TLorentzVector temp1;
		temp1.SetPtEtaPhiE(m_recoClusjet_pt->at(z), m_recoClusjet_eta->at(z), m_recoClusjet_phi->at(z), m_recoClusjet_p->at(z));
		TLorentzVector temp2;
		temp2.SetPtEtaPhiE(m_recoClusjet_pt->at(z2), m_recoClusjet_eta->at(z2), m_recoClusjet_phi->at(z2), m_recoClusjet_p->at(z2));
		TLorentzVector Merge; Merge=temp1+temp2;
			
		float avgZ=m_recoClusjet_avgZNum->at(z)+m_recoClusjet_avgZNum->at(z2);
		avgZ=avgZ/(m_recoClusjet_ntracks->at(z)+m_recoClusjet_ntracks->at(z2));
		m_recoPromptjet_eta->push_back(Merge.Eta());
		m_recoPromptjet_phi->push_back(Merge.Phi());
		m_recoPromptjet_pt->push_back(Merge.Pt());
		m_recoPromptjet_p->push_back(Merge.P());
		m_recoPromptjet_ntracks->push_back(m_recoClusjet_ntracks->at(z)+m_recoClusjet_ntracks->at(z2));
  		m_recoPromptjet_vz->push_back(avgZ);
		//std::cout<<"Old Jet "<<temp1.Eta()<<", "<<temp1.Phi()<<", "<<temp1.Pt()<<" z "<<m_recoClusjet_vz->at(z)<<std::endl;	
		//std::cout<<"Merged Jet "<<Merge.Eta()<<", "<<Merge.Phi()<<", "<<Merge.Pt()<<" z "<<m_recoClusjet_vz->at(z2)<<std::endl;	
		Merged=true;
		break;
		}
	}
for(unsigned int z3=0; z3<m_recoClusjet_vz->size(); ++z3){
	if(fabs(m_pv_HTz->at(0)-m_recoClusjet_vz->at(z3))>0.5 && fabs(m_pv_HTz->at(0)-m_recoClusjet_vz->at(z3))<1.0 ){
		float deta=m_recoClusjet_eta->at(z)-m_recoClusjet_eta->at(z3);		
		float dphi=m_recoClusjet_phi->at(z)-m_recoClusjet_phi->at(z3);
		float dz=fabs(m_recoClusjet_vz->at(z)-m_recoClusjet_vz->at(z3));	
		///if(sqrt(deta*deta+dphi*dphi)<0.3)std::cout<<"Merge these "<<dz<<std::endl;
		if(sqrt(deta*deta+dphi*dphi)<0.3 && dz<0.5){std::cout<<"Merge these "<<dz<<std::endl;
		TLorentzVector temp1;
		temp1.SetPtEtaPhiE(m_recoClusjet_pt->at(z), m_recoClusjet_eta->at(z), m_recoClusjet_phi->at(z), m_recoClusjet_p->at(z));
		TLorentzVector temp2;
		temp2.SetPtEtaPhiE(m_recoClusjet_pt->at(z3), m_recoClusjet_eta->at(z3), m_recoClusjet_phi->at(z3), m_recoClusjet_p->at(z3));
		TLorentzVector Merge; Merge=temp1+temp2;
		float avgZ=m_recoClusjet_avgZNum->at(z)+m_recoClusjet_avgZNum->at(z3);
		avgZ=avgZ/(m_recoClusjet_ntracks->at(z)+m_recoClusjet_ntracks->at(z3));
		m_recoPromptjet_eta->push_back(Merge.Eta());
		m_recoPromptjet_phi->push_back(Merge.Phi());
		m_recoPromptjet_pt->push_back(Merge.Pt());
		m_recoPromptjet_p->push_back(Merge.P());
		m_recoPromptjet_ntracks->push_back(m_recoClusjet_ntracks->at(z)+m_recoClusjet_ntracks->at(z3));
  		m_recoPromptjet_vz->push_back(avgZ);
		std::cout<<"Old Jet "<<temp1.Eta()<<", "<<temp1.Phi()<<", "<<temp1.Pt()<<" z "<<m_recoClusjet_vz->at(z)<<std::endl;	
		std::cout<<"Merged Jet "<<Merge.Eta()<<", "<<Merge.Phi()<<", "<<Merge.Pt()<<" z "<<m_recoClusjet_vz->at(z3)<<std::endl;	
                Merged=true;
		break;
		}	
	//std::cout<<"Neighbor Candidate Jets "<<m_recoClusjet_eta->at(z)<<", "<<m_recoClusjet_phi->at(z)<<std::endl;
	}
      }
	if(!Merged){
		m_recoPromptjet_eta->push_back(m_recoClusjet_eta->at(z));
                m_recoPromptjet_phi->push_back(m_recoClusjet_phi->at(z));
                m_recoPromptjet_pt->push_back(m_recoClusjet_pt->at(z));
                m_recoPromptjet_p->push_back(m_recoClusjet_p->at(z));
                m_recoPromptjet_ntracks->push_back(m_recoClusjet_ntracks->at(z));
                m_recoPromptjet_vz->push_back(m_recoClusjet_vz->at(z));
		}
    }
	std::cout<<"reco NJets "<<m_recoPromptjet_eta->size()<<std::endl;
	for(unsigned int j=0; j<m_recoPromptjet_eta->size(); ++j){
	std::cout<<" reco Jet "<<m_recoPromptjet_eta->at(j)<<", "<<m_recoPromptjet_phi->at(j)<<std::endl;
	}
}
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackNtupleMaker);
