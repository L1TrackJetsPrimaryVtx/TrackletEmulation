//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 23 12:17:05 2017 by ROOT version 5.34/37
// from TTree eventTree/Event tree
// found on file: TTBarFastJetsAllTracks.root
//////////////////////////////////////////////////////////

#ifndef PUAnalyzer_h
#define PUAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <TLorentzVector.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;
class PUAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
    double rho;
   vector<float>   *trk_p;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_chi2;
   vector<int>     *trk_nstub;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_consistency;
   vector<int>     *trk_genuine;
   vector<int>     *trk_loose;
   vector<int>     *trk_unknown;
   vector<int>     *trk_combinatoric;
   vector<int>     *trk_fake;
   vector<int>     *trk_matchtp_pdgid;
   vector<float>   *trk_matchtp_pt;
   vector<float>   *trk_matchtp_eta;
   vector<float>   *trk_matchtp_phi;
   vector<float>   *trk_matchtp_z0;
   vector<float>   *trk_matchtp_dxy;
   vector<int>     *trk_injet;
   vector<int>     *trk_injet_highpt;
   vector<float>   *tp_p;
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_momid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nloosematch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_eventid;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_ngenstublayer;
   vector<int>     *tp_nstublayer;
   vector<int>     *tp_injet;
   vector<int>     *tp_injet_highpt;
   vector<float>   *matchtrk_pt;
   vector<float>   *matchtrk_eta;
   vector<float>   *matchtrk_phi;
   vector<float>   *matchtrk_z0;
   vector<float>   *matchtrk_d0;
   vector<float>   *matchtrk_chi2;
   vector<int>     *matchtrk_nstub;
   vector<float>   *matchtrk_consistency;
   vector<int>     *matchtrk_injet;
   vector<int>     *matchtrk_injet_highpt;
   vector<float>   *loosematchtrk_pt;
   vector<float>   *loosematchtrk_eta;
   vector<float>   *loosematchtrk_phi;
   vector<float>   *loosematchtrk_z0;
   vector<float>   *loosematchtrk_d0;
   vector<float>   *loosematchtrk_chi2;
   vector<int>     *loosematchtrk_nstub;
   vector<float>   *loosematchtrk_consistency;
   vector<int>     *loosematchtrk_injet;
   vector<int>     *loosematchtrk_injet_highpt;
   vector<float>   *reliso;
   vector<float>   *absiso;
   vector<float>   *pv_L1reco;
   vector<float>   *pv_L1;
   vector<float>   *pv_MC;
    vector<int>   *MC_lep;

   vector<float>   *genjetak4_eta;
   vector<float>   *genjetak4_pt;
   vector<float>   *genjetak4_phi;
   vector<float>   *genjetak4_p;
    vector<float>   *genjetak4_metfrac;
    vector<float>   *genjetak4_chgfrac;
   vector<float>   *genjetak8_eta;
   vector<float>   *genjetak8_pt;
   vector<float>   *genjetak8_phi;
   vector<float>   *genjetak8_p;
   vector<float>   *jet_eta;
   vector<float>   *jet_vz;
   vector<float>   *jet_p;
   vector<float>   *jet_pt;
   vector<float>   *jet_phi;
   vector<int>     *jet_ntracks;
   vector<float>   *jet_tp_sumpt;
    vector<float>   *jet_truetp_sumpt;
    vector<float>   *recoClusjet_eta;
    vector<float>   *recoClusjet_vz;
    vector<float>   *recoClusjet_p;
    vector<float>   *recoClusjet_pt;
    vector<float>   *recoClusjet_phi;
    vector<int>     *recoClusjet_ntracks;
    
    vector<float>   *recojet_eta;
    vector<float>   *recojet_vz;
    vector<float>   *recojet_p;
    vector<float>   *recojet_pt;
    vector<float>   *recojet_phi;
    vector<int>     *recojet_ntracks;
   vector<float>   *jetPU_eta;
   vector<float>   *jetPU_vz;
   vector<float>   *jetPU_p;
   vector<float>   *jetPU_pt;
   vector<float>   *jetPU_phi;
   vector<int>     *jetPU_ntracks;
   vector<float>   *jetPU_tp_sumpt;
    vector<float>   *jetPU_truetp_sumpt;

   vector<float>   *jet_matchtrk_sumpt;
   vector<float>   *jet_loosematchtrk_sumpt;
   vector<float>   *jet_trk_sumpt;

   // List of branches
   TBranch        *b_trk_p;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_nstub;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_consistency;   //!
   TBranch        *b_trk_genuine;   //!
   TBranch        *b_trk_loose;   //!
   TBranch        *b_trk_unknown;   //!
   TBranch        *b_trk_combinatoric;   //!
   TBranch        *b_trk_fake;   //!
   TBranch        *b_trk_matchtp_pdgid;   //!
   TBranch        *b_trk_matchtp_pt;   //!
   TBranch        *b_trk_matchtp_eta;   //!
   TBranch        *b_trk_matchtp_phi;   //!
   TBranch        *b_trk_matchtp_z0;   //!
   TBranch        *b_trk_matchtp_dxy;   //!
   TBranch        *b_trk_injet;   //!
   TBranch        *b_trk_injet_highpt;   //!
   TBranch        *b_tp_p;   //!
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_momid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nloosematch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_ngenstublayer;   //!
   TBranch        *b_tp_nstublayer;   //!
   TBranch        *b_tp_injet;   //!
   TBranch        *b_tp_injet_highpt;   //!
   TBranch        *b_matchtrk_pt;   //!
   TBranch        *b_matchtrk_eta;   //!
   TBranch        *b_matchtrk_phi;   //!
   TBranch        *b_matchtrk_z0;   //!
   TBranch        *b_matchtrk_d0;   //!
   TBranch        *b_matchtrk_chi2;   //!
   TBranch        *b_matchtrk_nstub;   //!
   TBranch        *b_matchtrk_consistency;   //!
   TBranch        *b_matchtrk_injet;   //!
   TBranch        *b_matchtrk_injet_highpt;   //!
   TBranch        *b_loosematchtrk_pt;   //!
   TBranch        *b_loosematchtrk_eta;   //!
   TBranch        *b_loosematchtrk_phi;   //!
   TBranch        *b_loosematchtrk_z0;   //!
   TBranch        *b_loosematchtrk_d0;   //!
   TBranch        *b_loosematchtrk_chi2;   //!
   TBranch        *b_loosematchtrk_nstub;   //!
   TBranch        *b_loosematchtrk_consistency;   //!
   TBranch        *b_loosematchtrk_injet;   //!
   TBranch        *b_loosematchtrk_injet_highpt;   //!
   TBranch        *b_reliso;   //!
   TBranch        *b_absiso;   //!
   TBranch        *b_pv_L1reco;   //!
   TBranch        *b_pv_L1;   //!
   TBranch        *b_pv_MC;   //!
    TBranch        *b_MC_lep;   //!

   TBranch        *b_genjetak4_eta;   //!
   TBranch        *b_genjetak4_pt;   //!
   TBranch        *b_genjetak4_phi;   //!
   TBranch        *b_genjetak4_p;   //!
    TBranch        *b_genjetak4_metfrac;   //!
    TBranch        *b_genjetak4_chgfrac;   //!
   TBranch        *b_genjetak8_eta;   //!
   TBranch        *b_genjetak8_pt;   //!
   TBranch        *b_genjetak8_phi;   //!
   TBranch        *b_genjetak8_p;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_vz;   //!
   TBranch        *b_jet_p;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_ntracks;   //!
   TBranch        *b_jet_tp_sumpt;   //!
    TBranch        *b_jet_truetp_sumpt;   //!
    TBranch        *b_recoClusjet_eta;   //!
    TBranch        *b_recoClusjet_vz;   //!
    TBranch        *b_recoClusjet_p;   //!
    TBranch        *b_recoClusjet_pt;   //!
    TBranch        *b_recoClusjet_phi;   //!
    TBranch        *b_recoClusjet_ntracks;   //!
    TBranch        *b_recojet_eta;   //!
    TBranch        *b_recojet_vz;   //!
    TBranch        *b_recojet_p;   //!
    TBranch        *b_recojet_pt;   //!
    TBranch        *b_recojet_phi;   //!
    TBranch        *b_recojet_ntracks;   //!
   TBranch        *b_jetPU_eta;   //!
   TBranch        *b_jetPU_vz;   //!
   TBranch        *b_jetPU_p;   //!
   TBranch        *b_jetPU_pt;   //!
   TBranch        *b_jetPU_phi;   //!
   TBranch        *b_jetPU_ntracks;   //!
   TBranch        *b_jetPU_tp_sumpt;   //!
    TBranch        *b_jetPU_truetp_sumpt;   //!

   TBranch        *b_jet_matchtrk_sumpt;   //!
   TBranch        *b_jet_loosematchtrk_sumpt;   //!
   TBranch        *b_jet_trk_sumpt;   //!
    TBranch *b_rho;
   PUAnalyzer(TTree *tree=0);
   virtual ~PUAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PUAnalyzer_cxx
PUAnalyzer::PUAnalyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MinBias_D13_PU200New.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MinBias_D13_PU200New.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("MinBias_D13_PU200New.root:/L1TrackNtuple");
      dir->GetObject("eventTree",tree);

   }
   Init(tree);
}

PUAnalyzer::~PUAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PUAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PUAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PUAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_p = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_chi2 = 0;
   trk_nstub = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_consistency = 0;
   trk_genuine = 0;
   trk_loose = 0;
   trk_unknown = 0;
   trk_combinatoric = 0;
   trk_fake = 0;
   trk_matchtp_pdgid = 0;
   trk_matchtp_pt = 0;
   trk_matchtp_eta = 0;
   trk_matchtp_phi = 0;
   trk_matchtp_z0 = 0;
   trk_matchtp_dxy = 0;
   trk_injet = 0;
   trk_injet_highpt = 0;
   tp_p = 0;
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_pdgid = 0;
   tp_momid = 0;
   tp_nmatch = 0;
   tp_nloosematch = 0;
   tp_nstub = 0;
   tp_eventid = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_ngenstublayer = 0;
   tp_nstublayer = 0;
   tp_injet = 0;
   tp_injet_highpt = 0;
   matchtrk_pt = 0;
   matchtrk_eta = 0;
   matchtrk_phi = 0;
   matchtrk_z0 = 0;
   matchtrk_d0 = 0;
   matchtrk_chi2 = 0;
   matchtrk_nstub = 0;
   matchtrk_consistency = 0;
   matchtrk_injet = 0;
   matchtrk_injet_highpt = 0;
   loosematchtrk_pt = 0;
   loosematchtrk_eta = 0;
   loosematchtrk_phi = 0;
   loosematchtrk_z0 = 0;
   loosematchtrk_d0 = 0;
   loosematchtrk_chi2 = 0;
   loosematchtrk_nstub = 0;
   loosematchtrk_consistency = 0;
   loosematchtrk_injet = 0;
   loosematchtrk_injet_highpt = 0;
   reliso = 0;
   absiso = 0;
   pv_L1reco = 0;
   pv_L1 = 0;
    pv_MC = 0;
    MC_lep = 0;
    genjetak4_metfrac = 0;
    genjetak4_chgfrac = 0;

   genjetak4_eta = 0;
   genjetak4_pt = 0;
   genjetak4_phi = 0;
   genjetak4_p = 0;
   genjetak8_eta = 0;
   genjetak8_pt = 0;
   genjetak8_phi = 0;
   genjetak8_p = 0;
   jet_eta = 0;
   jet_vz = 0;
   jet_p = 0;
   jet_pt = 0;
   jet_phi = 0;
   jet_ntracks = 0;
   jet_tp_sumpt = 0;
    jet_truetp_sumpt = 0;

   jetPU_eta = 0;
   jetPU_vz = 0;
   jetPU_p = 0;
   jetPU_pt = 0;
   jetPU_phi = 0;
   jetPU_ntracks = 0;
   jetPU_tp_sumpt = 0;
    jetPU_truetp_sumpt = 0;

   jet_matchtrk_sumpt = 0;
   jet_loosematchtrk_sumpt = 0;
   jet_trk_sumpt = 0;
    recojet_eta = 0;
    recojet_vz = 0;
    recojet_p = 0;
    recojet_pt = 0;
    recojet_phi = 0;
    recojet_ntracks = 0;
    
    recoClusjet_eta = 0;
    recoClusjet_vz = 0;
    recoClusjet_p = 0;
    recoClusjet_pt = 0;
    recoClusjet_phi = 0;
    recoClusjet_ntracks = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
    fChain->SetBranchAddress("rho", &rho, &b_rho);

   fChain->SetBranchAddress("trk_p", &trk_p, &b_trk_p);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_consistency", &trk_consistency, &b_trk_consistency);
   fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
   fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
   fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
   fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
   fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
   fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
   fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
   fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
   fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
   fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
   fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
   fChain->SetBranchAddress("trk_injet", &trk_injet, &b_trk_injet);
   fChain->SetBranchAddress("trk_injet_highpt", &trk_injet_highpt, &b_trk_injet_highpt);
   fChain->SetBranchAddress("tp_p", &tp_p, &b_tp_p);
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_momid", &tp_momid, &b_tp_momid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nloosematch", &tp_nloosematch, &b_tp_nloosematch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_ngenstublayer", &tp_ngenstublayer, &b_tp_ngenstublayer);
   fChain->SetBranchAddress("tp_nstublayer", &tp_nstublayer, &b_tp_nstublayer);
   fChain->SetBranchAddress("tp_injet", &tp_injet, &b_tp_injet);
   fChain->SetBranchAddress("tp_injet_highpt", &tp_injet_highpt, &b_tp_injet_highpt);
   fChain->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
   fChain->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
   fChain->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
   fChain->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
   fChain->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
   fChain->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
   fChain->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
   fChain->SetBranchAddress("matchtrk_consistency", &matchtrk_consistency, &b_matchtrk_consistency);
   fChain->SetBranchAddress("matchtrk_injet", &matchtrk_injet, &b_matchtrk_injet);
   fChain->SetBranchAddress("matchtrk_injet_highpt", &matchtrk_injet_highpt, &b_matchtrk_injet_highpt);
   fChain->SetBranchAddress("loosematchtrk_pt", &loosematchtrk_pt, &b_loosematchtrk_pt);
   fChain->SetBranchAddress("loosematchtrk_eta", &loosematchtrk_eta, &b_loosematchtrk_eta);
   fChain->SetBranchAddress("loosematchtrk_phi", &loosematchtrk_phi, &b_loosematchtrk_phi);
   fChain->SetBranchAddress("loosematchtrk_z0", &loosematchtrk_z0, &b_loosematchtrk_z0);
   fChain->SetBranchAddress("loosematchtrk_d0", &loosematchtrk_d0, &b_loosematchtrk_d0);
   fChain->SetBranchAddress("loosematchtrk_chi2", &loosematchtrk_chi2, &b_loosematchtrk_chi2);
   fChain->SetBranchAddress("loosematchtrk_nstub", &loosematchtrk_nstub, &b_loosematchtrk_nstub);
   fChain->SetBranchAddress("loosematchtrk_consistency", &loosematchtrk_consistency, &b_loosematchtrk_consistency);
   fChain->SetBranchAddress("loosematchtrk_injet", &loosematchtrk_injet, &b_loosematchtrk_injet);
   fChain->SetBranchAddress("loosematchtrk_injet_highpt", &loosematchtrk_injet_highpt, &b_loosematchtrk_injet_highpt);
   fChain->SetBranchAddress("reliso", &reliso, &b_reliso);
   fChain->SetBranchAddress("absiso", &absiso, &b_absiso);
   fChain->SetBranchAddress("pv_L1reco", &pv_L1reco, &b_pv_L1reco);
   fChain->SetBranchAddress("pv_L1", &pv_L1, &b_pv_L1);
   fChain->SetBranchAddress("pv_MC", &pv_MC, &b_pv_MC);
    fChain->SetBranchAddress("MC_lep", &MC_lep, &b_MC_lep);

   fChain->SetBranchAddress("genjetak4_eta", &genjetak4_eta, &b_genjetak4_eta);
   fChain->SetBranchAddress("genjetak4_pt", &genjetak4_pt, &b_genjetak4_pt);
   fChain->SetBranchAddress("genjetak4_phi", &genjetak4_phi, &b_genjetak4_phi);
   fChain->SetBranchAddress("genjetak4_p", &genjetak4_p, &b_genjetak4_p);
    fChain->SetBranchAddress("genjetak4_metfrac", &genjetak4_metfrac, &b_genjetak4_metfrac);
    fChain->SetBranchAddress("genjetak4_chgfrac", &genjetak4_chgfrac, &b_genjetak4_chgfrac);
   fChain->SetBranchAddress("genjetak8_eta", &genjetak8_eta, &b_genjetak8_eta);
   fChain->SetBranchAddress("genjetak8_pt", &genjetak8_pt, &b_genjetak8_pt);
   fChain->SetBranchAddress("genjetak8_phi", &genjetak8_phi, &b_genjetak8_phi);
   fChain->SetBranchAddress("genjetak8_p", &genjetak8_p, &b_genjetak8_p);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_vz", &jet_vz, &b_jet_vz);
   fChain->SetBranchAddress("jet_p", &jet_p, &b_jet_p);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_ntracks", &jet_ntracks, &b_jet_ntracks);
   fChain->SetBranchAddress("jet_tp_sumpt", &jet_tp_sumpt, &b_jet_tp_sumpt);
   fChain->SetBranchAddress("jetPU_eta", &jetPU_eta, &b_jetPU_eta);
   fChain->SetBranchAddress("jetPU_vz", &jetPU_vz, &b_jetPU_vz);
   fChain->SetBranchAddress("jetPU_p", &jetPU_p, &b_jetPU_p);
   fChain->SetBranchAddress("jetPU_pt", &jetPU_pt, &b_jetPU_pt);
   fChain->SetBranchAddress("jetPU_phi", &jetPU_phi, &b_jetPU_phi);
   fChain->SetBranchAddress("jetPU_ntracks", &jetPU_ntracks, &b_jetPU_ntracks);
   fChain->SetBranchAddress("jetPU_tp_sumpt", &jetPU_tp_sumpt, &b_jetPU_tp_sumpt);
    fChain->SetBranchAddress("jet_truetp_sumpt", &jet_truetp_sumpt, &b_jet_truetp_sumpt);
    fChain->SetBranchAddress("jetPU_truetp_sumpt", &jetPU_truetp_sumpt, &b_jetPU_truetp_sumpt);

   fChain->SetBranchAddress("jet_matchtrk_sumpt", &jet_matchtrk_sumpt, &b_jet_matchtrk_sumpt);
   fChain->SetBranchAddress("jet_loosematchtrk_sumpt", &jet_loosematchtrk_sumpt, &b_jet_loosematchtrk_sumpt);
   fChain->SetBranchAddress("jet_trk_sumpt", &jet_trk_sumpt, &b_jet_trk_sumpt);
    fChain->SetBranchAddress("recojet_eta", &recojet_eta, &b_recojet_eta);
    fChain->SetBranchAddress("recojet_vz", &recojet_vz, &b_recojet_vz);
    fChain->SetBranchAddress("recojet_p", &recojet_p, &b_recojet_p);
    fChain->SetBranchAddress("recojet_pt", &recojet_pt, &b_recojet_pt);
    fChain->SetBranchAddress("recojet_phi", &recojet_phi, &b_recojet_phi);
    fChain->SetBranchAddress("recojet_ntracks", &recojet_ntracks, &b_recojet_ntracks);
    
    fChain->SetBranchAddress("recoClusjet_eta", &recoClusjet_eta, &b_recoClusjet_eta);
    fChain->SetBranchAddress("recoClusjet_vz", &recoClusjet_vz, &b_recoClusjet_vz);
    fChain->SetBranchAddress("recoClusjet_p", &recoClusjet_p, &b_recoClusjet_p);
    fChain->SetBranchAddress("recoClusjet_pt", &recoClusjet_pt, &b_recoClusjet_pt);
    fChain->SetBranchAddress("recoClusjet_phi", &recoClusjet_phi, &b_recoClusjet_phi);
    fChain->SetBranchAddress("recoClusjet_ntracks", &recoClusjet_ntracks, &b_recoClusjet_ntracks);
   Notify();
}

Bool_t PUAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PUAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PUAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PUAnalyzer_cxx
