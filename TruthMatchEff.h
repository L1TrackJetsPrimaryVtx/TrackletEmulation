/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  3 15:02:42 2017 by ROOT version 5.34/37
// from TTree eventTree/Event tree
// found on file: TTBar_D13_PU200HighStats.root
//////////////////////////////////////////////////////////

#ifndef TruthMatchEff_h
#define TruthMatchEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <TLorentzVector.h>
// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;
class TruthMatchEff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_chi2;
   vector<int>     *trk_nstub;
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
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_nstublayers;
   vector<int>     *tp_eventid;
   vector<int>     *tp_charge;

   vector<float>   *pv_L1;
   vector<float>   *pv_MC;
    vector<int>   *MC_lep;

   vector<float>   *truejet_eta;
   vector<float>   *truejet_vz;
   vector<float>   *truejet_p;
   vector<float>   *truejet_pt;
   vector<float>   *truejet_phi;
   vector<int>     *truejet_ntracks;
   vector<float>   *truejet_tp_sumpt;
   vector<float>   *truejet_truetp_sumpt;
    vector<float>   *recojet_eta;
    vector<float>   *recojet_vz;
    vector<float>   *recojet_p;
    vector<float>   *recojet_pt;
    vector<float>   *recojet_phi;
    vector<int>     *recojet_ntracks;
    
    
   vector<float>   *jet_eta;
   vector<float>   *jet_vz;
   vector<float>   *jet_p;
   vector<float>   *jet_pt;
   vector<float>   *jet_phi;
   vector<int>     *jet_ntracks;
   vector<float>   *jet_tp_sumpt;
   vector<float>   *jet_truetp_sumpt;
    vector<float>   *genjetak4_metfrac;
    vector<float>   *genjetak4_chgfrac;

   vector<float>   *genjetak4_eta;
   vector<float>   *genjetak4_phi;
   vector<float>   *genjetak4_p;
   vector<float>   *genjetak4_pt;

    vector<float>   *genjetchgak4_eta;
    vector<float>   *genjetchgak4_phi;
    vector<float>   *genjetchgak4_p;
    vector<float>   *genjetchgak4_pt;
    
   // List of branches
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_nstub;   //!
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
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_nstublayers;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_charge;   //!

   TBranch        *b_pv_L1;   //!
   TBranch        *b_pv_MC;   //!
    TBranch        *b_MC_lep;   //!

   TBranch        *b_truejet_eta;   //!
   TBranch        *b_truejet_vz;   //!
   TBranch        *b_truejet_p;   //!
   TBranch        *b_truejet_pt;   //!
   TBranch        *b_truejet_phi;   //!
   TBranch        *b_truejet_ntracks;   //!
   TBranch        *b_truejet_tp_sumpt;   //!
   TBranch        *b_truejet_truetp_sumpt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_vz;   //!
   TBranch        *b_jet_p;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_ntracks;   //!
   TBranch        *b_jet_tp_sumpt;   //!
   TBranch        *b_jet_truetp_sumpt;   //!
    TBranch        *b_recojet_eta;   //!
    TBranch        *b_recojet_vz;   //!
    TBranch        *b_recojet_p;   //!
    TBranch        *b_recojet_pt;   //!
    TBranch        *b_recojet_phi;   //!
    TBranch        *b_recojet_ntracks;   //!
    TBranch        *b_genjetak4_metfrac;   //!
    TBranch        *b_genjetak4_chgfrac;   //!

   TBranch        *b_genjetak4_eta;   //!
   TBranch        *b_genjetak4_phi;   //!
   TBranch        *b_genjetak4_p;   //!
   TBranch        *b_genjetak4_pt;   //!

    TBranch        *b_genjetchgak4_eta;   //!
    TBranch        *b_genjetchgak4_phi;   //!
    TBranch        *b_genjetchgak4_p;   //!
    TBranch        *b_genjetchgak4_pt;   //!
   TruthMatchEff(TTree *tree=0);
   virtual ~TruthMatchEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetNevents(){return fChain->GetEntries();};
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TruthMatchEff_cxx
TruthMatchEff::TruthMatchEff(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {					    
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTBar_D13_NoPU.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("TTBar_D13_NoPU.root","READ");
      }
      TDirectory * dir = (TDirectory*)f->Get("TTBar_D13_NoPU.root:/L1TrackNtuple");

      dir->GetObject("eventTree",tree);

   }
   Init(tree);
}

TruthMatchEff::~TruthMatchEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TruthMatchEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TruthMatchEff::LoadTree(Long64_t entry)
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

void TruthMatchEff::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_chi2 = 0;
   trk_nstub = 0;
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
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_pdgid = 0;
   tp_nmatch = 0;
   tp_nstub = 0;
   tp_nstublayers = 0;
   tp_eventid = 0;
   tp_charge = 0;
   pv_L1 = 0;

   pv_MC = 0;
MC_lep = 0;

   truejet_eta = 0;
   truejet_vz = 0;
   truejet_p = 0;
   truejet_pt = 0;
   truejet_phi = 0;
   truejet_ntracks = 0;
   truejet_tp_sumpt = 0;
   truejet_truetp_sumpt = 0;
   jet_eta = 0;
   jet_vz = 0;
   jet_p = 0;
   jet_pt = 0;
   jet_phi = 0;
   jet_ntracks = 0;
    recojet_eta = 0;
    recojet_vz = 0;
    recojet_p = 0;
    recojet_pt = 0;
    recojet_phi = 0;
    recojet_ntracks = 0;
    
   jet_tp_sumpt = 0;
   jet_truetp_sumpt = 0;
    genjetak4_metfrac = 0;
    genjetak4_chgfrac = 0;

   genjetak4_eta = 0;
   genjetak4_phi = 0;
   genjetak4_p = 0;
   genjetak4_pt = 0;
    
    genjetchgak4_eta = 0;
    genjetchgak4_phi = 0;
    genjetchgak4_p = 0;
    genjetchgak4_pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
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
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_nstublayers", &tp_nstublayers, &b_tp_nstublayers);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);

   fChain->SetBranchAddress("pv_L1", &pv_L1, &b_pv_L1);
   fChain->SetBranchAddress("pv_MC", &pv_MC, &b_pv_MC);
    fChain->SetBranchAddress("MC_lep", &MC_lep, &b_MC_lep);

   fChain->SetBranchAddress("truejet_eta", &truejet_eta, &b_truejet_eta);
   fChain->SetBranchAddress("truejet_vz", &truejet_vz, &b_truejet_vz);
   fChain->SetBranchAddress("truejet_p", &truejet_p, &b_truejet_p);
   fChain->SetBranchAddress("truejet_pt", &truejet_pt, &b_truejet_pt);
   fChain->SetBranchAddress("truejet_phi", &truejet_phi, &b_truejet_phi);
   fChain->SetBranchAddress("truejet_ntracks", &truejet_ntracks, &b_truejet_ntracks);
   fChain->SetBranchAddress("truejet_tp_sumpt", &truejet_tp_sumpt, &b_truejet_tp_sumpt);
   fChain->SetBranchAddress("truejet_truetp_sumpt", &truejet_truetp_sumpt, &b_truejet_truetp_sumpt);
    fChain->SetBranchAddress("recojet_eta", &recojet_eta, &b_recojet_eta);
    fChain->SetBranchAddress("recojet_vz", &recojet_vz, &b_recojet_vz);
    fChain->SetBranchAddress("recojet_p", &recojet_p, &b_recojet_p);
    fChain->SetBranchAddress("recojet_pt", &recojet_pt, &b_recojet_pt);
    fChain->SetBranchAddress("recojet_phi", &recojet_phi, &b_recojet_phi);
    fChain->SetBranchAddress("recojet_ntracks", &recojet_ntracks, &b_recojet_ntracks);
    
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_vz", &jet_vz, &b_jet_vz);
   fChain->SetBranchAddress("jet_p", &jet_p, &b_jet_p);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_ntracks", &jet_ntracks, &b_jet_ntracks);
   fChain->SetBranchAddress("jet_tp_sumpt", &jet_tp_sumpt, &b_jet_tp_sumpt);
   fChain->SetBranchAddress("jet_truetp_sumpt", &jet_truetp_sumpt, &b_jet_truetp_sumpt);
   fChain->SetBranchAddress("genjetak4_metfrac", &genjetak4_metfrac, &b_genjetak4_metfrac);
    fChain->SetBranchAddress("genjetak4_chgfrac", &genjetak4_chgfrac, &b_genjetak4_chgfrac);

    fChain->SetBranchAddress("genjetchgak4_eta", &genjetchgak4_eta, &b_genjetchgak4_eta);
    fChain->SetBranchAddress("genjetchgak4_phi", &genjetchgak4_phi, &b_genjetchgak4_phi);
    fChain->SetBranchAddress("genjetchgak4_p", &genjetchgak4_p, &b_genjetchgak4_p);
    fChain->SetBranchAddress("genjetchgak4_pt", &genjetchgak4_pt, &b_genjetchgak4_pt);
   fChain->SetBranchAddress("genjetak4_eta", &genjetak4_eta, &b_genjetak4_eta);
   fChain->SetBranchAddress("genjetak4_phi", &genjetak4_phi, &b_genjetak4_phi);
   fChain->SetBranchAddress("genjetak4_p", &genjetak4_p, &b_genjetak4_p);
   fChain->SetBranchAddress("genjetak4_pt", &genjetak4_pt, &b_genjetak4_pt);
   Notify();
}

Bool_t TruthMatchEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TruthMatchEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TruthMatchEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
    //if(fabs(pv_MC->at(0)-pv_L1->at(0))<1.0)return 1;
    //else return -1;
return 1;
}
#endif // #ifdef TruthMatchEff_cxx
