#define PUAnalyzer_cxx
#include "PUAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include<iostream>


void PUAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PUAnalyzer.C
//      Root > PUAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

    TH1D*hTPFastJetHT=new TH1D("hTPFastJetHT", "Fast Jet H_{T}",100, 0, 2000);
    TH1D*hTPFastJetMHT=new TH1D("hTPFastJetMHT", "Fast Jet H_{T}",100, 0, 2000);
    TH1D*hTPFastSingleJetPt=new TH1D("hTPFastSingleJetPt", "Gen H_{T}",100, 0, 500);
    TH1D*hTPFastDiJetPt=new TH1D("hTPFastDiJetPt", "Gen H_{T}",100, 0, 500);
    TH1D*hTPFastQuadJetPt=new TH1D("hTPFastQuadJetPt", "Gen H_{T}",100, 0, 500);

    TH1D*hFastJetNtracksHT=new TH1D("hFastJetNtracksHT", "Fast Jet H_{T}",100, 0, 2000);

    TH1D*hFastJetHT=new TH1D("hFastJetHT", "Fast Jet H_{T}",100, 0, 2000);
    TH1D*hFastJetMHT=new TH1D("hFastJetMHT", "Fast Jet H_{T}",100, 0, 2000);
    TH1D*hFastJetMHTNtracks=new TH1D("hFastJetMHTNtracks", "Fast Jet H_{T}",100, 0, 2000);

    TH1D*hFastSingleJetPt=new TH1D("hFastSingleJetPt", "Gen H_{T}",100, 0, 500);
    TH1D*hFastDiJetPt=new TH1D("hFastDiJetPt", "Gen H_{T}",100, 0, 500);
    TH1D*hFastQuadJetPt=new TH1D("hFastQuadJetPt", "Gen H_{T}",100, 0, 500);
    
    TH1D*hFastSingleJetPtNtracks=new TH1D("hFastSingleJetPtNtracks", "Gen H_{T}",100, 0, 500);
    TH1D*hFastDiJetPtNtracks=new TH1D("hFastDiJetPtNtracks", "Gen H_{T}",100, 0, 500);
    TH1D*hFastQuadJetPtNtracks=new TH1D("hFastQuadJetPtNtracks", "Gen H_{T}",100, 0, 500);
    
    TH1D*hGenJetMHT=new TH1D("hGenJetMHT", "Gen Missing H_{T}",30, 0, 300);
    TH1D*hGenJetHT=new TH1D("hGenJetHT", "Gen H_{T}",100, 0, 2000);

    TH1D*hGenJetHTPass200=new TH1D("hGenJetHTPass200", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetHTPass325=new TH1D("hGenJetHTPass325", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetHTPass250=new TH1D("hGenJetHTPass250", "Gen H_{T}",100, 0, 2000);
    
    TH1D*hGenJetRecoHTPass200=new TH1D("hGenJetRecoHTPass200", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass325=new TH1D("hGenJetRecoHTPass325", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass250=new TH1D("hGenJetRecoHTPass250", "Gen H_{T}",100, 0, 2000);

    TH1D*hGenJetRecoHTPass200Ntracks=new TH1D("hGenJetRecoHTPass200Ntracks", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass325Ntracks=new TH1D("hGenJetRecoHTPass325Ntracks", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass250Ntracks=new TH1D("hGenJetRecoHTPass250Ntracks", "Gen H_{T}",100, 0, 2000);
    
    TH1D*hGenJetMHTPass25=new TH1D("hGenJetMHTPass25", "Gen  Missing H_{T}",30, 0, 300);
    TH1D*hGenJetMHTPass50=new TH1D("hGenJetMHTPass50", "Gen Missing H_{T}",30, 0, 300);
    TH1D*hGenJetMHTPass100=new TH1D("hGenJetMHTPass100", "Gen Missing H_{T}",30, 0, 300);
    
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       if(MC_lep->at(0)>0)continue;

       float FastJetHT=0;
       float FastJetHTCorr=0;
       float LeadJetPt=0;
       float DiJetPt=0;
       float QuadJetPt=0;
       float LeadJetPtCorr=0;
       float DiJetPtCorr=0;
       float QuadJetPtCorr=0;
       if(jet_pt->size()>0){
           LeadJetPtCorr=jet_pt->at(0)*(jet_pt->at(0)*-0.00614846 + 2.58463);
LeadJetPt=jet_pt->at(0);
       }
       hTPFastSingleJetPt->Fill(LeadJetPt);

       if(jet_pt->size()>1){
         DiJetPt=jet_pt->at(1);
           DiJetPtCorr=jet_pt->at(1)*(jet_pt->at(1)*-0.00614846 + 2.58463);
       }
       hTPFastDiJetPt->Fill(DiJetPt);

       if(jet_pt->size()>3){
         QuadJetPt=jet_pt->at(3);
           QuadJetPtCorr=jet_pt->at(3)*(jet_pt->at(3)*-0.00614846 + 2.58463);
       }
       hTPFastQuadJetPt->Fill(QuadJetPt);
        FastJetHT=0;
        LeadJetPt=0;
        DiJetPt=0;
        QuadJetPt=0;
 
       if(recojet_pt->size()>0)LeadJetPt=recojet_pt->at(0);
       if(recojet_pt->size()>1)DiJetPt=recojet_pt->at(1);
       if(recojet_pt->size()>3)QuadJetPt=recojet_pt->at(3);
       hFastSingleJetPt->Fill(LeadJetPt);
       hFastDiJetPt->Fill(DiJetPt);
       hFastQuadJetPt->Fill(QuadJetPt);

       LeadJetPt=0;
       DiJetPt=0;
       QuadJetPt=0;
       std::vector<float>JetPtNtracks;
       for(unsigned int j=0; j<recojet_pt->size(); ++j){
           if(recojet_ntracks->at(j)>1)JetPtNtracks.push_back(recojet_pt->at(j));
    }
       if(JetPtNtracks.size()>0)LeadJetPt=JetPtNtracks[0];
       if(JetPtNtracks.size()>1)DiJetPt=JetPtNtracks[1];
       if(JetPtNtracks.size()>3)QuadJetPt=JetPtNtracks[3];
       hFastSingleJetPtNtracks->Fill(LeadJetPt);
       hFastDiJetPtNtracks->Fill(DiJetPt);
       hFastQuadJetPtNtracks->Fill(QuadJetPt);
       //std::cout<<"Lead Jet Pt"<<LeadJetPt<<std::endl;
  

       TLorentzVector tempMHT;
       TLorentzVector tempMHTCorr;

       float GenJetHT=0;

       
       
       TLorentzVector tempGenMHT;
       
       for(unsigned int j=0; j<genjetak4_eta->size(); ++j){
           if(genjetak4_pt->at(j)<10)continue;
           if(fabs(genjetak4_eta->at(j))>2.4)continue;
           GenJetHT=GenJetHT+genjetak4_pt->at(j);
           
           TLorentzVector temp;
           temp.SetPtEtaPhiE(genjetak4_pt->at(j), genjetak4_eta->at(j), genjetak4_phi->at(j), genjetak4_p->at(j));
           tempGenMHT=tempGenMHT+temp;
           
       }
       
       for(unsigned int j=0; j<jet_vz->size(); ++j){
           if(jet_pt->at(j)<10)continue;
           if(fabs(jet_eta->at(j))>2.4)continue;
           FastJetHT=FastJetHT+jet_pt->at(j);
           float jetPtCorr=1.59*jet_pt->at(j);
           if(jet_pt->at(j)>10 && jet_pt->at(j)<30)jetPtCorr=jet_pt->at(j)*(jet_pt->at(j)*-0.07307 + 4.17513);
           if(jet_pt->at(j)>30 && jet_pt->at(j)<100)jetPtCorr=jet_pt->at(j)*(jet_pt->at(j)*-0.00614846 + 2.58463);
           FastJetHTCorr=FastJetHTCorr+jetPtCorr;
           TLorentzVector temp;
           temp.SetPtEtaPhiE(jet_pt->at(j),jet_eta->at(j), jet_phi->at(j), jet_p->at(j));
           tempMHT=tempMHT+temp;
           temp.SetPtEtaPhiE(jetPtCorr,jet_eta->at(j), jet_phi->at(j), jet_p->at(j));
           tempMHTCorr=tempMHTCorr+temp;
       }
       hTPFastJetMHT->Fill(tempMHT.Pt());

       hTPFastJetHT->Fill(FastJetHT);
       if(FastJetHT>325)hGenJetHTPass325->Fill(GenJetHT);
       if(FastJetHT>250)hGenJetHTPass250->Fill(GenJetHT);
       if(FastJetHT>200)hGenJetHTPass200->Fill(GenJetHT);
       FastJetHT=0;
       float FastJetNtracks=0;
       int trackcount=0;
       TLorentzVector MHT;
       TLorentzVector MHTNtracks;

       for(unsigned int j=0; j<recojet_vz->size(); ++j){
           if(recojet_pt->at(j)<10)continue;
           if(fabs(recojet_eta->at(j))>2.4)continue;
           //if(trackcount>3)continue;
           //++trackcount;
           TLorentzVector temp;
           temp.SetPtEtaPhiE(recojet_pt->at(j),recojet_eta->at(j), recojet_phi->at(j),recojet_p->at(j));
           MHT=MHT+temp;
           
           FastJetHT=FastJetHT+recojet_pt->at(j);
           if(recojet_ntracks->at(j)>1){
               temp.SetPtEtaPhiE(recojet_pt->at(j),recojet_eta->at(j), recojet_phi->at(j),recojet_p->at(j));
               FastJetNtracks=FastJetNtracks+recojet_pt->at(j);
               MHTNtracks=MHTNtracks+temp;

           }
       }
       hFastJetHT->Fill(FastJetHT);
       hFastJetNtracksHT->Fill(FastJetNtracks);
       hFastJetMHT->Fill(MHT.Pt());
       hFastJetMHTNtracks->Fill(MHTNtracks.Pt());
       if(FastJetHT>360)hGenJetRecoHTPass325->Fill(GenJetHT);
       if(FastJetHT>280)hGenJetRecoHTPass250->Fill(GenJetHT);
       if(FastJetHT>200)hGenJetRecoHTPass200->Fill(GenJetHT);
       
       if(FastJetNtracks>360)hGenJetRecoHTPass325Ntracks->Fill(GenJetHT);
       if(FastJetNtracks>280)hGenJetRecoHTPass250Ntracks->Fill(GenJetHT);
       if(FastJetNtracks>200)hGenJetRecoHTPass200Ntracks->Fill(GenJetHT);
       
       if(tempMHT.Pt()>25)hGenJetMHTPass25->Fill(tempGenMHT.Pt());
       if(tempMHT.Pt()>50)hGenJetMHTPass50->Fill(tempGenMHT.Pt());
       if(tempMHT.Pt()>100)hGenJetMHTPass100->Fill(tempGenMHT.Pt());
       hGenJetMHT->Fill(tempGenMHT.Pt());
       hGenJetHT->Fill(GenJetHT);

   }
    TFile*fout=new TFile("PUJetsOutputMinBias.root", "RECREATE");
    fout->cd();
    
    
    TEfficiency*RecoThresh200=new TEfficiency(*hGenJetRecoHTPass200, *hGenJetHT);
    TEfficiency*RecoThresh250=new TEfficiency(*hGenJetRecoHTPass250, *hGenJetHT);
    TEfficiency*RecoThresh325=new TEfficiency(*hGenJetRecoHTPass325, *hGenJetHT);
    
    TEfficiency*Thresh200=new TEfficiency(*hGenJetHTPass200, *hGenJetHT);
    TEfficiency*Thresh250=new TEfficiency(*hGenJetHTPass250, *hGenJetHT);
    TEfficiency*Thresh325=new TEfficiency(*hGenJetHTPass325, *hGenJetHT);
    TEfficiency*ThreshMHT25=new TEfficiency(*hGenJetMHTPass25, *hGenJetMHT);
    TEfficiency*ThreshMHT50=new TEfficiency(*hGenJetMHTPass50, *hGenJetMHT);
    TEfficiency*ThreshMHT100=new TEfficiency(*hGenJetMHTPass100, *hGenJetMHT);
    ThreshMHT25->Write("Thresh25MHT");
    ThreshMHT50->Write("Thresh50MHT");
    ThreshMHT100->Write("Thresh100MHT");
    Thresh200->Write("TPThresh200HT");
    Thresh250->Write("TPThresh250HT");
    Thresh325->Write("TPThresh325HT");

    RecoThresh200->Write("RecoThresh200HT");
    RecoThresh250->Write("RecoThresh250HT");
    RecoThresh325->Write("RecoThresh325HT");
    
    TEfficiency*RecoThresh200Ntracks=new TEfficiency(*hGenJetRecoHTPass200Ntracks, *hGenJetHT);
    TEfficiency*RecoThresh250Ntracks=new TEfficiency(*hGenJetRecoHTPass250Ntracks, *hGenJetHT);
    TEfficiency*RecoThresh325Ntracks=new TEfficiency(*hGenJetRecoHTPass325Ntracks, *hGenJetHT);
    RecoThresh200Ntracks->Write("RecoThresh200HTNtracks");
    RecoThresh250Ntracks->Write("RecoThresh250HTNtracks");
    RecoThresh325Ntracks->Write("RecoThresh325HTNtracks");

    /*
    Thresh120->Write("Thresh120HT");
    Thresh300->Write("Thresh300HT");
    Thresh400->Write("Thresh400HT");
    */
    
    hTPFastSingleJetPt->Write("SingleJetPt");
    hTPFastDiJetPt->Write("DiJetPt");
    hTPFastQuadJetPt->Write("QuadJetPt");
    hTPFastJetHT->Write("hTPJetHT");
    hTPFastJetMHT->Write("hFastJetMHT");
    hGenJetHT->Write("hGenJetHT");
    hFastJetHT->Write("RecoJetHT");
    hFastJetNtracksHT->Write("RecoJetNtracksHT");
    hFastSingleJetPt->Write("RecoSingleJetPt");
    hFastSingleJetPtNtracks->Write("RecoSingleJetPtNtracks");
    hFastDiJetPt->Write("RecoDiJetPt");
    hFastDiJetPtNtracks->Write("RecoDiJetPtNtracks");
    hFastJetMHT->Write("RecoJetMHT");
    hFastJetMHTNtracks->Write("RecoJetMHTNtracks");

    hFastQuadJetPt->Write("RecoQuadJetPt");
    hFastQuadJetPtNtracks->Write("RecoQuadJetPtNtracks");

    fout->Close();

}
