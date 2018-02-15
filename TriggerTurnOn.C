#define TriggerTurnOn_cxx
#include "TriggerTurnOn.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include<iostream>
void TriggerTurnOn::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TriggerTurnOn.C
//      root> TriggerTurnOn t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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
    TH1D*hGenJetHTPass200=new TH1D("hGenJetHTPass200", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetHTPass325=new TH1D("hGenJetHTPass325", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetHTPass250=new TH1D("hGenJetHTPass250", "Gen H_{T}",100, 0, 2000);

    TH1D*hGenJetRecoHTPass200=new TH1D("hGenJetRecoHTPass200", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass325=new TH1D("hGenJetRecoHTPass325", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass250=new TH1D("hGenJetRecoHTPass250", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoMHTPass25=new TH1D("hGenJetRecoMHTPass25", "Gen  Missing H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoMHTPass50=new TH1D("hGenJetRecoMHTPass50", "Gen Missing H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoMHTPass100=new TH1D("hGenJetRecoMHTPass100", "Gen Missing H_{T}",100, 0, 2000);

    TH1D*hGenJetRecoHTPass200Ntracks=new TH1D("hGenJetRecoHTPass200Ntracks", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass325Ntracks=new TH1D("hGenJetRecoHTPass325Ntracks", "Gen H_{T}",100, 0, 2000);
    TH1D*hGenJetRecoHTPass250Ntracks=new TH1D("hGenJetRecoHTPass250Ntracks", "Gen H_{T}",100, 0, 2000);

    TH1D*hGenJetMHTPass25=new TH1D("hGenJetMHTPass25", "Gen  Missing H_{T}",100, 0, 2000);
    TH1D*hGenJetMHTPass50=new TH1D("hGenJetMHTPass50", "Gen Missing H_{T}",100, 0, 2000);
    TH1D*hGenJetMHTPass100=new TH1D("hGenJetMHTPass100", "Gen Missing H_{T}",100, 0, 2000);

    TH1D*hGenJetMHT=new TH1D("hGenJetMHT", "Gen Missing H_{T}",100, 0, 2000);
    TH1D*hGenJetHT=new TH1D("hGenJetHT", "Gen H_{T}",100, 0, 2000);

    TH1D*hGenSingleJetPt=new TH1D("hGenSingleJetPt", "Gen H_{T}",100, 0, 500);
    TH1D*hGenJetSingleJetPass100=new TH1D("hGenJetSingleJetPass100", "Gen Lead Jet p_{T}",100, 0, 2000);
    
    TH1D*hGenDiJetPt=new TH1D("hGenDiJetPt", "Gen H_{T}",100, 0, 500);
    
    TH1D*hGenQuadJetPt=new TH1D("hGenQuadJetPt", "Gen H_{T}",100, 0, 100);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       if(MC_lep->at(0)!=1)continue;
	//TP Reconstructed Quantities 
       float TPJetHT=0;
       float LeadJetPt=0;
       float DiJetPt=0;
       float QuadJetPt=0;
        if(jet_pt->size()>0)LeadJetPt=jet_pt->at(0);			
        if(jet_pt->size()>1)DiJetPt=jet_pt->at(1);			
        if(jet_pt->size()>3)QuadJetPt=jet_pt->at(3);
	       TLorentzVector tempMHT;
        for(unsigned int j=0; j<jet_vz->size(); ++j){
           if(jet_pt->at(j)<10)continue;
           if(fabs(jet_eta->at(j))>2.4)continue;
           TPJetHT=TPJetHT+jet_pt->at(j);
	   TLorentzVector temp;
           temp.SetPtEtaPhiE(jet_pt->at(j),jet_eta->at(j), jet_phi->at(j), jet_p->at(j));
           tempMHT=tempMHT+temp;		
	}
	//Tracklet Reconstructed Quantities 
       float FastJetHT=0;
       float FastJetLeadJetPt=0;
       float FastJetDiJetPt=0;
       float FastJetQuadJetPt=0;
	int TrackCount=0;
	/*
        if(recojet_pt->size()>0)FastJetLeadJetPt=recojet_pt->at(0);			
        if(recojet_pt->size()>1)FastJetDiJetPt=recojet_pt->at(1);			
        if(recojet_pt->size()>3)FastJetQuadJetPt=recojet_pt->at(3);
	*/
	       TLorentzVector tempRecoMHT;
        for(unsigned int j=0; j<recojet_vz->size(); ++j){
           if(recojet_pt->at(j)<10)continue;
           if(fabs(recojet_eta->at(j))>2.4)continue;
	   if(recojet_ntracks->at(j)<2)continue;
           FastJetHT=FastJetHT+recojet_pt->at(j);
	   TLorentzVector temp;
           temp.SetPtEtaPhiE(recojet_pt->at(j),recojet_eta->at(j), recojet_phi->at(j), recojet_p->at(j));
           tempRecoMHT=tempRecoMHT+temp;		
	++TrackCount;
	}
        if(TrackCount>0)FastJetLeadJetPt=recojet_pt->at(0);			
        if(TrackCount>1)FastJetDiJetPt=recojet_pt->at(1);			
        if(TrackCount>3)FastJetQuadJetPt=recojet_pt->at(3);

	//Compute Gen level quanties
	float GenJetHT=0;
               TLorentzVector tempGenMHT;	
	for(unsigned int j=0; j<genjetak4_eta->size(); ++j){
           if(genjetak4_pt->at(j)<30)continue;
           if(fabs(genjetak4_eta->at(j))>2.4)continue;
           GenJetHT=GenJetHT+genjetak4_pt->at(j);
	   TLorentzVector temp;
           temp.SetPtEtaPhiE(genjetak4_pt->at(j), genjetak4_eta->at(j), genjetak4_phi->at(j), genjetak4_p->at(j));
           tempGenMHT=tempGenMHT+temp;	
	
	}
	//now compute eff for turn on
       if(TPJetHT>325)hGenJetHTPass325->Fill(GenJetHT);
       if(TPJetHT>160)hGenJetHTPass250->Fill(GenJetHT);
       if(TPJetHT>200)hGenJetHTPass200->Fill(GenJetHT);			
       if(FastJetHT>200)hGenJetRecoHTPass200->Fill(GenJetHT);
       if(FastJetHT>160)hGenJetRecoHTPass250->Fill(GenJetHT);
       if(FastJetHT>325)hGenJetRecoHTPass325->Fill(GenJetHT);
       if(tempMHT.Pt()>25)hGenJetMHTPass25->Fill(tempGenMHT.Pt());
       if(tempMHT.Pt()>150)hGenJetMHTPass50->Fill(tempGenMHT.Pt());
       if(tempMHT.Pt()>100)hGenJetMHTPass100->Fill(tempGenMHT.Pt());
       if(tempRecoMHT.Pt()>25)hGenJetRecoMHTPass25->Fill(tempGenMHT.Pt());
       if(tempRecoMHT.Pt()>150)hGenJetRecoMHTPass50->Fill(tempGenMHT.Pt());
       if(tempRecoMHT.Pt()>100)hGenJetRecoMHTPass100->Fill(tempGenMHT.Pt());
    	 
	//Denominator for effciency
	hGenJetMHT->Fill(tempGenMHT.Pt());
	hGenJetHT->Fill(GenJetHT); 
	if(genjetak4_eta->size()>0)hGenSingleJetPt->Fill(genjetak4_pt->at(0));
	if(genjetak4_eta->size()>1)hGenDiJetPt->Fill(genjetak4_pt->at(1));
	if(genjetak4_eta->size()>3)hGenQuadJetPt->Fill(genjetak4_pt->at(3));
   }


    TEfficiency*RecoThresh200=new TEfficiency(*hGenJetRecoHTPass200, *hGenJetHT);
    TEfficiency*RecoThresh250=new TEfficiency(*hGenJetRecoHTPass250, *hGenJetHT);
    TEfficiency*RecoThresh325=new TEfficiency(*hGenJetRecoHTPass325, *hGenJetHT);

    TEfficiency*ThreshHT200=new TEfficiency(*hGenJetHTPass200, *hGenJetHT);
    TEfficiency*ThreshHT250=new TEfficiency(*hGenJetHTPass250, *hGenJetHT);
    TEfficiency*ThreshHT325=new TEfficiency(*hGenJetHTPass325, *hGenJetHT);
    TEfficiency*ThreshMHT25=new TEfficiency(*hGenJetMHTPass25, *hGenJetMHT);
    TEfficiency*ThreshMHT50=new TEfficiency(*hGenJetMHTPass50, *hGenJetMHT);
    TEfficiency*ThreshMHT100=new TEfficiency(*hGenJetMHTPass100, *hGenJetMHT);
    TEfficiency*RecoThreshMHT25=new TEfficiency(*hGenJetRecoMHTPass25, *hGenJetMHT);
    TEfficiency*RecoThreshMHT50=new TEfficiency(*hGenJetRecoMHTPass50, *hGenJetMHT);
    TEfficiency*RecoThreshMHT100=new TEfficiency(*hGenJetRecoMHTPass100, *hGenJetMHT);
    TFile*fout=new TFile("TriggerTurnOn.root", "RECREATE");
    fout->cd();
    ThreshMHT25->Write("TPThresh25MHT");
    ThreshMHT50->Write("TPThresh150MHT");
    ThreshMHT100->Write("TPThresh100MHT");
    ThreshHT200->Write("TPThresh200HT");
    ThreshHT250->Write("TPThresh160HT");
    ThreshHT325->Write("TPThresh325HT");

    RecoThresh200->Write("RecoThresh200HT");
    RecoThresh250->Write("RecoThresh160HT");
    RecoThresh325->Write("RecoThresh325HT");
    RecoThreshMHT25->Write("RecoThresh25MHT");
    RecoThreshMHT50->Write("RecoThresh150MHT");
    RecoThreshMHT100->Write("RecoThresh100MHT");
}
