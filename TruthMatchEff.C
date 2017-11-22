#define TruthMatchEff_cxx
#include "TruthMatchEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include<iostream>

void TruthMatchEff::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L TruthMatchEff.C
//      Root > TruthMatchEff t
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

    
    //Gen Level Efficiency Plots
    
    TH1D*TrueChgEnergyFracPt10Num=new  TH1D("TrueChgEnergyFracPt10Num", "", 100, 0.0,1.0);
    TH1D*TrueChgEnergyFracNum=new  TH1D("TrueChgEnergyFracNum", "", 100, 0.0,1.0);
    TH1D*TrueChgEnergyFracDen=new  TH1D("TrueChgEnergyFracDen", "", 100, 0.0,1.0);
    
    TH1D*TrueEtaPt10Num=new  TH1D("TrueEtaPt10Num", "", 50, -2.5, 2.5);
    TH1D*TrueEtaNum=new  TH1D("TrueEtaNum", "", 50, -2.5, 2.5);
    TH1D*TrueEtaDen=new  TH1D("TrueEtaDen", "", 50, -2.5, 2.5);
    
    TH1D*TruePtNumPt10=new  TH1D("TruePtNumPt10", "", 60, 0, 300);

    TH1D*TruePtNum=new  TH1D("TruePtNum", "", 60, 0, 300);
    TH1D*TruePtDen=new  TH1D("TruePtDen", "", 60, 0, 300);
    //Sim Level Efficiency Plots
    TH1D*TrkPartChgEnergyFracPt10Num=new  TH1D("TrkPartChgEnergyFracPt10Num", "", 100, 0.0,1.0);
    TH1D*TrkPartChgEnergyFracNum=new  TH1D("TrkPartChgEnergyFracNum", "", 100, 0.0,1.0);
    TH1D*TrkPartChgEnergyFracDen=new  TH1D("TrkPartChgEnergyFracDen", "", 100, 0.0,1.0);
    TH1D*TrkPartPtPt10Num=new  TH1D("TrkPartPtPt10Num", "", 50, -2.5, 2.5);

    TH1D*TrkPartPtNumPt10=new  TH1D("TrkPartPtNumPt10", "", 60, 0, 300);
    TH1D*TrkPartEtaPt10Num=new  TH1D("TrkPartEtaNumPt10", "", 50, -2.5, 2.5);

    TH1D*TrkPartEtaNum=new  TH1D("TrkPartEtaNum", "", 50, -2.5, 2.5);
    TH1D*TrkPartEtaDen=new  TH1D("TrkPartEtaDen", "", 50, -2.5, 2.5);
    
    TH1D*TrkPartPtNum=new  TH1D("TrkPartPtNum", "", 60, 0, 300);
    TH1D*TrkPartPtDen=new  TH1D("TrkPartPtDen", "", 60, 0, 300);
    
    //RECO Level Efficiency Plots
    TH1D*RecoChgEnergyFracPt10Num=new  TH1D("RecoChgEnergyFracPt10Num", "", 100, 0.0,1.0);
    TH1D*RecoChgEnergyFracNum=new  TH1D("RecoChgEnergyFracNum", "", 100, 0.0,1.0);
    TH1D*RecoChgEnergyFracDen=new  TH1D("RecoChgEnergyFracDen", "", 100, 0.0,1.0);
    TH1D*RecoEtaPt10Num=new  TH1D("RecoEtaPt10Num", "", 50, -2.5, 2.5);
    TH1D*RecoEtaNum=new  TH1D("RecoEtaNum", "", 50, -2.5, 2.5);
    TH1D*RecoEtaDen=new  TH1D("RecoEtaDen", "", 50, -2.5, 2.5);
    TH1D*RecoPtNum=new  TH1D("RecoPtNum", "", 60, 0, 300);
    TH1D*RecoPtDen=new  TH1D("RecoPtDen", "", 60, 0, 300);
    TH1D*RecoPtNtracksNum=new  TH1D("RecoPtNtracksNum", "", 60, 0, 300);
    TH1D*RecoEtaNtracksNum=new  TH1D("RecoEtaNtracksNum", "", 50, -2.5, 2.5);
    TH1D*RecoChgEnergyFracNtracksNum=new  TH1D("RecoChgEnergyFracNtracksNum", "", 100, 0.0,1.0);

    Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
       nbytes += nb;
       if(MC_lep->at(0)>0)continue;
       float genHT=0;
       for (int g=0; g<genjetak4_phi->size(); ++g){
           bool Match=false;
           if(fabs(genjetak4_eta->at(g))>2.2)continue;
           if(fabs(genjetak4_pt->at(g))<30)continue;
           genHT=genHT+genjetak4_pt->at(g);
           int jetindex=-1;
           
           for(unsigned int j=0; j<genjetchgak4_pt->size(); ++j){
               if(genjetchgak4_pt->at(j)<0)continue;
               if(fabs(genjetchgak4_eta->at(j))>2.4)continue;
               float deta=genjetchgak4_eta->at(j)-genjetak4_eta->at(g);
               float dphi=genjetchgak4_phi->at(j)-genjetak4_phi->at(g);
               float dR=sqrt((deta*deta)+(dphi*dphi));
               if(dR<0.4){Match=true;jetindex=j;break;}
           }
           if(Match){
               if( genjetchgak4_pt->at(jetindex)>10){
                   TrueChgEnergyFracPt10Num->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
                   TrueEtaPt10Num->Fill(genjetak4_eta->at(g));
                   TruePtNumPt10->Fill(genjetak4_pt->at(g));
               }
               TrueChgEnergyFracNum->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
               TruePtNum->Fill(genjetak4_pt->at(g));
               TrueEtaNum->Fill(genjetak4_eta->at(g));

           }
           TrueChgEnergyFracDen->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
           TrueEtaDen->Fill(genjetak4_eta->at(g));
           TruePtDen->Fill(genjetak4_pt->at(g));

           Match=false;
           jetindex=-1;
           for(unsigned int j=0; j<jet_pt->size(); ++j){
               // if(fabs(pv_MC->at(0)-jet_vz->at(j))>0.5)continue;
               if(jet_pt->at(j)<0)continue;
               //if(jet_truetp_sumpt->at(j)<0.5)continue;
               if(fabs(jet_eta->at(j))>2.4)continue;
               float deta=jet_eta->at(j)-genjetak4_eta->at(g);
               float dphi=jet_phi->at(j)-genjetak4_phi->at(g);
               float dR=sqrt((deta*deta)+(dphi*dphi));
               if(dR<0.4){Match=true;jetindex=j;break; }
           }
           if(Match){
               if( jet_pt->at(jetindex)>10){
                   TrkPartChgEnergyFracPt10Num->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
                   TrkPartEtaPt10Num->Fill(genjetak4_eta->at(g));
                   TrkPartPtNumPt10->Fill(genjetak4_pt->at(g));
               }
               TrkPartChgEnergyFracNum->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
               TrkPartPtNum->Fill(genjetak4_pt->at(g));
               TrkPartEtaNum->Fill(genjetak4_eta->at(g));
           }
           TrkPartChgEnergyFracDen->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
           TrkPartEtaDen->Fill(genjetak4_eta->at(g));
           TrkPartPtDen->Fill(genjetak4_pt->at(g));
           Match=false;
           jetindex=-1;
           for(int j=0; j<recojet_pt->size(); ++j){
               if(recojet_pt->at(j)<10)continue;
               float deta=recojet_eta->at(j)-genjetak4_eta->at(g);
               float dphi=recojet_phi->at(j)-genjetak4_phi->at(g);
               float dR=sqrt((deta*deta)+(dphi*dphi));
               if(dR<0.4){Match=true; jetindex=j;break;}
           }
           if(Match){
               RecoChgEnergyFracNum->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
               RecoPtNum->Fill(genjetak4_pt->at(g));
               RecoEtaNum->Fill(genjetak4_eta->at(g));
               if(recojet_ntracks->at(jetindex)>1){
                   RecoEtaNtracksNum->Fill(genjetak4_eta->at(g));
                   RecoPtNtracksNum->Fill(genjetak4_pt->at(g));
                   RecoChgEnergyFracNtracksNum->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
               }
		
           }

        //   if(fabs(genjetak4_eta->at(g))<1.0){
           RecoChgEnergyFracDen->Fill(genjetak4_chgfrac->at(g)/genjetak4_p->at(g));
           RecoPtDen->Fill(genjetak4_pt->at(g));
        //   }
           RecoEtaDen->Fill(genjetak4_eta->at(g));
           
       }


   }
    
    TFile*fout=new TFile("FastJetEffPlots.root", "RECREATE");
       //Gen Level
    TEfficiency*MatchEfficiencyChgFrac=new TEfficiency(*TrueChgEnergyFracNum, *TrueChgEnergyFracDen);
    TEfficiency*MatchEfficiencyEta=new TEfficiency(*TrueEtaNum, *TrueEtaDen);
    TEfficiency*MatchEfficiencyPt=new TEfficiency(*TruePtNum, *TruePtDen);
    TEfficiency*MatchEfficiencyChgFracPt10=new TEfficiency(*TrueChgEnergyFracPt10Num, *TrueChgEnergyFracDen);
    TEfficiency*MatchEfficiencyEtaPt10=new TEfficiency(*TrueEtaPt10Num, *TrueEtaDen);
    TEfficiency*MatchEfficiencyPtPt10=new TEfficiency(*TruePtNumPt10, *TruePtDen);

//Sim Level
    TEfficiency*TPMatchEfficiencyChgFrac=new TEfficiency(*TrkPartChgEnergyFracNum, *TrkPartChgEnergyFracDen);
    TEfficiency*TPMatchEfficiencyEta=new TEfficiency(*TrkPartEtaNum, *TrkPartEtaDen);
    TEfficiency*TPMatchEfficiencyPt=new TEfficiency(*TrkPartPtNum, *TrkPartPtDen);
    TEfficiency*TPMatchEfficiencyChgFracPt10=new TEfficiency(*TrkPartChgEnergyFracPt10Num, *TrkPartChgEnergyFracDen);
    TEfficiency*TPMatchEfficiencyEtaPt10=new TEfficiency(*TrkPartEtaPt10Num, *TrkPartEtaDen);
    TEfficiency*TPMatchEfficiencyPtPt10=new TEfficiency(*TrkPartPtNumPt10, *TrkPartPtDen);

    
    TEfficiency*RecoMatchEfficiencyChgFrac=new TEfficiency(*RecoChgEnergyFracNum, *RecoChgEnergyFracDen);

    TEfficiency*RecoMatchEfficiencyEta=new TEfficiency(*RecoEtaNum, *RecoEtaDen);

    TEfficiency*RecoMatchEfficiencyPt=new TEfficiency(*RecoPtNum, *RecoPtDen);

    TEfficiency*RecoMatchEfficiencyPtNtracks=new TEfficiency(*RecoPtNtracksNum, *RecoPtDen);

    TEfficiency*RecoMatchEfficiencyEtaNtracks=new TEfficiency(*RecoEtaNtracksNum, *RecoEtaDen);

    TEfficiency*RecoMatchEfficiencyChgFracNtracks=new TEfficiency(*RecoChgEnergyFracNtracksNum, *RecoChgEnergyFracDen);

    TEfficiency*RecoMatchEfficiencyChgFracPt10=new TEfficiency(*RecoChgEnergyFracPt10Num, *RecoChgEnergyFracDen);
    TEfficiency*RecoMatchEfficiencyEtaPt10=new TEfficiency(*RecoEtaPt10Num, *RecoEtaDen);
    
    MatchEfficiencyChgFrac->Write("GenTruthJetChgE");
    MatchEfficiencyEta->Write("GenTruthEta");
    MatchEfficiencyEtaPt10->Write("GenTruthJetEta10Pt");
    MatchEfficiencyPt->Write("GenTruthPt");
    MatchEfficiencyPtPt10->Write("GenTruthPt10Pt");
    MatchEfficiencyChgFracPt10->Write("GenTruthJetPt10ChgE");
    
    TPMatchEfficiencyChgFrac->Write("SIMTruthJetChgE");
    TPMatchEfficiencyEta->Write("SIMTruthEta");
    TPMatchEfficiencyEtaPt10->Write("SIMTruthJetPt10Eta");
    TPMatchEfficiencyPt->Write("SIMTruthPt");
    TPMatchEfficiencyPtPt10->Write("SIMTruthPt10Pt");
    TPMatchEfficiencyChgFracPt10->Write("SIMTruthJetPt10ChgE");

    RecoMatchEfficiencyChgFrac->Write("RecoTruthJetChgE");
    RecoMatchEfficiencyEta->Write("RecoTruthEta");
    RecoMatchEfficiencyEtaPt10->Write("RecoTruthJetPt10Eta");
    RecoMatchEfficiencyPt->Write("RecoTruthPt");
    RecoMatchEfficiencyPtNtracks->Write("RecoNtracksCut");
    RecoMatchEfficiencyEtaNtracks->Write("RecoNtracksCutEta");
    RecoMatchEfficiencyChgFracNtracks->Write("RecoNtracksChgFrac");
    RecoMatchEfficiencyChgFracPt10->Write("RecoTruthJetPt10ChgE");
    fout->Close();
}
