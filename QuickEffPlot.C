void QuickEffPlot(){
gROOT->ProcessLine(".L tdrStyle.C");
setTDRStyle();
TFile*fin=new TFile("FastJetEffQCDPtHat30Plots.root", "READ");
TEfficiency*EffGenChgE=(TEfficiency*)fin->Get("GenTruthJetChgE");
TEfficiency*EffGenChgEPt10=(TEfficiency*)fin->Get("GenTruthJetPt10ChgE");
TEfficiency*EffSimChgEPt10=(TEfficiency*)fin->Get("SIMTruthJetPt10ChgE");
TEfficiency*EffRecoChgEPt10=(TEfficiency*)fin->Get("RecoTruthJetChgE");
    TEfficiency*EffRecoChgENtracks=(TEfficiency*)fin->Get("RecoNtracksChgFrac");


TEfficiency*EffGenPt=(TEfficiency*)fin->Get("GenTruthPt");
TEfficiency*EffGenPt10=(TEfficiency*)fin->Get("GenTruthPt10Pt");
TEfficiency*EffSimPt10=(TEfficiency*)fin->Get("SIMTruthPt10Pt");
TEfficiency*EffRecoPt10=(TEfficiency*)fin->Get("RecoTruthPt");
    TEfficiency*EffRecoPt10=(TEfficiency*)fin->Get("RecoTruthPt");
    TEfficiency*EffRecoPtNtracks=(TEfficiency*)fin->Get("RecoNtracksCut");

    
    TEfficiency*EffGenEta=(TEfficiency*)fin->Get("GenTruthEta");
    TEfficiency*EffGenEta10=(TEfficiency*)fin->Get("GenTruthJetEta10Pt");
    TEfficiency*EffSimEta10=(TEfficiency*)fin->Get("SIMTruthJetPt10Eta");
    
    TEfficiency*EffRecoEta10=(TEfficiency*)fin->Get("RecoTruthEta");
    TEfficiency*EffRecoEtaNtracks=(TEfficiency*)fin->Get("RecoNtracksCutEta");
    //    TFile*fin2=new TFile("em_outQCDpt30_tp_16z.txt.root", "READ");
    //TEfficiency*EffRecoL1L2Pt=(TEfficiency*)fin2->Get("MatchingEfficiencyPt");
    //TEfficiency*EffRecoL1L2Eta=(TEfficiency*)fin2->Get("MatchingEfficiencyEta");
    
    TFile*fin3=new TFile("em_outQCDpt30Ntracks1_tp_16z.txt.root", "READ");
    TEfficiency*EffNtracks1RecoL1L2Pt=(TEfficiency*)fin3->Get("MatchingEfficiencyPt");
    TEfficiency*EffNtracks1RecoL1L2Eta=(TEfficiency*)fin3->Get("MatchingEfficiencyEta");
    fin->cd();
    
TCanvas*c1=new TCanvas("c1", "", 800, 800);
EffGenChgE->SetTitle("Track Jet Finding Efficiency; Charged Energy Fraction; Efficiency for Gen Jets (pT>30GeV)");
    
    leg=new TLegend(0.5363409,0.1587097,0.952381,0.3509677,NULL,"brNDC");
    leg->AddEntry(EffGenChgE, "Chg. Part. Gen Jet pT>0", "L");
    leg->AddEntry(EffGenChgEPt10, "Chg. Part. Gen Jet pT>10", "L");
    leg->AddEntry(SIMTruthJetPt10ChgE, "Trk. Part. Jets pT>10", "L");
    leg->AddEntry(EffRecoChgEPt10, "Fastjet Tracklet Jets pT>10 ", "L");
    leg->AddEntry(EffRecoChgENtracks, "Fastjet Tracklet Jets pT>10 N_{tracks}>1 ", "L");
    leg->AddEntry(EffNtracks1RecoL1L2Pt, "L1L2 Tracklet Jets pT>10 N_{tracks}>1 ", "L");


EffGenChgE->Draw();
    gPad->Update();
    TGraph* graph = (TGraph*)EffGenChgE->GetPaintedGraph();
    graph->GetYaxis()->SetTitleSize(0.05);

    graph->GetYaxis()->SetRangeUser(0.5, 1.05);
    graph->Draw("APE");
    EffGenChgE->Draw("same");


    EffGenChgE->SetLineWidth(2.0); EffGenChgE->SetLineColor(kBlack);
EffGenChgEPt10->Draw("same");
    EffGenChgEPt10->SetLineWidth(2.0); EffGenChgEPt10->SetLineColor(kRed);

SIMTruthJetPt10ChgE->Draw("same");
    SIMTruthJetPt10ChgE->SetLineWidth(2.0); SIMTruthJetPt10ChgE->SetLineColor(kMagenta);
    EffRecoChgEPt10->SetLineWidth(2.0);EffRecoChgEPt10->SetLineColor(kBlue);
EffRecoChgEPt10->Draw("same");
    EffRecoChgENtracks->SetLineWidth(2.0); EffRecoChgENtracks->SetLineColor(kCyan);
    EffRecoChgENtracks->Draw("same");


    leg->Draw();
    c1->Print("QCDPtHat30ChgE.pdf");
TCanvas*c2=new TCanvas("c2", "", 800, 800);

EffGenPt->SetTitle("Track Jet Finding Efficiency; Gen Jet p_{T} (GeV); Efficiency for Gen Jets (pT>30GeV)");
EffGenPt->Draw();
    EffGenPt->SetLineWidth(2.0); EffGenPt->SetLineColor(kBlack);
    gPad->Update();
    TGraph* graph = (TGraph*)EffGenPt->GetPaintedGraph();
    graph->GetYaxis()->SetTitleSize(0.05);

    graph->GetYaxis()->SetRangeUser(0.5, 1.05);
    graph->Draw("APE");
    EffGenPt->Draw("same");
EffGenPt10->Draw("same");
    EffGenPt10->SetLineWidth(2.0); EffGenPt10->SetLineColor(kRed);

SIMTruthPt10Pt->Draw("same");
    SIMTruthPt10Pt->SetLineWidth(2.0); SIMTruthPt10Pt->SetLineColor(kMagenta);

EffRecoPt10->Draw("same");
    EffRecoPt10->SetLineWidth(2.0); EffRecoPt10->SetLineColor(kBlue);
    EffRecoPtNtracks->SetLineWidth(2.0); EffRecoPtNtracks->SetLineColor(kCyan);
    EffRecoPtNtracks->Draw("same");
  //  fin2->cd();
   // EffRecoL1L2Pt->SetLineWidth(2.0);EffRecoL1L2Pt->SetLineColor(kGreen+2);
    //EffRecoL1L2Pt->Draw("same");
    leg->Draw();
    fin3->cd();
    EffNtracks1RecoL1L2Pt->SetLineWidth(2.0);EffNtracks1RecoL1L2Pt->SetLineColor(kGreen+2);
    EffNtracks1RecoL1L2Pt->Draw("same");
    c2->Print("QCDPtHat30Pt.pdf");

    fin->cd();
TCanvas*c3=new TCanvas("c3", "", 800, 800);

    EffGenEta->SetTitle("Track Jet Finding Efficiency; Gen Jet #eta ; Efficiency for Gen Jets (pT>30GeV)");
    EffGenEta->Draw();
    EffGenEta->SetLineWidth(2.0); EffGenEta->SetLineColor(kBlack);
    gPad->Update();
    TGraph* graph = (TGraph*)EffGenEta->GetPaintedGraph();
    graph->GetYaxis()->SetTitleSize(0.05);
    graph->GetYaxis()->SetRangeUser(0.5, 1.05);
    graph->Draw("APE");
    EffGenEta->Draw("same");
    EffGenEta10->SetLineWidth(2.0); EffGenEta10->SetLineColor(kRed);

    EffGenEta10->Draw("same");
    EffSimEta10->SetLineWidth(2.0); EffSimEta10->SetLineColor(kMagenta);

    EffSimEta10->Draw("same");
    EffRecoEta10->SetLineWidth(2.0); EffRecoEta10->SetLineColor(kBlue);

    EffRecoEta10->Draw("same");
    EffRecoEtaNtracks->SetLineWidth(2.0); EffRecoEtaNtracks->SetLineColor(kCyan);
    EffRecoEtaNtracks->Draw("same");
   // fin2->cd();
   // EffRecoL1L2Eta->SetLineWidth(2.0);EffRecoL1L2Eta->SetLineColor(kGreen+2);
   // EffRecoL1L2Eta->Draw("same");

    leg->Draw();
    EffNtracks1RecoL1L2Eta->SetLineWidth(2.0);EffNtracks1RecoL1L2Eta->SetLineColor(kGreen+2);
    EffNtracks1RecoL1L2Eta->Draw("same");
    c3->Print("QCDPtHat30Eta.pdf");
    

    
    
}
