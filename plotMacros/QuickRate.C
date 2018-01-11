void QuickRate(){
   gROOT->ProcessLine(".L tdrStyle.C");
   setTDRStyle();
   TFile*fin=new TFile("PUJetsOutputMinBiasPt10.root", "READ");
   MinBias=(TH1F*)fin->Get("hTPJetHT");
    MinBiasReco=(TH1F*)fin->Get("RecoJetHT");
    MinBiasRecoNtracks=(TH1F*)fin->Get("RecoJetNtracksHT");

    L1RateMinBias=(TH1F*)MinBias->Clone("L1RateMinBias");
    L1RateMinBiasReco=(TH1F*)MinBiasReco->Clone("L1RateMinBiasReco");
    L1RateMinBiasNtracksReco=(TH1F*)MinBiasReco->Clone("L1RateMinBiasNtracksReco");

    MinBias->Scale(1.0/MinBias->Integral());
    MinBiasReco->Scale(1.0/MinBiasReco->Integral());
    MinBiasRecoNtracks->Scale(1.0/MinBiasRecoNtracks->Integral());

    for(int i=1; i<=100;++i){
       L1RateMinBias->SetBinContent(i, 2760*11.246*MinBias->Integral(i,100)/MinBias->Integral());
        L1RateMinBiasReco->SetBinContent(i, 2760*11.246*MinBiasReco->Integral(i,100)/MinBiasReco->Integral());
        L1RateMinBiasNtracksReco->SetBinContent(i, 2760*11.246*MinBiasRecoNtracks->Integral(i,100)/MinBiasRecoNtracks->Integral());

    }
    TCanvas*c1=new TCanvas("c1", "", 800,800);
    c1->SetLogy();
    L1RateMinBias->SetTitle("MinBias Events PU 200; Track Jet H_{T}; L1 Rate (kHz) ");

   // L1RateMinBias->SetTitle("MinBias Events PU 200; Track Jet H_{T}; L1 Rate (kHz) ");
    L1RateMinBias->Draw();
    L1RateMinBias->SetLineWidth(3.0);L1RateMinBias->SetLineColor(kBlack);
    L1RateMinBias->GetXaxis()->SetRangeUser(0,600);
    TFile*fin2=new TFile("PUJetsOutputMinBiasPt30.root", "READ");
    MinBias30=(TH1F*)fin2->Get("hTPJetHT");
    MinBiasReco30=(TH1F*)fin2->Get("RecoJetHT");
    MinBiasRecoNtracks30=(TH1F*)fin2->Get("RecoJetNtracksHT");

    L1RateMinBias30=(TH1F*)MinBias->Clone("L1RateMinBias30");
    L1RateMinBiasReco30=(TH1F*)MinBiasReco30->Clone("L1RateMinBias");
    L1RateMinBiasRecoNtracks30=(TH1F*)MinBiasRecoNtracks30->Clone("L1RateMinBias");

    for(int i=1; i<=100;++i){
        L1RateMinBias30->SetBinContent(i, 2760*11.246*MinBias30->Integral(i,100)/MinBias30->Integral());
        L1RateMinBiasReco30->SetBinContent(i, 2760*11.246*MinBiasReco30->Integral(i,100)/MinBiasReco30->Integral());
        L1RateMinBiasRecoNtracks30->SetBinContent(i, 2760*11.246*MinBiasRecoNtracks30->Integral(i,100)/MinBiasRecoNtracks30->Integral());


    }
    L1RateMinBias30->SetLineWidth(3.0);L1RateMinBias30->SetLineColor(kBlue);

    L1RateMinBias30->Draw("same");
    //return;
    TCanvas*c2=new TCanvas("c2", "", 800,800);
    c2->SetLogy();
    
    L1RateMinBiasReco->SetTitle("MinBias Events PU 200; Fourth Track Jet p_{T}; L1 Rate (kHz) ");
    L1RateMinBiasReco->SetLineWidth(3.0);L1RateMinBiasReco->SetLineColor(kBlack);

    L1RateMinBiasReco->Draw("");
    L1RateMinBiasReco30->SetLineWidth(3.0);L1RateMinBiasReco30->SetLineColor(kBlue);

    L1RateMinBiasReco30->Draw("same");
    L1RateMinBiasRecoNtracks30->SetLineStyle(kDashed);
    L1RateMinBiasRecoNtracks30->SetLineWidth(3.0);L1RateMinBiasRecoNtracks30->SetLineColor(kBlue);
    L1RateMinBiasNtracksReco->SetLineStyle(kDashed);
    L1RateMinBiasNtracksReco->SetLineWidth(3.0);L1RateMinBiasNtracksReco->SetLineColor(kBlack);
    L1RateMinBiasRecoNtracks30->Draw("same");
    L1RateMinBiasNtracksReco->Draw("same");
    return;
    // L1RateMinBiasNtracksReco->Draw("same");

    TFile*fin3=new TFile("PUJetsOutputTTBarPt10.root", "READ");
    TTBar=(TH1F*)fin3->Get("hTPJetHT");
    TTBarEff=(TH1F*)TTBar->Clone("TTBarEff");
    RecoTTBar=(TH1F*)fin3->Get("RecoJetHT");
    RecoTTBarEff=(TH1F*)TTBar->Clone("RecoTTBarEff");
    RecoTTBarNtracks=(TH1F*)fin3->Get("RecoJetNtracksHT");
    RecoTTBarEffNtracks=(TH1F*)TTBar->Clone("RecoTTBarEffNtracks");
    //TTBarEff->Rebin(2);
    //L1RateMinBias->Rebin(2);
    TFile*fin4=new TFile("PUJetsOutputTTBarPt30.root", "READ");
    TTBar30=(TH1F*)fin4->Get("hTPJetHT");
    TTBarEff30=(TH1F*)TTBar30->Clone("TTBarEff30");
    
    TurnOnHT200TP=(TEfficiency*)fin3->Get("TPThresh200HT");
    TurnOnHT250TP=(TEfficiency*)fin3->Get("TPThresh250HT");
    TurnOnHT325TP=(TEfficiency*)fin3->Get("TPThresh325HT");

    TurnOnHT325TPPt30=(TEfficiency*)fin4->Get("TPThresh325HT");
    TurnOnHT325TPPt30->Draw();
    TurnOnHT325TP->Draw("same");
  //  return;
    for(int i=1; i<=100;++i){
        RecoTTBarEff->SetBinContent(i, 100*RecoTTBar->Integral(i, 100)/RecoTTBar->Integral());

        RecoTTBarEffNtracks->SetBinContent(i, 100*RecoTTBarNtracks->Integral(i, 100)/RecoTTBarNtracks->Integral());
    }

    TGraph*gr3=new TGraph();
    TGraph*gr4=new TGraph();

    for(int i=1; i<=100;++i){
        // gr->SetPoint(gr->GetN(),TTBarEffJet1Pt->GetBinContent(i),L1RateSingleJetMinBias->GetBinContent(i));
        gr3->SetPoint(gr3->GetN(),RecoTTBarEff->GetBinContent(i),L1RateMinBiasReco->GetBinContent(i));
        gr4->SetPoint(gr4->GetN(),RecoTTBarEffNtracks->GetBinContent(i),L1RateMinBiasNtracksReco->GetBinContent(i));
        
    }
    gr3->GetYaxis()->SetRangeUser(0,1000);
    TCanvas*c3=new TCanvas("c3", "", 800,800);
    c3->SetLogy();
    gr3->SetTitle("; Signal Efficiency in TTBar (%) ; L1 Rate (kHz);");
    gr3->Draw("ACP");
    gr4->SetMarkerColor(kRed);
    gr4->SetLineColor(kRed);

    gr4->Draw("CPSAME");
    return;
    
    for(int i=1; i<=100;++i){
        TTBarEff->SetBinContent(i, 100*TTBar->Integral(i, 100)/TTBar->Integral());
        TTBarEff30->SetBinContent(i, 100*TTBar30->Integral(i, 100)/TTBar30->Integral());
        
    }
    TTBarEff->Draw();
    TTBarEff30->Draw("same");
    TGraph*gr=new TGraph();
    TGraph*gr2=new TGraph();
    
    for(int i=1; i<=100;++i){
        // gr->SetPoint(gr->GetN(),TTBarEffJet1Pt->GetBinContent(i),L1RateSingleJetMinBias->GetBinContent(i));
        gr->SetPoint(gr->GetN(),TTBarEff->GetBinContent(i),L1RateMinBias->GetBinContent(i));
        gr2->SetPoint(gr2->GetN(),TTBarEff30->GetBinContent(i),L1RateMinBias30->GetBinContent(i));
        
    }
    gr->GetYaxis()->SetRangeUser(0,1000);
    TCanvas*c3=new TCanvas("c3", "", 800,800);
    c3->SetLogy();
    gr->SetTitle("; Signal Efficiency in TTBar (%) ; L1 Rate (kHz);");
    gr->Draw("ACP");
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineColor(kBlue);
    
    gr2->Draw("CPSAME");
    
    
    TCanvas*c4=new TCanvas("c4", "", 800,800);
    TurnOnHT200TP->SetTitle("; Gen HT (GeV); Efficiency to Pass Track HT Threshold;");
    TurnOnHT200TP->SetLineWidth(2.0); TurnOnHT200TP->SetLineColor(kBlack);
    TurnOnHT200TP->Draw();
    TurnOnHT250TP->SetLineWidth(2.0); TurnOnHT250TP->SetLineColor(kRed);

    TurnOnHT250TP->Draw("same");
    TurnOnHT325TP->SetLineWidth(2.0); TurnOnHT325TP->SetLineColor(kBlue);

    TurnOnHT325TP->Draw("same");


}
