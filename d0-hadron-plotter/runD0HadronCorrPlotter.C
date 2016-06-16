#include "d0_hadron_corr_plotter.cc"
#include "TLegend.h"
void runD0HadronCorrPlotter()
{
  int inputRebinFactor = 100;

  D0HadronCorrPlotter *plotter = new D0HadronCorrPlotter();
  plotter->init("root-files/dHCorr_P16.root");
  // plotter->init("root-files/corrTest13.root");
  // plotter->plotCorrelation();

  // D0HadronCorrPlotter *plotter_before = new D0HadronCorrPlotter();
  // plotter_before->init("root-files/dHadronCorPhiFixed.root");
  // plotter_before->getCorrelation(inputCentralityBin,inputRebinFactor);
  // plotter_before->plotCorrelation();
  //
  // TH1D *signalBefore = plotter_before->getSignalCorr();
  // TH1D *signal= plotter->getSignalCorr();
  //
  // TCanvas *comparePhiFix = new TCanvas();
  // comparePhiFix->cd();
  // signalBefore->Draw();
  // signal->Draw("same");
  // signalBefore->SetLineColor(1);
  // signal->SetLineColor(2);
  // TLegend *compareLeg = new TLegend(0.1,0.7,0.5,0.9);
  // compareLeg->AddEntry(signalBefore,"before phi fix");
  // compareLeg->AddEntry(signal,"after phi fix");
  // compareLeg->Draw("same");
  TCanvas *c = new TCanvas();
  c->Divide(3,3);
  TH1D *signal[9];
  for(int i=0;i<9;i++)
  {
    pair<int,int> inputCentralityBin(i+1,i+1);
    plotter->getCorrelation(inputCentralityBin,inputRebinFactor);
    signal[i] = plotter->getSignalCorr();
    c->cd(i+1);
    signal[i]->Draw();
  }

}

