/***********************************************************************************
 *  *
 *  * d0jet-corr-plotter
 *  *
 *  * Author: Leon He
 *  ***********************************************************************************
 *  *
 *  * Description: 
 *  *
 *  ***********************************************************************************
 *  *
 *  * Log:
 *  *
 *  ***********************************************************************************/
#include "fstream"
#include "iostream"
#include "TSTring.h"
class TFile;
class TH2F;

class D0HadronCorrPlotter
{
  public:
    D0HadronCorrPlotter() {}
    ~D0HadronCorrPlotter() {}
    void init(TString inputFileName);
    void finish();
    void getPxCut(double);
    void getCorrelation(std::pair<int,int> &, int);
    void plotSignificance(TH2D *);
    void plotCorrelation();
    double getSBRatio(TH1D *massHisto);
    TH1D *getCandCorr() {return candCorrelation;}
    TH1D *getBkgCorr() {return bkgCorrelation;}
    TH1D *getSignalCorr() {return signalCorrelation;}


  private:
    void fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3]);
    TFile *inputFile;
    ofstream mLog;
    TH1D *signalCorrelation;
    TH1D *candCorrelation;
    TH1D *bkgCorrelation;
    TH1D *candCorrelationClose;
    TH1D *bkgCorrelationClose;
    TH1D *candCorrelationFar;
    TH1D *bkgCorrelationFar;
};
