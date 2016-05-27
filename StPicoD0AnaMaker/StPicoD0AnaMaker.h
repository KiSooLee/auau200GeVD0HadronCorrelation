#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
//
#include "StThreeVectorF.hh"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
//// 

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StHFCuts;
class StPicoPrescales;
class StRefMultCorr;



class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);    
    ofstream fout1;

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();
    ofstream fout;

    bool isGoodPair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const*) const;
    int D0Reco(StThreeVectorF *);
    bool isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isGoodHadron(StPicoTrack const*) const;
    bool  isGoodGlobalHadron(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
    //StPicoDstMaker *
    StPicoDst *picoDst;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TFile* phiFix;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    // static const double     pi  = 3.14159265358979323846;
    TH2D *phiWeight2D;
    TH1D *phiWeight;
    TH1F *vtxz;
    TH1F *phiTrack;
    TH2F *phiRun;
    TH3F *massPt;
    TH3F *massPtPlus;
    TH3F *massPtMinus;
    TH2F *dCorPxPlus;
    TH2F *dCorPxMinus;

    TH2F *corClose;
    TH2F *corFar;
    TH2F *corCloseSB1;
    TH2F *corFarSB1;
    TH2F *corCloseSB2;
    TH2F *corFarSB2;

    // double fitsigma[6];
    // double fitmean[6];
    double mHadronV2[9];
    double efficiency[4][6];
    ClassDef(StPicoD0AnaMaker, 1)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
