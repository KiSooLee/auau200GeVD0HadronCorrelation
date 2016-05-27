#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <set>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
//
#include <stdio.h>
#include <time.h>
#include <algorithm>

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();
  fout.open("check.txt");
  fout1.open("check1.txt");

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  phiFix = new TFile("phiFix.root");
  phiWeight2D = (TH2D *)phiFix->Get("phiWeight")->Clone("phiWeight2D");

  mOutputFile->cd();
  vtxz = new TH1F("vtxz","",100,-10,10);
  phiRun = new TH2F("phiRun","",60018,15107000,15167018,1000,-1.*pi,1.*pi);
  phiTrack = new TH1F("phiTrack","phiTrack",1000,-1.*pi,1.*pi);
  dCorPxPlus = new TH2F("dCorPxPlus","",1000,-100,100,10,0,10);
  dCorPxMinus = new TH2F("dCorPxMinus","",1000,-100,100,10,0,10);
  float xbin[7] = {0,1,2,3,4,5,10};                                     
  float binMass[2001];
  float binPhi[2001];
  float binCent[10];
  for(int i=0;i<2001;i++)
    binPhi[i] = 0.005*i-5;                                              
  for(int i=0;i<2001;i++)                                               
    binMass[i] = 0.01*i; 
  for(int i=0;i<10;i++)
    binCent[i] = 1.0*i;
  
  massPt = new TH3F("massPt","",2000,binMass,6,xbin,9,binCent);
  massPtMinus = new TH3F("massPtMinus","",2000,binMass,6,xbin,9,binCent);
  massPtPlus = new TH3F("massPtPlus","",2000,binMass,6,xbin,9,binCent);
  massPt->Sumw2();
  massPtMinus->Sumw2();
  massPtPlus->Sumw2();

  corClose = new TH2F("corClose","",1000,-1.6,4.8,10,0,10);
  corFar = new TH2F("corFar","",1000,-1.6,4.8,10,0,10);
  corClose->Sumw2();
  corFar->Sumw2();

  corCloseSB1 = new TH2F("corCloseSB1","",1000,-1.6,4.8,10,0,10);
  corFarSB1 = new TH2F("corFarSB1","",1000,-1.6,4.8,10,0,10);
  corCloseSB1->Sumw2();
  corFarSB1->Sumw2();

  corCloseSB1->Sumw2();
  corFarSB1->Sumw2();

  corCloseSB2 = new TH2F("corCloseSB2","",1000,-1.6,4.8,10,0,10);
  corFarSB2 = new TH2F("corFarSB2","",1000,-1.6,4.8,10,0,10);
  corCloseSB2->Sumw2();
  corFarSB2->Sumw2();

  mOutputFile->cd();


  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  fout1.close();
  mOutputFile->cd();
  vtxz->Write();
  phiTrack->Write();
  phiRun->Write();
  dCorPxPlus->Write();
  dCorPxMinus->Write();
  massPt->Write();
  massPtPlus->Write();
  massPtMinus->Write();

  corClose->Write();
  corFar->Write();
  corCloseSB1->Write();
  corFarSB1->Write();
  corCloseSB2->Write();
  corFarSB2->Write();

     
  mOutputFile->Close();
  // phiFix->Close();
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  //StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();


  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  float field = event->bField();
  if(!(isGoodEvent()) || !event->isMinBias())//minBias trigger requires
  {
    LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  if(event) {
    pVtx = event->primaryVertex();
  }
  vtxz->Fill(pVtx.z());
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int centBin = 0;
  int runBin = phiWeight2D->GetXaxis()->FindBin(event->runId());
  phiWeight = (TH1D *)phiWeight2D->ProjectionY("phiWeight",runBin,runBin);
  cout<<"runid = "<<event->runId()<<"\t run bin = "<<runBin<<endl;
  if(centrality>=7) centBin=1;
  else if(centrality>=4)  centBin=2;
  else centBin=3;
  
    
  double reweight = mGRefMultCorrUtil->getWeight();
  double pi = TMath::Pi();
  //double pxCut[9] = {-3,-4,-8,-10,-15,-20,-25,-30,-35};
  //double pxPlusCut[9] = {-2.3,-3.9,-6.1,-9.5,-13.9,-19.3,-25.5,-29.9,-33.1};
  //double pxMinusCut[9] = {-2.3,-3.7,-6.1,-9.5,-13.7,-19.1,-24.7,-28.5,-30.9};10%cut for hadron 0.2
  double pxPlusCut[9] = {-0.7,-1.7,-3.3,-6.1,-9.5,-14.1,-19.5,-23.3,-26.1};
  double pxMinusCut[9] = {-0.9,-1.9,-3.3,-5.9,-9.5,-13.9,-18.7,-21.9,-24.1};//50%cut for hadron 0.2
//  double pxPlusCut[9] = {-0.7,-0.9,-1.1,-1.9,-3.1,-4.9,-7.1,-8.7,-9.9};
 // double pxMinusCut[9] = {-0.9,-0.9,-1.3,-1.9,-3.1,-4.7,-6.7,-8.1,-9.1};//50%cut for hadron 1
  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};

  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* mTrack = picoDst->track(i);
    double hadron_phi = mTrack->pMom().phi();
    int bin_phi = phiWeight->FindBin(hadron_phi);
    double mPhiReweight = phiWeight->GetBinContent(bin_phi);
    if(mTrack->pMom().perp()<0.2) continue;
    // phiRun->Fill(event->runId(),mTrack->pMom().phi());
    phiRun->Fill(event->runId(),mTrack->pMom().phi(),mPhiReweight);
    phiTrack->Fill(mTrack->pMom().phi());
  }
  
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
    if (!isTpcPion(pion)) continue;
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    int charge=0;
    float d0Pt = kp->pt();
    double d0Mass = kp->m();
    double d0Eta = kp->eta();
    if(d0Pt>10) continue;
    if(d0Pt<2)  continue;
    int fitindex = 5;
    if(d0Pt<5)
      fitindex = static_cast<int>(d0Pt);
    int ptIdx = 5;
    if(kp->pt()<5)
      ptIdx = static_cast<int>(kp->pt());
    double mean = fitmean[ptIdx];
    double sigma = fitsigma[ptIdx];

    double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    double mPxPlus = 0;
    double mPxMinus = 0;
    set<unsigned int> dDaughters;
    if(-1!=(charge=isD0Pair(kp))) continue;
    dDaughters.insert(kp->kaonIdx());
    dDaughters.insert(kp->pionIdx());

    for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
    {
      StPicoTrack const* mTrack = picoDst->track(i);
      if(!mTrack) continue;
      if(mTrack->pMom().perp()<0.2) continue;
      //if(mTrack->gPt()<0.2) continue;
      //if(!isGoodGlobalHadron(mTrack)) continue;
      if(!isGoodHadron(mTrack)) continue;
      if(dDaughters.find(i) != dDaughters.end())  continue;
      double hadron_phi = mTrack->pMom().phi();
      int bin_phi = phiWeight->FindBin(hadron_phi);
      double mPhiReweight = phiWeight->GetBinContent(bin_phi);
      cout<<"RUNNUMBER= "<<event->runId()<<"\tphi = "<<hadron_phi<<"\tphi reweight = "<<mPhiReweight<<endl;
      double deltaPhi = fabs(mTrack->pMom().phi()-kp->phi());
     // double deltaPhi = fabs(mTrack->gMom(pVtx,field).phi()-kp->phi());
      if(deltaPhi>pi)  deltaPhi = 2.*pi-deltaPhi;
      if(deltaPhi<0.5*pi) continue;
      if(mTrack->pMom().pseudoRapidity()>0.5) 
        mPxPlus += mTrack->pMom().perp()*cos(deltaPhi)*mPhiReweight;
      if(mTrack->pMom().pseudoRapidity()<-0.5) 
        mPxMinus += mTrack->pMom().perp()*cos(deltaPhi)*mPhiReweight;;
//      if(mTrack->gMom(pVtx,field).pseudoRapidity()>0.5) 
//        mPxPlus += mTrack->gPt()*cos(deltaPhi);
//      if(mTrack->gMom(pVtx,field).pseudoRapidity()<-0.5) 
//        mPxMinus += mTrack->gPt()*cos(deltaPhi);

    }
    //cout<<"PX = "<<mPx<<"\tidx = "<<idx<<"\tk index = "<<kp->kaonIdx()<<"\tp index = "<<kp->pionIdx()<<endl;
    dCorPxPlus->Fill(mPxPlus,centrality,reweight);//reweight*reweight_eff);
    dCorPxMinus->Fill(mPxMinus,centrality,reweight);//reweight*reweight_eff);
    massPt->Fill(kp->m(),d0Pt,centrality,reweight);
    if(mPxPlus<pxPlusCut[centrality]) massPtPlus->Fill(kp->m(),d0Pt,centrality,reweight);
    if(mPxMinus<pxMinusCut[centrality]) massPtMinus->Fill(kp->m(),d0Pt,centrality,reweight);
    if(mPxMinus>pxMinusCut[centrality] && mPxPlus>pxPlusCut[centrality]) continue; 

    for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
    {
      StPicoTrack const* mTrack = picoDst->track(i);
      if(!mTrack) continue;
      if(mTrack->pMom().perp()<1.0) continue;
      //if(mTrack->gPt()<1.0) continue;
      if(dDaughters.find(i) != dDaughters.end())  continue;
      if(!isGoodGlobalHadron(mTrack)) continue;
      double hadron_phi = mTrack->pMom().phi();
      int bin_phi = phiWeight->FindBin(hadron_phi);
      double mPhiReweight = phiWeight->GetBinContent(bin_phi);
      cout<<"phi reweight = "<<mPhiReweight<<endl;
      //if(!isGoodHadron(mTrack)) continue;
      double deltaPhi = (mTrack->pMom().phi()-kp->phi());
      double trackEta = mTrack->pMom().pseudoRapidity();
      //double deltaPhi = (mTrack->gMom(pVtx,field).phi()-kp->phi());
      //double trackEta = mTrack->gMom(pVtx,field).pseudoRapidity();

      if(deltaPhi<-0.5*pi)  deltaPhi += 2*pi;
      if(deltaPhi>1.5*pi)  deltaPhi -= 2*pi;
      if(d0Mass>mean-3*sigma && d0Mass<mean+3*sigma)
      {
        if(mPxPlus<pxPlusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corFar->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corClose->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
        if(mPxMinus<pxMinusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corClose->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corFar->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
      }
      if(d0Mass>mean-9*sigma && d0Mass<mean-4*sigma)
      {
        if(mPxPlus<pxPlusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corFarSB1->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corCloseSB1->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
        if(mPxMinus<pxMinusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corCloseSB1->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corFarSB1->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
      }
      if(d0Mass>mean+4*sigma && d0Mass<mean+9*sigma)
      {
        if(mPxPlus<pxPlusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corFarSB2->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corCloseSB2->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
        if(mPxMinus<pxMinusCut[centrality])
        {
          if(trackEta<0. && trackEta>-0.5)
            corCloseSB2->Fill(deltaPhi,centrality,reweight*mPhiReweight);
          if(trackEta>0. && trackEta<0.5)
            corFarSB2->Fill(deltaPhi,centrality,reweight*mPhiReweight);
        }
      }
    }
  }
  return kStOK;
}
//-----------------------------------------------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  int charge = kaon->charge() * pion->charge();
  bool pairCuts = kp->m()>1.6 && kp->m()<2.1 &&
    charge==-1;

  return (isTpcKaon(kaon,&pVtx) && isTpcPion(pion) && 
      pairCuts);
}


int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/

bool StPicoD0AnaMaker::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (event->triggerWord() & mycuts::triggerWord) &&
    fabs(event->primaryVertex().z()) < mycuts::vz &&
    fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
  //  return event->triggerWord() & mycuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StPhysicalHelixD helix = trk->helix();

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    double ptot = trk->dcaGeometry().momentum().mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }
  return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  return trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
bool StPicoD0AnaMaker::isGoodGlobalHadron(StPicoTrack const * const trk) const
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  float field = event->bField();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  return trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->gMom(pVtx,field).pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
