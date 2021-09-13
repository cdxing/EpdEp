//this code is for:
//conducting the v1 analysis, calculate v1_TPC relative to Psi_1^{EPD_full},
//calculate v1_EPD reative to Psi_1^{TPC_|eta|<0.8}. 11/18/2019 Xiaoyu
//Add one more dimension to the analysis, e.g low pT and high pT.11/21/2019
//I want to use the full statistics for the analysis, need a little bit modifications(because mNphibin=200 not 100):
//mTPCPhiWeightInput->RebinX();mTPCPhiAveragedInput->RebinX();
//mTPCPhiWeightInput->Divide(mTPCPhiAveragedInput);12/02/2019
//add two more "TrackTypes": pos/neg cahrged particles for v1TPC. 12/15/2019

#include "PicoAnalyzer.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBbcHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StThreeVectorF.hh"

#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StBbcGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"


#include "TChain.h"
//#include "TTree.h"
//#include "TLeaf.h"
#include "TMath.h"
#include "StBTofUtil/tofPathLength.hh"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>

#include "../run/badrun.h"

using namespace std;

ClassImp(PicoAnalyzer)                     //  Macro for CINT compatability
int Centrality(int gRefMult );

//=================================================
PicoAnalyzer::PicoAnalyzer(TString FileNameBase):mpTMin(0.15),mpTMax(2.0),mpTMink(0.2),mpTMaxk(2.0),mEtaMin(-1.2),mEtaMax(1.2),
mNPhibin(100),mVzMin(-70.0),mVzMax(70.0),
mNpTbin(20),mNVzbin(4),mNEtabin(60),
mVtxR(2.0),mDiffVzVPD(3.0),mNhitsfit(15),mNhitsfitratio(0.52),mDCAcut(3.0),mFourierOrder(8),
mNTPCSubEvents(2),mEPDMax(2.0),mNEPDSubEvents(10),mEPDthresh(0.3),mpTbound(0.425),mPionSigma(0.012),mKaonSigma(2.0),
mProtonSigma(0.012),d_KaonM2low(0.16),d_KaonM2high(0.36),mEtaMaxv1(2.0),mEtaMinv1(-2.0),dip_angle_cutLevel(0.04){
  mFileNameBase = FileNameBase;

  mPicoDst=0;
  mEpdHits=0;
  mBbcHits=0;
  mTracks=0;
  mEventClonesArray=0;

  mRunId=0;
  mRunEt=0;
  mRunCollisionSystem=0;

  mEpdGeom = new StEpdGeom;
  mBbcGeom = new StBbcGeom;
  mRan = new TRandom3;
  mRan->GetSeed();
}

//=================================================
PicoAnalyzer::~PicoAnalyzer(){
  /* no-op */
}


//=================================================
//void PicoAnalyzer::SetPicoDst(TTree* PicoDst){
void PicoAnalyzer::SetPicoDst(TChain* PicoDst){
  mPicoDst        = PicoDst;

  mEpdHits = new TClonesArray("StPicoEpdHit");
  mBbcHits = new TClonesArray("StPicoBbcHit");
  mTracks  = new TClonesArray("StPicoTrack");
  mEventClonesArray = new TClonesArray("StPicoEvent");
  mTraits = new TClonesArray("StPicoBTofPidTraits");

  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  unsigned int found;
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  cout << "EpdHit Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);
  mPicoDst->SetBranchStatus("Event*",1,&found);
  cout << "Event Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Event",&mEventClonesArray);


  mPicoDst->SetBranchStatus("Track*",1,&found);
  cout << "Track Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Track",&mTracks);

  mPicoDst->SetBranchStatus("BTofPidTraits*",1,&found);
  cout << "BTofPidTraits Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("BTofPidTraits",&mTraits);
}

//=================================================
short PicoAnalyzer::Init(char const* TPCWeightFile, char const* TPCShiftFile, char const* EPDPhiWeightFile){

  // ------------------- for the EP finder ------------------
  TString EpFinderOutputName = mFileNameBase;
  EpFinderOutputName += "EpFinderCorrectionsOUTPUT.root";

  mEpFinder = new StEpdEpFinder(9,EpFinderOutputName.Data(),"EPDcorrection.root");// EventType:centrality, read the INPUT file: EPDcorrection.root
  mEpFinder->SetnMipThreshold(0.3);
  mEpFinder->SetMaxTileWeight(2.0);
  mEpFinder->SetEpdHitFormat(2);     // 2=pico

  double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
  TH2D wt("Order1etaWeight","Order1etaWeight",100,1.5,6.5,9,0,9);
  for (int ix=1; ix<101; ix++){
    for (int iy=1; iy<10; iy++){
      double eta = wt.GetXaxis()->GetBinCenter(ix);
      wt.SetBinContent(ix,iy,lin[iy-1]*eta+cub[iy-1]*pow(eta,3));
    }
  }
  // mEpFinder->SetEtaWeights(1,wt);
  // cout<<"etaweight set"<<endl;
  cout<<"etaweight disabled"<<endl;
  // --------------------------------------------------------

  TString OutputRootFileName = mFileNameBase;
  OutputRootFileName += "_FlowHists.root";

  mHistoFile = new TFile(OutputRootFileName,"RECREATE");

  // Miscellaneous One-dimensional histograms
  //mHisto1D[0] = new TH1D("Vz","Vz",100,-80,80);
  mHisto1D[1] = new TH1D("RefMult","RefMult",100,-10,600);
  //mHisto1D[2] = new TH1D("dNdeta","dNdeta",100,-5.5,5.5);
  for(int cent=0;cent<9;cent++){
    mVz[cent]= new TH1D(Form("VzCent%d",cent),Form("VzCent%d",cent),16,-40,40);
  }

//--------------Prepare the constant for weighting the TPC tracks------------
  mNumberOfTrackTypes =mNpTbin*mNVzbin*mNEtabin*2;


//---------------- Make histograms for phi meson analysis --------------
  h2px = new TH2D("h2px","pMomX vs. P_vecX",600,-3.,3.,600,-3.,3.);
  h2py = new TH2D("h2py","pMomY vs. P_vecY",600,-3.,3.,600,-3.,3.);
  h2pz = new TH2D("h2pz","pMomZ vs. P_vecZ",600,-3.,3.,600,-3.,3.);
  h_dip_angle = new TH1D("h_dip_angle","h_dip_angle",1000,-1.,9.);
  h_Mass  = new TH1D("h_Mass","Same event invariant mass",200,0.98,1.08);
  h_Mass_rot  = new TH1D("h_Mass_rot","K+K- rotated invariant mass",200,0.98,1.08);
  hist_SE_PhiMeson_pT  = new TH1D("hist_SE_PhiMeson_pT","pT distribution of #phi",200,0.0,10);
  hist_SE_PhiMeson_mT  = new TH1D("hist_SE_PhiMeson_mT","mT distribution of #phi",200,0.0,10);
  hist_SE_PhiMeson_rap  = new TH1D("hist_SE_PhiMeson_rap","y distribution of #phi",200,-10.,10);
  hist_SE_PhiMeson_eta  = new TH1D("hist_SE_PhiMeson_eta","eta distribution of #phi",200,-10.,10);
  TString HistName = "Mass2_pt", HistName_rot = "Mass2_rot_pt";
  TString HistName_ptEta = "pT_eta", HistName_ptY = "pT_y";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.,5.0,200,0.98,1.08);
  h_Mass2_rot = new TH2F(HistName_rot.Data(),HistName_rot.Data(),20,0.,5.0,200,0.98,1.08);
  h2_pT_eta = new TH2F(HistName_ptEta.Data(),HistName_ptEta.Data(),200,-2.0,2.0,20,0.,5.0);
  h2_pT_y = new TH2F(HistName_ptY.Data(),HistName_ptY.Data(),200,-2.0,2.0,20,0.,5.0);


  hist_SE_pt_y_PhiMeson[0] = new TH2D("hist_SE_pt_y_PhiMeson_0","p_{T} [GeV/c] vs. y of #phi, 0-60% ",200,-2.0,2.0,20,0.,5.0);
  hist_SE_pt_y_Phi_tight_SigBkg[0] = new TH2D("hist_SE_pt_y_Phi_tight_SigBkg_0","p_{T} [GeV/c] vs. y of #phi, 0-60% ",40,-2.,2.,35,0.0,3.5);
  hist_SE_pt_y_Phi_tight_Bkg[0] = new TH2D("hist_SE_pt_y_Phi_tight_Bkg_0","p_{T} [GeV/c] vs. y of #phi^{Bkg}, 0-60% ",40,-2.,2.,35,0.0,3.5);
  hist_SE_pt_y_Phi_tight_Sig[0] = (TH2D*) hist_SE_pt_y_Phi_tight_SigBkg[0]->Clone("hist_SE_pt_y_Phi_tight_Sig_0");
  int centBES[4] = {0,10,40,80};
  for(int cent = 1; cent<4;cent++){
    hist_SE_pt_y_PhiMeson[cent] = new TH2D(Form("hist_SE_pt_y_PhiMeson_%d",cent),Form("p_{T} [GeV/c] vs. y of #phi, %d-%d%%",centBES[cent-1],centBES[cent]),200,-2.0,2.0,20,0.,5.0);
    hist_SE_pt_y_Phi_tight_SigBkg[cent] = new TH2D(Form("hist_SE_pt_y_Phi_tight_SigBkg_%d",cent),Form("p_{T} [GeV/c] vs. y of #phi, %d-%d%%",centBES[cent-1],centBES[cent]),40,-2.,2.,35,0.0,3.5);
    hist_SE_pt_y_Phi_tight_Bkg[cent] = new TH2D(Form("hist_SE_pt_y_Phi_tight_Bkg_%d",cent),Form("p_{T} [GeV/c] vs. y of #phi^{Bkg}, %d-%d%%",centBES[cent-1],centBES[cent]),40,-2.,2.,35,0.0,3.5);
    hist_SE_pt_y_Phi_tight_Sig[cent] = (TH2D*) hist_SE_pt_y_Phi_tight_SigBkg[cent]->Clone(Form("hist_SE_pt_y_Phi_tight_Sig_%d",cent));
  }
  TString Centrality_01[4] = {"0080","0010","1040","4080"};

  for(Int_t cent = 0; cent < Bin_Centrality_01; cent++)
  {
      for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++)
      {
        TString hist_name_SE = Form("InvMass_SE_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data());
        mHist_SE_InvM_rap_cent[rap_bin][cent] = new TH1F(hist_name_SE.Data() ,
        hist_name_SE.Data() ,
        200,0.98,1.08);
        mHist_SE_InvM_rap_cent[rap_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
        TString hist_name_rot = Form("InvMass_rot_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data());
        mHist_rotation_InvM_rap_cent[rap_bin][cent] = new TH1F(hist_name_rot.Data() ,
        hist_name_rot.Data() ,
        200,0.98,1.08);
        mHist_rotation_InvM_rap_cent[rap_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
        TString hist_name_profile = Form("flow_InvMass_rapbin%d_cent%s",rap_bin+1,Centrality_01[cent].Data());
        mProfile_flow_reso_rap_cent[rap_bin][cent] = new TProfile(hist_name_profile.Data(),
        hist_name_profile.Data(),
        100,0.98,1.08,
        0,0,"");
        mProfile_flow_reso_rap_cent[rap_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
        mProfile_flow_reso_rap_cent[rap_bin][cent]->GetYaxis()->SetTitle("<cos(1(#phi - #psi_{1}))>/R_{1}^{EPD}");
      }
  }

//----------------Make histograms for QA ----------------------------------
  href_vz = new TH1F("h_ref_vz","refmult_vz",1000,0.,1000.);
  hvz_b = new TH1F("h_vz_b","vz_dis_b",1000,-150,150);
  hvzvpdvz_b =new TH2F("h_vz_vpd_b","vz_vs_vpd_b",1000,-150,150,1000,-150,150);
  hvzvpdvzdiff_b =new TH1F("hvzvpdvzdiff_b","hvzvpdvzdiff_b",1000,-150,150);
  hvr_b = new TH2F("h_vr_b","vy_vs_vx_b",1000,-10,10,1000,-10,10);
  htofvsref_b = new TH2F("htofvsref_b","ref_vs_tof",2000,0.0,2000.0,2000,0.0,2000.0);
  htofmatchvsref_b=new TH2F("htofmatchvsref_b","",1500,0.0,1500.0,1500,0.0,1500.0);

  href = new TH1F("h_ref","refmult_dis",1000,0.,1000.);
  hvz = new TH1F("h_vz","vz_dis",1000,-100,100);
  hvzvpdvz =new TH2F("h_vz_vpd","vz_vs_vpd",1000,-150,150,1000,-150,150);
  hvzvpdvzdiff =new TH1F("hvzvpdvzdiff","hvzvpdvzdiff",1000,-150,150);
  htofvsref = new TH2F("htofvsref","ref_vs_tof",2000,0.0,2000.0,2000,0.0,2000.0);
  htofmatchvsref=new TH2F("htofmatchvsref","",1500,0.0,1500.0,1500,0.0,1500.0);
  hvr = new TH2F("h_vr","vy_vs_vx",1000,-10,10,1000,-10,10);
  hbtofYLocal = new TH1F("hbtofYLocal","",600,-6.,6.);
  hbtofYLocalvsMass2 = new TH2F("hbtofYLocalvsMass2","",800,-0.1,1.5,600,-6.,6.);

  hbetavsp =new TH2F("hbetavsp","beta_vs_p",1000,0,6,1000,0,10);
  hmassvsp =new TH2F("hmassvsp","m2_vs_p*q",3000,-6.,6.,2000,-0.2,15.8);
  hdedxvsp =new TH2F("hdedxvsp","dedx_vs_p*q",4000,-6.,6.,2000,0.,50.);

  h_eta_phi =new TH2F("h_etaphi","eta_vs_phi",1000,-6.3,6.3,300,-1.5,1.5);
  h_eta_phi_before =new TH2F("h_etaphi_before","eta_vs_phi w/o cut",1000,-6.3,6.3,300,-1.8,1.8);

  h_counter = new TH1F("h_counter","event_counter",58,-8,50.);

  h_runidvstofmult_b = new TProfile("runidvstofmult_b", "", 90000, 22031041, 22121041,"");
  h_runidvsrefmult_b = new TProfile("runidvsrefmult_b", "", 90000, 22031041, 22121041,"");

  h_runidvstofmult = new TProfile("runidvstofmult", "", 90000, 22031041, 22121041,"");
  h_runidvsrefmult = new TProfile("runidvsrefmult", "", 90000, 22031041, 22121041,"");

  h_pt = new TH1F("h_pt","",1000,0.,10.);
  h_eta_b= new TH1F("h_eta_b","",1000,-2.,2.);
  h_eta= new TH1F("h_eta","",1000,-2.,2.);
  h_nhitfit=new TH1F("h_nhitfit","",80,-0.5,79.5);
  h_nhitmax=new TH1F("h_nhitmax","",80,-0.5,79.5);
  h_nhitratio=new TH1F("h_nhitratio","",1000,-0.5,1.5);
  h_dca = new TH1F("h_dca","",1000,0.,5.);
  h_phi = new TH1F("h_phi","",1000,-6.28,6.28);
  //========================= Kaon PID ==========================================
  hist_pt_kaonPlus = new TH1D("hist_pt_kaonPlus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_eta_kaonPlus = new TH1D("hist_eta_kaonPlus","#eta",500,-1.5,1.5);
  hist_y_kaonPlus = new TH1D("hist_y_kaonPlus","y",500,-1.5,1.5);
  hist_phi_kaonPlus = new TH1D("hist_phi_kaonPlus","#phi [Radian]",1000,-1.5*TMath::Pi(),2.5*TMath::Pi());
  hist_rap_eta_kaonPlus = new TH2D("hist_rap_eta_kaonPlus","kaonPlus y versus #eta",500,-1.5,1.5,500,-1.5,1.5);
  hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y",200,-2.0,2.0,20,0.,5.0);
  hist_pt_eta_kaonPlus = new TH2D("hist_pt_eta_kaonPlus","p_{T} [GeV/c] vs. #eta",200,-2.0,2.0,20,0.,5.0);
  hist_dEdx_kaonPlus = new TH2D("hist_dEdx_kaonPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_beta_kaonPlus = new TH2D("hist_beta_kaonPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_mass_kaonPlus = new TH2D("hist_mass_kaonPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_pt_kaonMinus = new TH1D("hist_pt_kaonMinus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_eta_kaonMinus = new TH1D("hist_eta_kaonMinus","#eta",500,-1.5,1.5);
  hist_y_kaonMinus = new TH1D("hist_y_kaonMinus","y",500,-1.5,1.5);
  hist_phi_kaonMinus = new TH1D("hist_phi_kaonMinus","#phi [Radian]",1000,-1.5*TMath::Pi(),2.5*TMath::Pi());
  hist_rap_eta_kaonMinus = new TH2D("hist_rap_eta_kaonMinus","kaonMinus y versus #eta",500,-1.5,1.5,500,-1.5,1.5);
  hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y",200,-2.0,2.0,20,0.,5.0);
  hist_pt_eta_kaonMinus = new TH2D("hist_pt_eta_kaonMinus","p_{T} [GeV/c] vs. #eta",200,-2.0,2.0,20,0.,5.0);
  hist_dEdx_kaonMinus = new TH2D("hist_dEdx_kaonMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_beta_kaonMinus = new TH2D("hist_beta_kaonMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_mass_kaonMinus = new TH2D("hist_mass_kaonMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
//----------------Make histograms for shifting Psi_TPC--------------------
  mTPCCosShift = new TProfile3D("TPCCosShift","TPCCosShift",9,-0.5,8.5,mFourierOrder,0.5,0.5+mFourierOrder,_PsiOrderMax,0.5,(double)_PsiOrderMax+0.5);
  mTPCCosShift->Sumw2();
  mTPCSinShift = new TProfile3D("TPCSinShift","TPCSinShift",9,-0.5,8.5,mFourierOrder,0.5,0.5+mFourierOrder,_PsiOrderMax,0.5,(double)_PsiOrderMax+0.5);
  mTPCSinShift->Sumw2();
//---------------Make histograms for checking the flattening of EP---------
//---------------and calculate the resolution------------------------------
  mHisto2D[0] = new TH2D("EPDPsiEvsPsiW","EPDPsiEvsPsiW",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());

  for(int cent=0;cent<9;cent++){
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      mEPDFullPsiWeighted[cent][iorder] = new TH1D(Form("EPDFullPsi%dWeightedCent%d",iorder,cent),Form("EPDFullPsi%dWeightedCent%d",iorder,cent),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
      mEPDFullPsiShifted[cent][iorder] = new TH1D(Form("EPDFullPsi%dShiftedCent%d",iorder,cent),Form("EPDFullPsi%dShiftedCent%d",iorder,cent),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
      mResolution[cent][iorder] = new TProfile(Form("ResolutionPsi%dCent%d",iorder,cent),Form("ResolutionPsi%dCent%d",iorder,cent),10,0.5,10.5);//x axis corresponds to the <cos()> of different combinations of EP
      mResolution[cent][iorder]->Sumw2();
      mTPCPsiDisWeighted[cent][iorder] = new TH1D(Form("TPCPsi%dDisWeightedCent%d",iorder,cent),Form("TPCPsi%dDisWeightedCent%d",iorder,cent),100,-TMath::Pi()/(iorder+1.0),TMath::Pi()/(iorder+1.0));
      mTPCPsiDisShifted[cent][iorder] = new TH1D(Form("TPCPsi%dDisShiftedCent%d",iorder,cent),Form("TPCPsi%dDisShifteddCent%d",iorder,cent),100,-TMath::Pi()/(iorder+1.0),TMath::Pi()/(iorder+1.0));
    }
    /*
    for(int i=0;i<3;i++){
      mTPCPsiDisWeighted[cent][i] = new TH1D(Form("TPCPsiDisWeightedCent%dTT%d",cent,i),Form("TPCPsiDisWeightedCent%dTT%d",cent,i),100,-TMath::Pi(),TMath::Pi());
      mTPCPsiDisShifted[cent][i] = new TH1D(Form("TPCPsiDisShiftedCent%dTT%d",cent,i),Form("TPCPsiDisShiftedCent%dTT%d",cent,i),100,-TMath::Pi(),TMath::Pi());
    }//3 types of Psi_TPC:Full, pos, neg
    */
    for(int ew=0;ew<2;ew++){
      mAveEta[cent][ew] = new TProfile2D(Form("AveEtaCent%dEW%d",cent,ew),Form("AveEtaCent%dEW%d",cent,ew),16,0.5,16.5,16,0,16);
      mAveEta[cent][ew]->Sumw2();
    }
  }

//---------------Make histograms for measuring dNdphi--------------------
  /*
  for(int cent=0;cent<9;cent++){
    mThreeD[cent][0] = new TH3D(Form("dNdphidnMIPdetaCent%dEW0",cent),Form("dNdphidnMIPdetaCent%dEW0",cent),96,-1.0*TMath::Pi(),TMath::Pi(),80,0.0,8.0,10,-5.1,-2.1);//(phi-EP),nMIP,eta
    mThreeD[cent][1] = new TH3D(Form("dNdphidnMIPdetaCent%dEW1",cent),Form("dNdphidnMIPdetaCent%dEW1",cent),96,-1.0*TMath::Pi(),TMath::Pi(),80,0.0,8.0,10,2.1,5.1);//(phi-EP),nMIP,eta
    mThreeD[cent][0]->Sumw2();
    mThreeD[cent][1]->Sumw2();
  }
  */
  /*
  for(int cent=0;cent<9;cent++){
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      for(int rr=0;rr<2;rr++){
        for(int ivz=0;ivz<16;ivz++){
          //24 delta phi bins 11/15/20
          mThreeD[cent][0][iorder][rr][ivz] = new TH3D(Form("dNdphidnMIPdetaCent%dEW0Psi%dRR%dVz%d",cent,iorder,rr,ivz),Form("dNdphidnMIPdetaCent%dEW0Psi%dRR%dVz%d",cent,iorder,rr,ivz),24,-1.0*TMath::Pi()/(double)(iorder+1),TMath::Pi()/(double)(iorder+1),80,0.0,8.0,16,0.5,16.5);//(phi-EP),nMIP,Ring Id
          mThreeD[cent][1][iorder][rr][ivz] = new TH3D(Form("dNdphidnMIPdetaCent%dEW1Psi%dRR%dVz%d",cent,iorder,rr,ivz),Form("dNdphidnMIPdetaCent%dEW1Psi%dRR%dVz%d",cent,iorder,rr,ivz),24,-1.0*TMath::Pi()/(double)(iorder+1),TMath::Pi()/(double)(iorder+1),80,0.0,8.0,16,0.5,16.5);//(phi-EP),nMIP,Ring Id
          mThreeD[cent][0][iorder][rr][ivz]->Sumw2();
          mThreeD[cent][1][iorder][rr][ivz]->Sumw2();
        }
      }
    }
    //mTwoD[cent] = new TH2D(Form("nMipinRingCent%d",cent),Form("nMipinRingCent%d",cent),80,0.0,8.0,33,-16.5,16.5);//X:nMIP, Y:RingId
    //mTwoD[cent]->Sumw2();
  }
  */
  /*
  for(int ew=0;ew<2;ew++){
    mThreeD[ew] = new TH3D(Form("dNdphidnMIPdetaCent%dEW%dPsi%dRR%dVz%d",5,ew,0,0,7),Form("dNdphidnMIPdetaCent%dEW%dPsi%dRR%dVz%d",5,ew,0,0,7),24,-1.0*TMath::Pi(),TMath::Pi(),80,0.0,8.0,16,0.5,16.5);//(phi-EP),nMIP,Ring Id
    mThreeD[ew]->Sumw2();
    for(int tt=0;tt<24;tt++){
      mThreeDTile[ew][tt] = new TH3D(Form("dNdphidnMIPdetaTileCent%dEW%dPsi%dRR%dVz%dTT%d",5,ew,0,0,7,tt),Form("dNdphidnMIPdetaTileCent%dEW%dPsi%dRR%dVz%dTT%d",5,ew,0,0,7,tt),24,-1.0*TMath::Pi(),TMath::Pi(),80,0.0,8.0,16,0.5,16.5);//(phi-EP),nMIP,Ring Id
      mThreeDTile[ew][tt]->Sumw2();
    }
  }
  */
  //Making TProfiles for the vn analysis
  for(int cent=0;cent<9;cent++){
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      mTPCvn[cent][iorder] = new TProfile(Form("TPCvnCent%dPsi%d",cent,iorder),Form("TPCvnCent%dPsi%d",cent,iorder),10,mEtaMinv1,mEtaMaxv1);
      mTPCvn[cent][iorder]->Sumw2();
      mEPDvn[0][cent][iorder] = new TProfile(Form("EPDvnEW0Cent%dPsi%d",cent,iorder),Form("EPDvnEW0Cent%dPsi%d",cent,iorder),mNEPDSubEvents,-5.1,-2.1);
      mEPDvn[0][cent][iorder]->Sumw2();
      mEPDvn[1][cent][iorder] = new TProfile(Form("EPDvnEW1Cent%dPsi%d",cent,iorder),Form("EPDvnEW1Cent%dPsi%d",cent,iorder),mNEPDSubEvents,2.1,5.1);
      mEPDvn[1][cent][iorder]->Sumw2();
    }
  }


//---------------Open histograms for getting TPCPhiWeight-----------------
  mTPCWeightFile = new TFile(TPCWeightFile,"READ");

  if(mTPCWeightFile->IsZombie()){
    std::cout<<"TPCWeightFile doesn't exist, didn't apply any phi weight :("<<std::endl;
    mTPCPhiWeightInput = 0;
    mTPCPhiAveragedInput = 0;
  }
  else{
    mTPCPhiWeightInput = (TH2D*)mTPCWeightFile->Get("TPCPhiWeightOutput");
    mTPCPhiAveragedInput = (TH2D*)mTPCWeightFile->Get("TPCPhiAveraged");
    //mTPCPhiWeightInput->RebinX();
    //mTPCPhiAveragedInput->RebinX();
    mTPCPhiWeightInput->Divide(mTPCPhiAveragedInput);
  }
//----------------Open histograms for the TPC Psi-shifting----------------
  mTPCShiftFile = new TFile(TPCShiftFile,"READ");

  if(mTPCShiftFile->IsZombie()){
    std::cout<<"TPCShiftFile doesn't exist, didn't do Psi-shifting :("<<std::endl;
    mTPCCosShift_Input = 0;
    mTPCSinShift_Input = 0;
  }
  else{
    mTPCCosShift_Input = (TProfile3D*)mTPCShiftFile->Get("TPCCosShift");
    mTPCSinShift_Input = (TProfile3D*)mTPCShiftFile->Get("TPCSinShift");
  }

//----------opening histos for getting the EPD Phi weights-------------
mEPDPhiWeightFile = new TFile(EPDPhiWeightFile,"READ");

if(mEPDPhiWeightFile->IsZombie()){
  std::cout<<"EPDPhiWeightFile doesn't exist, only used unweighted TileWeight :("<<std::endl;
  for(int ew=0;ew<2;ew++){
    mEPDPhiWeights[ew]=0;
  }
}
else{
  for(int ew=0;ew<2;ew++){
    mEPDPhiWeights[ew] = (TH3D*)mEPDPhiWeightFile->Get(Form("PhiWeightEW%d",ew));
  }
}
  cout << "Init done\n";
  return 0;
}

//==============================================
// some things might need resetting when there is a new run
void PicoAnalyzer::NewRun(int runId){
  //mRunCollisionSystem = WhichSystem(runId); //tend to cause error.
  mRunId = runId;
  //mRunEt += 1;
}

//=================================================
short PicoAnalyzer::Make(int iEvent){
  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  StPicoEvent* mPicoEvent = (StPicoEvent*)((*mEventClonesArray)[0]);
  if (mPicoEvent->runId()!=mRunId){
    NewRun(mPicoEvent->runId());        // some things should be reset when there is a new run loaded
    cout << "New run detected: " << mRunId << " and it is collision system #" << mRunCollisionSystem << endl;
    //cout<<"RunEntry"<<mRunEt<<endl;
  }

  //----- done getting data; have fun! ------
  // if(!(mPicoEvent->isTrigger(610001)||mPicoEvent->isTrigger(610011)||mPicoEvent->isTrigger(610021)||mPicoEvent->isTrigger(610031)||mPicoEvent->isTrigger(610041)||mPicoEvent->isTrigger(610051))) return 0;
  cout << "test 0 ! "  << endl;
  if(!(mPicoEvent->isTrigger(640002)||mPicoEvent->isTrigger(640012)||mPicoEvent->isTrigger(640022)||mPicoEvent->isTrigger(640032)||mPicoEvent->isTrigger(640001)||mPicoEvent->isTrigger(640011)||mPicoEvent->isTrigger(640021)||mPicoEvent->isTrigger(640031))) return 0;
  cout << "test 1 ! "  << endl;
  //  StThreeVectorF primaryVertex = mPicoEvent->primaryVertex();
  //  TVector3 vertexPos(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
  TVector3 vertexPos = mPicoEvent->primaryVertex();
  int mRunId = mPicoEvent->runId();
  // int RUNYear = 2018;
  // float RUNEnergy = 27.;
  // int RunDay = floor( (mRunId - (RUNYear-2000)*pow(10,6))/pow(10,3) );
  // int DayBinId = RunDay-89;

  double pi = TMath::Pi();

//-----------------Get the multiplicity by the StRefMulCorr class--------------
  int CentId=-1;
  int tofMult=mPicoEvent->btofTrayMultiplicity();
  int tofmatch=mPicoEvent->nBTOFMatch();
  int refMult=mPicoEvent->refMult();
  // bool ISRefMultCorrBadRun=false;
  // mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  //mRefMultCorr = new StRefMultCorr();
  // mRefMultCorr->init(mRunId);
  // if(mRefMultCorr->getBeginRun(RUNEnergy,RUNYear)==-1) return 0;
  // ISRefMultCorrBadRun=mRefMultCorr->isBadRun(mRunId);
  // if(ISRefMultCorrBadRun) return 0;
  // for(int ii=0;ii<44;ii++)
  // {
  //   if(mRunId == badrun[ii]) return 0;
  // }
  //mRefMultCorr->initEvent(refMult,vertexPos.Z(),mPicoEvent->ZDCx());
  //refMult = mRefMultCorr->getRefMultCorr();
  CentId = Centrality(mPicoEvent->refMult());//An integer between 0 (70-80%) and 8 (0-5%)
  //int CentIdmy = FindCent(mPicoEvent->refMult());   // returns an integer between 0 (70-80%) and 8 (0-5%)
  //-------------remove the pile-up events-----------------
  //if(mRefMultCorr->passnTofMatchRefmultCut(1.*mPicoEvent->refMult(), 1.*mPicoEvent->nBTOFMatch())!=1) return 0;

  hvz_b->Fill(vertexPos.Z());
  cout << "vertex z = "<<vertexPos.Z()<<endl;
  hvr_b->Fill(vertexPos.X(),vertexPos.Y());
  hvzvpdvz_b->Fill(vertexPos.Z(),mPicoEvent->vzVpd());
  hvzvpdvzdiff_b->Fill(vertexPos.Z()-mPicoEvent->vzVpd());
  htofvsref_b->Fill(refMult,tofMult);
  htofmatchvsref_b->Fill(refMult,tofmatch);
  h_runidvstofmult_b->Fill(mRunId,tofMult);
  h_runidvsrefmult_b->Fill(mRunId,refMult);
  if (fabs(vertexPos.Z())>=mVzMax) return 0;
  if (sqrt(pow(vertexPos.X(),2)+pow(vertexPos.Y(),2))>mVtxR) return 0;
  //if (abs(vertexPos.Z()-mPicoEvent->vzVpd())>mDiffVzVPD) return 0;//Get rid of the VPD cut, Prithwish said it does no good at low energy 06/18/2020
  if(TMath::Abs(vertexPos.Z()) < 10 ){
      href_vz->Fill(mPicoEvent->refMult());
  }
  href->Fill(mPicoEvent->refMult());
  hvz->Fill(vertexPos.Z());
  hvr->Fill(vertexPos.X(),vertexPos.Y());
  hvzvpdvz->Fill(vertexPos.Z(),mPicoEvent->vzVpd());
  hvzvpdvzdiff->Fill(vertexPos.Z()-mPicoEvent->vzVpd());
  htofvsref->Fill(refMult,tofMult);
  htofmatchvsref->Fill(refMult,tofmatch);

  h_runidvstofmult->Fill(mRunId,tofMult);
  h_runidvsrefmult->Fill(mRunId,refMult);
  if (CentId<0) return 0;            // 80-100% - very peripheral

  mVz[CentId]->Fill(vertexPos.Z());

  double d_resolution_EPD[9] = {0.26618012, 0.39441821, 0.53429421, 0.63668343, 0.68304687,
       0.67352165, 0.59120378, 0.44391744, 0.27105964};

  int VzBin;
  double VzArr[17]={-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40};
  for(int i=0;i<17;i++){
    if(vertexPos.Z()>VzArr[i]&&vertexPos.Z()<=VzArr[i+1]){
      VzBin=i;
      break;
    }
  }

  // if(VzBin!=7) return 0;

  StEpdEpInfo result = mEpFinder->Results(mEpdHits,vertexPos,CentId);  // and now you have all the EP info you could ever want :-)
  //StEpdEpInfo result = mEpFinder->Results(mEpdHits,vertexPos,1);  // testing. and now you have all the EP info you could ever want :-)
  double EpAngle[_PsiOrderMax][3];//EP order,east/west/full
  for(int iorder=0;iorder<_PsiOrderMax;iorder++){
    EpAngle[iorder][0]=result.EastPhiWeightedAndShiftedPsi(iorder+1);
    EpAngle[iorder][1]=result.WestPhiWeightedAndShiftedPsi(iorder+1);
    EpAngle[iorder][2]=result.FullPhiWeightedAndShiftedPsi(iorder+1);
  }


  mHisto2D[0]->Fill((double)EpAngle[0][1],(double)EpAngle[0][0]);//East vs. West
  //-----------Check the flattening of Psi_EPDFull----------
  for(int iorder=0;iorder<_PsiOrderMax;iorder++){
    mEPDFullPsiWeighted[CentId][iorder]->Fill(result.FullPhiWeightedPsi(iorder+1));
    mEPDFullPsiShifted[CentId][iorder]->Fill(result.FullPhiWeightedAndShiftedPsi(iorder+1));
  }

  //------------Prepare Qvectors for TPC EP (used in v1EPD)--------------
  TH1D* PCosPhi;
  TH1D* PSinPhi;
  //PCosPhi = new TH1D(Form("PCosPhi"),Form("PCosPhi"),3,0.5,3.5);
  //PSinPhi = new TH1D(Form("PSinPhi"),Form("PSinPhi"),3,0.5,3.5);//x is three types of Psi_TPC
  PCosPhi = new TH1D(Form("PCosPhi"),Form("PCosPhi"),_PsiOrderMax,0.5,_PsiOrderMax+0.5);
  PSinPhi = new TH1D(Form("PSinPhi"),Form("PSinPhi"),_PsiOrderMax,0.5,_PsiOrderMax+0.5);
  const Float_t   mField = mPicoEvent->bField(); // Magnetic field
  std::vector<StPicoTrack *> v_KaonPlus_tracks;
  std::vector<StPicoTrack *> v_KaonMinus_tracks;
  //------------Begin loop over TPC tracks--------------------------
  for(int itrk=0; itrk<mTracks->GetEntries(); itrk++){
    StPicoTrack* track = (StPicoTrack*)((*mTracks)[itrk]);
    if (!track->isPrimary()) continue; // I should check DCA too, but am confused how
    double nHitsFitRatio = track->nHitsFit()*1.0/track->nHitsMax();
    //cout<<nHitsFitRatio<<endl;
    TVector3 pMom = track->pMom();//track->gMom() if I want to look at the global tracks.
    StPicoPhysicalHelix helix = track->helix(mField);
    double Tphi = pMom.Phi();
    double Teta = pMom.Eta();
    double TPt = pMom.Pt();//
    double dca = track->gDCA(vertexPos).Mag();

    h_pt->Fill(pMom.Perp());
    h_eta_b->Fill(pMom.PseudoRapidity());
    h_nhitfit->Fill(track->nHitsFit());
    h_nhitmax->Fill(track->nHitsMax());
    h_nhitratio->Fill((float)track->nHitsFit()/(float)track->nHitsMax());
    h_dca->Fill(dca);
    h_phi->Fill(pMom.Phi());
    h_eta_phi_before->Fill(pMom.Phi(),pMom.PseudoRapidity());

    if(pMom.Perp() < 0.15) continue;
    if(pMom.Mag() > 10.0) continue;
    if (dca>mDCAcut) continue;
    if (track->nHitsFit()<mNhitsfit) continue;
    if (nHitsFitRatio<mNhitsfitratio) continue;//get rid of the nhitsfitratio cut, Prithwish said it is a very old cut used by STAR. 06/18/2020
    h_eta_phi->Fill(pMom.Phi(),pMom.PseudoRapidity());
    h_eta->Fill(pMom.PseudoRapidity());

    /* from Shaowei lan */
    int tofIndex = track->bTofPidTraitsIndex();
    Int_t   btofMatchFlag =  0;
    Float_t btofYLocal    =  -999;
    float tof = 0, L=0, beta=0.0, mass2= 0.0;
    if(tofIndex>=0) {
        StPicoBTofPidTraits *tofPid = (StPicoBTofPidTraits*)((*mTraits)[tofIndex]);
        btofMatchFlag = tofPid->btofMatchFlag();
        btofYLocal    = tofPid->btofYLocal();
        if(tofPid) {
            beta = tofPid->btofBeta();
            tof = tofPid->btof();
            if(beta<1e-4) {
              TVector3 btofHitPos_ = tofPid->btofHitPos();
              const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());
              const StThreeVectorF *vertexPos_ = new StThreeVectorF(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
              L = tofPathLength(vertexPos_, btofHitPos, helix.curvature());
                if(tof>0) beta = L/(tof*(TMath::C()/1.e7));
                else beta = -1;
            }
            // cout << "L = " << L << endl;
            // cout << "time of flight = " <<  tof << endl;
            // cout << "tof xpostion = " <<  tofPid->btofHitPos().X() << endl;
            // cout << "tof*(TMath::C()/1.e7 = " << (tof*(TMath::C()/1.e7)) << endl;
            // cout << "beta = " << beta << endl;
        }
    }
    bool isGoodTof = btofMatchFlag >0 && beta > 0 /*&& fabs(btofYLocal) < 1.8*/;
    if(isGoodTof){
      hbtofYLocal -> Fill(btofYLocal);
      mass2 = pMom.Mag()*pMom.Mag()*(1./pow(beta,2)-1);
      hbtofYLocalvsMass2 -> Fill(mass2,btofYLocal);
    }  else mass2 = -999;

    if(TMath::Abs(beta)>1e-5)  hbetavsp->Fill(pMom.Mag(), 1/beta);
    else hbetavsp->Fill(pMom.Mag(), 0);
    hmassvsp->Fill(pMom.Mag()/track->charge(), mass2);
    hdedxvsp->Fill(pMom.Mag()/track->charge(),track->dEdx());
    /* above from Shaowei lan */
    int Tch=track->charge();
    double rig=Tch*pMom.Mag();
    // double dEdx=track->dEdx();

    // Kaons PID: require both TPC and TOF
    TLorentzVector ltrackk;
    if(
      TMath::Abs(track->nSigmaKaon()) < mKaonSigma &&
      isGoodTof && mass2 > d_KaonM2low && mass2 < d_KaonM2high
      && TPt >= mpTMink
      && TPt <= mpTMaxk
    ){
      ltrackk.SetXYZM(pMom.X(),pMom.Y(),pMom.Z(),_massKaon);
      if(Tch > 0){
        v_KaonPlus_tracks.push_back(track); // push back kp tracks
        // Fill histograms
        hist_pt_kaonPlus->Fill(TPt);
        hist_eta_kaonPlus->Fill(Teta);
        hist_y_kaonPlus->Fill(ltrackk.Rapidity());
        hist_phi_kaonPlus->Fill(Tphi);
        hist_rap_eta_kaonPlus->Fill(Teta,ltrackk.Rapidity());
        hist_pt_y_kaonPlus->Fill(ltrackk.Rapidity(),TPt,1);
        hist_pt_eta_kaonPlus->Fill(Teta,TPt,1);
        hist_dEdx_kaonPlus->Fill(rig,track->dEdx());
        hist_beta_kaonPlus->Fill(rig,1.0/beta);
        hist_mass_kaonPlus->Fill(rig,mass2);
      } else { // charge < 0
        v_KaonMinus_tracks.push_back(track); // push back km tracks
        // Fill histograms
        hist_pt_kaonMinus->Fill(TPt);
        hist_eta_kaonMinus->Fill(Teta);
        hist_y_kaonMinus->Fill(ltrackk.Rapidity());
        hist_phi_kaonMinus->Fill(Tphi);
        hist_rap_eta_kaonMinus->Fill(Teta,ltrackk.Rapidity());
        hist_pt_y_kaonMinus->Fill(ltrackk.Rapidity(),TPt,1);
        hist_pt_eta_kaonMinus->Fill(Teta,TPt,1);
        hist_dEdx_kaonMinus->Fill(rig,track->dEdx());
        hist_beta_kaonMinus->Fill(rig,1.0/beta);
        hist_mass_kaonMinus->Fill(rig,mass2);
      }
    }
    if(TMath::Abs(Teta)>=mEtaMaxv1) continue;
    //mHisto1D[3]->Fill(TPt);

    if(TPt<=mpTMin||TPt>=mpTMax) continue;//pT range [0.15,2.0] for Psi_TPC
    //----------Prepare Q to calculate the Psi_TPC-----------------
    int TtrId=FindTrackId(Tch,vertexPos.Z(),Teta,TPt);
    if(TtrId<=0||TtrId>mNumberOfTrackTypes){
      cout<<"Ah-oh, invalid TrackId :("<<endl;
    }
    double TrackPhiWeight=1.0;
    if(mTPCPhiWeightInput!=0) TrackPhiWeight=1.0/mTPCPhiWeightInput->GetBinContent(mTPCPhiWeightInput->FindBin(Tphi,TtrId));
    if(!TMath::Finite(TrackPhiWeight)) continue;

    double TPCWeight[_PsiOrderMax];
    TPCWeight[0] = TrackPhiWeight*(-1.0*Teta);//weight the TPC tracks with (-eta) for Psi1
    TPCWeight[1] = TrackPhiWeight*TPt;//weight the TPC tracks with pT for Psi2

    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      PCosPhi->Fill(iorder+1,TPCWeight[iorder]*cos((double)(iorder+1)*Tphi));
      PSinPhi->Fill(iorder+1,TPCWeight[iorder]*sin((double)(iorder+1)*Tphi));
    }

    //Fill the TProfiles for the vn(TPC)
    for(int iorder=0;iorder<2;iorder++){
      mTPCvn[CentId][iorder]->Fill(Teta,cos((double)(iorder+1)*(Tphi-EpAngle[iorder][2])));
    }

  }
//-------------End loop over TPC tracks and get kaon track vectors-----------------------

//-------------Now let's play with the EP---------------------
  //double TPCPsi1[3]={0.0};// Three types of Psi_TPC: 0_Full, 1_pos, 2_neg
  double TPCPsi[_PsiOrderMax]={0.0};
  for(int i=1; i<=_PsiOrderMax;i++){
    TPCPsi[i-1] = atan2(PSinPhi->GetBinContent(i),PCosPhi->GetBinContent(i))/(double)i;
    //atan2 returns angle between -pi and pi, remember to divide by i!!!! 06/16/2020
  }

//--------Fill the histos for Psi-shifting the TPC_EP---------
  for(int k=1;k<=mFourierOrder;k++){
    for(int i=1;i<=_PsiOrderMax;i++){
      double tmp=(double)(i*k);
      mTPCCosShift->Fill(CentId,k,i,cos(tmp*TPCPsi[i-1]));
      mTPCSinShift->Fill(CentId,k,i,sin(tmp*TPCPsi[i-1]));
    }
  }
//-------------Shift the TPC EP if the shifting histos exist------------
  double TPCPsiShifted[_PsiOrderMax];
  for(int i=0;i<_PsiOrderMax;i++){
    TPCPsiShifted[i]=TPCPsi[i];
  }
  if(mTPCCosShift_Input!=0&&mTPCSinShift_Input!=0){
    double deltapsi[_PsiOrderMax]={0.0};
    for(int i=1;i<=_PsiOrderMax;i++){
      for(int k=1;k<=mFourierOrder;k++){
        double tmp=(double)(i*k);
        deltapsi[i-1]+=2.0*(-mTPCSinShift_Input->GetBinContent(CentId+1,k,i)*cos(tmp*TPCPsi[i-1])+mTPCCosShift_Input->GetBinContent(CentId+1,k,i)*sin(tmp*TPCPsi[i-1]))/tmp;
      }
      TPCPsiShifted[i-1]+=deltapsi[i-1];
      if(TPCPsiShifted[i-1]<-pi/(double)i) TPCPsiShifted[i-1]+=2*pi/(double)i;
      if(TPCPsiShifted[i-1]>pi/(double)i) TPCPsiShifted[i-1]-=2*pi/(double)i;
    }
  }

//-------------Fill the TProfiles for calculating the Resolution---------
  for(int i=0;i<_PsiOrderMax;i++){
    mResolution[CentId][i]->Fill(1,cos((double)(i+1)*(EpAngle[i][0]-TPCPsiShifted[i])));
    mResolution[CentId][i]->Fill(2,cos((double)(i+1)*(EpAngle[i][1]-TPCPsiShifted[i])));
    mResolution[CentId][i]->Fill(3,cos((double)(i+1)*(EpAngle[i][0]-EpAngle[i][1])));
  }

//----Fill the histograms for checking the flattening of the TPC EPs-----
  for(int i=0;i<_PsiOrderMax;i++){
    mTPCPsiDisWeighted[CentId][i]->Fill(TPCPsi[i]);
    mTPCPsiDisShifted[CentId][i]->Fill(TPCPsiShifted[i]);
  }
//--------------Begin loop over EPD hits---------------------------------
  for(int hit=0;hit<mEpdHits->GetEntries();hit++){
    int tileId, ring, TT, PP, EW, ADC;
    float nMip;
    StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[hit]);
    tileId=epdHit->id();
    EW = (tileId<0)?0:1;
    ring = epdHit->row();//[1,16]
    TT = epdHit->tile();
    PP = epdHit->position();
    ADC = epdHit->adc();
    nMip = epdHit->nMIP();

    int TTxyId;
    if(EW==0){
      int OE=(TT%2)?2:1;//odd 2, even 1
      int PPId=(PP>9)?(21-PP):(9-PP);
      TTxyId=PPId*2+OE;
    }
    else{
      int OE=(TT%2)?1:2;//odd 1, even 2
      int PPId=(PP<4)?(PP+8):(PP-4);
      TTxyId=PPId*2+OE;
    }

    //if (nMip<mEPDthresh) continue;
    double TileWeight = (nMip<4.0)?nMip:4.0;//Note: here a different nMIPMax from the EpFinder was used
    TVector3 StraightLine = mEpdGeom->RandomPointOnTile(tileId) - vertexPos;
    double Hphi = StraightLine.Phi();
    double Heta = StraightLine.Eta();

    mAveEta[CentId][EW]->Fill(ring,VzBin,Heta);

    //Fill the TProfiles for EPDvn
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      mEPDvn[EW][CentId][iorder]->Fill(Heta,cos((double)(iorder+1)*(Hphi-TPCPsiShifted[iorder])),TileWeight);
    }

    double deltaphi[_PsiOrderMax][2];//order, reference: 0-TPC,1-the other side of EPD
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      deltaphi[iorder][0]=Hphi-TPCPsiShifted[iorder];
      deltaphi[iorder][1]=Hphi-EpAngle[iorder][1-EW];//Note EPD has different EP ranges as TPC or phi
    }

    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      for(int rr=0;rr<2;rr++){
        double tmpn=pi/(double)(iorder+1);
        while(deltaphi[iorder][rr]<-tmpn||deltaphi[iorder][rr]>tmpn){
          if(deltaphi[iorder][rr]<-tmpn) deltaphi[iorder][rr]+=2.0*tmpn;
          else if(deltaphi[iorder][rr]>tmpn) deltaphi[iorder][rr]-=2.0*tmpn;
        }
      }
    }
    /*
    for(int iorder=0;iorder<2;iorder++){
      for(int rr=0;rr<2;rr++){
        mThreeD[CentId][EW][iorder][rr][VzBin]->Fill(deltaphi[iorder][rr],nMip,ring);
      }
    }
    */
    /*
    mThreeD[EW]->Fill(deltaphi[0][0],nMip,ring);
    mThreeDTile[EW][TTxyId-1]->Fill(deltaphi[0][0],nMip,ring);
    */

    //double RingId=0.0;
    //if(EW==0) RingId=-1.0*(double)ring;
    //else RingId=(double)ring;
    //mTwoD[CentId]->Fill(nMip,RingId);

  }
  //-----------------End looping over EPD hits------------------------

  // ------------------ Phi meson flow analysis ----------------------
  TLorentzVector ltrackA, ltrackB;
  for(unsigned int i = 0; i < v_KaonPlus_tracks.size(); i++){
    StPicoTrack * picoTrackA = v_KaonPlus_tracks.at(i); // i-th K+ track
    if(!picoTrackA) continue;
    // K+ Variables
    StPicoPhysicalHelix   trackhelixA = picoTrackA->helix(mField);
    TVector3 p_vecA = picoTrackA->pMom();  // primary momentum
    // TVector3 p_vecA = trackhelixA.cat(trackhelixA.pathLength(vertexPos));  // primary momentum
    // p_vecA *= (double)picoTrackA->pMom().Mag();  // primary momentum
    ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),_massKaon);
    double d_chargeA  = picoTrackA->charge();
    h2px  -> Fill(picoTrackA->pMom().X(),p_vecA.X());
    h2py  -> Fill(picoTrackA->pMom().Y(),p_vecA.Y());
    h2pz  -> Fill(picoTrackA->pMom().Z(),p_vecA.Z());
    Double_t d_ptA = ltrackA.Perp(), d_pzA = ltrackA.Pz(), d_momA = ltrackA.P();
    // StPicoBTofPidTraits *traitA = NULL;
    // double d_tofBeta0    = -999.;
    // double d_inv_tofBeta0    = -999.;
    // if(picoTrackA->isTofTrack()) traitA = (StPicoBTofPidTraits*)((*mTraits)[picoTrackA->bTofPidTraitsIndex()]);
    // if(traitA)        d_tofBeta0 = traitA->btofBeta();
    // double d_M0   = _massKaon;
    // double d_E0   = sqrt((d_pxA*d_pxA+d_pyA*d_pyA+d_pzA*d_pzA)+_massKaon*_massKaon);
    // double d_y0   = ((d_E0-d_pzA) != 0.0) ? 0.5*TMath::Log( (d_E0 + d_pzA) / (d_E0 - d_pzA) ) : -999.0;
    // double eta0   = ((d_momA - d_pzA) != 0.0) ? 0.5*TMath::Log( (d_momA + d_pzA) / (d_momA - d_pzA) ) : -999.0;
    // double d_mT0  = sqrt(d_ptA*d_ptA + d_M0*d_M0);
    // double d_pq0   = fabs(d_momA) * d_chargeA;
    for(unsigned int j = 0; j < v_KaonMinus_tracks.size(); j++){
      StPicoTrack * picoTrackB = v_KaonMinus_tracks.at(j); // j-th K- track
      if(!picoTrackB) continue;
      // K- Variables
      double d_chargeB  = picoTrackB->charge();
      if(d_chargeA == d_chargeB) continue; // same charge cut
      TVector3 p_vecB = picoTrackB->pMom();  // primary momentum
      ltrackB.SetXYZM(p_vecB.X(),p_vecB.Y(),p_vecB.Z(),_massKaon);
      Double_t d_ptB = ltrackB.Perp(), d_pzB = ltrackB.Pz(), d_momB = ltrackB.P();

      // StPicoBTofPidTraits *traitB = NULL;
      // double d_tofBeta1    = -999.;
      // double d_inv_tofBeta1    = -999.;
      // if(picoTrackB->isTofTrack()) traitB = (StPicoBTofPidTraits*)((*mTraits)[picoTrackB->bTofPidTraitsIndex()]);
      // if(traitB)        d_tofBeta1 = traitB->btofBeta();
      // double d_M1   = _massKaon;
      // double d_E1   = sqrt((d_px1*d_px1+d_py1*d_py1+d_pzB*d_pzB)+_massKaon*_massKaon);
      // double d_y1   = ((d_E1-d_pzB) != 0.0) ? 0.5*TMath::Log( (d_E1 + d_pzB) / (d_E1 - d_pzB) ) : -999.0;
      // double eta1   = ((d_momB - d_pzB) != 0.0) ? 0.5*TMath::Log( (d_momB + d_pzB) / (d_momB - d_pzB) ) : -999.0;
      // double d_mT1  = sqrt(d_ptB*d_ptB + d_M1*d_M1);
      // double d_pq1   = fabs(d_momA) * d_chargeA;
      // phi Variables
      TLorentzVector trackAB      = ltrackA+ltrackB;
      Double_t InvMassAB          = trackAB.M();
      Double_t d_dip_angle = TMath::ACos((d_ptA*d_ptB+d_pzA*d_pzB) / (d_momA*d_momB) );

      Double_t pt = trackAB.Perp();
      Double_t eta = trackAB.Eta();
      Double_t rap = trackAB.Rapidity();
      Double_t d_mT_phi = sqrt(pt*pt + _massPhi*_massPhi );

      Double_t randomNumber = gRandom->Uniform(1);
      // std::cout << "randomNumber " << randomNumber  << std::endl;
      double d_randAngle = TMath::Pi()*randomNumber;
      // std::cout << "randomAngle " << d_randAngle  << std::endl;
      TLorentzVector ltrackB_rot = ltrackB;
      ltrackB_rot.RotateZ(d_randAngle);
      TLorentzVector trackAB_rot      = ltrackA+ltrackB_rot;
      Double_t InvMassAB_rot          = trackAB_rot.M();
      Double_t pt_rot = trackAB_rot.Perp();
      Double_t rap_rot =  trackAB_rot.Rapidity(); // Rotation don't  influence y

      h_Mass    ->Fill(InvMassAB);
      h_Mass_rot    ->Fill(InvMassAB_rot);
      h_Mass2    ->Fill(pt,InvMassAB);
      h_Mass2_rot    ->Fill(pt_rot,InvMassAB_rot);
      hist_SE_PhiMeson_pT ->Fill(pt);
      h2_pT_eta->Fill(eta,pt);
      h2_pT_y->Fill(rap,pt);
      hist_SE_PhiMeson_mT ->Fill(d_mT_phi);
      hist_SE_PhiMeson_rap ->Fill(rap);
      hist_SE_PhiMeson_eta ->Fill(eta);
      if(CentId >= 7 && CentId <= 8){ // 0-10%
        hist_SE_pt_y_PhiMeson[0] ->Fill(rap,pt);
        hist_SE_pt_y_PhiMeson[1] ->Fill(rap,pt);
        if(InvMassAB >= 1.005 && InvMassAB <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_SigBkg[0] -> Fill(rap,pt);
          hist_SE_pt_y_Phi_tight_SigBkg[1] -> Fill(rap,pt);
        }
        if(InvMassAB_rot >= 1.005 && InvMassAB_rot <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_Bkg[0] -> Fill(rap_rot,pt_rot);
          hist_SE_pt_y_Phi_tight_Bkg[1] -> Fill(rap_rot,pt_rot);
        }
      }
      if(CentId >= 4 && CentId <= 6){ // 10-40%
        hist_SE_pt_y_PhiMeson[0] ->Fill(rap,pt);
        hist_SE_pt_y_PhiMeson[2] ->Fill(rap,pt);
        if(InvMassAB >= 1.005 && InvMassAB <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_SigBkg[0] -> Fill(rap,pt);
          hist_SE_pt_y_Phi_tight_SigBkg[2] -> Fill(rap,pt);
        }
        if(InvMassAB_rot >= 1.005 && InvMassAB_rot <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_Bkg[0] -> Fill(rap_rot,pt_rot);
          hist_SE_pt_y_Phi_tight_Bkg[2] -> Fill(rap_rot,pt_rot);
        }
      }
      if(CentId >= 0 && CentId <= 3){ // 40-80%
        hist_SE_pt_y_PhiMeson[0] ->Fill(rap,pt);
        hist_SE_pt_y_PhiMeson[3] ->Fill(rap,pt);
        if(InvMassAB >= 1.005 && InvMassAB <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_SigBkg[0] -> Fill(rap,pt);
          hist_SE_pt_y_Phi_tight_SigBkg[3] -> Fill(rap,pt);
        }
        if(InvMassAB_rot >= 1.005 && InvMassAB_rot <= 1.033){ // tight phi-mass cut
          hist_SE_pt_y_Phi_tight_Bkg[0] -> Fill(rap_rot,pt_rot);
          hist_SE_pt_y_Phi_tight_Bkg[3] -> Fill(rap_rot,pt_rot);
        }
      }
      // ---------------- phi-meson cuts: decay length, dip angle ------------
      StPicoPhysicalHelix    trackhelixB = picoTrackB->helix(mField);
      pair<double,double> pairLengths = trackhelixA.pathLengths(trackhelixB);
      TVector3 v3D_x_daughterA = trackhelixA.at(pairLengths.first);
      TVector3 v3D_x_daughterB = trackhelixB.at(pairLengths.second);
      TVector3 v3D_x_AB    = (v3D_x_daughterA+v3D_x_daughterB)*0.5;
      TVector3 v3D_xvec_decayl = v3D_x_AB - vertexPos;
      // double d_AB_decay_length =  v3D_xvec_decayl.Mag();
      // hist_AB_decay_length->Fill(d_AB_decay_length);
      // if(d_AB_decay_length > d_cut_AB_decay_length_PHI) continue; //decay length cut

      h_dip_angle         ->Fill(d_dip_angle);
      if(d_dip_angle <= dip_angle_cutLevel) continue; // dip-angle cut
      // --------------------- phi-meson flows -------------------------------
      TVector3 v3D_p_daughterA = trackhelixA.momentumAt(pairLengths.first, mField*kilogauss);
      TVector3 v3D_p_daughterB = trackhelixB.momentumAt(pairLengths.second, mField*kilogauss);
      TVector3 v3D_p_AB    = v3D_p_daughterA+v3D_p_daughterB;
      // double d_pmom = v3D_xvec_decayl.Dot(v3D_p_AB);
      // double d_dca_AB = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_AB.Mag2()) );
      double d_phi_azimuth = v3D_p_AB.Phi();
      if(d_phi_azimuth < 0.0            ) d_phi_azimuth += 2.0*TMath::Pi();
      if(d_phi_azimuth > 2.0*TMath::Pi()) d_phi_azimuth -= 2.0*TMath::Pi();
      double d_flow_PHI_raw[2] = {-999.0,-999.0}; // v1, v2 raw flow
      double d_flow_PHI_resolution[2] = {-999.0,-999.0}; // v1, v2 flow corrected by resolution
      if(EpAngle[0][2]!=-999.0){// Using EPD-full
        for(int km=0;km<2;km++){ // km - flow order
          d_flow_PHI_raw[km]        = TMath::Cos((double)(km+1.) * (d_phi_azimuth - EpAngle[0][2]));
          d_flow_PHI_resolution[km] = TMath::Cos((double)(km+1.) * (d_phi_azimuth - EpAngle[0][2]))/(d_resolution_EPD[CentId]); // km {0,1}, centrality [1,9]
          // d_flow_PHI_resolution[km] = d_flow_PHI_raw[km];
        }
      }
      int cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
      int cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
      double rap_low_phi[8] = {-1.0, -0.6, -0.3, -0.1, 0., 0.1, 0.3, 0.6};
      double rap_up_phi[8]  = {-0.6, -0.3, -0.1, 0.,  0.1, 0.3, 0.6, 1.0};


      for(Int_t cent = 0; cent < Bin_Centrality_01; cent++)
      {
          for(Int_t rap_bin = 0; rap_bin < Bin_rap; rap_bin++)
          {
            if(cent_low[cent]<= CentId && CentId <= cent_up[cent] &&
               rap_low_phi[rap_bin] <= rap && rap < rap_up_phi[rap_bin])
               {
                 mHist_SE_InvM_rap_cent[rap_bin][cent]->Fill(InvMassAB);
                 // std::cout << "invM = " << InvMassAB << std::endl;
                 if(!(EpAngle[0][2] == -999.0 || d_flow_PHI_resolution[0] == -999.0))
                 {
                   // if(rap_bin==0)std::cout << "EpAngle  = " << EpAngle << std::endl;
                   // if(rap_bin==0)std::cout << "Res_EP  = " << Res_EP << std::endl;
                   // if(rap_bin==0)std::cout << "ptbin1 d_flow_PHI_resolution[0] = " << d_flow_PHI_resolution[0] << std::endl;
                   mProfile_flow_reso_rap_cent[rap_bin][cent]->Fill(InvMassAB,d_flow_PHI_resolution[0]);
                 }
               }
            if(cent_low[cent]<= CentId && CentId <= cent_up[cent] &&
               rap_low_phi[rap_bin] <= rap_rot && rap_rot <= rap_up_phi[rap_bin])
               {
                 // std::cout << "invM rot = " << InvMassAB_rot << std::endl;
                 mHist_rotation_InvM_rap_cent[rap_bin][cent]->Fill(InvMassAB_rot);
               }

          }
      }

    }
  }


  delete PCosPhi;
  delete PSinPhi;
  v_KaonPlus_tracks.clear();
  v_KaonMinus_tracks.clear();
  return 0;
}
//=================================================
short PicoAnalyzer::Finish(){

  // subtraction
  for(int cent=0;cent<4;cent++){
    hist_SE_pt_y_Phi_tight_Sig[cent] = (TH2D*) hist_SE_pt_y_Phi_tight_SigBkg[cent]->Clone(Form("hist_SE_pt_y_Phi_tight_Sig_%d",cent));
    hist_SE_pt_y_Phi_tight_Sig[cent]->Add(hist_SE_pt_y_Phi_tight_Bkg[cent],-1.);
  }
  cout << "In PicoAnalyzer::Finish - calling StEpdEpFinder::Finish()\n";
  mEpFinder->Finish();

  cout << "I have called it\n";

  mTPCWeightFile->Close();
  mTPCShiftFile->Close();
  mEPDPhiWeightFile->Close();
  cout<<"Hello, I just want some weighting hisots"<<endl;

  mHistoFile->Write();
  mHistoFile->Close();

  cout << "Finish!!\n\n";

  return 0;
}
/*
bool PicoAnalyzer::Runlist(int runId){
    bool runpass=0;
    int runlist[]={19130071,19130077,19130078,19130079,19130084,19130085,19130086,19131001,19131003,19131004,19131005,19131006,19131007,19131009,19131010,19131012,19131013,19131014,19131015,19131016,19131019,19131020,19131026,19131027,19131030,19131031,19131037,19131039,19131040,19131041,19131042,19131044,19131045,19131046,19131048,19131049,19131050,19131051,19131052,19131054,19131055,19131056,19131057,19131060,19131061,19131062,19132001,19132004,19132005,19132006,19132007,19132008,19132011,19132013,19132014,19132016,19132017,19132019,19132020,19132021,19132022,19132026,19132027,19132029,19132030,19132031,19132032,19132034,19132035,19132036,19132037,19132038,19132039,19132043,19132044,19132045,19132046,19132047,19132048,19132063,19132064,19132065,19132070,19132071,19132073,19132074,19132075,19132076,19132078,19132079,19132080,19132081,19132082,19132083,19133002,19133003,19133004,19133005,19133008,19133009,19133010,19133012,19133013,19133014,19133016,19133017,19133018,19133021,19133022,19133023,19133025,19133026,19133027,19133028,19133030,19133031,19133032,19133033,19133038,19133039,19133040,19133041,19133043,19133044,19133045,19133047,19133048,19133049,19133050,19133053,19133054,19133055,19133056,19133058,19133059,19133060,19133061,19134002,19134003,19134004,19134005,19134007,19134008,19134009,19134010,19134011,19134013,19134014,19134016,19134017,19134018,19134019,19134023,19134024,19134025,19134027,19134028,19134029,19134030,19134035,19134036,19134037,19134038,19134040,19134041,19134042,19134044,19134045,19134046,19134047,19134049,19134050,19135001,19135003,19135004,19135005,19135012,19135013,19135014,19135016,19135017,19135018,19135020,19135021,19135022,19135027,19135028,19135029,19135037,19135038,19135039,19135040,19135042,19135043,19136001,19136003,19136004,19136005,19136007,19136008,19136009,19136011,19136012,19136013,19136014,19136016,19136017,19136018,19136040,19136041,19136042,19136044,19136045,19136046,19136047,19136049,19137001,19137002,19137003,19137004,19137006,19137007,19137008,19137009,19137010,19137011,19137013,19137014,19137015,19137016,19137019,19137020,19137022,19137023,19137024,19137025,19137026,19137027,19137028,19137029,19137036,19137037,19137038,19137040,19137041,19137045,19137047,19137048,19137050,19137051,19137052,19137053,19137054,19137056,19137057,19137058,19137059,19138002,19138003,19138004,19138005,19138006,19138008,19138009,19138010,19138011,19138012,19138014,19138015,19138016,19138019,19138020,19138021,19138025,19138026,19138027,19138028,19139023,19139024,19139026,19139027,19139028,19139032,19139033,19139034,19139037,19139038,19139039,19139041,19139042,19139043,19139044,19139050,19139058,19139063,19139064,19139066,19139067,19139068,19139069,19139071,19139072,19139073,19140002,19140003,19140004,19140006,19140007,19140008,19140010,19140011,19140012,19140014,19140015,19140016,19140017,19140020,19140021,19140022,19140025,19140030,19140031,19140035,19140036,19140037,19140040,19140041,19140042,19140043,19140045,19140046,19140047,19140051,19140052,19140053,19140055,19140056,19141001,19141003,19141004,19141005,19141006,19141008,19141009,19141010,19141011,19141013,19141014,19141015,19141018,19141019,19141020,19141021,19141022,19141024,19141025,19141026,19141028,19141029,19141030,19141047,19141048,19141049,19141051,19141052,19141053,19142001,19142002,19142003,19142005,19142006,19142007,19142008,19142010,19142011,19142012,19142014,19142015,19142016,19142018,19142019,19142020,19142021,19142022,19142027,19142034,19142035,19142038,19142039,19142041,19142045,19142048,19142049,19142053,19142054,19142055,19142057,19142058,19142059,19142062,19142064,19142065,19142068,19143001,19143003,19143006,19143007,19143008,19143009,19143010,19143011,19143012,19143013,19143014,19143015,19143016,19143017,19144012,19144013,19144014,19144018,19144019,19144020,19144024,19144025,19144026,19144031,19144032,19144033,19144036,19144037,19144038,19144042,19144043,19144044,19144046,19144047,19145001,19145004,19145005,19145006,19145008,19145009,19145010,19145011,19145013,19145014,19145015,19145017,19145019,19145020,19145028,19145031,19145034,19145035,19145036,19145038,19145039,19145040,19145042,19145043,19145044,19145047,19145048,19145050,19146002,19146003,19146004,19146006,19146007,19146008,19146009,19146012,19146013,19146014,19146016,19146017,19146019,19146020,19146024,19146025,19146026,19147007,19147008,19147009,19147014,19147015,19147016,19147021,19147022,19147023,19147025,19147026,19147027,19147029,19147030,19147031,19147033,19147034,19147035,19147038,19147039,19147040,19147042,19147043,19147044,19147046,19147047,19147048,19148002,19148003,19148004,19148007,19148008,19148009,19148011,19148012,19148013,19148015,19148016,19148017,19148021,19148022,19148023,19148024,19148050,19148051,19148052,19149002,19149003,19149007,19149008,19149009,19149012,19149013,19149014,19149015,19149017,19149018,19149023,19149024,19149025,19149030,19149031,19149032,19149035,19149037,19149038,19149039,19149042,19149043,19149044,19149046,19149047,19149048,19149050,19149051,19149052,19149054,19149055,19150001,19150005,19150006,19150007,19150009,19150010,19150011,19150013,19150014,19150016,19150017,19150018,19155057,19155058,19156001,19156002,19156004,19156005,19156006,19156008,19156009,19156011,19156013,19156014,19156015,19156018,19156019,19156030,19156031,19156032,19156033,19156042,19156043,19156044,19156045,19156046,19156047,19157002,19157003,19157004,19157006,19157007,19157008,19157012,19157013,19157015,19157017,19157018,19158020,19158059,19158060,19158062,19159001,19159003,19159004,19159006,19159007,19159008,19159010,19159011,19159012,19159014,19159015,19159016,19159018,19159019,19159021,19159023,19159024,19159025,19159039,19159040,19159041,19160002,19160003,19160004,19160006,19160007,19160008,19160010,19160011,19160012,19160014,19160015,19160016,19160019,19160020,19160021,19160023,19160024,19160025,19160027,19160028,19160029,19161003,19161004,19161005,19161007,19161008,19161009,19161011,19161012,19161013,19161015,19161016,19161017,19161048,19161049,19161050,19161051,19161053,19162002,19162006,19162008,19162009,19162010,19162013,19162014,19162015,19162017,19162018,19162019,19162021,19162022,19162023,19162030,19162031,19162032,19162036,19162037,19162038,19162040,19162041,19162042,19163006,19163007,19163008,19163010,19163011,19163012,19163014,19163015,19163016,19163018,19163019,19163020,19163039,19163040,19163041,19164004,19164005,19164006,19164008,19164009,19164010,19164011,19164012,19164014,19164015,19164016,19164017,19165002,19165003,19165004,19165006,19165007,19165008,19165009,19165011,19165012,19165013,19165015,19165016,19165017,19165018,19165020,19165021,19165027,19165028,19165029,19166006,19166007,19166008,19166010,19166011,19166012,19166014,19166015,19167003,19167004,19167005,19167007,19167008,19167009,19167011,19167012,19167013,19167015,19167016,19167017,19167021,19167022,19167023,19167025,19167026,19167027,19167028,19167030,19167031,19167032,19167034,19167035,19167042,19167044,19167045,19167046,19167049,19168025,19168026,19168028,19168029,19168030,19168033,19168034,19168036,19168038,19168039,19168040};
    for(int i=0;i<788;i++){
      if(runId==runlist[i]){
        runpass=1;
        break;
      }
    }
    return runpass;
}
*/

//TrackId starts at 1!
int PicoAnalyzer::FindTrackId(int Trkch,double TrkVz,double TrkEta,double TrkpT){
  int TrackId=0;
  int chId=(Trkch==-1)?0:1;
  int VzId=-1;
  if(TrkVz<0.0) VzId=-(int)(mNVzbin*abs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2-1;
  else VzId=(int)(mNVzbin*abs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2;
  int EtaId=-1;
  if(TrkEta<0.0) EtaId=-(int)(mNEtabin*abs(TrkEta)/(mEtaMax-mEtaMin))+mNEtabin/2-1;
  else EtaId=(int)(mNEtabin*abs(TrkEta)/(mEtaMax-mEtaMin))+mNEtabin/2;
  int pTId=-1;
  pTId=(int)(mNpTbin*abs(TrkpT-mpTMin)/(mpTMax-mpTMin));

  TrackId=chId*mNpTbin*mNVzbin*mNEtabin+VzId*mNpTbin*mNEtabin+EtaId*mNpTbin+pTId+1;

  return TrackId;
}

int Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={4, 9,17,30,50,78, 116,170,205};
    if      (gRefMult>=centFull[8]) centrality=8; // 0 - 5%
    else if (gRefMult>=centFull[7]) centrality=7; // 5 - 10%
    else if (gRefMult>=centFull[6]) centrality=6; // 10 - 20%
    else if (gRefMult>=centFull[5]) centrality=5; // 20 - 30%
    else if (gRefMult>=centFull[4]) centrality=4; // 30 - 40%
    else if (gRefMult>=centFull[3]) centrality=3; // 40 - 50%
    else if (gRefMult>=centFull[2]) centrality=2; // 50 - 60%
    else if (gRefMult>=centFull[1]) centrality=1; // 60 - 70%
    else if (gRefMult>=centFull[0]) centrality=0; // 70 - 80%
    else centrality = -1;
    //else centrality = 9;
    return centrality;
}
