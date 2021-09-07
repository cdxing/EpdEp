#ifndef PicoAnalyzer__
#define PicoAnalyzer__

#include "TObject.h"

class TChain;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
class TH3D;
class TH1F;
class TH2F;
class TProfile;
class TProfile2D;
class TProfile3D;
class TNtuple;
class TFile;
class StEpdGeom;
class StBbcGeom;
class TRandom3;
class StPicoEpdHit;
class StEpdEpFinder;
class TLorentzVector;
class StRefMultCorr;
#include "StRefMultCorr/CentralityMaker.h"
#include "TString.h"
#include "TMath.h"

#define _PsiOrderMax 1 //Maximum order of EP to worry about

const Double_t _massKaon     = 0.493677;
const Double_t _massPhi = 1.019461;

class PicoAnalyzer : public TObject {
 public:
  PicoAnalyzer(TString FileNameBase="MikesStuff");
  ~PicoAnalyzer();

  void SetPicoDst(TChain*);
  short Init(char const* TPCweightFileName="TPCWeightFile.root",char const*  TPCShiftFileName="TPCShiftFile.root",char const* EPDPhiWeightFile="EPDcorrection.root");
  short Make(int iEvent);
  short Finish();

 private:

  TString mFileNameBase;

  // parameters relevant to my analysis
  double mNmipQtB;  // ADC value of MIP peak on rings 6-16 (read out thru QT32Bs)
  double mNmipQtC;  // ADC value of MIP peak on rings 1-5  (read out thru QT32Cs)
  double mnMipThreshold;  // low-signal threshold, to cut out noise basically.

  // the data objects
  TChain*   mPicoDst;
  TClonesArray* mEpdHits;
  TClonesArray* mBbcHits;
  TClonesArray* mTracks;
  TClonesArray* mEventClonesArray;  // kind of hilarious that the StPicoEvent is stored as a one-element TClonesArray :-)
  TClonesArray* mTraits;

  // ntuples
  TNtuple* mQ1vectorNtuple;      // Q1 vectors ring-by-ring. For offline weight optimization
  TNtuple* mQ2vectorNtuple;      // Q2 vectors ring-by-ring. For offline weight optimization

  int mRunId;                         // when this changes, refresh some information.
  int mRunEt;                    //Run entry
  short mRunCollisionSystem;



  // internal methods
  int FindTrackId(int Trkch,double TrkVz,double Trketa,double TrkpT );//this method only works when bin numbers are even and non-zero!!!!
  int FindCent(int RefMult);   // utility class just giving centrality bin.  Copied directly from Isaac 1 May 2018
  double GetBbcPmtPhi(short PmtId);
  void ReadInSystems();    // reads a text file that idenfies the collision system for every run
  short WhichSystem(int runId);
  void NewRun(int runId);    // invoked when Make() detects that a new run has been loaded
  //void FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight);     // fills the histograms used for (a later job's) phi weighting
  bool Runlist(int runId);
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double v1Weight(int CentId, double eta);

  // useful objects kept by PicoAnalyzer
  StEpdGeom* mEpdGeom;
  StBbcGeom* mBbcGeom;
  StEpdEpFinder* mEpFinder;
  StRefMultCorr* mRefMultCorr;
  TRandom3* mRan;//seems like RCF only likes TRandom3
  char mCollidingSystem[365][500];    // index1=day of year;  index2=run of day

  static const int mTPCphibin = 80;
  double mTPCphibinwidth = TMath::Pi()/40;

  int Bin_Centrality_01 = 4;
  int Bin_rap = 8;
  // Centrality bin

  // int cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
  // int cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
  // double rap_low_phi[8] = {-1.0, -0.6, -0.3, -0.1, 0., 0.1, 0.3, 0.6};
  // double rap_up_phi[8]  = {-0.6, -0.3, -0.1, 0.,  0.1, 0.3, 0.6, 1.0};
  // TString Centrality_01[4] = {"0080","0010","1040","4080"};

  /*
  //  static const int mNumberOfEpdSubEvents = 6;
  //  double mEpdEtaBbounds[mNumberOfEpdSubEvents+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  static const int mNumberOfEpdSubEvents = 3;
  double mEpdEtaBbounds[mNumberOfEpdSubEvents+1] = {2.0, 3.0, 4.0, 5.0};
  static const int mNumberOfTpcSubEvents = 16;
  TProfile* mEastLongDecorProfile[mNumberOfEpdSubEvents][4];    // second index is order of event plane
  TProfile* mWestLongDecorProfile[mNumberOfEpdSubEvents][4];
  */

  // ---------- Now, my histograms, ntuples, TFiles, etc.  All stuff particular to my analysis
  // 1D histograms
  TH1D* mHisto1D[40];            // miscellaneous 1D histograms

  // 2D histograms
  TH2D* mHisto2D[40];           // miscellaneous 2D histograms

  TH2D* mTPCPhiWeightOutput;    // this is used for "phi weighting" TPC tracks
  TH2D* mTPCPhiAveraged;        // this is just used for the normalization of the above. - not ever saved.  Only internal use
  TH2D* mTPCPhiWeightInput;     // "phi weighting" correction factors that were calculated and saved in a PREVIOUS run
  TH2D* mTPCPhiAveragedInput;

double mpTMin;
double mpTMax;
double mpTMink;
double mpTMaxk;
double mEtaMin;//for phi-weighting, used in FindTrackId() and when MakeWeight
double mEtaMax;////for phi-weighting
int mNPhibin;
double mVzMin;
double mVzMax;
double mEPDMax;
double mEPDthresh;
int mNpTbin;
int mNVzbin;
int mNEtabin;
int mNumberOfTrackTypes;
double mVtxR;
double mDiffVzVPD;
int mNhitsfit;
double mNhitsfitratio;
double mDCAcut;
int mFourierOrder;
int mNTPCSubEvents;
int mNEPDSubEvents;
double mpTbound;//For looking at low/high pT tracks
double mPionSigma;//1/beta sigma
double mKaonSigma;//1/beta sigma
double mProtonSigma;//1/beta sigma
double d_KaonM2low;//1/beta sigma
double d_KaonM2high;//1/beta sigma
double mEtaMaxv1;//for v1 analysis
double mEtaMinv1;//for v1 analysis

//histos for QA
TH1F *href_vz, *hvz_b;
TH1F *href, *hvz, *hbtofYLocal;
TH2F *hvzvpdvz_b, *hvr_b, *hvzvpdvz, *hvr, *hmassvsp, *hdedxvsp, *htofvsref_b, *htofvsref, *hbtofYLocalvsMass2;
TH2F *htofmatchvsref, *htofmatchvsref_b;
TH2F *hbetavsp;
TH2F *h_eta_phi;
TH2F *h_eta_phi_before;
TH1F *h_counter;

TProfile *h_runidvstofmult_b, *h_runidvsrefmult_b;
TProfile *h_runidvstofmult, *h_runidvsrefmult;

TH1F *h_pt, *h_eta_b, *h_eta, *h_nhitfit, *h_nhitmax, *h_nhitratio, *h_dca, *h_phi;

// Kaon PID
TH1D *hist_pt_kaonPlus;
TH1D *hist_eta_kaonPlus;
TH1D *hist_y_kaonPlus;
TH1D *hist_phi_kaonPlus;
TH2D *hist_rap_eta_kaonPlus;
TH2D *hist_pt_y_kaonPlus;
TH2D *hist_pt_eta_kaonPlus;
TH2D *hist_dEdx_kaonPlus;
TH2D *hist_beta_kaonPlus;
TH2D *hist_mass_kaonPlus;

TH1D *hist_pt_kaonMinus;
TH1D *hist_eta_kaonMinus;
TH1D *hist_y_kaonMinus;
TH1D *hist_phi_kaonMinus;
TH2D *hist_rap_eta_kaonMinus;
TH2D *hist_pt_y_kaonMinus;
TH2D *hist_pt_eta_kaonMinus;
TH2D *hist_dEdx_kaonMinus;
TH2D *hist_beta_kaonMinus;
TH2D *hist_mass_kaonMinus;

//histos for checking the phi-weighted eta-phi distribution for the TPC tracks
  TH2D* mEtaPhiDisPhiWeighted[9];//phi-weighted eta phi distribution for 9 centralitites.
  TH2D* mEtaPhiDisRaw[9];
//histos for checking the pT distribution
  TH1F* mPtRaw[9];//9 centralities
  //TH1D* mpTPhiWeighted[9][10];
  //TH1F* mPtRaw;
  TH1D* mVz[9];//Vz dis for nine centralites
  TH1D* mdNdeta[9];//dNdeta for nine centralities

//histos for checking the flattening of the TPC Psi
  TH1D* mTPCPsiDisWeighted[9][3];//Entry is the centrality and "TrackType", e.g sometimes I want to look at low pT or high pT tracks.
  TH1D* mTPCPsiDisShifted[9][3];
//histos for checking the flattening of the Psi_EODfull
  TH1D* mEPDFullPsiWeighted[9][_PsiOrderMax];
  TH1D* mEPDFullPsiShifted[9][_PsiOrderMax];

//TProfile for calculating the resolutions.
  TProfile* mResolution[9][_PsiOrderMax];//Indicies are cent.

//TProfile for vn(TPC) w/ respect to Psi_EPDfull and vn(EPD) w/ respect to Psi_TPC using trucated nMIP(0.3,3.0)
  TProfile* mTPCvn[9][_PsiOrderMax];//cent, Psi order
  TProfile* mEPDvn[2][9][_PsiOrderMax];//EPD vn of E/W EPD, cent, Psi order

//TH3D for measuring dNdphi
  //TH3D* mThreeD[9][2];//centralities, ew
  //TH3D* mThreeD[9][2][_PsiOrderMax][2][16];//centralities, ew, order of harmonics, reference: 0-TPC, 1-the other side of EPD, 16 Vz bins in [-60,60]
  TH3D* mThreeD[2];//centralities, ew, order of harmonics, reference: 0-TPC, 1-the other side of EPD, 16 Vz bins in [-60,60]
  TH3D* mThreeDTile[2][24];//centralities, ew, order of harmonics, reference: 0-TPC, 1-the other side of EPD, 16 Vz bins in [-60,60]

  //Looking at the percentage of differnt nMIPs
  //TH2D* mTwoD[9];//centralities, X:nMIp, Y:RingId (West-1, East-0)
/*
//TProfiles for flow decorrelation analysis
  TProfile* mRefMultvsVz;
  TProfile* mPhiPsievent[2][9];//indicies are e/w and 9 centrality windows
  TProfile* mTPC1TPC2[9];
  TProfile* mPhiPsievent3[2][9];//indicies are e/w and 9 centrality windows
*/
  TProfile3D* mTPCCosShift;
  TProfile3D* mTPCSinShift;

  TProfile3D* mTPCCosShift_Input;
  TProfile3D* mTPCSinShift_Input;

  TProfile2D* mAveEta[9][2]; //nine centralities, EW(ring#, ivz,eta)

//cp TH3D EPDPhiWeight histos from the StEpdEpFinder OUTPUT.
  TH3D* mEPDPhiWeights[2];//e/w

  // TFiles (to store histograms and data)
  TFile* mHistoFile;
  TFile* mTPCWeightFile;
  TFile* mTPCShiftFile;
  TFile* mEPDPhiWeightFile;
  //TFile* mRecenQFile;

  // Phi meson analysis plots
  TH2D *h2px;
  TH2D *h2py;
  TH2D *h2pz;

  TH1D * h_dip_angle;
  TH1D * h_Mass;
  TH1D * h_Mass_rot;
  TH1D *hist_SE_PhiMeson_pT;
  TH1D *hist_SE_PhiMeson_mT;
  TH1D *hist_SE_PhiMeson_rap;
  TH1D *hist_SE_PhiMeson_eta;
  TH2F *h_Mass2;
  TH2F *h_Mass2_rot;
  TH2F *h2_pT_eta;
  TH2F *h2_pT_y;


  TH2D *hist_SE_pt_y_PhiMeson[8];
  TH2D *hist_SE_pt_y_Phi_tight_SigBkg[8];
  TH2D *hist_SE_pt_y_Phi_tight_Bkg[8];
  TH2D *hist_SE_pt_y_Phi_tight_Sig[8];

  TH1F *mHist_SE_InvM_rap_cent[8][4];
  TH1F *mHist_rotation_InvM_rap_cent[8][4];
  TProfile *mProfile_flow_reso_rap_cent[8][4];

  static const int mEPTPCMaxTerm = 6;
  ClassDef(PicoAnalyzer, 1)                     //  Macro for CINT compatability

};


#endif
