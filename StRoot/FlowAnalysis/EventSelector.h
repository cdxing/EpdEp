#ifndef PicoAnalyzer__
#define PicoAnalyzer__

#include "TObject.h"

class TChain;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
class TProfile;
class TProfile2D;
class TNtuple;
class TFile;
class StEpdGeom;
class StBbcGeom;
class TRandom3;
class StPicoEpdHit; 
class StEpdEpFinder;

class PicoAnalyzer : public TObject {
 public:
  PicoAnalyzer();
  ~PicoAnalyzer();

  void SetPicoDst(TChain*);
  short Init();
  short Make(int iEvent);
  short Finish();

 private:

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

  // ntuples
  TNtuple* mQ1vectorNtuple;      // Q1 vectors ring-by-ring. For offline weight optimization
  TNtuple* mQ2vectorNtuple;      // Q2 vectors ring-by-ring. For offline weight optimization

  int mRunId;                         // when this changes, refresh some information.
  short mRunCollisionSystem;



  // internal methods
  int FindCent(int RefMult);   // utility class just giving centrality bin.  Copied directly from Isaac 1 May 2018
  double GetBbcPmtPhi(short PmtId);
  void ReadInSystems();    // reads a text file that idenfies the collision system for every run
  short WhichSystem(int runId);
  void NewRun(int runId);    // invoked when Make() detects that a new run has been loaded
  void FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight);     // fills the histograms used for (a later job's) phi weighting

  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double v1Weight(int CentId, double eta);

  // useful objects kept by PicoAnalyzer
  StEpdGeom* mEpdGeom;
  StBbcGeom* mBbcGeom;
  StEpdEpFinder* mEpFinder;
  TRandom3* mRan;
  char mCollidingSystem[365][500];    // index1=day of year;  index2=run of day




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
  TH1D* mNmipDists[2][12][31];   // nMIP distributions for all tiles
  TH1D* mRefMultHist[9];         // convenient to have refMult distributions for the 9 centrality bins

  TH1D* mAdcDists[2][12][31];   // ADC distributions for all tiles
  TH1D* mTdcDists[2][12][31];   // TDC distributions for all tiles
  TH1D* mTacDists[2][12][31];   // TAC distributions for all tiles


  // 2D histograms
  TH2D* mHisto2D[40];           // miscellaneous 2D histograms
  TH2D* mEpdEwPsi[3][3];        // correlation between East and West EPD event planes : [order n-1][correction level].
  TH2D* mBbcEwPsi[3][3];        // correlation between East and West BBC event planes : [order n-1][correction level].

  TH2D* mEpdEwPsi_midCentral[3][3];        // [order n-1][correction level].  Just to focus on 10-40%, for example
  TH2D* mBbcEwPsi_midCentral[3][3];        // [order n-1][correction level].  Just to focus on 10-40%, for example
  TH2D* mBbcXY[2][16];    // shows locations of EPD hits when a BBC PMT does *not* fire

  TH2D* mPhiWeightOutput[2];    // this is used for "phi weighting"
  TH2D* mPhiAveraged[2];        // this is just used for the normalization of the above. - not ever saved.  Only internal use
  TH2D* mPhiWeightInput[2];     // "phi weighting" correction factors that were calculated and saved in a PREVIOUS run


  TH2D* mTpcEtaPhi;                 // used for weights when calculating Q-vector in TPC (output)
  TH2D* mTpcEtaPhiInput;            // used for weights when calculating Q-vector in TPC (input)

  // 1D profiles

  TProfile* mScalarProduct[10][3][3];   // cent, ewFull, order

  TProfile* mEpdAveCos[3][3];  // [order n-1][correction level].
  TProfile* mBbcAveCos[3][3];  // [order n-1][correction level].

  // flow!
  //  TProfile* mEpdRawVn[9][3][3];                // <cos(phiTile-Psi1)> where Psi1=East/West/Full if the second index is 0/1/2   it is <cos(EastTile-WestEP[1])>. Axis is eta.  One profile for each centrality ID
  //  TProfile* mEpdWeight[9][3];                  // this is only used for a normalization of the above at the end.  It is deleted, so not even written out.

  TH1D* mEpdRawVn[9][3][3];                // <cos(phiTile-Psi1)> where Psi1=East/West/Full if the second index is 0/1/2   it is <cos(EastTile-WestEP[1])>. Axis is eta.  One profile for each centrality ID
  TH1D* mEpdWeight[9][3];                  // this is only used for a normalization of the above at the end.  It is deleted, so not even written out.

  TProfile* mTpcRawV1_FullPsiRotate[9];        // the TPC is also filled in the mEpdRawV1, but that is using just ONE of the EPDs.  This one uses "full event"

  // 2D profiles
  TProfile2D* mProfile2D[40];    // miscellaneous 2D profiles

  //  these are shift correction factors that we MAKE now and write out
  TProfile2D* mEpdShiftOutput_sin[2][3];    // [ew][order-1]
  TProfile2D* mEpdShiftOutput_cos[2][3];    // [ew][order-1]
  TProfile2D* mBbcShiftOutput_sin[2][3];    // [ew][order-1]
  TProfile2D* mBbcShiftOutput_cos[2][3];    // [ew][order-1]

  //  these are shift correction factors that we made before, and USE now
  TProfile2D* mEpdShiftInput_sin[2][3];    // [ew][order-1]
  TProfile2D* mEpdShiftInput_cos[2][3];    // [ew][order-1]
  TProfile2D* mBbcShiftInput_sin[2][3];    // [ew][order-1]
  TProfile2D* mBbcShiftInput_cos[2][3];    // [ew][order-1]
  


  // TFiles (to store histograms and data)
  TFile* mHistoFile;
  TFile* mQ1NtupleFile;
  TFile* mQ2NtupleFile;
  TFile* mCorrectionInputFile;
  TFile* mCorrectionOutputFile;

  TFile* mTpcWeightFile;      // just holds eta-phi 2D histogram
  TFile* mTpcWeightInputFile;
  
  ClassDef(PicoAnalyzer, 1)                     //  Macro for CINT compatability

};


#endif
