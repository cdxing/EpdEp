#include <TSystem>
#include <string>
#include <iostream>
#include <fstream>


void Analysis(Int_t nRootFilesToProcess=1, TString FileListName, TString JobIdName)
{

  if (nRootFilesToProcess==0) nRootFilesToProcess=10000000;  // i.e. all of them

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StEpdUtil.so");  // malisa
  gSystem->Load("FlowAnalysis.so");  // note that the flow code requires StEpdUtil library to be loaded first.
  gSystem->Load("StRefMultCorr");

  PicoAnalyzer* pA = new PicoAnalyzer(JobIdName);
  pA->Init();                                  // this also calls init of EP finder

  TString rootFileName;
  TChain* picoDst = new TChain("PicoDst");
  std::ifstream ifs(FileListName,std::ifstream::in);
  if (!ifs){std::cout << "File list " << FileListName << " does not exist!  I quit.\n ";  return;}
  for (int iFile=0; iFile<nRootFilesToProcess; iFile++){
    ifs >> rootFileName;
    if (!(ifs.good())){
      std::cout << "I ran out of files to add after " << iFile << " files\n";
      break;
    }
    std::cout << " - Adding file " << rootFileName.Data() << endl;
    picoDst->Add(rootFileName.Data());
  }
  pA->SetPicoDst(picoDst);
  std::cout << "test 0 "<< std::endl;
  int NeventsToAnalyze = picoDst->GetEntries();
  std::cout << "test 1 "<< std::endl;
  std::cout << "Preparing to analyze " << NeventsToAnalyze << " events...\n";
  for (int ievent=0; ievent<picoDst->GetEntries(); ievent++){
    if (ievent%1000==0) std::cout << "Processed " << ievent << " events - " << 100.0*(double)ievent/(double)NeventsToAnalyze << "% \n";
    pA->Make(ievent);
  }

  pA->Finish();

  cout << "Well hello there!\n";
  cout << "Well, that was fun -- goodbye!!\n";


  delete pA;
  delete picoDst;



  //  chain                            =  new StChain();
  //  StPicoDstMaker*    picoMaker     =  new StPicoDstMaker(2,InputFileList,"picoDst");
  //  MyAnalysisMaker*   AnalysisCode  =  new MyAnalysisMaker(picoMaker) ;  malisa - MUST PASS JOBNAME IN CONSTRUCTOR!!
  //  MyAnalysisMaker*   AnalysisCode  =  new MyAnalysisMaker(picoMaker,JobIdName) ;


  /* no.  this is all bad - malisa
     get the unique identifier in at the very start, then open your files then book your histograms/trees/ntuples/whatever
     that is the most robust way, rather than relying on the user to invoke methods in a particular order.  also, if another
     code (e.g. EventPlane finder) opens TFiles on its own, every ttree/ntuple/whatever MUST be already associated with
     an open TFile, or else there is confusion, and the TTree gets flushed into the most-recently-opened file when it grows
     large enough that root would like to flush it.  this causes .root file corruption that is actually really hard to
     debug!

     TString Name = JobIdName ;  Name.Append(".histograms.root") ;
     TString FemtoDSTName = JobIdName ;  FemtoDSTName.Append(".FemtoDst.root") ;
     TString Name_Broken = JobIdName ;  Name_Broken.Append(".histograms.broken.root") ;
     TString FemtoDSTName_Broken = JobIdName ;  FemtoDSTName_Broken.Append(".FemtoDst.broken.root") ;
     AnalysisCode -> SetOutputFileName(Name) ;       // Name the output file for histograms
     AnalysisCode -> SetBrokenOutputFileName(Name_Broken);
     AnalysisCode -> SetFemtoDstFileName(FemtoDSTName) ;       // Name the output file for histograms
     AnalysisCode -> SetBrokenFemtoDstFileName(FemtoDSTName_Broken) ;       // Name the output file for histograms
  */

  /*
  if( chain->Init()==kStErr ){
    cout<<"chain->Init();"<<endl;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++){

    if(i%1000==0)
    cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) { cout << "Bad return code!" << iret << endl; break;}

    total++;
  }
  */
  /*
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  delete AnalysisCode;
  delete picoMaker;
  delete chain;
  */
}
