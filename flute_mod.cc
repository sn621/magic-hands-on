/* ======================================================================== * \
!
! *
! * This file is part of MARS, the MAGIC Analysis and Reconstruction
! * Software. It is distributed to you in the hope that it can be a useful
! * and timesaving tool in analysing Data of imaging Cerenkov telescopes.
! * It is distributed WITHOUT ANY WARRANTY.
! *
! * Permission to use, copy, modify and distribute this software and its
! * documentation for any purpose is hereby granted without fee,
! * provided that the above copyright notice appears in all copies and
! * that both that copyright notice and this permission notice appear
! * in supporting documentation. It is provided "as is" without expressed
! * or implied warranty.
! *
!
!
!   Authors: Abelardo Moralejo 04/2011 <moralejo@ifae.es>
!            Stefan Klepser    02/2012 <klepser@ifae.es>
!            Markus Gaug       11/2014 <markus.gaug@uab.cat>
!
!   Copyright: MAGIC Software Development, 2000-2014
!
!
\* ======================================================================== */

////////////////////////////////////////////////////////////////////////////////
//
//    flute.cc
//
//    FLUTE: FLUx vs. Time and Energy
//
//           Code for calculating energy spectra and light curves
//           It works only on stereoscopic wobble data (and using theta2 as
//           angular variable)
//
//    OUTPUT:
//           In order to keep compatibility with the unfolding code, we use
//           MHExcessEnergyTheta as the output container for the histograms
//           of event excess vs. energy and vs. zenith angle (actually, only
//           one bin is allowed in zenith angle). But we do not use any of
//           the machinery in that class to obtain the excesses. Instead, the
//           class MCalcExcess is used.
//
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>

#include <TApplication.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TGTab.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNtupleD.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeIndex.h>

#include "MAnalysisProblems.h"
#include "MAnalysisWizard.h"
#include "MArgs.h"
#include "MAstro.h"
#include "MAverageEnergy.h"
#include "MBinning.h"
#include "MCalcDeltaT.h"
#include "MCalcExcess.h"
#include "MContinue.h"
#include "MFillH.h"
#include "MGeomCamMagic.h"
#include "MGraphicsWizard.h"
#include "MEnergyEst.h"
#include "MEnv.h"
#include "MEvtLoop.h"
#include "MF.h"
#include "MFHadAlpha.h"
#include "MH3.h"
#include "MHadAlphaCut.h"
#include "MHadronness.h"
#include "MHAlphaEnergyTheta.h"
#include "MHCTC.h"
#include "MHEffectiveOnTime.h"
#include "MHExcessEnergyTheta.h"
#include "MHFlux.h"
#include "MHillas.h"
#include "MHMcCollectionArea.h"
#include "MHMcEnergyMigration.h"
#include "MLidarInterpolate.h"
#include "MLog.h"
#include "MLogManip.h"
#include "MMarsVersion.h"
#include "MMcCollectionAreaCalc.h"
#include "MParList.h"
#include "MPointingPos.h"
#include "MRawRunHeader.h"
#include "MReadEBLtau.h"
#include "MReadMarsFile.h"
#include "MReadTree.h"
#include "MSrcPosCalc.h"
#include "MStatusDisplay.h"
#include "MStereoPar.h"
#include "MTaskList.h"
#include "MTheta2vsEest.h"
#include "MTransmissionEnergyCorr.h"
#include "MSrcPosCam.h"

using namespace std;

void AddSubsetOfFiles(MReadMarsFile *mrmf, TString InputFiles, TString substring = "");
void ReadRADec(TString fullfilename, Double_t* pointingRA, Double_t* pointingDec);
void RoundTheta2cut (Double_t* theta2Cut, Int_t nBinsTheta2, Double_t maxTheta2);
Bool_t findCuts(TString mcdata, Float_t sizemin,
		Double_t lowE, Double_t upE, Int_t NbinsE,
		Double_t* frac1, Double_t* hadcuts,
		Double_t* frac2, Double_t* theta2cuts,
		Int_t nBinsTheta2, Double_t maxTheta2,
	        vector<Double_t> allowedHadRange, vector<Double_t> allowedTheta2Range,
		Double_t minZd, Double_t maxZd,
		TString nameEnergyEst, TString nameHadronness,
		TString nameStereoPar, MStatusDisplay* display = NULL);

void GetMjdBinLimits(const TH1D& hist, Int_t bin, Double_t* minmjd, Double_t* maxmjd);
void GetNightlyLCbinEdges(MHEffectiveOnTime* hEffTime, TArrayD* lowedge, TArrayD* upedge);
void GetSingleBinLCbinEdges(MHEffectiveOnTime* hEffTime, TArrayD* LCbinlowedge, TArrayD* LCbinupedge);
void FillMatrix(Double_t *x, Int_t nmax, MEnv* environment, const Char_t* flag);
Bool_t Check_offset_source(TString InputFilesRealData, double wobble);
Bool_t Compute_IRF(MHMcCollectionArea *collarea,MBinning *BinningDE, MBinning *BinningImpact, Double_t minZd, Double_t maxZd, MBinning *zdbinning, MBinning *eestbinning, MBinning *etruebinning,
		TString InputFilesMC,TString statusOutName, TString posContainerName,  MStatusDisplay* display, MFHadAlpha *mfhadtheta2cuts, MHadAlphaCut HadTheta2Cuts,MContinue *mfconthadtheta2,
		Bool_t AeffAtmCorr,TString nameEnergyEst,TF1 *assumedspectrum, Bool_t EnergyReCalc, MContinue *filterSize, MContinue*  usercuts, MHMcEnergyMigration *energyMigration,
		MH3 *hdeltaetot_eest,MH3 *hdeltaetot_etrue, TH2D EffOnTimevsAzvsZd, Bool_t IsStandardAnalysis,MHMcCollectionArea *collareaVsEest);

TGraph* tau_vs_log10e;
TF1* intrinsicspectrum;
// Functions below needed for case source redshift > 0
Double_t AbsorbedSpectrum(Double_t *x, Double_t *par) // x = E/GeV;  Provided spectral shape times EBL transmission factor
{
  return par[0]*intrinsicspectrum->Eval(x[0])*exp(-tau_vs_log10e->Eval(log10(x[0])-3., 0, "S"));
}

Double_t AbsorbedSED(Double_t *x, Double_t *par) // x = E/GeV;  Provided spectral shape times EBL transmission factor
{
  return par[0]*x[0]*x[0]*1.e-6*intrinsicspectrum->Eval(x[0])*exp(-tau_vs_log10e->Eval(log10(x[0])-3., 0, "S")); // TeV cm-2 s-1
}


// Explain the usage of flute
static void Usage()
{
  gLog << all << endl;
  gLog << "Usage of flute:" << endl;
  gLog << "   modified flute [options]" << endl << endl;
  gLog << "Options:" << endl;
  gLog << "   -h, --help: show this help" << endl;
  gLog << "   --config=flute.rc : use .rc file for configuration" << endl;
  gLog << "   --debug[=n] : set debug level to 1 [or n]" << endl;
  gLog << "   -b: Batch mode (no graphical output to screen)" << endl;
  gLog << "   -q: quit after finishing" << endl;
  gLog << "   --onlyWobblePos=xxx: Use only wobble position xxx (0, 1,...) to get ON-source data" << endl;
  gLog << endl;
}

MGraphicsWizard * gWizard = new MGraphicsWizard(0.05, 42);
MAnalysisWizard * aWizard = new MAnalysisWizard();
//////////////////////////////////////////////////////////////////////////////////
//
//  The main code starts here
//
int main(int argc, char **argv)
{

  // Read in the command-line arguments:

  MArgs arg(argc, argv, kTRUE);
  if (arg.HasOnly("-h") || arg.HasOnly("--help"))
    {
      Usage();
      return 0;
    }
  TApplication app("flute", &argc, argv);


  // Read in name of configuration .rc file:
  TString defaultConfigFile = gSystem->Getenv("MARSSYS")+TString("/mrcfiles/")+TString("flute.rc");
  TString ConfigFile = arg.GetStringAndRemove("--config=", defaultConfigFile);
  ConfigFile = gSystem->ExpandPathName(ConfigFile.Data());
  if (gSystem->AccessPathName(ConfigFile, kFileExists))
    {
      gLog << err << "Sorry, " << ConfigFile << " doesn't exist." << endl;
      return MAnalysisProblems::kWrongConfigFile;
    }
  else
    gLog << inf << "Using configuration file " << ConfigFile.Data() << ")." << endl;
  MEnv* env = new MEnv(ConfigFile);




  // Get the name of the log file to be written
  if(!arg.HasOption("--log="))
    {
      TString defaultLogName = ConfigFile;
      defaultLogName.ReplaceAll(".rc", ".log");
      defaultLogName.Insert(defaultLogName.Last('/')+1, "Log_");
      gLog.SetOutputFile(defaultLogName.Data(), 1);
      gLog << inf << "Setting output log file to default (" << defaultLogName.Data() << ")." << endl;
    }
  else
    {
      gLog << inf << "Setting output log file to " << arg.GetString("--log=") << endl;
      gLog.Setup(arg);
    }

  TString statusOutName = ConfigFile;
  statusOutName.ReplaceAll(".rc", ".root");
  statusOutName.Insert(statusOutName.Last('/')+1, "Status_");
  TString matrixOutName = statusOutName;
  matrixOutName.ReplaceAll("Status_", "Output_");

  // Should we obtain the spectrum using all data, from all wobble pointings, or from one selected pointing?
  Int_t processWobblePos = -1; // This default (-1) means all wobble pointings will be used to calculate the spectrum
  // One can use the command-line option below to select one of the wobble pointings to calculate the spectrum (the bckg will
  // be obtained from the other wobble pointings):
  if(arg.HasOption("--onlyWobblePos="))
    processWobblePos = arg.GetIntAndRemove("--onlyWobblePos=");


  // Setup the global root debugging
  gDebug = arg.HasOption("--debug=") ? arg.GetIntAndRemove("--debug=") : 0;
  if (gDebug == 0 && arg.HasOnlyAndRemove("--debug"))
    gDebug=1;
  // Search for the "quit" and "batch" options:
  const Bool_t  kQuit  = arg.HasOnlyAndRemove("-q");
  const Bool_t  kBatch = arg.HasOnlyAndRemove("-b");
  if ((!gROOT->IsBatch() && !gClient) || (gROOT->IsBatch() && !kBatch))
    {
      gLog << err << "Cannot open display... maybe your DISPLAY variable is not set correctly!" << endl;
      return -1;
    }

  const Bool_t kShowCTC = arg.HasOnlyAndRemove("--showctc");

  if (arg.GetNumArguments() > 0 || arg.GetNumOptions() > 0)
    {
      gLog << err << "Unknown argument!" << endl;
      arg.Print();
      Usage();
      return MAnalysisProblems::kWrongArguments;
    }
  // Create the status display to plot all results:
  MStatusDisplay *display = new MStatusDisplay(1024, 768);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleAlign(13);
  gWizard->FixStyle(gStyle);


  // Read the names of the input files (real data):
  TString InputFilesRealData;
  if ( ! env->GetRequiredValue("flute.data", InputFilesRealData) )
    return MAnalysisProblems::kWrongConfigFile;

  // Read the names of the input files (Monte Carlo):
  TString InputFilesMC;
  if ( ! env->GetRequiredValue("flute.mcdata", InputFilesMC) )
    return MAnalysisProblems::kWrongConfigFile;



  // Read the sky position at which gamma-ray flux has to be calculated:
  Double_t SourceRA = env->GetValue("flute.sourceRa", NAN);
  Double_t SourceDec = env->GetValue("flute.sourceDec", NAN);

  if (TMath::IsNaN(SourceRA) || TMath::IsNaN(SourceDec))
    {
      vector<Double_t> sourceRaDec = aWizard->GetSourceRaDec(InputFilesRealData);
      if(sourceRaDec.size() == 2)
	{
	  SourceRA  = sourceRaDec.at(0);
	  SourceDec = sourceRaDec.at(1);
	}
    }

  gLog << inf << endl << "Source coordinates:  RA = " << SourceRA << " hours,   Dec = " << SourceDec << " degree" << endl;


  // Check what is the requested background calculation mode. Available now:
  //    0 -> off from wobble partner;  3 -> pulsar mode
  Int_t bckgMode = env->GetValue("flute.bckgMode", 0);

  Bool_t DoLightCurve = kTRUE;

  Int_t numSimultaneousBgPositions = env->GetValue("flute.numSimultaneousBgPositions", 3);

  // Fixed factor by which the # of off events will be multiplied before subtracting them from the on events.
  // (e.g. in pulsar mode, or when onSelection and offSelection are used).
  Float_t fixedBckgNormalization = env->GetValue("flute.fixedBckgNormalization", -1.);

  TString onSelection  = env->GetValue("flute.onSelection", "");
  TString offSelection = env->GetValue("flute.offSelection", "");
  // Cuts to be applied (respectively) to select on and off REAL data events. Use sparingly!!
  // Compulsory in pulsar mode (bckgMode=3), in which all events are from the on-region around sky position
  // SourceRA, SourceDec

  TString ebltemplate = env->GetValue("flute.EBLModel", "D11");
  gLog << inf << endl << "EBL model selected: " << ebltemplate.Data() << endl;

  if (onSelection != "" || offSelection != "")
    {
      if (fixedBckgNormalization < 0.)  // e.g. if setting is missing in .rc file
	{
	  gLog << err <<  "ERROR: rc file settings flute.onSelection and flute.offSelection MUST be accompanied"
	       << endl << "by a user-calculated off-normalization factor set through flute.fixedBckgNormalization: xx"
	       << endl << "Exiting..." << endl << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
    }

  Float_t fluxCorrectionFactor = env->GetValue("flute.fluxCorrectionFactor", 1.);
  // Overall flux correction factor. Use sparingly! It is only to obtain the right flux when using onSelection to remove part
  // of the ON data (like Crab data around the pulsar peaks via phase cuts).
  if (fluxCorrectionFactor != 1.)
    gLog << warn << endl << "NOTE: a user-selected correction factor " << Form("%.4f", fluxCorrectionFactor) << " will be applied to all calculated fluxes!" << endl;


  if (bckgMode == 3  && (onSelection == "" || offSelection == "") )
    {
      gLog << err << "ERROR: Background mode = 3 requires that onSelection and offSelection filters are specified in the .rc file!" << endl;
      return MAnalysisProblems::kWrongConfigFile;
    }
  else if (bckgMode == 0  && (onSelection != "" || offSelection != ""))
    {
      gLog << err << "ERROR: Background mode = 0 is incompatible with onSelection and offSelection filters! Fix your .rc file!" << endl;
      return MAnalysisProblems::kWrongConfigFile;
    }

  if (fixedBckgNormalization  != -1. && onSelection == "" && offSelection == "")
    {
      gLog << err << "ERROR: fixedBckgNormalization is only for use with onSelection and/or offSelection filters! Fix your .rc file!" << endl;
      return MAnalysisProblems::kWrongConfigFile;
    }

  // Read in the name of the container which should be used to get the event direction:
  TString posContainerName;
  if ( ! env->GetRequiredValue("flute.posContainer", posContainerName) )
    return MAnalysisProblems::kWrongConfigFile;
  gLog << inf << endl << posContainerName << " will be used to obtain the events' directions" << endl;


  Bool_t ReDoEnergyAverage = env->GetValue("flute.ReDoEnergyAverage", kFALSE);;
  // If kTRUE, it will re-do the energy averaging, both for MC and data (to solve issue with wrong energy average
  // in melibea files created with  Mars <= V2-17-0)
  // NOTE: this will not be applied in the findCuts routine - differences in the final energy are not large (except for few events),
  // and the small changes that would result in the cut values are irrelevant for the final spectrum.

  Bool_t EnergyReCalc = ReDoEnergyAverage;  // Force energy recalculation if user sets ReDoEnergyAverage to TRUE

  // AM, 20180725: took out the light/energy scaling options below. They can now
  // be done in fold (or fitebl), in a more reliable way.
  //
  // Read scaling factor for Eest of real data (to test absolute energy scale):
  //  Float_t eest_scaling_factor = env->GetValue("flute.EestScaling", 1.);

  // Read scaling factor for M1's Eest of real data (to correct for possible miscalibrations)
  //  Float_t eest1_scaling_factor = env->GetValue("flute.M1EestScaling", 1.);
  // Read scaling factor for M2's Eest of real data (to correct for possible miscalibrations)
  //  Float_t eest2_scaling_factor = env->GetValue("flute.M2EestScaling", 1.);
  Float_t eest1_scaling_factor = 1.;
  Float_t eest2_scaling_factor = 1.;
  Float_t eest_scaling_factor = 1.;

  // Now set EnergyReCalc to TRUE also if any of the user-configurable scaling factors is not 1:
  if (fabs(1.-eest_scaling_factor)>1.e-6 || fabs(1.-eest1_scaling_factor)>1.e-6 || fabs(1.-eest2_scaling_factor)>1.e-6)
    EnergyReCalc = kTRUE; // will recalculate the estimated energy of real data using the above scaling factors

  if (EnergyReCalc)
    gLog << inf << "I will re-calculate the average of M1 and M2 energies for each event, both for real data and MC" << endl;


  if (fabs(eest1_scaling_factor-1.)>1.e-6 || fabs(eest_scaling_factor-1.)>1.e-6)
    {
      gLog << inf << endl << "NOTE: M1's Eest of real data will be scaled by a factor "
	   << eest1_scaling_factor << "*" << eest_scaling_factor << endl << endl;
    }
  if (fabs(eest2_scaling_factor-1.)>1.e-6 || fabs(eest_scaling_factor-1.)>1.e-6)
    {
      gLog << inf << endl << "NOTE: M2's Eest of real data will be scaled by a factor "
	   << eest2_scaling_factor << "*" << eest_scaling_factor << endl << endl;
    }

  Bool_t AeffAtmCorr    = env->GetValue("flute.AeffAtmCorr",kFALSE);
  Bool_t AtmCorrDisplay = env->GetValue("flute.AtmCorrDisplay",kFALSE);
  Bool_t LIDARRecalibrate = env->GetValue("flute.LIDARRecalib",kFALSE);
  Bool_t LIDARUseGDAS     = env->GetValue("flute.LIDARUseGDAS",kFALSE);

  if (LIDARUseGDAS && !LIDARRecalibrate)
    {
      gLog << warn << "LIDAR with use of GDAS has been chosen, but not LIDAR re-calibration. Will re-calibrate the LIDAR anyhow!" << endl;
      LIDARRecalibrate = kTRUE;
    }

  Bool_t StoreMcEvents  = env->GetValue("flute.StoreMcEvents",kFALSE);
  Bool_t McRingAnalysis  = env->GetValue("flute.McRingAnalysis",kFALSE);
  vector<Double_t> RingMcEdges;
  vector<MHMcCollectionArea*> collarea_per_ring;
  int NRing=0;
  if(McRingAnalysis){
    RingMcEdges = env->GetNumbers("flute.RingMcEdges");
		if (!RingMcEdges.size())
	    {
	      gLog << err << "Sorry, you didn't give any ring to compute the IRF from the diffuse MC" << endl;
	      return MAnalysisProblems::kWrongConfigFile;
	    }
    NRing = RingMcEdges.size() - 1;
  }

  //MHMcCollectionArea collarea_per_ring[NRing];
  // AM 20180725: commented out the size scaling options. The modification of
  // the absolute light scale can be done now in fold, in a more precise way
  // Read scaling factor for M1's and M2's Size of real data (to correct for possible miscalibrations):
  // Float_t size1_scaling_factor = env->GetValue("flute.M1SizeScaling", 1.);
  // Float_t size2_scaling_factor = env->GetValue("flute.M2SizeScaling", 1.);

  // Global rescaling factor of Size of real data (of both telescopes - just multiplies the telescope-specific factors defined above)

  // Float_t size_scaling_factor = env->GetValue("flute.SizeScaling", 1.);
  // size1_scaling_factor *= size_scaling_factor;
  // size2_scaling_factor *= size_scaling_factor;

  Float_t size1_scaling_factor = 1.;
  Float_t size2_scaling_factor = 1.;

  // Read in the binning in estimated energy:
  Int_t    nBinsEnergyEst = -1;
  Double_t minEnergyEst   = NAN;
  Double_t maxEnergyEst   = NAN;
  if ( ! env->GetRequiredValue("flute.nBinsEnergyEst", nBinsEnergyEst) )
    return MAnalysisProblems::kWrongConfigFile;
  if ( ! env->GetRequiredValue("flute.minEnergyEst", minEnergyEst) )
    return MAnalysisProblems::kWrongConfigFile;
  if ( ! env->GetRequiredValue("flute.maxEnergyEst", maxEnergyEst) )
    return MAnalysisProblems::kWrongConfigFile;
  // Do not allow too low maximum energy, or too high minimum energy; it would cut out events, later affecting the unfolding!
  if (minEnergyEst > 10.)
    {
      gLog << err << endl << "Sorry, the maximum allowed minEnergyEst is 10 GeV!" << endl;
      return -1;
    }
  if (maxEnergyEst < 3.e4)
    {
      gLog << err << endl << "Sorry, the minimum allowed maxEnergyEst is 3e4 GeV!" << endl;
      return -1;
    }
  gLog << inf << endl << "(Log10) binning in estimated energy (nbins, Emin, Emax): "
       << nBinsEnergyEst << ", " << minEnergyEst << ", " << maxEnergyEst << "  GeV" << endl;

  MBinning eestbinning("BinningEest"); // The naming is needed so that task MHMcCollectionArea can find it in the task list
  eestbinning.SetEdgesLog(nBinsEnergyEst, minEnergyEst, maxEnergyEst);

  TString  nameEnergyEst = env->GetValue("flute.NameEnergyEst","MEnergyEst");

  if (AeffAtmCorr && nameEnergyEst != "MEnergyEstAtmCorr")
    {
      TString checkVal = env->GetValue("flute.NameEnergyEst", "");
      if (!checkVal.IsNull() && nameEnergyEst != "MEnergyEst")
	{
	  gLog << warn << endl << "WARNING: You chose atmospheric corrections, but set the name of the MEnergyEst container explicitly to \"" << nameEnergyEst << "\", instead of MEnergyEstAtmCorr or MEnergyEst" << endl;
	  gLog << err << "This feature is not yet implemented!" << endl;
	  gLog << err << "If you like to use a different energy estimator than the default one, you cannot use atmospheric corrections yet. " << endl;
	  gLog << err << "In case you think, this feature is absolutely necessary, please contact Julian Sitarek and Markus Gaug in order to provide a fix" << endl;
	  exit(2);
	}
      gLog << warn << endl << "You chose atmospheric corrections, will modify the name of MEnergyEst to MEnergyEstAtmCorr" << endl;
      nameEnergyEst = "MEnergyEstAtmCorr";
    }

  gLog << inf << "Using NameEnergyEst: " << nameEnergyEst << "END" << endl;

  // Now set a binning in true energy, which will be used in the energy migration and in the collection area object "collarea"
  // for later use by the unfolding. The number of bins has to be smaller than that of the estimated energy binning.
  Double_t estTrueFactor = env->GetValue("flute.estTrueFactor", sqrt(2.));
  Int_t nBinsEnergyTrue = (Int_t)(nBinsEnergyEst/estTrueFactor);
  MBinning etruebinning("BinningE"); // The naming is needed so that task MHMcCollectionArea can find it in the task list
  etruebinning.SetEdgesLog(nBinsEnergyTrue, minEnergyEst, maxEnergyEst);

  // Read in the energy range for the light curve, and find the corresponding Eest bins
  Int_t EminLCbin, EmaxLCbin;
  Float_t EminLC = env->GetValue("flute.EminLC", 300.);
  Float_t EmaxLC;
  TString EmaxLCstr = env->GetValue("flute.EmaxLC", "inf");

  if (EmaxLCstr == "inf")
    {
      // We will modify the user-selected Eest binning and the corresponding Etrue binning
      // (also for the spectrum in order to have the light curve exactly above the requested
      // EminLC energy, instead of a modified value.
      // We will keep the same number of bins, and the same bin widths, and just shift everything
      // so that EminLC corresponds to the border between two Eest bins:

      Double_t* ed = eestbinning.GetEdges();
      Double_t factor = EminLC/ed[eestbinning.FindRoundEdge(EminLC)];
      minEnergyEst*= factor;
      maxEnergyEst *= factor;  // With these new limits, EminLC will be right between two Eest bins
      // New binnings:
      eestbinning.SetEdgesLog(nBinsEnergyEst, minEnergyEst, maxEnergyEst);
      etruebinning.SetEdgesLog(nBinsEnergyTrue, minEnergyEst, maxEnergyEst);
      if ( fabs(factor-1.) > 1.e-5)
	{
	  gLog << endl;
	  gLog << inf << "Energy binning for the spectrum was modified to match the chosen energy cut for the LC" << endl;
	  gLog << inf << "New minEnergyEst: " << minEnergyEst << " GeV;  New maxEnergyEst: " << maxEnergyEst << endl << endl;
	}

      EminLCbin = eestbinning.FindRoundEdge(EminLC)+1;
      EmaxLCbin = nBinsEnergyEst;
    }
  else  //
    {
      EminLCbin = eestbinning.FindRoundEdge(EminLC)+1;
      Float_t EmaxLC = EmaxLCstr.Atof();
      EmaxLCbin = eestbinning.FindRoundEdge(EmaxLC);
    }

  // Updated values:
  EminLC = eestbinning.GetEdges()[EminLCbin-1];
  EmaxLC = eestbinning.GetEdges()[EmaxLCbin];

  gLog << inf << endl << "Light curve will be produced in range "
       << EminLC << " to " << EmaxLC << " GeV, with Eest bins "
       << EminLCbin << " to " << EmaxLCbin << endl << endl;


  // Binning in Azimuth
  Int_t nBinsAz = env->GetValue("flute.nBinsAz", 1);


  Double_t StartMJD = env->GetValue("flute.StartMJD", NAN);
  Double_t EndMJD = env->GetValue("flute.EndMJD", NAN);

  MContinue* MJDrangeSelection = 0;
  if (!TMath::IsNaN(StartMJD) && !TMath::IsNaN(EndMJD))
    {
      MJDrangeSelection = new MContinue(Form("(MTime_1.fMjd<%.6lf) || (MTime_1.fMjd>%.6lf)", StartMJD, EndMJD));
      MJDrangeSelection->SetName("MJDrangeSelection");

      if (bckgMode != 1 && bckgMode != 3)
	{
	  gLog << err << "Error: MJD range selection option is currently only available with Bacground mode = 1 or 3!  Exiting..." << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
    }


  // Read in the binning in theta^2
  Int_t nBinsTheta2;
  Double_t maxTheta2;
  if ( ! env->GetRequiredValue("flute.nBinsTheta2", nBinsTheta2) )
    return MAnalysisProblems::kWrongConfigFile;
  if ( ! env->GetRequiredValue("flute.maxTheta2", maxTheta2) )
    return MAnalysisProblems::kWrongConfigFile;


  // Read in the theta^2 range in which we want to normalize the OFF histograms to the ON ones:
  Double_t normRangeTheta2Min, normRangeTheta2Max;
  if ( ! env->GetRequiredValue("flute.normRangeTheta2Min", normRangeTheta2Min) )
    return MAnalysisProblems::kWrongConfigFile;
  if ( ! env->GetRequiredValue("flute.normRangeTheta2Max", normRangeTheta2Max) )
    return MAnalysisProblems::kWrongConfigFile;


  // Read in the zenith angle range:
  Double_t minZd, maxZd;
  if ( ! env->GetRequiredValue("flute.minZd", minZd) )
    return MAnalysisProblems::kWrongConfigFile;
  if ( ! env->GetRequiredValue("flute.maxZd", maxZd) )
    return MAnalysisProblems::kWrongConfigFile;
  gLog << inf << endl << "Zenith distance range: " << minZd << " to " << maxZd << "  degree" << endl;

  // Read in the size cut (to be applied to both CTs):
  Double_t minSize;
  if ( ! env->GetRequiredValue("flute.minSize", minSize) )
    return MAnalysisProblems::kWrongConfigFile;
  gLog << inf << endl << "Minimum size (in both telescopes): " << minSize << " phe-" << endl;
  if (fabs(size1_scaling_factor-1.)>1.e-6)
    {
      gLog << inf << endl << "NOTE: M1's size of real data will be scaled by a factor " << size1_scaling_factor << endl;
      gLog << inf << endl << "      (just for the size cuts!)"<< endl << endl;
    }
  if (fabs(size2_scaling_factor-1.)>1.e-6)
    {
      gLog << inf << endl << "NOTE: M2's size of real data will be scaled by a factor " << size2_scaling_factor << endl;
      gLog << inf << endl << "      (just for the size cuts!)"<< endl << endl;
    }

  // Create the filters for tasklists:
  MContinue filterecorr("MEnergyEstAtmCorr.fEnergy<0");
  MContinue filterSize(Form("(MHillas_1.fSize<%f) || (MHillas_2.fSize<%f)", minSize, minSize), "filterSize");
  MContinue filterSizeRealData(Form("(MHillas_1.fSize<%f) || (MHillas_2.fSize<%f)",
				    minSize/size1_scaling_factor, minSize/size2_scaling_factor),
			       (fabs(size1_scaling_factor-1.)>1.e-6 || fabs(size2_scaling_factor-1.)>1.e-6?
				"filterSizeScaled" : "filterSize"));

  TString LCbinningOption = env->GetValue("flute.LCbinning", "night-wise");

  TArrayD LCbinlowedge;  // Will contain lower ends of the MJD "bins" for the light curve
  TArrayD LCbinupedge;   // Will contain upper ends of the MJD "bins" for the light curve
  // If the user has selected a custom LC binning, read in bin limits:
  if (LCbinningOption == "custom")
    {
      TString str;
      if ( ! env->GetRequiredValue("flute.LCbinlowedge", str) )
	{
	  gLog << err << "Custom LC binning chosen, but no bin low edges provided! Exiting." << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
      Int_t numvalues = str.Tokenize(",")->GetEntries();
      Double_t* val = new Double_t[numvalues];
      FillMatrix(val, numvalues, env, "flute.LCbinlowedge");
      LCbinlowedge.Set(numvalues, val);

      if ( ! env->GetRequiredValue("flute.LCbinupedge", str) )
	{
	  gLog << err << "Custom LC binning chosen, but no bin up edges provided! Exiting." << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
      if (numvalues != str.Tokenize(",")->GetEntries())
	{
	  gLog << err << "Custom LC binning: number of bin low edges does not match number of bin up edges!" << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
      FillMatrix(val, numvalues, env, "flute.LCbinupedge");
      LCbinupedge.Set(numvalues, val);
    }


  // User - defined cuts. To be used later in task lists:
  TString UserCuts = env->GetValue("flute.UserCuts", "");

  MContinue* usercuts = 0;
  if (UserCuts != "")
    {
      usercuts = new MContinue(UserCuts);
      usercuts->SetName("UserCuts");
      usercuts->SetInverted();
    }



  // Hadronness and theta2 cuts.
  // There are three possibilities:

  Bool_t FindCutsFromEfficiency = kFALSE;
  if(!env->GetRequiredValue("flute.FindCutsFromEfficiency", FindCutsFromEfficiency))
    return MAnalysisProblems::kWrongConfigFile;

  // True if you want to compute the hadronness cut from efficiency but having the theta2 cut fix.
  Bool_t FindHadCutsFromEfficiencyTheta2Fix = env->GetValue("flute.FindHadCutsFromEfficiencyTheta2Fix", kFALSE);
  vector<Double_t> allowedHadRange, allowedTheta2Range;
  if (FindCutsFromEfficiency || FindHadCutsFromEfficiencyTheta2Fix){
    allowedHadRange   = env->GetNumbers("flute.allowedHadRange");
    if(allowedHadRange.empty()){
      allowedHadRange.push_back(0.15);
      allowedHadRange.push_back(0.95);
    }
    if (! FindHadCutsFromEfficiencyTheta2Fix)
    {
    	allowedTheta2Range = env->GetNumbers("flute.allowedTheta2Range");
    	if(allowedTheta2Range.empty()){
    		allowedTheta2Range.push_back(0.01);
    		allowedTheta2Range.push_back(0.2);
    	}
    	gLog << inf << endl << Form("Allowed range in Theta2: [%.4f, %.4f]\nAllowed range in Had: [%.4f, %.4f]",
    		allowedTheta2Range.at(0),  allowedTheta2Range.at(1),  allowedHadRange.at(0), allowedHadRange.at(1)) << endl;
    }
    else{
    	gLog << inf << endl << Form("Allowed range in Had: [%.4f, %.4f]", allowedHadRange.at(0), allowedHadRange.at(1)) << endl;
    }
  }

  vector<Double_t> hadCuts;
  vector<Double_t> theta2Cuts;

  vector<Double_t> thetaMin2Cuts;
  // ^^ Lower theta2 cut value, to study "ring" region of interest, e.g. dark matter halo excluding central astrophysical source

  if (FindCutsFromEfficiency || FindHadCutsFromEfficiencyTheta2Fix)
    {
      Double_t HadEffi, Theta2Effi;
      if ( ! env->GetRequiredValue("flute.HadEffi", HadEffi) )
	return MAnalysisProblems::kWrongConfigFile;
      if (FindHadCutsFromEfficiencyTheta2Fix)
      {
    	  // Extract the fix Theta2 cut
    	  theta2Cuts = env->GetNumbers("flute.theta2Cuts");
    	  if(theta2Cuts.size() > 0)
    	  {
    		  if(theta2Cuts.size() != UInt_t(nBinsEnergyEst))
    		  {
    			  gLog << err << "ERROR: Number of theta2Cuts must be the same as number of energy bins!" << endl;
    			  return MAnalysisProblems::kWrongConfigFile;
    		  }
    		  for(vector<Double_t>::iterator it = theta2Cuts.begin(); it != theta2Cuts.end(); it++)
    			  RoundTheta2cut (&(*it), nBinsTheta2, maxTheta2);
    	  }
    	  else  // The same theta2 cut for all energy bins
    	  {
    		  Double_t theta2Cut = NAN;
    		  if ( ! env->GetRequiredValue("flute.theta2Cut", theta2Cut) )
    			  return MAnalysisProblems::kWrongConfigFile;
    		  RoundTheta2cut (&theta2Cut, nBinsTheta2, maxTheta2); // Round the cut so that it matches a bin edge in the histograms
    		  theta2Cuts.assign(nBinsEnergyEst, theta2Cut);
    		  gLog << inf << endl << "Global theta2 cut: < " << theta2Cut << endl;
    	  }
          gLog << inf << "NOTE that the theta^2 cut(s) may have been rounded to match the closest bin edge in the theta^2 histograms!" << endl << endl;
          gLog << inf << endl << "User-selected Hadronness cut efficiencies: " << HadEffi << endl;
      }
      else
      {
    	  if ( ! env->GetRequiredValue("flute.Theta2Effi", Theta2Effi) )
    		  return MAnalysisProblems::kWrongConfigFile;
    	  gLog << inf << endl << "User-selected cut efficiencies: " << HadEffi << " (Hadronness cut)  and "
    	  	   << Theta2Effi << " (Theta2 cut)" << endl;
      }

      Double_t* frac1 = new Double_t[nBinsEnergyEst];
      Double_t* frac2 = new Double_t[nBinsEnergyEst];
      for (Int_t ebin = 0; ebin < nBinsEnergyEst; ebin++)
	{
	  frac1[ebin] = HadEffi;
	  // If True will calculate only the hadronness cut from efficiency
	  if (FindHadCutsFromEfficiencyTheta2Fix)
		  frac2[ebin] = -1;
	  else
		  frac2[ebin] = Theta2Effi;
	}

      Double_t* had = new Double_t[nBinsEnergyEst];
      Double_t* th2 = new Double_t[nBinsEnergyEst];

      // Find the hadronness and theta2 cuts which, for every Eest bin, result in the requested efficiency:
      // Here the nameEnergyEst string must not be used, since findCuts loops only over MC events
      if (! findCuts(InputFilesMC, minSize,
		     minEnergyEst, maxEnergyEst, nBinsEnergyEst, frac1, had,
		     frac2, th2, nBinsTheta2, maxTheta2, allowedHadRange, allowedTheta2Range, minZd, maxZd,
		     AeffAtmCorr ? TString("MEnergyEst") : nameEnergyEst, TString("MHadronness"),
		     posContainerName, display))
	return MAnalysisProblems::kWrongConfigFile;

      hadCuts.resize(nBinsEnergyEst);
      if (! FindHadCutsFromEfficiencyTheta2Fix)
    	  theta2Cuts.resize(nBinsEnergyEst);

      for (Int_t ebin = 0; ebin < nBinsEnergyEst; ebin++)
	{
	  hadCuts[ebin] = had[ebin];
	  if (! FindHadCutsFromEfficiencyTheta2Fix)
		  theta2Cuts[ebin] = th2[ebin];
	}
    }
  else  // User-determined explicit hadronness and theta2 cuts (either a different value per Eest bin, or the same for all of them)
    {
      // Read in hadronness and theta2 cuts, for the case in which they are fixed energy-wise by the user:
      hadCuts = env->GetNumbers("flute.hadCuts");
      if(hadCuts.size() > 0)
	{
	  if(hadCuts.size() != UInt_t(nBinsEnergyEst))
	    {
	      gLog << err << "ERROR: Number of hadCuts must be the same as number of energy bins!" << endl;
	      return MAnalysisProblems::kWrongConfigFile;
	    }
	}
      else  // The same hadronness cut for all energy bins
	{
	  Double_t hadCut = NAN;
	  if ( ! env->GetRequiredValue("flute.hadCut", hadCut) )
	    return MAnalysisProblems::kWrongConfigFile;
	  hadCuts.assign(nBinsEnergyEst, hadCut);
	  gLog << inf << endl << "Global hadronnes cut: < " << hadCut << endl;
	}

      theta2Cuts = env->GetNumbers("flute.theta2Cuts");
      if(theta2Cuts.size() > 0)
	{
	  if(theta2Cuts.size() != UInt_t(nBinsEnergyEst))
	    {
	      gLog << err << "ERROR: Number of theta2Cuts must be the same as number of energy bins!" << endl;
	      return MAnalysisProblems::kWrongConfigFile;
	    }
	  for(vector<Double_t>::iterator it = theta2Cuts.begin(); it != theta2Cuts.end(); it++)
	    RoundTheta2cut (&(*it), nBinsTheta2, maxTheta2);
	}
      else  // The same theta2 cut for all energy bins
	{
	  Double_t theta2Cut = NAN;
	  if ( ! env->GetRequiredValue("flute.theta2Cut", theta2Cut) )
	    return MAnalysisProblems::kWrongConfigFile;
	  RoundTheta2cut (&theta2Cut, nBinsTheta2, maxTheta2); // Round the cut so that it matches a bin edge in the histograms
	  theta2Cuts.assign(nBinsEnergyEst, theta2Cut);
	  gLog << inf << endl << "Global theta2 cut: < " << theta2Cut << endl;
	}

      thetaMin2Cuts = env->GetNumbers("flute.thetaMin2Cuts");
      Double_t thetaMin2Cut = NAN;
      if(thetaMin2Cuts.size() > 0)
	{
	  if(thetaMin2Cuts.size() != UInt_t(nBinsEnergyEst))
	    {
	      gLog << err << "ERROR: Number of thetaMin2Cuts must be the same as number of energy bins!" << endl;
	      return MAnalysisProblems::kWrongConfigFile;
	    }
	  for(vector<Double_t>::iterator it = thetaMin2Cuts.begin(); it != thetaMin2Cuts.end(); it++)
	    RoundTheta2cut (&(*it), nBinsTheta2, maxTheta2);
	}
      else // The same lower theta2 cut for all energy bins
	{
	  thetaMin2Cut = env->GetValue("flute.thetaMin2Cut", -1.);
	  if (thetaMin2Cut > 0.)
	    RoundTheta2cut (&thetaMin2Cut, nBinsTheta2, maxTheta2); // Round the cut so that it matches a bin edge in the histograms

	  thetaMin2Cuts.assign(nBinsEnergyEst, thetaMin2Cut);
	  gLog << inf << endl << "Global thetaMin2 cut: > " << thetaMin2Cut << endl;
	}

      gLog << inf << "NOTE that the theta^2 cut(s) may have been rounded to match the closest bin edge in the theta^2 histograms!" << endl << endl;
    }


  tau_vs_log10e = new TGraph;
  tau_vs_log10e->SetName("tau_vs_log10e");
  Double_t SourceRedshift = 0.;
  if (env->Defined("flute.SourceRedshift"))
    SourceRedshift = env->GetValue("flute.SourceRedshift", 0.);

  if (SourceRedshift > 0.)
    {
      MReadEBLtau eblreader;
      TGraph2D tau_E_z;
      if (!eblreader.ReadModel(ebltemplate, &tau_E_z))
	{
	  gLog << "Could not obtain the EBL optical depth. Exiting!" << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
      // Optical depths will be stored in global TGraph tau_vs_log10e

      eblreader.GetTauVsE(tau_vs_log10e, SourceRedshift, &tau_E_z);  // NOTE!  Energy in TeV!

      gLog <<"Read the EBL optical depth vs energy for redshift " << SourceRedshift << " from model " << ebltemplate.Data() << endl << endl;
    }


  // Assumed differential g-ray photon index (intrinsic) for effective area calculations
  Double_t spectralIndex = env->GetValue("flute.spectralIndex", 2.6);

  TString spec;
  if (!env->Defined("flute.AssumedSpectrum"))
    spec = Form("pow(x/500.,-%f)", spectralIndex);
  else // Read formula describing spectral shape, as provided by the user:
    spec = env->GetValue("flute.AssumedSpectrum", "");

  TF1* assumedspectrum;

  if (SourceRedshift > 0.)  // Need to apply the EBL absorption to obtain the observed spectrum
    {
      // Intrinsic spectral SHAPE:
      intrinsicspectrum  = new TF1("IntrinsicSpectralShape", spec.Data(), 10., 1.e5);
      if (intrinsicspectrum->IsZombie())
	{
	  gLog << err << "Wrong function provided in flute.AssumedSpectrum !  Exiting." << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
      // AbsorbedSpectrum is a C-style function (see beginning of code) which uses intrinsicspectrum and the EBL optical depth.

      assumedspectrum = new TF1("Spectrum", AbsorbedSpectrum, 10., 1.e5, 1);
      assumedspectrum->SetParameter(0, 1.);

      // NOTE!! : tau NOT reliable beyond 30 TeV!! (just extrapolated!) Might be an issue only for a very very bright source very close to us...
      assumedspectrum->SetNpx(100000);
      // Important! This function will cannot be stored as a formula, precision depends on # of points! In this way we have 1 point/GeV.
    }

  else  // Galactic or nearby source. Just use the spectral shape as provided by user
    {
      assumedspectrum = new TF1("Spectrum", spec.Data());
      if (assumedspectrum->IsZombie())
	{
	  gLog << err << "Wrong function provided in flute.AssumedSpectrum !  Exiting." << endl;
	  return MAnalysisProblems::kWrongConfigFile;
	}
    }
  assumedspectrum->SetRange(10., 1.e5); // GeV
  assumedspectrum->SetTitle("Assumed Spectral Shape");

  // NOTE: EBL optical depth not reliable beyond 30 TeV (just extrapolated!) MIGHT be relevant for a very very bright and close source!

  // Maximum relative error delta_flux/flux in an energy bin: above this value
  // an upper limit will be calculated. Default is 50%:
  Double_t maximumRelError = env->GetValue("flute.maximumRelError", 0.5);
  // The same but for Light Curve bins:
  Double_t maximumRelErrorLC = env->GetValue("flute.maximumRelErrorLC", 0.5);

  // Upper limit calculation details:
  // Do we allow the use of "negative excesses" from fluctuations? (otherwise Non will be increased to match the background estimate):
  Bool_t AllowNegativeExcessInUpperLimit = env->GetValue("flute.AllowNegativeExcessInUpperLimit", kTRUE);
  // What is the minimum fraction of background that we allow as an upper limit? (safeguard against background systematics)
  Double_t MinAllowedUpperLimitAsFractionOfBackground = env->GetValue("flute.MinAllowedUpperLimitAsFractionOfBackground", 0.03);


  // Put hadronness and Theta2 cuts into HadTheta2Cut (of type MHadAlphaCut):
  MHadAlphaCut HadTheta2Cuts;
  HadTheta2Cuts.SetName("HadTheta2Cuts");
  MBinning zdbinning("BinningTheta"); // The naming is needed so that tasks MHEffectiveOnTime and MHMcCollectionArea can find it in the task list
  zdbinning.SetEdges(1, minZd, maxZd);
  HadTheta2Cuts.SetBinnings(eestbinning, zdbinning);
  TArrayD hadCutVsEest;
  TArrayD theta2CutVsEest;
  TArrayD thetaMin2CutVsEest;

  for (Int_t bineest = 0; bineest < nBinsEnergyEst; bineest++)
    {
      hadCutVsEest.Set(bineest+1);
      hadCutVsEest.AddAt(hadCuts.at(bineest), bineest);
      if(thetaMin2Cuts.size() > 0)  // Only for donought-like searches (e.g. ring to exclude central source in Perseus cluster)
	{
	  thetaMin2CutVsEest.Set(bineest+1);
	  thetaMin2CutVsEest.AddAt(thetaMin2Cuts.at(bineest), bineest);
	}
      theta2CutVsEest.Set(bineest+1);
      theta2CutVsEest.AddAt(theta2Cuts.at(bineest), bineest);
    }

  if (thetaMin2Cuts.size() > 0)
    HadTheta2Cuts.SetCutValues(1, hadCutVsEest, thetaMin2CutVsEest, theta2CutVsEest);
  else
    HadTheta2Cuts.SetCutValues(1, hadCutVsEest, theta2CutVsEest);

  Double_t fixedDeadTime;
  if ( ! env->GetRequiredValue("flute.deadTimePerEvent", fixedDeadTime) )
    return MAnalysisProblems::kWrongConfigFile;

  // Read time format for time axes:
  TString TimeFormatString;
  TimeFormatString = "%d/%m/%y %F1995-01-01 00:00:00 GMT";
  if (strcmp(env->GetValue("flute.TimeFormat", ""),"HMS") == 0)
    TimeFormatString = "%H:%M:%S %F1995-01-01 00:00:00 GMT";


  Bool_t PropagateEffAreaError = kTRUE;  // Propagate error in effective area (vs. Eest) to spectral points?

  Int_t check = env->GetValue("flute.PropagateEffAreaError", -1);
  if (check != -1)
    env->GetRequiredValue("flute.PropagateEffAreaError", PropagateEffAreaError);


  // Now the same but for the LC rather than for the spectrum:
  Bool_t PropagateEffAreaErrorToLightCurve = env->GetValue("flute.PropagateEffAreaErrorToLightCurve", kTRUE);


  gLog << endl << "Please check below if there are any unused (untouched) commands in the rc file!  :" << endl;
  env->PrintUntouched();
  gLog << endl;

  MBinning runedgetimes(0, 0, 0);

  // Now determine how many different telescope pointings are there in the data. This is done using the file names, i.e.
  // looking for the "source name" in the filename. File names are assumed to be of the type:
  // YYYYMMDD_RRRRRRRR_Q_SOURCENAME-W0.40+150.root  We use TChain to get the list of files in an easy way.
  // We find how many different "source names" are among the files, store those source names in a TList (ListOfPointingNames) for
  // later use, and read from the run header of one file of each set the coordinates (RA and Dec) of the telescope pointing,
  // which are stored in a TList of MPointingPos (ListOfPointings)
  TChain ch("RunHeaders");
  if (ch.Add(InputFilesRealData) == 0)
    {
      gLog << err << "ERROR: input files " << InputFilesRealData.Data() << " not found! Exiting." << endl;
      return MAnalysisProblems::kWrongPath;
    }
  TList ListOfPointings; // List of MPointingPos's containing the pointings of the telescopes in the different wobble positions
  TList ListOfPointingNames; // List of TObjStrings containing the different telescope pointings present in the input data files.
  gLog << inf << endl << "Input data files: " << endl;
  for (Int_t ifile = 0; ifile < ch.GetListOfFiles()->GetEntries(); ifile++)
    {
      TString fullfilename = ch.GetListOfFiles()->At(ifile)->GetTitle();
      if (!fullfilename.Contains("_Q_")  || !fullfilename.EndsWith(".root"))
	{
	  gLog << err << "ERROR: found an input file ("<< fullfilename << ") which is not a MARS ROOT file of _Q_ (melibea) type!" << endl;
	  return -1;
	}
      gLog << inf << fullfilename << endl;

      TChain chevt("Events");
      chevt.Add(fullfilename);
      MTime* time = 0;
      chevt.SetBranchStatus("*", 0);
      chevt.SetBranchStatus("MTime_1.*", 1);
      chevt.SetBranchAddress("MTime_1.", &time);

      // Get MJD of first and last event in run
      chevt.GetEntry(0);
      Double_t mjd1 = time->GetMjd();
      chevt.GetEntry(chevt.GetEntries()-1);
      Double_t mjd2 = time->GetMjd();
      if (runedgetimes.GetNumBins() == 0)
	runedgetimes.SetEdges(1, mjd1, mjd2);
      else
	{
	  runedgetimes.AddEdge((Axis_t) mjd1);
	  runedgetimes.AddEdge((Axis_t) mjd2);
	}

      if (LCbinningOption == "run-wise") // Set MJD limits of light curve bin as start and end of run:
	{
	  LCbinlowedge.Set(LCbinlowedge.GetSize()+1);
	  LCbinlowedge[LCbinlowedge.GetSize()-1] = mjd1;
	  LCbinupedge.Set(LCbinupedge.GetSize()+1);
	  LCbinupedge[LCbinupedge.GetSize()-1] = mjd2;
	}

      // Strip path, leave just the file name:
      TString filename = fullfilename(1+fullfilename.Last('/'), fullfilename.Length()-fullfilename.Last('/')-1);
      // Now leave just the "source" field, including the wobble orientation indications:
      TString sourcename = filename(filename.First('Q')+2, filename.Length()-filename.First('Q')-2-5);
      if (ifile == 0)
	{
	  TObjString* dummy = new TObjString(sourcename);
	  ListOfPointingNames.Add(dummy);
	  Double_t RA, Dec;
	  ReadRADec(fullfilename, &RA, &Dec); // Get the RA and Dec of the telescope pointing from the run header
	  MPointingPos* ppos = new MPointingPos("WobblePointing_0");
	  ppos->SetSkyPosition(RA, Dec);
	  ListOfPointings.Add(ppos);
	  continue;
	}

      TIter Next(&ListOfPointingNames);
      Bool_t knownpointing = kFALSE;
      for (;;)  // Check all previusly identified telescope pointings
	{
	  TObjString* ostr = (TObjString*) Next();
	  if (!ostr)
	    break;

	  if (ostr->GetString() == sourcename)  // => This is a known pointing, previously identified. No need to check further.
	    {
	      knownpointing = kTRUE;
	      break;
	    }
	}
      if (!knownpointing) // => this is a new pointing, not previously identified. Add it to the list.
	{
	  TObjString* dummy = new TObjString(sourcename);
	  ListOfPointingNames.Add(dummy); // Add new pointing.
	  Double_t RA, Dec;
	  ReadRADec(fullfilename, &RA, &Dec); // Get the RA and Dec of the telescope pointing from the run header
	  MPointingPos* ppos = new MPointingPos(Form("WobblePointing_%d", ListOfPointings.GetSize()));
	  ppos->SetSkyPosition(RA, Dec);
	  ListOfPointings.Add(ppos);
	}
    }
  gLog << inf << endl << "Identified telescope pointings: " << endl;
  TIter Next(&ListOfPointingNames);
  TObjString* pointingname = 0;
  Int_t ipointing = -1;

  // We will try to determine which combinations of Wobble pointings are valid for doing "off from wobble partner"
  // If two pointings are too close, the Off obtained with one of them (to be used as backg estimation for the other)
  // may contain gamma-rays from the source! We keep this info in the boolean array below:
  Bool_t* isOffValid = new Bool_t[ListOfPointings.GetSize()*ListOfPointings.GetSize()];

  while ( (pointingname = (TObjString*)Next()) )
    {
      ipointing++;
      Double_t thisRa  = ((MPointingPos*)ListOfPointings.At(ipointing))->GetRa();
      Double_t thisDec = ((MPointingPos*)ListOfPointings.At(ipointing))->GetDec();

      gLog << inf << pointingname->GetString() << "    RA: " << Form("%.3f", thisRa )
	   << " hours;  Dec = " << Form ("%.2f", thisDec) << " deg" << endl;

      // Check if this pointing is too close to any of the others (in that case they cannot be used by each other for the background calculation!)
      for (Int_t jpointing = 0; jpointing < ListOfPointings.GetSize(); jpointing++)
	{
	  if (ipointing == jpointing)
	    {
	      isOffValid[ipointing*ListOfPointings.GetSize() + ipointing] = kFALSE; // Obviously a file cannot be used as off (from wobble partner) estimator for itself
	      continue;
	    }
	  Double_t jRa  = ((MPointingPos*)ListOfPointings.At(jpointing))->GetRa();
	  Double_t jDec = ((MPointingPos*)ListOfPointings.At(jpointing))->GetDec();

	  Double_t angulardistance = TMath::RadToDeg() *
	    MAstro::AngularDistance((90.-thisDec)*TMath::DegToRad(), thisRa/12.*TMath::Pi(),
				    (90.-jDec)*TMath::DegToRad(),    jRa/12.*TMath::Pi());

	  if (angulardistance < 0.35)  // In degrees: this will make 6 evenly spaced wobble pointings (0.4 deg off-axis) still possible
	    isOffValid[ipointing*ListOfPointings.GetSize() + jpointing] = kFALSE;
	  else
	    isOffValid[ipointing*ListOfPointings.GetSize() + jpointing] = kTRUE;
	}
    }

  gLog << endl << inf << " Matrix of valid on - off combinations: " << endl << endl;
  for (Int_t ipointing = 0; ipointing < ListOfPointings.GetSize(); ipointing++)
    {
      gLog << "       ";
      for (Int_t jpointing = 0; jpointing < ListOfPointings.GetSize(); jpointing++)
	gLog << inf << (Int_t) isOffValid[ipointing*ListOfPointings.GetSize() + jpointing] << "  ";
      gLog << endl;
    }
  gLog << endl;


  if (processWobblePos >= 0)
    gLog << endl << inf
	 << "The spectrum will be calculated only from the signal in " << ((TObjString*)ListOfPointingNames.At(processWobblePos))->GetString()
	 << " files!" << endl << endl;


  // Create the object of type MTheta2vsEest where the theta2 histograms (on-source and off-source) for each energy bin will be stored
  MTheta2vsEest Theta2vsEest;
  Theta2vsEest.SetNormRangeTheta2(normRangeTheta2Min, normRangeTheta2Max);
  Theta2vsEest.SetBckgMode(bckgMode);


  Theta2vsEest.SetAllowNegativeExcessInUpperLimit(AllowNegativeExcessInUpperLimit);
  Theta2vsEest.SetMinAllowedUpperLimitAsFractionOfBackground(MinAllowedUpperLimitAsFractionOfBackground);

  // The Object Theta2vsEest2 is only needed for atm. corrections, but needs to be initialized here
  MTheta2vsEest Theta2vsEest2("MTheta2vsEest2");
  Theta2vsEest2.SetNormRangeTheta2(normRangeTheta2Min, normRangeTheta2Max);
  Theta2vsEest2.SetBckgMode(bckgMode);


  Theta2vsEest2.SetAllowNegativeExcessInUpperLimit(AllowNegativeExcessInUpperLimit);
  Theta2vsEest2.SetMinAllowedUpperLimitAsFractionOfBackground(MinAllowedUpperLimitAsFractionOfBackground);

  if (bckgMode == 3) // pulsar mode (i.e. on and off events from different pulsar phases, but same sky position)
    {
      Theta2vsEest. Init(2, eestbinning, nBinsTheta2, maxTheta2); // 2 is the minimum dimension needed to have on and off histograms inside MTheta2vsEest.
      Theta2vsEest2.Init(2, eestbinning, nBinsTheta2, maxTheta2); // 2 is the minimum dimension needed to have on and off histograms inside MTheta2vsEest.
    }
  else if (bckgMode == 0) // Off from wobble partner
    {
      Theta2vsEest.SetNormalizeEventNumbers(kTRUE);
      Theta2vsEest.Init(ListOfPointings.GetSize(), eestbinning, nBinsTheta2, maxTheta2); // Dimension determined by # of wobble positions
      Theta2vsEest.SetValidOff(isOffValid);
      Theta2vsEest2.SetNormalizeEventNumbers(kTRUE);
      Theta2vsEest2.Init(ListOfPointings.GetSize(), eestbinning, nBinsTheta2, maxTheta2); // Dimension determined by # of wobble positions
      Theta2vsEest2.SetValidOff(isOffValid);
    }
  else if (bckgMode == 1) // Off from simultaneous data (off positions in FoV)
    {
      Theta2vsEest.SetNormalizeEventNumbers(kFALSE); // We will use as normalization just: 1 / #of off positions
      Theta2vsEest.SetNumSimultaneousBgPositions(numSimultaneousBgPositions);
      Theta2vsEest.Init(ListOfPointings.GetSize(), eestbinning, nBinsTheta2, maxTheta2); // Dimension determined by # of wobble positions
      Theta2vsEest.SetValidOff(isOffValid);
      // The latter is just for the OffWP approach, but we set it because in this way we can check the histograms with OffWP background
      // that are stored in the MTheta2vsEest object in flute's output file.
      Theta2vsEest2.SetNormalizeEventNumbers(kFALSE); // We will use as normalization just: 1 / #of off positions
      Theta2vsEest2.SetNumSimultaneousBgPositions(numSimultaneousBgPositions);
      Theta2vsEest2.Init(ListOfPointings.GetSize(), eestbinning, nBinsTheta2, maxTheta2); // Dimension determined by # of wobble positions
      Theta2vsEest2.SetValidOff(isOffValid);
    }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // First loop over the real data:
  // Calculation of effective on-time. In making the calculation, we obviously have to apply the requested zenith angle cuts. This is
  // done by MHEffectiveOnTime, which reads in the zenith angle binning stored in zdbinning.
  // But it is better to omit other cuts, like those on Size, which change the event rate (and the average energy of the events in the
  // resulting sample) but not the observation time. Other such cuts are theta2 and hadronness cuts; we want to have as many events
  // as possible for the calculation of effective on-time.
  // We will loop over all the real data. The calculations are done by MHEffectiveOnTime.
  //
  // If processWobblePos >= 0, then only the wobble pointing indicated by processWobblePos will be used in spectrum calculations, and so
  // we will make the loop to calculate the eff. on-time only over those files.
  //

  MReadMarsFile readDataTcalc("Events");
  if (processWobblePos < 0) // Loop over all of the real data
    readDataTcalc.AddFile(InputFilesRealData);
  else // Loop only over pointing processWobblePos
    AddSubsetOfFiles(&readDataTcalc, InputFilesRealData, ((TObjString*)ListOfPointingNames.At(processWobblePos))->GetString());

  readDataTcalc.EnableBranch("MRawEvtHeader_2.*");
  readDataTcalc.EnableBranch("MPointingPos_1.*");
  readDataTcalc.EnableBranch("MTime_1.*");
  readDataTcalc.EnableBranch("MTime_2.*");
  readDataTcalc.EnableBranch("MStereoParDisp.*");
  readDataTcalc.VetoBranch("MBadPixelsCam_1");
  readDataTcalc.VetoBranch("MBadPixelsCam_2");
  MParList dataParListTcalc;
  MHMcCollectionArea collarea("MHMcCollectionAreaEtrue", "Collection area vs. true energy", nBinsAz); // This coll. area is for use by the unfolding program!
  collarea.SetStoreMCEvents(StoreMcEvents);
  collarea.SetMcRingAnalysis(McRingAnalysis);
  // Now we create another MHMcCollectionArea object, in which we will calculate the coll. area using as "efficiency" in a given
  // bin of energy (E1, E2) the ratio N_gammaMC,final (E1<Eest<E2) / N_gammaMC,generated (E1<Etrue<E2)  When we divide the excess
  // of real events in a given E-bin (E1<Eest<E2) by that coll. area, we should get back the number of g-rays in the corresponding bin
  // of tre energy (E1<Etrue<E2) assuming that the true spectrum matches the one assumed in the calculation of the collection area.
  // This object collareaVsEest will be used within flute, to calculate the diff. energy spectrum.
  MHMcCollectionArea collareaVsEest("collareaVsEest", "Collection area vs. estimated energy", nBinsAz);
  //collareaVsEest.SetStoreMCEvents(StoreMcEvents);
  collareaVsEest.SetMcRingAnalysis(McRingAnalysis);

  // For the computation of the migration matrix
  // energyMigration_IRFcalc used in the function compute_IRF, otherwise the binning are deleted at the end of the function so we need a copy.
  MHMcEnergyMigration energyMigration;
  //MHMcEnergyMigration *energyMigration_bis;
  if (AeffAtmCorr)
    {
      //
      // retrieve the fine binning, needed for the atm. corrections
      // the second time, the class MMcCollectionAreaCalc will call these functions
      // and retrieve exactly the same binning. However, then it will be too late
      // for the hdeltaetot histograms
      //
      collarea.SetCoarseBinnings(etruebinning,zdbinning);
      collarea.RedoFineEbinning();

      collareaVsEest.SetCoarseBinnings(eestbinning,zdbinning);
      collareaVsEest.RedoFineEbinning();
    }

  MBinning binsThetaFine("BinningThetaFine");
  binsThetaFine.SetEdges(*(collarea.GetHistAll()->GetZaxis()));
  dataParListTcalc.AddToList(&binsThetaFine);
  MBinning binsPhiFine("BinningPhiFine");
  binsPhiFine.SetEdges(*(collarea.GetHistAll()->GetYaxis()));
  dataParListTcalc.AddToList(&binsPhiFine);

  MBinning binsEFine("BinningEFine");    // needed for the atm. corrections
  binsEFine.SetEdges(*(collarea.GetHistAll()->GetXaxis()));
  MBinning binsEestFine("BinningEestFine");    // needed for the atm. corrections
  binsEestFine.SetEdges(*(collareaVsEest.GetHistAll()->GetXaxis()));

  MBinning binsDeltaT("BinningDeltaT"); // DeltaT binning for the histograms of MHEffectiveOnTime
  binsDeltaT.SetEdges(500, 0., 0.2);
  dataParListTcalc.AddToList(&binsDeltaT);
  dataParListTcalc.AddToList(&zdbinning);
  MHEffectiveOnTime* hEffTime = new MHEffectiveOnTime; // Task to calculate the effective on-time
  hEffTime->SetUseFixedDeadTimeMethod(kTRUE);          // Use the method which relies on the knowledge of the dead time per event.
  hEffTime->SetFixedDeadTime(fixedDeadTime);

  hEffTime->SetTimeDisplay(kTRUE);
  TAxis *xaxe = (TAxis*) hEffTime->GetHEffOnTime() .GetXaxis();
  xaxe->SetTimeFormat(TimeFormatString);
  xaxe = (TAxis*)hEffTime->GetHTimeProb()  .GetXaxis();
  xaxe->SetTimeFormat(TimeFormatString);
  xaxe = (TAxis*)hEffTime->GetHTimeLambda().GetXaxis();
  xaxe->SetTimeFormat(TimeFormatString);

  dataParListTcalc.AddToList(hEffTime);
  MFillH fillEffTime(hEffTime, "MTime_1");
  fillEffTime.SetBit(MFillH::kDoNotDisplay);
  MTaskList dataTaskListTcalc;
  dataTaskListTcalc.AddToList(&readDataTcalc);
  // Zenith angle cuts:
  MContinue *MinZdCut = new MContinue(Form("MPointingPos_1.fZd<%.2f", minZd), "MinZdCut");
  MContinue *MaxZdCut = new MContinue(Form("MPointingPos_1.fZd>%.2f", maxZd), "MaxZdCut");
  dataTaskListTcalc.AddToList(MinZdCut);
  dataTaskListTcalc.AddToList(MaxZdCut);

  // MJD cuts, if any:
  if (MJDrangeSelection)
    dataTaskListTcalc.AddToList(MJDrangeSelection);

  if (AeffAtmCorr)
    dataTaskListTcalc.AddToList(&filterecorr);

  dataTaskListTcalc.AddToList(&fillEffTime);
  dataParListTcalc.AddToList(&dataTaskListTcalc);
  MEvtLoop dataTcalcEvtLoop;
  dataTcalcEvtLoop.SetDisplay(display);
  dataTcalcEvtLoop.SetParList(&dataParListTcalc);


  // Calculation of (something like) the "Cherenkov Transparency coefficient"
  MHCTC ChTranspCoeff("ChTranspCoeff");
  ChTranspCoeff.SetTimeBinning(&runedgetimes);
  MFillH fillCTC(&ChTranspCoeff, "MEnergyEst");
  fillCTC.SetName("fillCTC");
  fillCTC.SetBit(MFillH::kDoNotDisplay); // do not show it automatically!
  dataTaskListTcalc.AddToList(&fillCTC);


  // Run the loop:
  if (!dataTcalcEvtLoop.Eventloop())
    return -1;
  dataTaskListTcalc.PrintStatistics();

  //
  // END of calculation of the effective on-time
  //
  // Display results of the effective on-time calculation:
  TString efftimetabname = "Eff. Time";
  TCanvas &canvtime = display->AddTab(efftimetabname);
  gWizard->WhiteBack(canvtime);
  hEffTime->Draw();
  gPad->Modified();
  gPad->Update();
  // Modified & Update needed, otherwise the tab won't save properly!

  if (kShowCTC)
    {
      TCanvas &canvCTC = display->AddTab("CTC");
      gWizard->WhiteBack(canvCTC);
      ChTranspCoeff.SetEffTime(hEffTime);
      ChTranspCoeff.Calc();
      ChTranspCoeff.Draw();
    }

  // Create filter for the hadronness cuts in data and MC (in MC it will also cut in theta2 for coll. area calculations)
  MFHadAlpha mfhadtheta2cuts;
  mfhadtheta2cuts.SetHadAlphaCutName("HadTheta2Cuts");
  mfhadtheta2cuts.SetIsStereo();
  mfhadtheta2cuts.SetNameHadronness("MHadronness");
  //  mfhadtheta2cuts.SetNameEnergyEst(nameEnergyEst);
  if (AeffAtmCorr)
    // this ensures that the cuts are applied in reconstructed energy,
    // and the effective area is that of reconstructed energy.
    // The correction on effective area due to the shift from MEnergyEst to MEnergyEstAtmCorr is
    // applied later on
    mfhadtheta2cuts.SetNameEnergyEst("MEnergyEst");
  else
    mfhadtheta2cuts.SetNameEnergyEst(nameEnergyEst);
  mfhadtheta2cuts.SetNameStereoPar(posContainerName);
  // We will now use an MContinue:
  MContinue mfconthadtheta2;
  mfconthadtheta2.SetFilter(&mfhadtheta2cuts);
  mfconthadtheta2.SetName("mfconthadtheta2");
  mfconthadtheta2.SetInverted();  // Note: this inverts the filter mfhadtheta2cuts itself!


  // Atmospheric Coll area correction stuff
  MH3      *hdeltaetot_eest  = 0;
  MH3      *hdeltaetot_etrue = 0;     // Etrue refers to the binning here only!
  MH3      *hdeltae_eest [ListOfPointings.GetSize()];
  MH3      *hdeltae_etrue[ListOfPointings.GetSize()];
  MFillH   *filldeltae_eest [ListOfPointings.GetSize()];
  MFillH   *filldeltae_etrue[ListOfPointings.GetSize()];

  MBinning *bdeltaex_eest [ListOfPointings.GetSize()];
  MBinning *bdeltaey_eest [ListOfPointings.GetSize()];
  MBinning *bdeltaez_eest [ListOfPointings.GetSize()];
  MBinning *bdeltaex_etrue[ListOfPointings.GetSize()];
  MBinning *bdeltaey_etrue[ListOfPointings.GetSize()];
  MBinning *bdeltaez_etrue[ListOfPointings.GetSize()];

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Second "loop" over the real data: in reality, it will be split in several loops, separating the data in the different pointings
  // of the telescopes.
  // We will fill the theta2 histograms from which the # of excess events (vs. estimated energy) will be obtained.
  // We treat separatedly the different pointings (FoVs) observed by the telescopes, since we want to use the method sometimes
  // called "off from wobble partner": the background under the signal "contained" in a given set of files is estimated from different
  // data files, by looking at the events in a position on the sky which follows the same path on the camera as the candidate source.
  //
  // We will make one loop for every different pointing. In each loop we will fill the "ON" histogram (theta2 calculated w.r.t. the
  // candidate source position) and the OFF histograms (a total of N-1, N being the number of different pointings in the sample,
  // i.e. N = ListOfPointings.GetSize() ). Note that those N-1 OFF histograms are the OFF for the "ON" filled in the other loops!
  //

  Next.Reset();
  ipointing = -1;
  while ( (pointingname = (TObjString*)Next()) )
    {
      ipointing++;
      gLog << inf << endl << endl << "Looping over " << pointingname->GetString() << " data files..." << endl << endl;

      // Loop over the files:
      MReadMarsFile readData("Events");
      AddSubsetOfFiles(&readData, InputFilesRealData, pointingname->GetString());
      // readData.DisableAutoScheme(); // This would enable reading of all branches, making things slower!
      // Enable branches needed in the loop:

      readData.EnableBranch(Form("%s.fEnergy",nameEnergyEst.Data()));
      readData.EnableBranch(Form("MEnergyEst.fEnergy"));

      if (EnergyReCalc)
	{
	  readData.EnableBranch("MEnergyEst_1.fEnergy");
	  readData.EnableBranch("MEnergyEst_2.fEnergy");
	  readData.EnableBranch("MEnergyEst_1.fUncertainty");
	  readData.EnableBranch("MEnergyEst_2.fUncertainty");
	}

      readData.EnableBranch("MHadronness.fHadronness");
      readData.EnableBranch("MPointingPos_1.*");
      readData.EnableBranch("MTime_1.*");
      readData.EnableBranch("MHillas_1.fSize");
      readData.EnableBranch("MHillas_2.fSize");
      readData.EnableBranch(posContainerName+TString(".*"));
      // Disabling the MBadPixelsCam* branches below is _very_ important for speed; even if they are in
      // the run headers, the "Reset" function will be called for every event in the loop over the Events tree!
      readData.VetoBranch("MBadPixelsCam_1");
      readData.VetoBranch("MBadPixelsCam_2");

      MTaskList dataTaskList;
      MParList  dataParList;
      dataParList.AddToList(&HadTheta2Cuts);
      dataParList.AddToList(&Theta2vsEest);
      if (AeffAtmCorr)
	dataParList.AddToList(&Theta2vsEest2);
      dataTaskList.AddToList(&readData);

      // Zd cuts:
      dataTaskList.AddToList(MinZdCut);
      dataTaskList.AddToList(MaxZdCut);

      // MJD cuts if any:
      if (MJDrangeSelection)
	dataTaskList.AddToList(MJDrangeSelection);

      MAverageEnergy averageEnergies; // Task to re-calculate the average of M1 & M2 energies
      if (EnergyReCalc)
	{
	  averageEnergies.SetE1scaling(eest_scaling_factor*eest1_scaling_factor);
	  averageEnergies.SetE2scaling(eest_scaling_factor*eest2_scaling_factor);
	  dataTaskList.AddToList(&averageEnergies);
	}


      // Size cuts:
      dataTaskList.AddToList(&filterSizeRealData);

      // User cuts:
      if (usercuts)
	dataTaskList.AddToList(usercuts);

      // Now the hadronness cut:
      mfhadtheta2cuts.SetCutType(MFHadAlpha::kHadCut);
      dataTaskList.AddToList(&mfconthadtheta2);

      // LIDAR re-calibration:
      MLidarInterpolate lextr;
      lextr.SetOutname(statusOutName);
      lextr.SetSerialNumber(1);
      lextr.AddFiles(InputFilesRealData);

      // LIDAR corrections
      //    lextr.SetUseGFS();   // not yet, depends on the GFS reports written to files
      if (LIDARUseGDAS)
	lextr.SetUseGDAS();
      lextr.SetUseWeather();
      lextr.SetUseTimeCorr();

      MBinning lbins("BinningHeight");
      lbins.SetEdges(62, 0.5, 16.);    // km

      MBinning rbins("BinningHeightRep");
      rbins.SetEdges(62,0.5, 16.);    // km

      MTransmissionEnergyCorr tcorr;
      if (LIDARRecalibrate)
	{
	  dataParList.AddToList(&lbins);
	  dataParList.AddToList(&rbins);
	  dataTaskList.AddToList(&lextr);
	  dataTaskList.AddToList(&tcorr);
	}

      if (AeffAtmCorr)
	dataTaskList.AddToList(&filterecorr);

      MSrcPosCalc srcposcalc;
      srcposcalc.SetSourcePos(SourceRA, SourceDec);
      srcposcalc.SetSerialNumber(1);
      dataTaskList.AddToList(&srcposcalc);


      if (bckgMode == 0)  // "Off from Wobble Partner"
	{
	  // Create the tasks needed to find where the off positions lie in the camera.
	  // Loop over the other telescope pointings; for each of them (different from the current pointing "ipointing") we have to
	  // calculate one "off" position on the camera for the current set of files. In the most common case of two wobble pointings
	  // the number of off positions will be one. If there are 4 wobble positions, each has 3 offs, and so on.
	  // We Use a task MSrcPosCalc to calculate each off position; in that task we set the pointing position to the corresponding
	  // wobble pointing. In the common case of two wobble pointings, we set, for the off position calculation, the other wobble
	  // pointing -> the resultin off position is the "anti-source", i.e. the point opposite to the source with respect to the
	  // camera center.
	  MSrcPosCalc* offposcalc = new MSrcPosCalc[ListOfPointings.GetSize()-1];
	  Int_t ioff = 0;
	  for (Int_t ipoint = 0; ipoint < ListOfPointings.GetSize(); ipoint++)
	    {
	      if ( ipoint == ipointing)
		continue;
	      offposcalc[ioff].SetName(Form("MSrcPosCalc_off_%d_%d", ipointing, ipoint));
	      // Note that the task will add to the name above the telescope tag (serial number): "_1"
	      offposcalc[ioff].SetPointPos((MPointingPos*)ListOfPointings.At(ipoint));
	      offposcalc[ioff].SetSourcePos(SourceRA, SourceDec);
	      offposcalc[ioff].SetExactNameSrcPosCam(Form("OffPosCam_%d_%d", ipointing, ipoint));
	      offposcalc[ioff].SetSerialNumber(1); // Use M1's MTime, RunHeader, etc.
	      dataTaskList.AddToList(&offposcalc[ioff]);
	      ioff++;
	    }
	}
      else  // "Off for pulsars mode, or in special cases of simultaneous Off mode"
	    // (most often FilterOn and FilterOff are filters on the pulsar phase to select ON and OFF)
	{
	  if (onSelection != "")
	    {
	      MF* FilterOn  = new MF(onSelection.Data(),  "FilterOn");
	      dataTaskList.AddToList(FilterOn);
	      // Enable branches needed for the cuts:
	      TObjArray* vars =  FilterOn->GetDataMember().Tokenize(",");
	      for (Int_t ivar = 0; ivar < vars->GetEntries(); ivar++)
		readData.EnableBranch( ((TObjString*)vars->At(ivar))->GetString().Data() );
	    }
	  if (offSelection != "")
	    {
	      MF* FilterOff = new MF(offSelection.Data(), "FilterOff");
	      dataTaskList.AddToList(FilterOff);
	      // Enable branches needed for the cuts:
	      TObjArray* vars = FilterOff->GetDataMember().Tokenize(",");
	      for (Int_t ivar = 0; ivar < vars->GetEntries(); ivar++)
		readData.EnableBranch( ((TObjString*)vars->At(ivar))->GetString().Data() );
	    }
	}

      MCalcExcess calcexcess;
      MCalcExcess calcexcess2("MCalcExcess2");

      calcexcess.SetStereoParName(posContainerName); // Container from which stereo information has to be obtained
      calcexcess.SetBckgMode(bckgMode);
      calcexcess2.SetStereoParName(posContainerName); // Container from which stereo information has to be obtained
      calcexcess2.SetBckgMode(bckgMode);

      calcexcess. SetPointingId(ipointing); // Indicate which of the wobble pointings we are processing
      calcexcess2.SetPointingId(ipointing); // Indicate which of the wobble pointings we are processing
      calcexcess. SetNameEnergyEst(nameEnergyEst);  // in case of atm. corrections, nameEnergyEst will be "EnergyEstAtmCorr"
      calcexcess2.SetNameEnergyEst("MEnergyEst");   // in case of atm. corrections, will be different from nameEnergyEst
      calcexcess2.SetNameTheta2vsEest("MTheta2vsEest2");

      dataTaskList.AddToList(&calcexcess);

      MF* efilt = 0;
      MCalcDeltaT* calcdelta = 0;

      if (AeffAtmCorr)
	{
	  dataTaskList.AddToList(&calcexcess2);

	  readData.EnableBranch("MEnergyEst.fEnergy");       // just to be sure...
	  readData.EnableBranch("MEnergyEstAtmCorr.fEnergy");   //
	  readData.EnableBranch("MPointingPos_1.fZd");

	  hdeltae_eest[ipointing] = new MH3("MEnergyEst.fEnergy","MEnergyEstAtmCorr.fEnergy","MPointingPos_1.fZd");
	  hdeltae_eest[ipointing]->SetName(Form("HDeltaE_eest%d",ipointing));

	  // Here reconstructed energy available only (and hence required)
	  // BUT the binning is done in true energy bins!
	  hdeltae_etrue[ipointing] = new MH3("MEnergyEst.fEnergy","MEnergyEstAtmCorr.fEnergy","MPointingPos_1.fZd");
	  hdeltae_etrue[ipointing]->SetName(Form("HDeltaE_etrue%d",ipointing));

	  bdeltaex_eest[ipointing] = new MBinning(Form("BinningHDeltaE_eest%dX",ipointing));
	  bdeltaey_eest[ipointing] = new MBinning(Form("BinningHDeltaE_eest%dY",ipointing));
	  bdeltaez_eest[ipointing] = new MBinning(Form("BinningHDeltaE_eest%dZ",ipointing));

	  bdeltaex_etrue[ipointing] = new MBinning(Form("BinningHDeltaE_etrue%dX",ipointing));
	  bdeltaey_etrue[ipointing] = new MBinning(Form("BinningHDeltaE_etrue%dY",ipointing));
	  bdeltaez_etrue[ipointing] = new MBinning(Form("BinningHDeltaE_etrue%dZ",ipointing));

	  eestbinning.Copy(*(bdeltaex_eest[ipointing]));
	  eestbinning.Copy(*(bdeltaey_eest[ipointing]));
	  //	  zdbinning  .Copy(*(bdeltaez_eest[ipointing]));  // This would be just one bin!
	  binsThetaFine.Copy(*(bdeltaez_eest[ipointing]));

	  binsEFine.Copy(*(bdeltaex_etrue[ipointing]));
	  binsEFine.Copy(*(bdeltaey_etrue[ipointing]));
	  binsThetaFine.Copy(*(bdeltaez_etrue[ipointing]));

	  efilt = new MF("(MEnergyEst.fEnergy > 0.) && (MEnergyEstAtmCorr.fEnergy > 0.)", "FilterValidE");
	  calcdelta = new MCalcDeltaT;
	  // calculates time between current and previously processed event.
	  // We use this as weight to fill the deltae histograms, so that
	  // each transmission condition is weighted by the time spent in
	  // such condition
	  dataTaskList.AddToList(efilt);
	  calcdelta->SetFilter(efilt);
	  dataTaskList.AddToList(calcdelta);

	  filldeltae_eest[ipointing]  = new MFillH(hdeltae_eest[ipointing],"", Form("FillDeltaE_eest%d", ipointing));
	  filldeltae_eest[ipointing]->SetFilter(efilt); // the filter is just to avoid misleading event stats in TH3 projections, due to underflows
	  filldeltae_eest[ipointing]->SetWeight("DeltaTime"); // Weight is calculated by calcdelta
	  filldeltae_eest[ipointing]->SetBit(MFillH::kDoNotDisplay);
	  dataTaskList.AddToList(filldeltae_eest[ipointing]);

	  filldeltae_etrue[ipointing] = new MFillH(hdeltae_etrue[ipointing],"",Form("FillDeltaE_etrue%d",ipointing));
	  filldeltae_etrue[ipointing]->SetFilter(efilt);
	  filldeltae_etrue[ipointing]->SetWeight("DeltaTime");
	  filldeltae_etrue[ipointing]->SetBit(MFillH::kDoNotDisplay);
	  dataTaskList.AddToList(filldeltae_etrue[ipointing]);

	  dataParList.AddToList(bdeltaex_eest[ipointing]);
	  dataParList.AddToList(bdeltaey_eest[ipointing]);
	  dataParList.AddToList(bdeltaez_eest[ipointing]);
	  dataParList.AddToList(bdeltaex_etrue[ipointing]);
	  dataParList.AddToList(bdeltaey_etrue[ipointing]);
	  dataParList.AddToList(bdeltaez_etrue[ipointing]);
	}

      dataParList.AddToList(&dataTaskList);
      MEvtLoop dataEvtLoop;
      dataEvtLoop.SetDisplay(display);
      dataEvtLoop.SetParList(&dataParList);

      // Run the loop:
      if (!dataEvtLoop.Eventloop())
	return -1;
      dataTaskList.PrintStatistics();
      if (calcdelta)
	delete calcdelta;
      if (efilt)
	delete efilt;
    }
  gLog << endl;

  Int_t filledbins_eest  = 0;
  Int_t filledbins_etrue = 0;

  if (fixedBckgNormalization>0) // use the fixed normalization stated in the .rc file (e.g. in pulsar mode):
    {
      Theta2vsEest. SetDefaultNormFactor(fixedBckgNormalization);
      Theta2vsEest. SetForceUseDefaultNormFactor(kTRUE);
      Theta2vsEest2.SetDefaultNormFactor(fixedBckgNormalization);
      Theta2vsEest2. SetForceUseDefaultNormFactor(kTRUE);
    }

  if (bckgMode != 3)
    {
      // Display Theta^2 plots for each wobble pointing (and their corresponding normalized OFFs):
      for (Int_t ipoint = 0; ipoint < ListOfPointings.GetSize(); ipoint++)
	{
	  if (processWobblePos >= 0 && ipoint != processWobblePos)  // Just display the pointing that has been searched for signal!
	    continue;

	  TCanvas &canv = display->AddTab(Form("WOB%d", ipoint));
	  gWizard->WhiteBack(canv);
	  canv.SetGridy();

	  Theta2vsEest.DrawNormalizedHists(Theta2vsEest.GetHist(), ipoint, 1, nBinsEnergyEst);

	  if (AeffAtmCorr)
	    {
	      Int_t nbins = ((TH3&)hdeltae_etrue[ipoint]->GetHist()).GetNbinsZ();
	      filledbins_etrue  = 0;
	      for (Int_t b=1;b<=nbins;b++)
		{
		  TH3 &h = (TH3&)hdeltae_etrue[ipoint]->GetHist();
		  h.GetZaxis()->SetRange(b,b);
		  TH2 *hist = (TH2*)h.Project3D("yx");
		  if (hist->GetEntries()>0)
		    filledbins_etrue++;
		  delete hist;
		}

	      nbins = ((TH3&)hdeltae_eest[ipoint]->GetHist()).GetNbinsZ();

	      filledbins_eest = 0;
	      TH3 &h = (TH3&)hdeltae_eest[ipoint]->GetHist();
	      TH1 *hist = (TH1*)h.Project3D("z");
	      for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++)
		if (hist->GetBinContent(ibin) > 0)
		  filledbins_eest++;
	      delete hist;

	      if (AtmCorrDisplay)
		{
		  TCanvas &cana = display->AddTab(Form("ATM%d_EEST", ipoint));
		  const Int_t npads = filledbins_eest==1 ? 1 : (Int_t) (sqrt((Float_t)filledbins_eest)+1.);
		  Int_t mpads;
		  for (mpads = 1; mpads <= npads; mpads++)
		    if (npads*mpads >= filledbins_eest)
		      break;
		  cana.Divide(npads, mpads);
		  gWizard->WhiteBack(cana);

		  filledbins_eest = 0;
		  for (Int_t b = 1; b <= nbins; b++)
		    {
		      TH3 &h = (TH3&)hdeltae_eest[ipoint]->GetHist();
		      h.GetZaxis()->SetRange(b,b);
		      TH2 *hist = (TH2*)h.Project3D("yx");
		      hist->SetName(Form("HDELTAE_W%1d_EST%02d", ipoint, b));
		      hist->SetTitle(Form("HDELTAE_EST, Zd = %.1f - %.1f deg",
					  h.GetZaxis()->GetBinLowEdge(b),
					  h.GetZaxis()->GetBinUpEdge(b)));

		      //  gLog << dbg << Form("HIST%d: ",ipoint) << hist->GetEntries() << endl;
		      if (hist->GetEntries()>0)
			{
			  cana.cd(++filledbins_eest);
			  hist->Draw("colz");
			  gPad->SetLogx();
			  gPad->SetLogy();
			  gPad->SetGridx();
			  gPad->SetGridy();
			  gPad->Modified();
			  gPad->Update();
			}
		      // delete hist;
		    }
		}
	    }
	}
    }

  if (AeffAtmCorr)
    {
      hdeltaetot_eest  = (MH3*)hdeltae_eest[0]->New();
      hdeltaetot_eest ->SetName("HDeltaEtot_eest");
      hdeltaetot_etrue = (MH3*)hdeltae_etrue[0]->New();
      hdeltaetot_etrue->SetName("HDeltaEtot_etrue");

      TH3 &htot_eest  = (TH3&)hdeltaetot_eest ->GetHist();
      TH3 &htot_etrue = (TH3&)hdeltaetot_etrue->GetHist();

      TH3 &h0_eest    = (TH3&)hdeltae_eest [0]->GetHist();
      TH3 &h0_etrue   = (TH3&)hdeltae_etrue[0]->GetHist();

      h0_eest .GetZaxis()->SetRange(1,zdbinning.GetNumBins());
      h0_etrue.GetZaxis()->SetRange(1,zdbinning.GetNumBins());

      //      MH::SetBinning(&htot_eest,  &eestbinning, &eestbinning, &zdbinning);
      //      MH::SetBinning(&htot_etrue, &eestbinning, &eestbinning, &zdbinning);
      MH::SetBinning(&htot_eest, &eestbinning, &eestbinning, &binsThetaFine);
      MH::SetBinning(&htot_etrue,&binsEFine,&binsEFine,&binsThetaFine);

      for (Int_t ipoint = 0; ipoint < ListOfPointings.GetSize(); ipoint++)
	{
	  TH3 &h_eest    = (TH3&)hdeltae_eest [ipoint]->GetHist();
	  TH3 &h_etrue   = (TH3&)hdeltae_etrue[ipoint]->GetHist();

	  h_eest.    GetZaxis()->SetRange(1, h_eest.GetNbinsZ());
	  h_etrue.   GetZaxis()->SetRange(1, h_etrue.GetNbinsZ());
	  htot_eest. GetZaxis()->SetRange(1, htot_eest.GetNbinsZ());
	  htot_etrue.GetZaxis()->SetRange(1, htot_etrue.GetNbinsZ());

	  for (Int_t ibin=0;ibin<=htot_eest.GetNbinsX();ibin++)
	    for (Int_t jbin=0;jbin<=htot_eest.GetNbinsY();jbin++)
	      for (Int_t kbin=0;kbin<=htot_eest.GetNbinsZ();kbin++)
		htot_eest.SetBinContent(ibin,jbin,kbin,htot_eest.GetBinContent(ibin,jbin,kbin)+h_eest.GetBinContent(ibin,jbin,kbin));

	  for (Int_t ibin=0;ibin<=htot_etrue.GetNbinsX();ibin++)
	    for (Int_t jbin=0;jbin<=htot_etrue.GetNbinsY();jbin++)
	      for (Int_t kbin=0;kbin<=htot_etrue.GetNbinsZ();kbin++)
	  	htot_etrue.SetBinContent(ibin,jbin,kbin,htot_etrue.GetBinContent(ibin,jbin,kbin)+h_etrue.GetBinContent(ibin,jbin,kbin));
	}

      TCanvas *cana = NULL;
      Int_t npads = 0;
      if (AtmCorrDisplay)
	{
	  filledbins_eest = 0;
	  TH1 *hist = (TH1*)htot_eest.Project3D("z");
	  for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++)
	    if (hist->GetBinContent(ibin) > 0)
	      filledbins_eest++;
	  delete hist;

	  cana = &(display->AddTab("ATMTOT_EEST"));
	  npads = filledbins_eest==1 ? 1 : (filledbins_eest==2 ? 2 : (Int_t) (sqrt((Float_t)filledbins_eest)+1.));
	  Int_t mpads;
	  for (mpads = 1; mpads <= npads; mpads++)
	    if (npads*mpads >= filledbins_eest)
	      break;

	  cana->Divide(npads, mpads);
	  gWizard->WhiteBack(*cana);

	  filledbins_eest = 0;
	  for (Int_t b=1; b<=htot_eest.GetNbinsZ(); b++)
	    {
	      htot_eest.GetZaxis()->SetRange(b,b);
	      TH2 *hist = (TH2*)htot_eest.Project3D("yx");
	      hist->SetName(Form("HDELTAE_EST%02d",b));
	      hist->SetTitle(Form("HDELTAE_EST, Zd = %.1f - %.1f deg",
				  htot_eest.GetZaxis()->GetBinLowEdge(b),
				  htot_eest.GetZaxis()->GetBinUpEdge(b)));
	      //	  gLog << dbg << "HIST: " << hist->GetEntries() << endl;
	      if (hist->GetEntries()>0)
		{
		  // There used to be a normalization to 1 of  hist's rows here, but it is not really needed
		  // (and it assumed histogram contents were events - now they have units of time - seconds)
		  cana->cd(++filledbins_eest);
		  hist->Draw("colz");
		  gPad->SetLogx();
		  gPad->SetLogy();
		  gPad->SetGridx();
		  gPad->SetGridy();
		  gPad->Modified();
		  gPad->Update();
		}
	      // delete hist;
	    }
	}

      //Put back full range of zenith axes:
      htot_eest. GetZaxis()->SetRange(1, htot_eest.GetNbinsZ());

      TCanvas *canb = NULL;
      if (AtmCorrDisplay)
	{
	  filledbins_etrue = 0;
	  TH1 *hist = (TH1*)htot_etrue.Project3D("z");
	  for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++)
	    if (hist->GetBinContent(ibin) > 0)
	      filledbins_etrue++;
	  delete hist;

	  canb = &(display->AddTab("ATMTOT_ETRUE"));
	  npads = filledbins_etrue==1 ? 1 : (Int_t) (sqrt((Float_t)filledbins_etrue)+1.);
	  Int_t mpads;
	  for (mpads = 1; mpads <= npads; mpads++)
	    if (npads*mpads >= filledbins_etrue)
	      break;
	  canb->Divide(npads, mpads);
	  gWizard->WhiteBack(*canb);

	  filledbins_etrue = 0;
	  for (Int_t b=1; b<=htot_etrue.GetNbinsZ(); b++)
	    {
	      htot_etrue.GetZaxis()->SetRange(b,b);
	      TH2 *hist = (TH2*)htot_etrue.Project3D("yx");
	      hist->SetName(Form("HDELTAE_TRUE%02d",b));
	      hist->SetTitle(Form("HDELTAE_TRUE, Zd = %.1f - %.1f deg",
				  htot_etrue.GetZaxis()->GetBinLowEdge(b),
				  htot_etrue.GetZaxis()->GetBinUpEdge(b)));

	      //	  gLog << dbg << "THETA: " << b << ": " << endl;
	      if (hist->GetEntries()>0)
		{
		  // There used to be a normalization to 1 of  hist's rows here, but it is not really needed
		  // (and it assumed histogram contents were events - now they have units of time - seconds)
		  //	      gLog << dbg << "HISTN" << b << ": " << hist->GetEntries() << endl;
		  canb->cd(++filledbins_etrue);
		  hist->Draw("colz");
		  gPad->SetLogx();
		  gPad->SetLogy();
		  gPad->SetGridx();
		  gPad->SetGridy();
		  gPad->Modified();
		  gPad->Update();
		}
	      //	  delete hist;
	    }
	  htot_etrue.GetZaxis()->SetRange(1, htot_etrue.GetNbinsZ());
	}
    }

  delete MinZdCut;
  delete MaxZdCut;

  // Now the overall Theta2 plot, adding up all pointings (or only the chosen one if processWobblePos >= 0):
  TCanvas &canvtheta2 = display->AddTab("Theta2, ALL");
  gWizard->WhiteBack(canvtheta2);
  canvtheta2.SetGridy();

  Theta2vsEest.DrawNormalizedHists(Theta2vsEest.GetHist(), processWobblePos, 1,nBinsEnergyEst);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Start IRF calculation done in two parts: one first loop over the original and one second loop over the reconstructed one
  /// If the diffuse MCs are used, then the collective area is calculated in the different predefined ring.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // For the weighting of the effectivearea and Migration matrix depending on the EffectiveOntime distribution of the events
  const TH2D& EffOnTimevsAzvsZd = hEffTime->GetHEffOnPhiTheta();

  // The binning below are needed by MHMcEnergyMigration:
  MBinning BinningDE("BinningDE");
  BinningDE.SetEdges(60, -1.2, 1.2);    // Binning for (Eest-Etrue)/Eest
  MBinning BinningImpact("BinningImpact");
  BinningImpact.SetEdges(60, 0., 600.); // Binning for impact parameter

  // standard flute analysis: if single offset MC simulation used the 0.4 degree simulation, if the diffuse MCs are used, determine the ring
  //to compute the IRFs
  if(McRingAnalysis){
      	  //LeaJouvin: If diffuse Analysis, have to determine the diffuse MC ring matching the source position
      	  // Ring defined as the center on the source position and a width of 0.2 degree around the center
      	  vector<TString> wobble_name;
      	  for (Int_t ifile = 0; ifile < ch.GetListOfFiles()->GetEntries(); ifile++){
      		  TString fullfilename = ch.GetListOfFiles()->At(ifile)->GetTitle();
      		  wobble_name.push_back(fullfilename(fullfilename.Last('W')+1, 4));
      	  }
      	  for (unsigned int i=0 ; i < wobble_name.size(); i++){
      		  if(wobble_name[0]!=wobble_name[i]){
      			  gLog << err << "The wobble offset value are not compatible between the data files!  Exiting..." << endl;
      			  return MAnalysisProblems:: kUnknownProblem;
      		  }
      	  }
      	  double wobble;
      	  wobble=wobble_name[0].Atof();
      	  double off_min, off_max;
      	  off_min= wobble - 0.1;
      	  off_max= wobble + 0.1;
      	  gLog << "IRF calculated from diffuse MC in the following ring:"<< endl;
      	  gLog << off_min <<  "-" << off_max << " degree"<< endl;
      	  collarea.SetMcRingAnalysis(McRingAnalysis);
      	  collarea.SetOffsetRingMin(off_min);
      	  collarea.SetOffsetRingMax(off_max);
      	  collarea.SetOffsetRingMed((off_min + off_max)/2.);
        }
    Bool_t IsComputeIRFOk;
    IsComputeIRFOk = Compute_IRF(&collarea, &BinningDE, &BinningImpact, minZd, maxZd, &zdbinning, &eestbinning, &etruebinning, InputFilesMC, statusOutName, posContainerName, display,&mfhadtheta2cuts,
    		HadTheta2Cuts,&mfconthadtheta2, AeffAtmCorr, nameEnergyEst,	assumedspectrum, EnergyReCalc, &filterSize,usercuts, &energyMigration ,hdeltaetot_eest, hdeltaetot_etrue,EffOnTimevsAzvsZd,
			kTRUE, &collareaVsEest);

    if(!IsComputeIRFOk){
    		return -1;
    	}

    // for each M ring: compute the true effective area and fill the TNtuple MC to compute a posteriori the IRFs  each ring for the conversion to DL3
    if((McRingAnalysis) & (StoreMcEvents)){
    Bool_t IsComputeIRFOk_ringMC;
  	  for (int i=0; i< NRing; i++){
  		MHMcCollectionArea *collarea_temp = new MHMcCollectionArea("MHMcCollectionAreaEtrue", "Collection area vs. true energy", nBinsAz);
  		collarea_temp->SetMcRingAnalysis(McRingAnalysis);
  		collarea_temp->SetStoreMCEvents(StoreMcEvents);
		if (AeffAtmCorr)
		{
			//
			// retrieve the fine binning, needed for the atm. corrections
			// the second time, the class MMcCollectionAreaCalc will call these functions
			// and retrieve exactly the same binning. However, then it will be too late
			// for the hdeltaetot histograms
			//
			collarea_temp->SetCoarseBinnings(etruebinning,zdbinning);
			collarea_temp->RedoFineEbinning();
		 }


		collarea_temp->SetMcRingAnalysis(McRingAnalysis);
		collarea_temp->SetOffsetRingMin(RingMcEdges[i]);
		collarea_temp->SetOffsetRingMax(RingMcEdges[i+1]);
		collarea_temp->SetOffsetRingMed((RingMcEdges[i]+RingMcEdges[i+1])/2.);

  		IsComputeIRFOk_ringMC = Compute_IRF(collarea_temp, &BinningDE, &BinningImpact, minZd, maxZd, &zdbinning, &eestbinning, &etruebinning, InputFilesMC, statusOutName, posContainerName, display,&mfhadtheta2cuts,
  	    		HadTheta2Cuts,&mfconthadtheta2, AeffAtmCorr, nameEnergyEst,	assumedspectrum, EnergyReCalc, &filterSize,usercuts, &energyMigration ,hdeltaetot_eest, hdeltaetot_etrue,EffOnTimevsAzvsZd,
  				kFALSE, &collareaVsEest);
  		if(!IsComputeIRFOk_ringMC){
  			return -1;
  		}

  		  MHMcCollectionArea* test_area;
  		  test_area = (MHMcCollectionArea*) collarea_temp->Clone();
  		  collarea_per_ring.push_back(test_area);
  	  }

	  }

  //
  // Display the collection areas:
  //
  TCanvas &canvCollAreaRaw = display->AddTab("Coll. Area Raw");
  gWizard->WhiteBack(canvCollAreaRaw);
  canvCollAreaRaw.SetRightMargin(0.12);
  canvCollAreaRaw.SetTopMargin(0.13);
  TH2D* AeffvsZdvsE = (TH2D*)collarea.GetHist()->Project3D("zx");
  // Average in all azimuths, assumes all bins filled &
  // flat distribution in azimuth:
  AeffvsZdvsE->GetYaxis()->SetTitle("Zd (deg)");
  AeffvsZdvsE->Scale(1./(Float_t)collarea.GetHist()->GetNbinsY());
  AeffvsZdvsE->SetName("AeffvsZdvsE");
  AeffvsZdvsE->SetTitle("A_{eff}(m^{2}) vs. E_{true} & Zd");
  canvCollAreaRaw.SetLogx();
  canvCollAreaRaw.SetLogz();
  canvCollAreaRaw.SetGridx();
  canvCollAreaRaw.SetGridy();
  AeffvsZdvsE->SetStats(kFALSE);
  AeffvsZdvsE->Draw("zcol");
  TCanvas &canvCollArea = display->AddTab("Coll. Area");
  gWizard->WhiteBack(canvCollArea);
  canvCollArea.SetLogx();
  canvCollArea.SetLogy();
  canvCollArea.SetGridx();
  canvCollArea.SetGridy();
  TH1D* AeffVsEtrue = collarea.GetHistCoarse()->ProjectionX("AeffVsEtrue",1,1,"e");
  AeffVsEtrue->SetStats(kFALSE);
  AeffVsEtrue->SetMarkerColor(4);
  AeffVsEtrue->SetFillColor(4);
  AeffVsEtrue->SetMinimum(1.);
  AeffVsEtrue->SetMaximum(3.e6);
  AeffVsEtrue->Draw("E2");
  // We choose here and below for AeffVsEest the option E2 and same color for maker and filling so that the bin center
  // (in linear scale) is not shown. Otherwise it may seem that the Aeff value refers to that specific energy, whereas
  // actually it is an average in the energy bin (an average weighted with the assumed gamma-ray spectrum)

  TCanvas &canvCollAreaEest = display->AddTab("Coll. Area Eest");
  gWizard->WhiteBack(canvCollAreaEest);
  canvCollAreaEest.SetLogx();
  canvCollAreaEest.SetLogy();
  canvCollAreaEest.SetGridx();
  canvCollAreaEest.SetGridy();
  TH1D* AeffVsEest = collareaVsEest.GetHistCoarse()->ProjectionX("AeffVsEest",1,1,"e");
  AeffVsEest->SetStats(kFALSE);
  AeffVsEest->SetMarkerColor(4);
  AeffVsEest->SetMinimum(1.);
  AeffVsEest->SetMaximum(3.e6);
  AeffVsEest->SetFillColor(4);
  AeffVsEest->Draw("E2");

  if (AeffAtmCorr)
    {
      AeffVsEest = collareaVsEest.GetHistCoarse()->ProjectionX("AeffVsEest",1,1,"e");
      AeffVsEest->SetStats(kFALSE);
      AeffVsEest->SetMarkerColor(3);
      AeffVsEest->SetMinimum(1.);
      AeffVsEest->SetMaximum(3.e6);
      AeffVsEest->SetFillColor(4);
      AeffVsEest->Draw("E2 same");
    }


  // Calculate the migration matrix for the assumed source spectrum, and weighting it according to the distribution of
  // the data (in terms of obs. time) in zenith and (possibly) azimuth bins.
  energyMigration.SetSpectrum(assumedspectrum);
  energyMigration.SetPhiThetaDistr(&EffOnTimevsAzvsZd);

  if (AeffAtmCorr)
    {
      TH3F* hdetrue = (TH3F*) &((TH3F&)hdeltaetot_etrue->GetHist());
      if (!energyMigration.Calc(hdetrue))
	{
	  gLog << err << "energyMigration.Calc() failed! Exiting..." << endl;
	  return -1;
	}
    }
  else if (!energyMigration.Calc())
    {
      gLog << err << "energyMigration.Calc() failed! Exiting..." << endl;
      return -1;
    }
  TH3D* migrmatrix = (TH3D*) energyMigration.GetHistMigCoarse()->Clone("migrmatrix");
  migrmatrix->SetTitle("Migration matrix");
  // Name needed like this by the unfolding program.


  // Display Migration Matrix
  TCanvas &canvMigMatrix = display->AddTab("Mig. Matrix");
  gWizard->WhiteBack(canvMigMatrix);
  canvMigMatrix.SetRightMargin(0.12);
  canvMigMatrix.SetLogx();
  canvMigMatrix.SetLogy();
  canvMigMatrix.SetLogz();
  canvMigMatrix.SetGridx();
  canvMigMatrix.SetGridy();

  TH2D* migrmatrix2d = (TH2D*) migrmatrix->Project3D("eyx"); // z-axis is just zenith angle (here, one bin only!)
  xaxe = (TAxis*) migrmatrix2d->GetXaxis();
  xaxe->SetTitle("E_{est} (GeV)");
  TAxis *yaxe = (TAxis*) migrmatrix2d->GetYaxis();
  yaxe->SetTitle("E_{true} (GeV)");
  migrmatrix2d->SetTitle("Migration matrix");
  migrmatrix2d->SetStats(kFALSE);
  migrmatrix2d->SetMinimum(1.e-5);
  migrmatrix2d->Draw("zcol");

  // Plot the theta2 plots for the energy bins:
  TCanvas &canvTheta2vsE = display->AddTab("Theta2 plots vs. E");

  gWizard->WhiteBack(canvTheta2vsE);
  Theta2vsEest.DrawAllEbins(processWobblePos, theta2CutVsEest);

  // Create the MHExcessEnergyTheta container where the main output (excess vs. estimated energy and vs. zenith angle)
  // will be stored. Only one zenith angle bin allowed!
  // Create also another object of type MHExcessEnergyTheta where the UPPER LIMITs to the number of events will be stored,
  // to be plotted in energy bins where signal is not significant.


  // First we will create and fill TH2D objects which will be used later in the creation of the MHExcessEnergyTheta objects:
  //
  TH2D* hExcEZd  = new TH2D("MHExcessEnergyTheta", "# Excess events vs. E and Zd", nBinsEnergyEst, eestbinning.GetEdges(), 1, minZd, maxZd);
  TH2D* hExcEZd2 = new TH2D("MHExcessEnergyTheta2","# Excess events vs. E and Zd", nBinsEnergyEst, eestbinning.GetEdges(), 1, minZd, maxZd);
  // NOTE: the label "2" refers to the use of UNcorrected energy, in the case atmospheric corrections are switched ON. The "default" objects
  // like hExcEZd, hEstBckgE or Theta2vsEest will be filled with the (LIDAR-) corrected energies!


  TH2D* hULEZd = new TH2D("NULvsEnergyTheta", "Upper limit of # events vs. E and Zd", nBinsEnergyEst, eestbinning.GetEdges(), 1, minZd, maxZd);

  // Also TH1D's to keep the number of estimated background events, and the significance, vs. Eest:
  TH1D* hEstBckgE = new TH1D("hEstBckgE", "# bckg events vs. E", nBinsEnergyEst, eestbinning.GetEdges());
  hEstBckgE->GetXaxis()->SetTitle("E_{est} (GeV)");
  hEstBckgE->GetYaxis()->SetTitle("Normalized background events");
  TH1D* hEstBckgE2 = (TH1D*) hEstBckgE->Clone("hEstBckgE2");
  // ^^^^ hEstBckgE2 will contain estimated background vs. UN-corrected Eest, in case AeffAtmCorr is kTRUE



  TH1D* hSigniE = new TH1D("hSigniE", "Significance vs. E", nBinsEnergyEst, eestbinning.GetEdges());
  hSigniE->GetXaxis()->SetTitle("E_{est} (GeV)");
  hSigniE->GetYaxis()->SetTitle("Statistical Significance (Li&Ma)");

  TH1D* hOffNormRegionE = new TH1D("hOffNormRegionE", "# OFF events in theta2 normalization region vs. Eest", nBinsEnergyEst, eestbinning.GetEdges());
  hOffNormRegionE->GetXaxis()->SetTitle("E_{est} (GeV)");
  hOffNormRegionE->GetYaxis()->SetTitle("Normalized OFF events in th2 norm. region");
  TH1D* hOnNormRegionE = (TH1D*) hOffNormRegionE->Clone("hOnNormRegionE");
  hOnNormRegionE->SetTitle("# ON events in theta2 normalization region vs. Eest");
  hOnNormRegionE->GetYaxis()->SetTitle("Normalized ON events in th2 norm. region");


  for (Int_t ebin = 1; ebin <= nBinsEnergyEst; ebin++)
    {
      TH1D* hdiff  = new TH1D;
      TH1D* hdiff2 = new TH1D;
      TH1D* hon    = 0;
      TH1D* hoff   = 0;
      TH1D* hon2   = 0;
      TH1D* hoff2  = 0;
      Double_t non,  noff,  delta_noff,  alpha,  significance;
      Double_t non_normregion, noff_normregion, delta_noff_normregion;
      Double_t non2, noff2, delta_noff2, alpha2, significance2;

      Float_t energy = hExcEZd->GetXaxis()->GetBinCenter(ebin);
      Double_t th2cut = (Double_t) HadTheta2Cuts.GetAlphaCut(energy);

      // J. Palacio: 12/09/2017
      // MTheta2vsEest::GetNormalizedHists modified to add thetaMin2 cuts
      Double_t thMin2cut = (Double_t) HadTheta2Cuts.GetAlphaMinCut(energy);

      if (!Theta2vsEest.GetNormalizedHists(Theta2vsEest.GetHist(), processWobblePos, ebin, ebin, hdiff, hon, hoff,
					   th2cut, thMin2cut, &non, &noff, &delta_noff, &significance, &alpha,
					   &non_normregion, &noff_normregion, &delta_noff_normregion, kFALSE))
	{
	  delete hdiff;
	  continue;
	}

      // 20161004, AM: now we calculate the excess from the returned non, noff values (instead of from hdiff), and its uncertainty from
      // propagating sqrt(non) and delta_noff. This is not always the same as getting the uncertainty from the hdiff integral, because the uncertainty
      // delta_noff includes the contribution from the uncertainty in the background normalization factor (in case the Off was normalized o the On,
      // e.g. if you choose Off from Wobble partner background), whereas the hdiff error bars does not include that, because it is an error which is
      // totally correlated among theta2 bins).

      hExcEZd->SetBinContent(ebin, 1, non-noff);
      hExcEZd->SetBinError  (ebin, 1, sqrt(non + delta_noff*delta_noff));

      hEstBckgE->SetBinContent(ebin, noff);
      hEstBckgE->SetBinError(ebin, delta_noff);

      hOnNormRegionE->SetBinContent(ebin, non_normregion);
      hOnNormRegionE->SetBinError(ebin, sqrt(non_normregion));
      hOffNormRegionE->SetBinContent(ebin, noff_normregion);
      hOffNormRegionE->SetBinError(ebin, delta_noff_normregion);

      hSigniE->SetBinContent(ebin, significance);

      delete hdiff;

      // Now the upper limits:
      Double_t nul;
      if (!Theta2vsEest.GetUpperLimit(Theta2vsEest.GetHist(), processWobblePos, ebin, ebin, &HadTheta2Cuts, &nul, kTRUE))
	continue;

      if (TMath::IsNaN(nul))  // <--- Unsuccessful calculation of UL (e.g. 0 event statistics somewhere.
	                      //  TO BE DONE: some estimate is possible, with Aeff and Teff known even if no events are observed
	continue;

      hULEZd->SetBinContent(ebin, 1, nul);
      hULEZd->SetBinError(ebin, 1, 0.);

      // Now the Eest hists, in case aeff atm. corr is used

      if (!AeffAtmCorr)
	continue;

      // If atmospheric corrections are switched on, all of the above has been done with *LIDAR-corrected* energies. We need however some
      // outputs made with the UNcorrected energy, for later steps in the analysis:

      energy = hExcEZd2->GetXaxis()->GetBinCenter(ebin);
      //      th2cut = (Double_t) HadTheta2Cuts.GetAlphaCut(energy);
      //      thMin2cut = (Double_t) HadTheta2Cuts.GetAlphaMinCut(energy);

      // J. Palacio: 12/09/2017
      // MTheta2vsEest::GetNormalizedHists modified to add thetaMin2 cuts
      if (!Theta2vsEest2.GetNormalizedHists(Theta2vsEest2.GetHist(), processWobblePos, ebin, ebin, hdiff2, hon2, hoff2,
					    th2cut, thMin2cut, &non2, &noff2, &delta_noff2, &significance2, &alpha2, 0, 0, 0, kFALSE))
	{
	  delete hdiff2;
	  continue;
	}

      hExcEZd2->SetBinContent(ebin, 1, non2-noff2);
      hExcEZd2->SetBinError  (ebin, 1, sqrt(non2+delta_noff2*delta_noff2));

      hEstBckgE2->SetBinContent(ebin, noff2);
      hEstBckgE2->SetBinError(ebin, delta_noff2);

      delete hdiff2;
    }

  // Put the excesses into an object of type MHExcessEnergyTheta. Note: often, as in this case, "Theta" is used confusingly in
  // MARS class or variable names to refer to Zenith angle! (here there is only one zenith angle bin allowed, anyway).
  MHExcessEnergyTheta ExcessVsEest(hExcEZd);
  ExcessVsEest.SetAngleType(MHAlphaEnergyTheta::kTheta2);

  MHExcessEnergyTheta ExcessVsEest2(hExcEZd2);
  ExcessVsEest2.SetAngleType(MHAlphaEnergyTheta::kTheta2);

  // Put also the upper limits into an object of type MHExcessEnergyTheta.
  MHExcessEnergyTheta UpperLimitVsEest(hULEZd);
  UpperLimitVsEest.SetName("UpperLimitVsEest");
  UpperLimitVsEest.SetAngleType(MHAlphaEnergyTheta::kTheta2);


  TCanvas &canvExcess = display->AddTab("Excess vs. Eest");
  gWizard->WhiteBack(canvExcess);
  TH1D* excessvse = ExcessVsEest.GetHist()->ProjectionX("excessvse", 1, 1, "e");
  excessvse->GetXaxis()->SetTitle("E_{est} (GeV)");
  excessvse->GetYaxis()->SetTitle("Excess events");
  canvExcess.SetLogx();
  canvExcess.SetLogy();
  canvExcess.SetGridx();
  canvExcess.SetGridy();
  excessvse->SetStats(kFALSE);
  excessvse->SetMinimum(0.1);
  excessvse->Draw();
  // Write statistical significance of each bin with a >0 excess of events:
  for (Int_t ebin = 1; ebin <= excessvse->GetNbinsX(); ebin++)
    {
      if (excessvse->GetBinContent(ebin) <= 0.)
	continue;
      if (TMath::IsNaN(hSigniE->GetBinContent(ebin)))
	  continue;

      TLatex* siglabel = new TLatex(excessvse->GetXaxis()->GetBinCenter(ebin),
				    2.*excessvse->GetBinContent(ebin),Form("%.1f \\sigma",
									   hSigniE->GetBinContent(ebin) ));
      siglabel->SetTextAngle(80.);
      siglabel->SetTextSize(0.02);
      siglabel->Draw();
    }


  TCanvas &canvBckg = display->AddTab("Background vs. Eest");
  gWizard->WhiteBack(canvBckg);
  canvBckg.SetLogx();
  canvBckg.SetLogy();
  canvBckg.SetGridx();
  canvBckg.SetGridy();
  hEstBckgE->SetStats(kFALSE);
  hEstBckgE->SetMinimum(0.1);
  hEstBckgE->Draw("e");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate and display the differential flux:
  //
  MHFlux* hFlux = new MHFlux();
  if (SourceRedshift > 0)
    hFlux->SetOpticalDepthEinTeV(tau_vs_log10e);
  hFlux->SetMaximumRelError(maximumRelError); // if error in excess events is larger, no flux point will be calculated
  hFlux->SetAssumedSpectrum(assumedspectrum); // Needed for placing the points at the right energies
  hFlux->SetPropagateEffAreaError(PropagateEffAreaError);
  hFlux->SetFluxCorrectionFactor(fluxCorrectionFactor);

  hFlux->Calc(&ExcessVsEest, &collareaVsEest, hEffTime, kFALSE);
  gWizard->WhiteBack(display->AddTab("dF/dE"));
  hFlux->Draw("gz");

  TGraphAsymmErrors* dfde = hFlux->GetFlux();

  TString formula = "[0]*(" + assumedspectrum->GetExpFormula() + ")";
  // Fit from 50 GeV to 50 TeV
  TF1* fittedassumedspectrum;

  // We do the fit in the range of the light curve:
  Float_t e1 = EminLC;
  Float_t e2 =  (EmaxLCstr=="inf")? 5.e4 : EmaxLC;

  if (SourceRedshift > 0)
    {
      fittedassumedspectrum = new TF1("fittedassumedspectrum", AbsorbedSpectrum, e1, e2, 1);
      fittedassumedspectrum->SetNpx(5e4); // Needed for decent interpolation!! Will be stored as a histogram (then interpolated in lin-lin scale!)
    }
  else
    fittedassumedspectrum = new TF1("fittedassumedspectrum", formula, e1, e2);


  fittedassumedspectrum->SetParameter(0, dfde->Eval(300.) / assumedspectrum->Eval(300.));
  fittedassumedspectrum->SetTitle("Assumed Spectral Shape (normalized to data)");

  dfde->Fit(fittedassumedspectrum, "0QN");

  fittedassumedspectrum->SetLineColor(8);
  fittedassumedspectrum->SetLineWidth(2);
  fittedassumedspectrum->Draw("same");
  TLatex* fitlabel = new TLatex;
  fitlabel->SetTextColor(8);
  fitlabel->SetNDC();
  fitlabel->DrawLatex(0.2, 0.2, "Assumed spectral shape");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate and display the spectral energy distribution:

  gWizard->WhiteBack(display->AddTab("SED"));
  hFlux->DrawSED("z");

  TF1* fittedassumedSED;
  if (SourceRedshift > 0)
    {
      fittedassumedSED = new TF1("fittedassumedSED", AbsorbedSED, 50., 5.e4, 1);
      fittedassumedSED->SetNpx(5e4); // Needed for decent interpolation!! Will be stored as a histogram (then interpolated in lin-lin scale!)
    }
  else
    {
      TString formula2 = "x*x*1.e-6*[0]*(" + assumedspectrum->GetExpFormula() + ")";
      fittedassumedSED = new TF1("fittedassumedSED", formula2, 50., 5.e4);
    }
  fittedassumedSED->SetTitle("Assumed SED shape (normalized to data)");
  fittedassumedSED->SetParameter(0, fittedassumedspectrum->GetParameter(0));

  fittedassumedSED->SetLineColor(8);
  fittedassumedSED->SetLineWidth(2);
  fittedassumedSED->Draw("same");
  fitlabel->DrawLatex(0.2, 0.2, "Assumed spectral shape");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upper Limits:
  //
  MHFlux* FluxUL = new MHFlux("FluxUL", "Flux UL vs. E and Theta");
  if (SourceRedshift > 0)
    FluxUL->SetOpticalDepthEinTeV(tau_vs_log10e);
  FluxUL->SetMaximumRelError(1.e4); // --> do calculations for all spectral points! (see MHFlux::Calc)
  FluxUL->SetAssumedSpectrum(assumedspectrum);
  FluxUL->GetHist()->SetName("FluxUL");
  FluxUL->GetHist()->SetTitle("Flux UL vs. E");
  FluxUL->SetPropagateEffAreaError(PropagateEffAreaError);
  FluxUL->SetFluxCorrectionFactor(fluxCorrectionFactor);
  FluxUL->Calc(&UpperLimitVsEest, &collareaVsEest, hEffTime, kFALSE);

  TH1D* histFlux = hFlux->GetHist()->ProjectionX("histFlux", 1, 1, "e");
  TH1D* histFluxUL = FluxUL->GetHist()->ProjectionX("histFluxUL", 1, 1, "e");

  for (Int_t ebin = 1; ebin <= histFluxUL->GetNbinsX(); ebin++)
    {
      if ( histFlux->GetBinContent(ebin) / histFlux->GetBinError(ebin) > 1./maximumRelError )  // No need to show the upper limit
	continue;

      Double_t edge_low  = histFluxUL->GetXaxis()->GetBinLowEdge(ebin);
      Double_t edge_high = histFluxUL->GetXaxis()->GetBinUpEdge(ebin);
      // Lafferty-Wyatt prescription for point position (place where dPhi/dE equals its average value in the bin):
      Double_t x = assumedspectrum->GetX( assumedspectrum->Integral(edge_low, edge_high) / (edge_high-edge_low),
					  edge_low, edge_high);

      Double_t y = histFluxUL->GetBinContent(ebin);
      if (! (y > 0.) )
	continue;

      display->CdCanvas("dF/dE");

      TArrow* arr = new TArrow(x, y, x, y/5., 0.02);
      arr->Draw();

      ((TH1D*)gPad->FindObject("DifferentialFlux1"))->GetListOfFunctions()->Add(arr);

      Double_t norm = y / assumedspectrum->Eval(x);
      TF1* specline;
      if (SourceRedshift > 0)
	{
	  specline = new TF1(Form("specline_%d", ebin),
			     AbsorbedSpectrum,
			     histFluxUL->GetXaxis()->GetBinLowEdge(ebin),
			     histFluxUL->GetXaxis()->GetBinUpEdge(ebin), 1);
	  specline->SetParameter(0, norm);
	}
      else
	specline = new TF1(Form("specline_%d", ebin),
			   Form("%e*%s", norm, assumedspectrum->GetExpFormula().Data()),
			   histFluxUL->GetXaxis()->GetBinLowEdge(ebin),
			   histFluxUL->GetXaxis()->GetBinUpEdge(ebin));

      specline->SetLineWidth(1);
      specline->Draw("same");


      // Now the UL's in the SED:
      display->CdCanvas("SED");
      TArrow* arr2 = new TArrow(x, y*x*x*1.e-6, x, y*x*x*1.e-6/2., 0.02);  // Converting GeV to TeV since SED is in TeV cm-1 s-1
      arr2->Draw();
      ((TH1D*)gPad->FindObject("SEDframe"))->GetListOfFunctions()->Add(arr2);

      TF1* sedline;
      if (SourceRedshift > 0)
	{
	  sedline = new TF1(Form("sedline_%d", ebin),
			    AbsorbedSED,
			    histFluxUL->GetXaxis()->GetBinLowEdge(ebin),
			    histFluxUL->GetXaxis()->GetBinUpEdge(ebin), 1);
	  sedline->SetParameter(0, norm);
	}
      else
	sedline = new TF1(Form("sedline_%d", ebin),
			  Form("%e*x*x*1.e-6*%s", norm, assumedspectrum->GetExpFormula().Data()),
			  histFluxUL->GetXaxis()->GetBinLowEdge(ebin),
			  histFluxUL->GetXaxis()->GetBinUpEdge(ebin));

      sedline->SetLineWidth(1);
      sedline->Draw("same");
    }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Now we go for the lightcurve calculations. We will apply the same hadronness cuts (vs. E) as in the
  // spectral calculations.
  //

  TGraphErrors* LightCurve = new TGraphErrors;
  LightCurve->SetName("LightCurve");
  TGraphErrors* BackgroundCurve = new TGraphErrors;
  BackgroundCurve->SetName("BackgroundCurve");
  TGraphErrors* AeffLC = new TGraphErrors;
  AeffLC->SetName("AeffLC");
  AeffLC->SetTitle("Total Aeff for LC");
  TGraphErrors* teffLC = new TGraphErrors;
  teffLC->SetName("teffLC");
  teffLC->SetTitle("Effective time");
  TGraphErrors* excessLC = new TGraphErrors;
  excessLC->SetName("excessLC");
  excessLC->SetTitle("Excess events");
  TGraphErrors* offLC = new TGraphErrors;
  offLC->SetName("offLC");
  offLC->SetTitle("Estimated off events");
  TGraph* UpperLimLC = new TGraph;
  UpperLimLC->SetName("UpperLimLC");
  UpperLimLC->SetTitle("Upper Limits in LC bins");
  UpperLimLC->SetMarkerStyle(23);

  if (DoLightCurve)
    {
      // The time bin limits are defined in MJD.
      // If the selected option is "one LC bin = one night", then we obtain from the MHEffectiveOnTime the MJD limits:

      if (LCbinningOption == "night-wise")
	GetNightlyLCbinEdges(hEffTime, &LCbinlowedge, &LCbinupedge);

      if (LCbinningOption == "single-bin")
	GetSingleBinLCbinEdges(hEffTime, &LCbinlowedge, &LCbinupedge);

      Int_t NumLCbins = LCbinlowedge.GetSize();
      gLog << inf << endl << NumLCbins << " light curve bins (in MJD) :" << endl;
      for (Int_t i= 0; i < LCbinlowedge.GetSize(); i++)
	gLog << inf << Form("#%d   %.6f   to  %.6f", i, LCbinlowedge[i], LCbinupedge[i]) << endl;


      TCanvas &canvLCtheta2 = display->AddTab("Theta2, LC");
      gWizard->WhiteBack(canvLCtheta2);
      Int_t npadsLC = (Int_t) (sqrt(LCbinlowedge.GetSize())+1.);
      canvLCtheta2.Divide(npadsLC,npadsLC);

      TH1D** AeffVsEest_LC = new TH1D*[NumLCbins];

      Int_t LCpoint = 0;
      Int_t LC_UL_point = 0;

      if (fluxCorrectionFactor != 1.)
	gLog << warn << "NOTE: Using user-selected flux correction factor: "
	      << Form("%.4f", fluxCorrectionFactor) << endl;

      for (Int_t LCbin = 0; LCbin < NumLCbins; LCbin++)
	{
	  TH1D* hon = new TH1D();
	  TH1D* hoff = new TH1D();
	  TH1D* hdiff = new TH1D();

	  TH2D* efftimevsAzvsZd = (TH2D*) hEffTime->GetHEffOnPhiTheta().Clone("efftimevsAzvsZd");
	  efftimevsAzvsZd->Reset();

	  // Fill the theta2 histos for this LC bin, as well as the eff time vs. Az and Zd which is needed
	  // for the collection area calculation:
	  Theta2vsEest.FillTheta2HistsLC(LCbinlowedge[LCbin], LCbinupedge[LCbin], hEffTime, efftimevsAzvsZd);

	  // Calculate coll. area to be used for this LC bin:
	  if (AeffAtmCorr)
	    {
	      TH3F* hde = (TH3F*) &((TH3F&)hdeltaetot_eest->GetHist());
	      collareaVsEest.Calc(assumedspectrum, efftimevsAzvsZd, hde, kFALSE);
	    }
	  else
	    collareaVsEest.Calc(assumedspectrum, efftimevsAzvsZd, 0, kFALSE);


	  AeffVsEest_LC[LCbin] = collareaVsEest.GetHistCoarse()->ProjectionX(Form("AeffVsEest_LCbin%d", LCbin),1,1,"e");
	  AeffVsEest_LC[LCbin]->SetTitle(Form("Coll. area for LC bin #%d", LCbin));

	  delete efftimevsAzvsZd;


	  // Obtain the overall theta2 histogram in this LC time bin, for all bins between EminLCbin and EmaxLCbin:
	  Double_t non, noff, delta_noff, significance, alpha;
      // J. Palacio: 12/09/2017
      // MTheta2vsEest::GetNormalizedHists modified to add thetaMin2 cuts
	  Theta2vsEest.GetNormalizedHists(Theta2vsEest.GetHistLCbin(), processWobblePos, EminLCbin, EmaxLCbin, hdiff, hon, hoff,
					  -1.,-1, &non, &noff, &delta_noff, &significance, &alpha, 0, 0, 0, kFALSE);
	  // Now plot it:
	  hon->SetName(Form("Theta2_on_LC_bin_#%d", LCbin));
	  hoff->SetName(Form("Theta2_off_LC_bin_#%d", LCbin));
	  hdiff->SetName(Form("Theta2_diff_LC_bin_#%d", LCbin));
	  hon->SetDirectory(0);
	  hoff->SetDirectory(0);
	  hdiff->SetDirectory(0);

	  canvLCtheta2.cd(LCbin+1);
	  gPad->SetGridx();
	  gPad->SetGridy();

	  hon->SetStats(kFALSE);

	  if (EmaxLCstr == "inf")
	    hon->SetTitle(Form("E_{est}> %d GeV, LC bin #%d", (Int_t)TMath::Floor(EminLC+0.5), LCbin));

	  else
	    hon->SetTitle(Form("%d < E_{est} < %d GeV, LC bin #%d",
			       (Int_t)TMath::Floor(EminLC+0.5),
			       (Int_t)TMath::Floor(EmaxLC+0.5), LCbin));

	  hon->SetMinimum(1.1*(hdiff->GetMinimum()-hdiff->GetBinError(hdiff->GetMinimumBin())));
	  hon->Draw("e");
	  hon->GetYaxis()->SetTitle("Events");
	  hon->GetXaxis()->SetTitle("\\theta^{2} [deg^{2}]");
	  hoff->SetLineColor(2);
	  hoff->Draw("same,e");
	  hdiff->SetLineColor(4);
	  hdiff->Draw("same,e");

	  // Now go on with flux calculations:

	  // Total effective time in this LC bin:
	  Double_t teff = hEffTime->GetEffTime(LCbinlowedge[LCbin], LCbinupedge[LCbin]);

	  Double_t total_non = 0.;
	  Double_t total_noff = 0.;
	  Double_t nul = 0.;

	  Double_t integralflux = 0.;
	  Double_t delta_integralflux = 0.;

	  Double_t integralbckg = 0.;
	  Double_t delta_integralbckg = 0.;


	  if (! Theta2vsEest.GetEventStatistics(Theta2vsEest.GetHistLCbin(), processWobblePos, EminLCbin, EmaxLCbin, &HadTheta2Cuts,
						&total_non, &total_noff, &delta_noff, &nul, 0, 0, 0, kFALSE))
	    {
	      gLog << "LC bin #" << LCbin << " : could not calculate a flux point! Theta2vsEest::GetEventStatistics failed."
		   << endl;
	      continue;
	    }
	  Double_t sum_delta_noff_squared = delta_noff * delta_noff;


	  Double_t delta_aeff;  // NOTE!!! : Will be propagated to final LC flux error, BUT the Aeff's are highly correlated for the LC bins!
	  // That means we will have overestimated LC error bars (relevant mostly if flux is high, otherwise real data even statistics dominates!
	  Double_t total_aeff = collareaVsEest.GetTotalAeffInCoarseBinRange(assumedspectrum, EminLCbin, EmaxLCbin, &delta_aeff);

	  // In case of very short-timescale light curves, or light curves (e.g. daily ones on consecutive nights) in which, due to similar ZD-Az
	  // range of the data, each bin gets its Aeff from about the same MC events, it may be better NOT to propagate the uncertainties in Aeff
	  // to the LC points (because the potential error would affect LC points equally, whereas the propagation assumes they can fluctuate
	  // independently):
	  if (PropagateEffAreaErrorToLightCurve == kFALSE)
	    delta_aeff = 0.;

	  if (total_aeff < 1.e-6)    // safety check
	    continue;

	  // In case eff. time calculation failed for some reason!
	  if (teff < 0.)
	    {
	      gLog << warn << "LC bin #" << LCbin << ", eff. time calculation failed! (MJD range out of data range?)" << endl;
	      continue;
	    }

	  // In case teff is 0 or very small, aeff can be 0 and then one gets inf, or a NaN if the excess is 0.
	  // At least for the latter case we set the flux to 0 to avoid trouble:
	  if ( fabs(total_non-total_noff) > 1.e-6)
	    integralflux = (total_non-total_noff)/total_aeff/teff;
	  else
	    integralflux = 0.;

	  // Same as above:
	  if (total_noff > 0)
	    integralbckg = total_noff/total_aeff/teff;
	  else
	    integralbckg = 0.;

	  //	  delta_integralflux = sqrt(total_non+sum_delta_noff_squared)/total_aeff/teff;  // OLD,before 20160602, no Aeff uncertainty
	  //
	  //	  Double_t rel_err = sqrt( (total_non+sum_delta_noff_squared)/pow(total_non-total_noff,2.)  + pow(delta_aeff/total_aeff,2.)) ;
	  //      ^^^^^^^^ OLD, before 20170213 - Relative error. But was wrong! : integralflux = TMath::Abs(integralflux) * rel_err  is 0 when excess is 0!
	  //

	  // 20170213, fixed bug that gave delta_integralflux = 0 when the excess in a LC point was 0:
	  delta_integralflux = sqrt( (total_non+sum_delta_noff_squared)/pow(total_aeff,2.) +
				     pow(delta_aeff*(total_non-total_noff)/pow(total_aeff,2.), 2.))
				/ teff;

	  Double_t total_delta_noff = sqrt(sum_delta_noff_squared);
	  Double_t rel_err = sqrt( sum_delta_noff_squared/pow(total_noff,2.) + pow(delta_aeff/total_aeff,2.));
	  delta_integralbckg = integralbckg * rel_err;  // Here it is ok to use rel_err; if integralbckg is 0 there is no uncertainty we can calculate!


	  integralflux *= 1.e-4; // Convert from m-2 to cm-2
	  delta_integralflux *= 1.e-4;

	  integralbckg *= 1.e-4; // Convert from m-2 to cm-2
	  delta_integralbckg *= 1.e-4;


	  Double_t fluxupperlimit = 0.;
	  if (!TMath::IsNaN(nul))
	    fluxupperlimit = nul/total_aeff/teff * 1.e-4;

	  gLog << "LC bin #" << LCbin << ",  non: " << total_non << ",  noff: " << total_noff << " +/- " << total_delta_noff;
	  if (integralflux/delta_integralflux < 1./maximumRelErrorLC)
	    gLog << ", N_UL = " << nul << endl;
	  else
	    gLog << endl;
	  gLog << "             Coll. Area: " << total_aeff << " m2,  teff = " << teff << "s." << endl;

	  if (fluxCorrectionFactor == 1.)
	    gLog << "             Flux: " << integralflux << " +/- " << delta_integralflux;
	  else
	    gLog << "             Flux: " << Form("%.4f", fluxCorrectionFactor) << " * (" << integralflux << " +/- " << delta_integralflux << ")";

	  if (integralflux/delta_integralflux < 1./maximumRelErrorLC)
	    {
	      if (fluxCorrectionFactor == 1.)
		gLog << "  (UL = " << fluxupperlimit << ") ";
	      else
		gLog << "  (UL = " << Form("%.4f", fluxCorrectionFactor) << "*" << fluxupperlimit << ") ";
	    }

	  gLog << " cm-2 s-1" << endl;

	  /* Listing ON events (by M. Takahashi) */
          Long64_t nOnEvt = 0;
          //Int_t NRequiredEvt = 9;
          Int_t MjdOffset = 55000;
          ofstream file_arrivaltime("EventArrivalTime.csv");
          file_arrivaltime << "MJD,Eest,theta2" << endl;
          for(Long64_t jpoint=0;jpoint<ListOfPointings.GetSize();jpoint++)
            {
              //gLog << "Pointing No." << jpoint << endl;
	      TNtupleD * nt = (TNtupleD*)Theta2vsEest.GetNtuple(jpoint,jpoint);
              Double_t zenith, azimuth, eest, time, dirx, diry, srcx, srcy;

              nt->SetBranchAddress("zenith", &zenith);
              nt->SetBranchAddress("azimuth", &azimuth);
              nt->SetBranchAddress("eest", &eest);
              nt->SetBranchAddress("time", &time);
              nt->SetBranchAddress("srcx", &srcx);
              nt->SetBranchAddress("srcy", &srcy);
              nt->SetBranchAddress("dirx", &dirx);
              nt->SetBranchAddress("diry", &diry);

              Float_t theta2 = -1;

              nt->BuildIndex("time");
              TTreeIndex* index = (TTreeIndex*)nt->GetTreeIndex();
              for( Long64_t ievent = 0; ievent<index->GetN(); ievent++ )
                {
                  //Long64_t local = nt->LoadTree( index->GetIndex()[ievent] );
                  nt->GetEntry(index->GetIndex()[ievent]);
                  //nt->GetEvent(ievent);
                  if (time < LCbinlowedge[LCbin]-MjdOffset || time > LCbinupedge[LCbin]-MjdOffset)
                    continue;
                  if(eest>=EminLC && eest<EmaxLC)
                    {
                      if (TMath::IsNaN(dirx) || TMath::IsNaN(diry))
                        theta2 = -1;
                      Float_t delta_x = dirx - srcx;
                      Float_t delta_y = diry - srcy;
                      theta2 = delta_x*delta_x + delta_y*delta_y;
                      if(theta2<HadTheta2Cuts.GetAlphaCut(eest))
                        {
                          nOnEvt++;
                          //if(nOnEvt%NRequiredEvt==0 || nOnEvt%NRequiredEvt==NRequiredEvt-1)
                          file_arrivaltime << setprecision(11) << MjdOffset+time << "," << setprecision(4) << eest << "," << setprecision(6) << theta2 << endl;
                        }
                    }
                }
            }
          gLog << "Total number of ON event: " << nOnEvt << endl;
          file_arrivaltime.close();
          /* End of modification */



	  if (total_aeff < 1.e-6) // Avoid null values (when too few statistics in bin to calculate weighted average!)
	    continue;

	  integralflux *= fluxCorrectionFactor;
	  delta_integralflux *= fluxCorrectionFactor;
	  fluxupperlimit *= fluxCorrectionFactor;

	  LightCurve->SetPoint(LCpoint,
			       0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
			       integralflux);

	  LightCurve->SetPointError(LCpoint,
				    0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
				    delta_integralflux);

	  BackgroundCurve->SetPoint(LCpoint,
				    0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
				    integralbckg);

	  BackgroundCurve->SetPointError(LCpoint,
					 0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
					 delta_integralbckg);

	  AeffLC->SetPoint(LCpoint,
			   0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
			   total_aeff);
	  AeffLC->SetPointError(LCpoint,
				0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
				delta_aeff);

	  teffLC->SetPoint(LCpoint,
			   0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
			   teff);
	  teffLC->SetPointError(LCpoint,
				0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
				0.);

	  excessLC->SetPoint(LCpoint,
			     0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
			     total_non-total_noff);
	  excessLC->SetPointError(LCpoint,
				  0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
				  sqrt(total_non+sum_delta_noff_squared));

	  offLC->SetPoint(LCpoint,
			  0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
			  total_noff);
	  offLC->SetPointError(LCpoint,
			       0.5*(LCbinupedge[LCbin]-LCbinlowedge[LCbin]),
			       total_delta_noff);

	  LCpoint++;  // Increase the counter of valid LC points.

	  if (integralflux/delta_integralflux < 1./maximumRelErrorLC)  // Set upper limit in graph
	    {
	      if (!TMath::IsNaN(nul))  // <--- Avoid unsuccessful calculations of UL (e.g. 0 event statistics somewhere).
		{
		  UpperLimLC->SetPoint(LC_UL_point++,
				       0.5*(LCbinlowedge[LCbin]+LCbinupedge[LCbin]),
				       fluxupperlimit);
		}
	    }
	}


      TCanvas &canvLCaeff = display->AddTab("Coll area Eest, LC");
      gWizard->WhiteBack(canvLCaeff);
      npadsLC = (Int_t) (sqrt((Float_t)NumLCbins)+1.);
      canvLCaeff.Divide(npadsLC,npadsLC);

      for (Int_t LCbin = 0; LCbin < NumLCbins; LCbin++)
	{
	  canvLCaeff.cd(LCbin+1);

	  gPad->SetGridx();
	  gPad->SetGridy();
	  gPad->SetLogx();
	  gPad->SetLogy();

	  AeffVsEest_LC[LCbin]->SetStats(kFALSE);
	  AeffVsEest_LC[LCbin]->SetMarkerColor(4);
	  AeffVsEest_LC[LCbin]->SetFillColor(4);
	  AeffVsEest_LC[LCbin]->SetMinimum(1.);
	  AeffVsEest_LC[LCbin]->SetMaximum(3.e6);
	  AeffVsEest_LC[LCbin]->Draw("E2");

	  Float_t x1 = AeffVsEest_LC[LCbin]->GetXaxis()->GetBinLowEdge(EminLCbin);
	  TLine* Eminline = new TLine(x1, 1., x1, 3.e6);
	  Eminline->Draw();
	}

      TCanvas &canvLCaux = display->AddTab("LC aux");
      gWizard->WhiteBack(canvLCaux);
      canvLCaux.Divide(2,2);

      canvLCaux.cd(1);
      gPad->SetGridx();
      gPad->SetGridy();
      AeffLC->SetMarkerStyle(20);
      AeffLC->SetMarkerSize(0.5);
      AeffLC->Draw("ap");
      AeffLC->GetXaxis()->SetNdivisions(505);
      AeffLC->GetYaxis()->SetTitle("Aeff (m^{2})");
      canvLCaux.cd(2);
      gPad->SetGridx();
      gPad->SetGridy();
      teffLC->SetMarkerStyle(20);
      teffLC->SetMarkerSize(0.5);
      teffLC->Draw("ap");
      teffLC->GetXaxis()->SetNdivisions(505);
      teffLC->GetYaxis()->SetTitle("teff (s)");
      canvLCaux.cd(3);
      gPad->SetGridx();
      gPad->SetGridy();
      excessLC->SetMarkerStyle(20);
      excessLC->SetMarkerSize(0.5);
      excessLC->Draw("ap");
      excessLC->GetXaxis()->SetNdivisions(505);
      excessLC->GetYaxis()->SetTitle("N_{excess}");
      canvLCaux.cd(4);
      gPad->SetGridx();
      gPad->SetGridy();
      offLC->SetMarkerStyle(20);
      offLC->SetMarkerSize(0.5);
      offLC->Draw("ap");
      offLC->GetXaxis()->SetNdivisions(505);
      offLC->GetYaxis()->SetTitle("N_{off}");

      TCanvas &canvLightCurve = display->AddTab("Light Curve");
      gWizard->WhiteBack(canvLightCurve);
      canvLightCurve.SetGridx();
      canvLightCurve.SetGridy();

      // Crab Nebula flux in same energy range (spectrum from MAGIC-1, ApJ 674):
      // TF1* CrabSpec = new TF1("CrabFlux","6.e-10*pow(x/300.,-2.31-0.26*log10(x/300.))",40,1.e5); // TeV-1 s-1 cm-2

      // Crab Nebula flux in same energy range (spectrum from MAGIC stereo, arXiv:1406.6892   ):
      TF1* CrabSpec2 = new TF1("CrabSpec2","3.23e-11*pow(x/1000.,-2.47-0.24*log10(x/1000.))",40,1.e5); // TeV-1 s-1 cm-2

      // Crab Nebula flux in same energy range (spectrum from MAGIC stereo, arXiv:1409.5594   ):
      TF1* CrabSpec3 = new TF1("CrabSpec3","[0]*pow(x/1000.,[1]+[2]*log10(x/1000.))",40,1.e5); // TeV-1 s-1 cm-2
      CrabSpec3->SetParameters(3.394e-11, -2.511, -0.2143);

      e1 = EminLC;
      e2 = (EmaxLCstr=="inf")? 1.e5 : EmaxLC;

      // NOTE: since V2-16-1, integral flux is from MAGIC stereo paper arXiv:1406.6892:
      Float_t crabflux2 = CrabSpec2->Integral(e1, e2);
      crabflux2 *= 1.e-3; // Correct for the fact that y-axis was in TeV-1, while x-axis in GeV


      // Another reference (with uncertainty): Crab flux from stereo post-upgrade paper:
      Double_t arrMu[3] = {CrabSpec3->GetParameter(0), CrabSpec3->GetParameter(1), CrabSpec3->GetParameter(2)};
      Double_t arrCov[9] = {7.65909e-25 , 2.00659e-15 , -9.69137e-15, 2.00659e-15 , 0.00049751 ,
			    0.000434075, -9.69137e-15 , 0.000434075 , 0.000753692 };

      Float_t crabflux3 = CrabSpec3->Integral(e1, e2) * 1.e-3;
      // 1.e-3  corrects for the fact that y-axis was in TeV-1, while x-axis in GeV
      Float_t errcrabflux3 = CrabSpec3->IntegralError(e1, e2, arrMu , arrCov) * 1.e-3;


      Float_t ymax = TMath::Max(crabflux2, crabflux3+errcrabflux3);
      Float_t ymin = 0.;

      if (LightCurve->GetN() > 0)
	{
	  Float_t lcmax = TMath::MaxElement(LightCurve->GetN(), LightCurve->GetY());
	  if (lcmax > ymax)
	    ymax = lcmax;
	  ymin = TMath::MinElement(LightCurve->GetN(), LightCurve->GetY());
	}

      if (UpperLimLC->GetN() > 0)
	{
	  Float_t lcmax = TMath::MaxElement(UpperLimLC->GetN(), UpperLimLC->GetY());
	  if (lcmax > ymax)
	    ymax = lcmax;
	}

      ymax *= 1.1;
      ymin *= 0.9;

      // Pulsar fluxes are much smaller than 1 Crab.
      if(bckgMode==3 && ymin>0.){
        canvLightCurve.SetLogy();
      }


      TH1* lcframe = canvLightCurve.DrawFrame(LCbinlowedge[0] -
					      0.15*(LCbinupedge[NumLCbins-1]-LCbinlowedge[0]),
					      ymin,
					      LCbinupedge[NumLCbins-1] +
					      0.15*(LCbinupedge[NumLCbins-1]-LCbinlowedge[0]),
					      ymax);
      lcframe->GetXaxis()->SetNdivisions(505);
      lcframe->GetXaxis()->SetTitle("MJD");

      LightCurve->GetXaxis()->SetTitle("MJD");

      if (EmaxLCbin == nBinsEnergyEst)
	lcframe->GetYaxis()->SetTitle(Form("Flux (cm^{-2} s^{-1}) for E > %d GeV",
					   (Int_t)TMath::Floor(EminLC+0.5)));
      else
	lcframe->GetYaxis()->SetTitle(Form("Flux (cm^{-2} s^{-1}) for %d < E < %d GeV",
					   (Int_t)TMath::Floor(EminLC+0.5),
					   (Int_t)TMath::Floor(EmaxLC+0.5) ));
      LightCurve->GetYaxis()->SetTitle(lcframe->GetYaxis()->GetTitle()); // so that axis title is saved...

      LightCurve->SetMarkerStyle(20);
      LightCurve->SetMarkerSize(0.7);
      LightCurve->Draw("p");
      LightCurve->GetXaxis()->SetNdivisions(505);

      // The settings of max and min below should not be necessary, but otherwise sometimes the range is totally wrong:
      UpperLimLC->GetHistogram()->SetMaximum(lcframe->GetYaxis()->GetXmax());
      UpperLimLC->GetHistogram()->SetMinimum(lcframe->GetYaxis()->GetXmin());
      UpperLimLC->GetHistogram()->GetXaxis()->SetLimits(lcframe->GetXaxis()->GetXmin(), lcframe->GetXaxis()->GetXmax());
      UpperLimLC->GetXaxis()->SetNdivisions(505);

      // Draw upper limits:
      if (UpperLimLC->GetN() > 0)
	{
	  Double_t *xul = UpperLimLC->GetX();
	  Double_t *yul = UpperLimLC->GetY();
	  Double_t arrsize = 0.1 * (lcframe->GetYaxis()->GetXmax() - lcframe->GetYaxis()->GetXmin());
	  for (Int_t iul = 0; iul < UpperLimLC->GetN(); iul++)
	    {
	      TArrow* arr = new TArrow(xul[iul], yul[iul], xul[iul], yul[iul]-arrsize, 0.02);
	      Int_t k;
	      for (k = 0; k < NumLCbins; k++)
		if (LCbinupedge[k] > xul[iul])
		  break;
	      TLine* lin = new TLine(LCbinlowedge[k], yul[iul], LCbinupedge[k], yul[iul]);
	      arr->SetLineColor(4);
	      lin->SetLineColor(4);
	      lin->Draw();
	      arr->Draw();
	      lcframe->GetListOfFunctions()->Add(arr);
	    }
	}

      // Plot Crab Nebula flux in same energy range (spectrum from arXiv:1406.6892):

      TLine *crabref2 = new TLine(lcframe->GetXaxis()->GetXmin(), crabflux2, lcframe->GetXaxis()->GetXmax(), crabflux2);
      crabref2->SetLineStyle(2);
      crabref2->SetLineColor(4);
      crabref2->Draw();
      lcframe->GetListOfFunctions()->Add(crabref2);

      TGraph *crabref3 = new TGraph;
      crabref3->SetPoint(0, lcframe->GetXaxis()->GetXmin(), crabflux3-errcrabflux3);
      crabref3->SetPoint(1, lcframe->GetXaxis()->GetXmax(), crabflux3-errcrabflux3);
      crabref3->SetPoint(2, lcframe->GetXaxis()->GetXmax(), crabflux3+errcrabflux3);
      crabref3->SetPoint(3, lcframe->GetXaxis()->GetXmin(), crabflux3+errcrabflux3);

      crabref3->SetFillColor(6);
      crabref3->SetFillStyle(3003);
      crabref3->Draw("f");

      // We now integrate the *tentative* spectral shape (with normalization factor fitted to the data) in
      // the same range of the light curve, for comparison:
      Float_t integrated_dFdE = (EmaxLCstr=="inf")? fittedassumedspectrum->Integral(EminLC, 1.e5) :
	  fittedassumedspectrum->Integral(EminLC, EmaxLC);
      integrated_dFdE *= 1.e-3; // Correct for the fact that y-axis was in TeV-1, while x-axis in GeV
      TLine *ref_integrated_dFdE = new TLine(lcframe->GetXaxis()->GetXmin(), integrated_dFdE,
					     lcframe->GetXaxis()->GetXmax(), integrated_dFdE);
      ref_integrated_dFdE->SetLineStyle(2);
      ref_integrated_dFdE->SetLineColor(2);
      ref_integrated_dFdE->Draw();

      TLegend * leg = new TLegend(0.65,0.85,0.99,0.99);
      leg->SetName("leg");
      leg->AddEntry(crabref2, "Crab, MAGIC, arXiv:1406.6892", "l");
      leg->AddEntry(crabref3, "Crab, MAGIC, arXiv:1409.5594", "f");
      leg->AddEntry(ref_integrated_dFdE, "Integral flux from tentative spectrum (for crosscheck)", "l");
      leg->Draw();

      // Now plot the "light curve" of the residual background:

      TCanvas &canvBackgroundCurve = display->AddTab("Backg Curve");
      gWizard->WhiteBack(canvBackgroundCurve);
      canvBackgroundCurve.SetGridx();
      canvBackgroundCurve.SetGridy();

      ymin = 0.;
      ymax = 1.;
      if (BackgroundCurve->GetN() > 0)
	ymax = 1.1*TMath::MaxElement(BackgroundCurve->GetN(), BackgroundCurve->GetY());

      TH1* backgframe = canvBackgroundCurve.DrawFrame(LCbinlowedge[0] -
						      0.15*(LCbinupedge[NumLCbins-1]-LCbinlowedge[0]),
						      ymin,
						      LCbinupedge[NumLCbins-1] +
						      0.15*(LCbinupedge[NumLCbins-1]-LCbinlowedge[0]),
						      ymax);

      backgframe->GetXaxis()->SetNdivisions(505);
      backgframe->GetXaxis()->SetTitle("MJD");
      backgframe->GetYaxis()->SetTitle(Form("Background %s", lcframe->GetYaxis()->GetTitle()));

      BackgroundCurve->GetXaxis()->SetNdivisions(505);
      BackgroundCurve->SetMarkerStyle(20);
      BackgroundCurve->SetMarkerSize(0.7);
      BackgroundCurve->Draw("p");


      TCanvas &canvSourcePos = display->AddTab("Source position, LC");
      gWizard->WhiteBack(canvSourcePos);
      canvSourcePos.Divide(npadsLC,npadsLC);

      for (Int_t LCbin = 0; LCbin < NumLCbins; LCbin++)
	{
	  TString hname = Form("histpsi_LCbin#%d", LCbin);
	  TH1D* histpsi = new TH1D(hname, hname, 24, -180., 180.);
	  TH1D* hpsi = (TH1D*) histpsi->Clone("hpsi");
          Int_t lop = ListOfPointings.GetSize();
          if(bckgMode==3) lop = 2;
	  for (Int_t ipoint = 0; ipoint < lop; ipoint++)
	    {
	      if (processWobblePos >= 0 && ipoint != processWobblePos)
		continue;

	      TString cut = Form("time>%lf && time<%lf",
				 LCbinlowedge[LCbin]-55000., LCbinupedge[LCbin]-55000.);
	      hpsi->Reset();
	      TNtupleD* nt = (TNtupleD*) Theta2vsEest.GetNtuple(ipoint, ipoint);
	      nt->Project("hpsi", "atan2(srcy,srcx)*57.2957795", cut.Data());
	      histpsi->Add(hpsi);
	    }
	  delete hpsi;
	  canvSourcePos.cd(LCbin+1);
	  gPad->SetGridx();
	  histpsi->SetStats(0);
	  histpsi->Draw("histo");
	}

      // Call again the Calc function of the Aeff vs. Eest with the global AzvsZd distribution, otherwise it will keep the values
      // for the last processed LC bin, and it will be saved as that in the Output file...

      if (AeffAtmCorr)
	{
	  TH3F* hde = (TH3F*) &((TH3F&)hdeltaetot_eest->GetHist());
	  collareaVsEest.Calc(assumedspectrum, &EffOnTimevsAzvsZd, hde, kFALSE);
	}
      else
	collareaVsEest.Calc(assumedspectrum, &EffOnTimevsAzvsZd, 0, kFALSE);
    }  // end of if (DoLightCurve)

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Here we save the outputs right before finishing Flute:

  // Save the status display:
  display->SaveAsRoot(statusOutName.Data());



  // Write the output file (which will be the input of the unfolding programs: fold or CombUnfold)

  TFile* outputfile = new TFile(matrixOutName.Data(), "RECREATE");

  MMarsVersion marsvers("MMarsVersion_flute");
  marsvers.Write();

  hEffTime->Write();
  if (AeffAtmCorr)
    {
      // NOTE: when running with atmospheric corrections ON, the *variable* names using the label "2" through the code refer always to objects
      // which use the UN-corrected estimated energy MEnergyEst.fEnergy. However, when *writing* them to disk we rename them (their "ROOT" name)
      // with their default names, i.e. without the "2". This is because CombUnfold needs some of these, and searches them with thir usual name.

      ExcessVsEest2.Write("MHExcessEnergyTheta"); // use here the *UN-corrected* energy, needed for CombUnfold.C (the energy-correction probability
                                                  // goes into the migration matrix)
      hEstBckgE2->Write("hEstBckgE"); // ditto. This is estimated background vs. UN-corrected estimated energy, for the record
      Theta2vsEest2.Write("Theta2vsEest");   // This MTheta2vsEest object is made with the UN-corrected energies (in case AeffAtmCorr == kTRUE)

      // We also write out the objects which use the corrected energy (will be used by foam when merging flute outputs), with adequate names:
      ExcessVsEest.Write("MHExcessEnergyTheta_AtmCorr");
      hEstBckgE->Write("hEstBckgE_AtmCorr");
      Theta2vsEest.Write("Theta2vsEest_AtmCorr");

      hOnNormRegionE->Write(TString(hOnNormRegionE->GetName())+"_AtmCorr");    // These two are only needed for foam, therefore
      hOffNormRegionE->Write(TString(hOffNormRegionE->GetName())+"_AtmCorr");  // we do not write out the un-corrected version
    }
  else
    {
      ExcessVsEest.Write();
      Theta2vsEest.Write();
      hEstBckgE->Write();
      hOnNormRegionE->Write();
      hOffNormRegionE->Write();
    }


  HadTheta2Cuts.Write();
  assumedspectrum->Write();
  collarea.Write();
  if((McRingAnalysis) & (StoreMcEvents)){
	  ostringstream converttostring;
	  string ring_offset_min, ring_offset_max;
	  string name_collarea;
	  for (int i=0; i< NRing; i++){
		  converttostring.str("");
		  converttostring << RingMcEdges[i];
		  ring_offset_min = converttostring.str();
		  converttostring.str("");
		  converttostring << RingMcEdges[i+1];
		  ring_offset_max = converttostring.str();
		  name_collarea="MHMcCollectionAreaEtrue_RingMC_offset_"+ring_offset_min+"_"+ring_offset_max;
		  collarea_per_ring[i]->Write(name_collarea.c_str());
	  }
  }

  collareaVsEest.Write();
  // NOTE: in case of applied atmospheric correction (AeffAtmCorr==kTRUE) collareaVsEest is the one already corrected!! To be used together with
  // excess events in corrected energy!

  if (UserCuts != "")
      {
	TObjString cuts(UserCuts.Data());
	cuts.Write("UserCuts");
      }

  // J. Palacio: 12/09/2017
  // Write proccesWobble as TObjString in the Output_flute.root file
  if (processWobblePos >= 0)
      {
        TObjString oProcessWobblePos(Form("%i",processWobblePos));
        oProcessWobblePos.Write("processWobblePos");
      }

  energyMigration.Write();
  migrmatrix->Write();


  hFlux->Write();
  hFlux->GetFlux()->Write(); // For easy access to dF/dE points
  hFlux->GetSED()->Write();  // For easy access to SED points

  if (SourceRedshift > 0)
    {
      hFlux->GetDeabsorbedSED()->Write();
      intrinsicspectrum->Write();
    }

  FluxUL->Write();

  if (DoLightCurve)
    {
      LightCurve->Write();
      BackgroundCurve->Write();
      AeffLC->Write();
      teffLC->Write();
      excessLC->Write();
      offLC->Write();
      UpperLimLC->Write();
    }

  if (AeffAtmCorr)
    {
      hdeltaetot_eest->Write();
      hdeltaetot_etrue->Write();
    }

  outputfile->Close();

  if (kQuit || kBatch)
    delete display;
  else
    {
      // From now on each 'Close' means: Terminate the application
      display->SetBit(MStatusDisplay::kExitLoopOnClose);

      // Wait until the user decides to exit the application:
      app.Run(kFALSE);
    }

  gLog << inf << endl << "Flute finished!" << endl;
  TCollection::StartGarbageCollection();  // to avoid crashes

  return 0;
}



//////////////////////////////////////////////////////////////////////////////////
//
// Function to add to an MReadMarsFile (mrmf) all the files in InputFiles whose
// names contain substring.
//
void AddSubsetOfFiles(MReadMarsFile *mrmf, TString InputFiles, TString substring)
{
  TChain ch("RunHeaders");
  ch.Add(InputFiles);
  for (Int_t ifile = 0; ifile < ch.GetListOfFiles()->GetEntries(); ifile++)
    {
      TString filename = ch.GetListOfFiles()->At(ifile)->GetTitle();
      if ( ! filename.Contains(substring) )
	continue;
      mrmf->AddFile(filename);
    }
  return;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Function to read the pointing RA and Dec from the RunHeader of the file called
// fullfilename
//
void ReadRADec(TString fullfilename, Double_t* pointingRA, Double_t* pointingDec)
{
  TChain ch("RunHeaders");
  ch.Add(fullfilename);
  ch.SetBranchStatus("MRawRunHeader_1.*", 1);
  MRawRunHeader* rawrh = 0;
  ch.SetBranchAddress("MRawRunHeader_1.", &rawrh);
  ch.GetEvent(0);

  *pointingRA  = rawrh->GetTelescopeRA()/3600.;  // RA converted to hours
  *pointingDec = rawrh->GetTelescopeDEC()/3600.; // Dec converted to degrees

  ch.ResetBranchAddresses();

  return;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Function to round the requested theta2 cut so that it matches the closest bin
// edge in the theta2 histograms!
//
void RoundTheta2cut (Double_t* theta2Cut, Int_t nBinsTheta2, Double_t maxTheta2)
{
  *theta2Cut = TMath::Max(TMath::Nint(*theta2Cut*Double_t(nBinsTheta2)/maxTheta2), 1)*(maxTheta2/Double_t(nBinsTheta2));
  return;
}



/////////////////////////////////////////////////////////////////////////////////
//
// Find the hadronness cuts which keep, at each bin of estimated energy,
// a given fraction of MC gammas. The fractions are contained in the array
// frac1, and the cut values are returned in the array hadcuts. The size and
// zenith angle cuts are applied in advance, and hence their effect is not
// included in the computation of the efficiency.
//
// Find also the theta2 cuts (after all other cuts, including hadronness) which
// keep a certain fraction of MC gammas. The fractions are set in the array frac2, and
// the cut values are returned  in the theta2cuts array. The cuts are also limited
// to ranges specified in allowed{Had,Theta2}Range.
//
// This function works using non-Mars loops over the Events trees.
//
// The theta2 cuts must match the edges of the theta2 binning chosen by the user,
// so that the efficiency cannot always be the requested one. The theta2 cut is rounded
// such that the efficiency is equal or larger than the requested one. The hadronness
// is much more finely discretized, in steps of 2.5e-4
//
// NOTE: The efficiency of the cuts as a function of TRUE MC energy is not expected
// to be constant. This calculation is made for bins of estimated energy,
// since the cuts in real data can be applied only in bins of estimated energy.
// Therefore even if one sets inside the "frac1" & "frac2" arrays constant
// efficiencies for all Eest bins, the efficiencies vs. Etrue will not be constant.
// (note also that for some events there is no successful energy reconstruction, and
// the efficiency vs. Eest can be tuned only for the events with reconstructed Eest).
//
// Besides, the efficiencies are those of the Theta2 and Hadronness cuts
// applied to the MC without any of the additional user cuts one might set in flute.rc.
// These other cuts will further reduce the overall gamma efficiency of the analysis.
//


Bool_t findCuts(TString mcdata, Float_t sizemin,
		Double_t lowE, Double_t upE, Int_t NbinsE,
		Double_t* frac1, Double_t* hadcuts,
		Double_t* frac2, Double_t* theta2cuts,
		Int_t nBinsTheta2, Double_t maxTheta2,
	        vector<Double_t> allowedHadRange, vector<Double_t> allowedTheta2Range,
		Double_t minZd, Double_t maxZd,
		TString nameEnergyEst, TString nameHadronness,
		TString nameStereoPar, MStatusDisplay* display)
{
  if (frac2[0]< 0)
	  gLog << inf << endl << "Finding hadronness for the selected efficiencies..." << endl << endl;
  else
	  gLog << inf << endl << "Finding hadronness and theta2 cuts for the selected efficiencies..." << endl << endl;
  MGeomCamMagic mgeom;
  Float_t ConvMm2Deg = mgeom.GetConvMm2Deg();

  // Histogram of hadronness vs log10(Eest) :
  TH2F* hadvsloge = new TH2F("hadvsloge", "Hadronness vs. log10 Eest, MC gammas",
			     NbinsE, log10(lowE), log10(upE),
			     4080, -0.01, 1.01);
  hadvsloge->GetXaxis()->SetTitle("log_{10}(E_{est}/GeV)");

  // Histogram of theta2 vs log10(Eest) :
  TH2F* theta2vsloge = new TH2F("theta2vsloge", "Theta2 vs. log10 Eest, MC gammas",
				NbinsE, log10(lowE), log10(upE),
				nBinsTheta2, 0., maxTheta2);
  theta2vsloge->GetXaxis()->SetTitle("log_{10}(E_{est}/GeV)");
  theta2vsloge->GetYaxis()->SetTitle("\\theta^{2} (deg^{2})");

  TChain* chain = new TChain("Events");
  if (!chain->Add(mcdata))
    {
      gLog << err << endl << "Could not load MC files: " << mcdata << endl << endl;
      return kFALSE;
    }
  chain->SetBranchStatus("*", 0);

  chain->SetBranchStatus(Form("%s.fHadronness", nameHadronness.Data()), 1);
  chain->SetBranchStatus(Form("%s.fEnergy",     nameEnergyEst.Data()),  1);
  chain->SetBranchStatus(Form("%s.*",           nameStereoPar.Data()),  1);

  chain->SetBranchStatus("MHillas_1.fSize",    1);
  chain->SetBranchStatus("MHillas_2.fSize",    1);
  chain->SetBranchStatus("MPointingPos_1.fZd", 1);
  chain->SetBranchStatus("MSrcPosCam_1.fX",    1);
  chain->SetBranchStatus("MSrcPosCam_1.fY",    1);


  MHadronness   *mhad    = 0;
  MEnergyEst    *meest   = 0;
  MHillas       *mhil1   = 0;
  MHillas       *mhil2   = 0;
  MPointingPos  *mpoint  = 0;
  MSrcPosCam    *msrcpos = 0;
  MStereoPar    *mstereo = 0;

  chain->SetBranchAddress(Form("%s.", nameHadronness.Data()), &mhad);
  chain->SetBranchAddress(Form("%s.", nameEnergyEst.Data()),  &meest);
  chain->SetBranchAddress(Form("%s.", nameStereoPar.Data()),  &mstereo);
  chain->SetBranchAddress("MHillas_1.", &mhil1);
  chain->SetBranchAddress("MHillas_2.", &mhil2);
  chain->SetBranchAddress("MPointingPos_1.", &mpoint);
  chain->SetBranchAddress("MSrcPosCam_1.", &msrcpos);

  Float_t nentries = chain->GetEntries();

  TCanvas &canvCuts = display->AddTab("Cuts");
  gWizard->WhiteBack(canvCuts);
  canvCuts.Divide(2,1);

  display->StartUpdate();

  gLog << inf << endl << "Determining hadronness cuts..." << endl;

  for (Int_t ievent = 0; ievent < nentries; ievent++)
    {
      if (ievent % 10 == 0)
	display->SetProgressBarPosition((Float_t)ievent/(Float_t)nentries);
      gSystem->ProcessEvents();

      chain->GetEvent(ievent);

      if (mpoint->GetZd() < minZd || mpoint->GetZd() > maxZd)
	continue;
      if (mhil1->GetSize() < sizemin)
	continue;
      if (mhil2->GetSize() <sizemin)
	continue;

      Float_t energy = meest->GetEnergy();
      if (energy < lowE || energy > upE)
	continue;
      Float_t log10energy = log10(energy);

      hadvsloge->Fill(log10energy, mhad->GetHadronness());
    }

  canvCuts.cd(1);
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetRightMargin(0.15);
  hadvsloge->SetStats(kFALSE);
  hadvsloge->Draw("zcol");


  TGraph *grhadcut = new TGraph;
  grhadcut->SetMarkerStyle(23);
  grhadcut->SetMarkerSize(1);
  Int_t graphpoint = 0;

  for (Int_t ebin = 1; ebin <= NbinsE; ebin++)
    {
      hadcuts[ebin-1] = -1.;
      // If no cut can be determined, e.g. because there is no MC in a certain energy bin,
      // the cut value is set to -1 so that no event survives the cut on the real data.

      TH1D* slice = hadvsloge->ProjectionY("hadslice", ebin, ebin);

      if (slice->Integral() < 1)
        continue;

      // If MC is too sparse to judge, take highest possible had cut (at lowest energy bin),
      // or whatever was the outcome of the previous bin
      if(slice->Integral() < 10){
	if(ebin > 1 && hadcuts[ebin-2]>0){
	  hadcuts[ebin-1] = hadcuts[ebin-2];
	}
	else{
	  hadcuts[ebin-1] = allowedHadRange.at(1);
	}
	grhadcut->SetPoint(graphpoint++, hadvsloge->GetXaxis()->GetBinCenter(ebin), hadcuts[ebin-1]);
	continue;
      }


      // Loop in hadronness bins:
      for (Int_t hbin = 1; hbin <= hadvsloge->GetNbinsY(); hbin++)
	{

	  // Check if fraction of events after hadronness cut (integral
	  // from 1 to nbin), over the total number (before hadronness cut)
	  // is larger than the required value stored in frac1:

	  if (slice->Integral(1, hbin) > frac1[ebin-1]*slice->Integral())
	    {
	      hadcuts[ebin-1] = TMath::Min(TMath::Max(slice->GetXaxis()->GetBinUpEdge(hbin), allowedHadRange.at(0)), allowedHadRange.at(1));
	      grhadcut->SetPoint(graphpoint++, hadvsloge->GetXaxis()->GetBinCenter(ebin), hadcuts[ebin-1]);
	      break;
	    }
	}
      delete slice;
    }

  grhadcut->Draw("p");
  // If frac2 <0, only the hadronness cut i calculated from efficiency.
  if (frac2[0]< 0)
  {
	  canvCuts.Update();
	  display->StopUpdate();
	  return kTRUE;
  }
  /////////////////////////////////////////////////////////////////////////////
  // Now fill the theta2 histogram after hadronness cuts,
  // to determine the cut values:
  //

  gLog << inf << "Determining theta2 cuts..." << endl << endl;

  for (Int_t ievent = 0; ievent < nentries; ievent++)
    {
      if (ievent % 10 == 0)
	display->SetProgressBarPosition((Float_t)ievent/(Float_t)nentries);
      gSystem->ProcessEvents();

      chain->GetEvent(ievent);

      if (mhil1->GetSize() < sizemin)
	continue;
      if (mhil2->GetSize() < sizemin)
	continue;

      Float_t energy = meest->GetEnergy();
      if (energy < lowE || energy > upE)
	continue;
      Float_t log10energy = log10(energy);

      Int_t ebin = TMath::Max(1, hadvsloge->GetXaxis()->FindBin(log10energy));

      if (mhad->GetHadronness() > hadcuts[ebin-1])
	continue;

      Double_t theta2 = pow(mstereo->GetDirectionX()-msrcpos->GetX()*ConvMm2Deg, 2) +
	                pow(mstereo->GetDirectionY()-msrcpos->GetY()*ConvMm2Deg, 2);

      theta2vsloge->Fill(log10energy, theta2);
    }

  chain->ResetBranchAddresses();

  canvCuts.cd(2);
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  theta2vsloge->SetStats(kFALSE);
  theta2vsloge->Draw("zcol");

  TGraph *grth2cut = new TGraph;
  grth2cut->SetMarkerStyle(23);
  grth2cut->SetMarkerSize(1);
  graphpoint = 0;

  for (Int_t ebin = 1; ebin <= NbinsE; ebin++)
    {
      // Default theta2 cut value:
      theta2cuts[ebin-1] = -1.;
      // If no cut can be determined, e.g. because there is no MC in a certain energy bin,
      // the cut value is set to -1 so that no event survives the cut on the real data.

      // Slice of theta2 histogram for a given energy:
      //
      TH1D* slice = theta2vsloge->ProjectionY("py",ebin,ebin);

      if (slice->Integral() < 1)
        continue;

      // If MC is too sparse to judge, take highest possible had cut (at lowest energy bin),
      // or whatever was the outcome of the previous bin
      if(slice->Integral() < 10)
	{
	  if(ebin > 1 && theta2cuts[ebin-2]>0)
	    theta2cuts[ebin-1] = theta2cuts[ebin-2];
	  else
	    theta2cuts[ebin-1] = slice->GetXaxis()->GetBinUpEdge(slice->GetXaxis()->FindBin(allowedTheta2Range.at(1)-1.e-6));
	  // (round to theta2 histogram bin edges: the -1.e-6 is for the case in which allowedTheta2Range.at(1)
	  // is at the edge of a bin)

	  grth2cut->SetPoint(graphpoint++, theta2vsloge->GetXaxis()->GetBinCenter(ebin), theta2cuts[ebin-1]);
	  continue;
	}


      for (Int_t th2bin = 1; th2bin <= theta2vsloge->GetNbinsY(); th2bin++)
	{
	  // Check if fraction integrated from 1 to th2bin is larger than
	  // the required value stored in frac2:

	  if (slice->Integral(1, th2bin) > frac2[ebin-1] * slice->Integral())
	    {
  	      Float_t theta2cut = TMath::Min(TMath::Max(slice->GetXaxis()->GetBinUpEdge(th2bin), allowedTheta2Range.at(0)), allowedTheta2Range.at(1));
	      // Make sure it matches the edge of one of the theta2 histogram bins:
	      theta2cuts[ebin-1] = slice->GetXaxis()->GetBinUpEdge(slice->GetXaxis()->FindBin(theta2cut-1.e-6));
	      grth2cut->SetPoint(graphpoint++, theta2vsloge->GetXaxis()->GetBinCenter(ebin), theta2cuts[ebin-1]);
	      break;
	    }
	}

      delete slice;
    }

  grth2cut->Draw("p");
  canvCuts.Update();
  display->StopUpdate();

  return kTRUE;
}

//
// Return MJD bin limits of bin "bin" in histogram "hist" which has time in the x-axis
//
void GetMjdBinLimits(const TH1D& hist, Int_t bin, Double_t* minmjd, Double_t* maxmjd)
{
  if (hist.GetXaxis()->GetTimeDisplay()) // units in the axis are ROOT's "axis time"
    {
      MTime dummyt;
      dummyt.SetAxisTime(hist.GetXaxis()->GetBinLowEdge(bin));
      *minmjd = dummyt.GetMjd();

      dummyt.SetAxisTime(hist.GetXaxis()->GetBinUpEdge(bin));
      *maxmjd = dummyt.GetMjd();
    }
  else
    {
      *minmjd = hist.GetXaxis()->GetBinLowEdge(bin);
      *maxmjd = hist.GetXaxis()->GetBinUpEdge(bin);
    }

  return;
}

//
// Calculate LC bin edges (in MJD) for the case of night-wise ligth curve:
//
void GetNightlyLCbinEdges(MHEffectiveOnTime* hEffTime, TArrayD* lowedge, TArrayD* upedge)
{
  lowedge->Reset();
  upedge->Reset();

  const TH1D& timehist = hEffTime->GetHEffOnTime();
  Double_t minmjd, maxmjd;

  // Note: first and last bin of timehist are empty (they are there just for display purposes):
  GetMjdBinLimits(timehist, 2, &minmjd, &maxmjd);

  lowedge->Set(1);
  upedge->Set(1);
  (*lowedge)[0] = minmjd;
  (*upedge)[0]  = maxmjd;

  Int_t prevmjd = (Int_t)(maxmjd+0.5);  // rounded to midnight's MJD value

  Int_t LCbin = 0;
  for (Int_t ibin = 2; ibin < timehist.GetNbinsX(); ibin++)
    {
      GetMjdBinLimits(timehist, ibin, &minmjd, &maxmjd);

      // Is this maxmjd already in the next night?
      if ( ((Int_t)(maxmjd+0.5)) == prevmjd )  // Not yet!
	{
	  if (timehist.GetBinContent(ibin) < 0.1) // : empty bin - ignore
	    continue;

	  (*upedge)[LCbin] = maxmjd;  // set upper end of LC bin

	  continue;
	}

      // If we arrived here, we are in a new night. Move forward until you find a non-empty bin,
      // then set low edge of new LC bin, and go on:

      while(timehist.GetBinContent(ibin) < 0.1)
	{
	  ibin++;
	  if (ibin >= timehist.GetNbinsX())
	    break;
	}

      if (ibin >= timehist.GetNbinsX())
	break;

      GetMjdBinLimits(timehist, ibin, &minmjd, &maxmjd);

      LCbin++;
      lowedge->Set(LCbin+1);
      upedge->Set(LCbin+1);
      (*lowedge)[LCbin] = minmjd;
      (*upedge)[LCbin] = maxmjd;

      prevmjd = (Int_t)(maxmjd+0.5);
    }

  return;
}

//
// Calculate LC bin edges (in MJD) for the case of a single bin ligth curve:
//
void GetSingleBinLCbinEdges(MHEffectiveOnTime* hEffTime, TArrayD* lowedge, TArrayD* upedge)
{
  lowedge->Reset();
  upedge->Reset();

  const TH1D& timehist = hEffTime->GetHEffOnTime();
  Double_t minmjd, maxmjd;

  // Note: first and last bin of timehist are empty (they are there just for display purposes):
  GetMjdBinLimits(timehist, 2, &minmjd, &maxmjd);
  lowedge->Set(1);
  (*lowedge)[0] = minmjd;

  GetMjdBinLimits(timehist, timehist.GetNbinsX()-1, &minmjd, &maxmjd);
  upedge->Set(1);
  (*upedge)[0] = maxmjd;

  return;
}

////////////////////////////////////////////////////////////////////////////
//
// Read in values from configuration file and fill in an array
//
void FillMatrix(Double_t *x, Int_t nmax, MEnv* environment, const Char_t* flag)
{
  TString dummy = environment->GetValue(flag, "");

  if (dummy.IsNull())
    {
      for (Int_t i = 0; i < nmax; i++)
	x[i] = 0.;
      return;
    }
  Int_t icount = 0;
  while(1)
    {
      if (icount > nmax-1)
	{
	  gLog << inf << "Note: more than " << nmax << " initializers for " << flag << endl;
	  gLog << inf << "The extra ones will be ignored!" << endl;
	      break;
	}

      Ssiz_t siz = dummy.First(",");
      if (siz <= 0)
      siz = dummy.Length();
      TSubString substr = dummy(0, siz); //DD fixing the bug with the reading of two values in same line of inputcard (e.g.: za_edges=1.,29.)

      Char_t *endptr;
      x[icount] = strtod(substr.Data(), &endptr);

      if (endptr > substr.Data())  // Something was read; Otherwise the two pointers are equal.
	icount++;

      if (dummy.First(",") < 0)
	break;

      dummy.Remove(0, dummy.First(",")+1);
    }
  if (icount != nmax)
    {
      gLog << err << "ERROR, less than " << nmax << " initializers for " << flag << ". Exiting!" << endl;
      exit(-1);
    }

  return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// IRF calculation done in two parts: first loop over the original MC and second loop over the reconstructed one
  /// If the diffuse MCs are used, then the collective area is calculated in the different predefined ring.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t Compute_IRF(MHMcCollectionArea *collarea,MBinning *BinningDE, MBinning *BinningImpact, Double_t minZd, Double_t maxZd, MBinning *zdbinning, MBinning *eestbinning, MBinning *etruebinning,TString InputFilesMC,TString statusOutName,
		TString posContainerName,  MStatusDisplay* display, MFHadAlpha *mfhadtheta2cuts, MHadAlphaCut HadTheta2Cuts,MContinue *mfconthadtheta2,
		Bool_t AeffAtmCorr,TString nameEnergyEst,TF1 *assumedspectrum, Bool_t EnergyReCalc, MContinue *filterSize, MContinue*  usercuts, MHMcEnergyMigration *energyMigration,
		MH3 *hdeltaetot_eest,MH3 *hdeltaetot_etrue, TH2D EffOnTimevsAzvsZd, Bool_t IsStandardAnalysis,MHMcCollectionArea *collareaVsEest){

	MReadMarsFile readMC("OriginalMC");
	if (!readMC.AddFile(InputFilesMC)){
		display->SaveAsRoot(statusOutName.Data());
		return -1;
	}
	readMC.EnableBranch("MMcEvtBasic_1.*");
	readMC.EnableBranch("MSrcPosCam_1.*");
	readMC.EnableBranch("MSrcPosCam_2.*");
	readMC.VetoBranch("MBadPixelsCam_1");
	readMC.VetoBranch("MBadPixelsCam_2");

	MReadMarsFile readMC2("Events", InputFilesMC);
	readMC2.EnableBranch("MMcEvt_1.*");
	readMC2.EnableBranch("MHillas_1.fSize");
	readMC2.EnableBranch("MHillas_2.fSize");
	readMC2.EnableBranch(TString(posContainerName)+".*");
	readMC2.EnableBranch("MEnergyEst.fEnergy");
	readMC2.EnableBranch("MEnergyEst_1.fEnergy");
	readMC2.EnableBranch("MEnergyEst_2.fEnergy");
	readMC2.EnableBranch("MEnergyEst_1.fUncertainty");
	readMC2.EnableBranch("MEnergyEst_2.fUncertainty");
	readMC2.VetoBranch("MBadPixelsCam_1");
	readMC2.VetoBranch("MBadPixelsCam_2");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Now the loops over the Monte Carlo files to calculate the collection area. We weight the MC according to the distribution of the
	// data in zenith angle, i.e. the weight is proportional to the effective on-time spent in a given bin of zenith angle in the
	// analyzed data sample.
	// We need two loops: one over the "OriginalMC" tree, which contains all the generated MC events (even those which did not even
	// trigger the telescope), and one over the usual "Events" tree in which the analysis cuts will be applied.
	//
	// First MC loop: go through the OriginalMC tree, to obtain the number of generated MC events in each bin of (true) energy and
	// zenith angle.
	//

	MContinue *MinZdCut = new MContinue(Form("MMcEvtBasic_1.fTelescopeTheta*kRad2Deg<%.2f", minZd), "MinZdCut");
	MContinue *MaxZdCut = new MContinue(Form("MMcEvtBasic_1.fTelescopeTheta*kRad2Deg>%.2f", maxZd), "MaxZdCut");

	MParList *MCParList = new MParList();
	MTaskList *MCTaskList = new MTaskList();
    MCParList->AddToList(etruebinning);
    MCParList->AddToList(zdbinning);
    MCParList->AddToList(collarea);
    if(IsStandardAnalysis) {
    		MCParList->AddToList(collareaVsEest);
    		MCParList->AddToList(eestbinning);
    	}

    MCTaskList->AddToList(&readMC);


    // zenith angle cuts:
    MCTaskList->AddToList(MinZdCut);
    MCTaskList->AddToList(MaxZdCut);

    if(collarea->IsMcRingAnalysis()){
    	mfhadtheta2cuts->SetCutType(MFHadAlpha::kOffset);
    	MCParList->AddToList(&HadTheta2Cuts);
    	mfconthadtheta2->SetFilter(mfhadtheta2cuts);
    	MCTaskList->AddToList(mfconthadtheta2);
    }

	MMcCollectionAreaCalc collareaCalc("collareaCalc", "", "BinningE");
	collareaCalc.SetCallCalcInPostProcess(kFALSE);
	collareaCalc.SetCollectionAreaName("MHMcCollectionAreaEtrue");
	MCTaskList->AddToList(&collareaCalc);
	collareaCalc.SetSpectrum(assumedspectrum);
	MMcCollectionAreaCalc collareaCalcVsEest("collareaCalcVsEest", "", "BinningEest");
	if(IsStandardAnalysis) {
		collareaCalcVsEest.SetCallCalcInPostProcess(kFALSE);
		collareaCalcVsEest.SetCollectionAreaName("collareaVsEest");
		if (AeffAtmCorr)
			collareaCalcVsEest.SetNameEnergyEst("MEnergyEst");
		else
			collareaCalcVsEest.SetNameEnergyEst(nameEnergyEst);

		MCTaskList->AddToList(&collareaCalcVsEest);
		collareaCalcVsEest.SetSpectrum(assumedspectrum);
	}
	MEvtLoop MCEvtLoop;
    MCEvtLoop.SetDisplay(display);
    MCParList->AddToList(MCTaskList);
    MCEvtLoop.SetParList(MCParList);
    // Run the loop:
    if (!MCEvtLoop.Eventloop())
  	  return -1;
    MCTaskList->PrintStatistics();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Second MC loop: go through the Events tree, to obtain the number of MC events after all cuts, in each bin of (true and estimated)
	// energy and zenith angle, in order to calculate the collection area. We also fill the energy migration matrix (essentially, a histogram
	// of estimated energy vs. true energy, to be used later in the unfolding).
	//
	delete MinZdCut;
	delete MaxZdCut;
	MinZdCut = new MContinue(Form("MMcEvt_1.fTelescopeTheta*kRad2Deg<%.2f", minZd), "MinZdCut");
	MaxZdCut = new MContinue(Form("MMcEvt_1.fTelescopeTheta*kRad2Deg>%.2f", maxZd), "MaxZdCut");

	MParList *MCParList2 = new MParList();
	MTaskList *MCTaskList2 = new MTaskList();
	MCParList2->AddToList(collarea);
	if(IsStandardAnalysis) MCParList2->AddToList(collareaVsEest);
	MCParList2->AddToList(zdbinning);
	if(IsStandardAnalysis) MCParList2->AddToList(eestbinning);
	MCParList2->AddToList(etruebinning); // The binnings needed for the migration matrix!


    MCTaskList2->AddToList(&readMC2);
    // Zd cuts:

    MCTaskList2->AddToList(MinZdCut);
    MCTaskList2->AddToList(MaxZdCut);

	MAverageEnergy averageEnergies; // Task to re-calculate the average of M1 & M2 energies
	if (EnergyReCalc)
	{
		averageEnergies.SetE1scaling(1.);
		averageEnergies.SetE2scaling(1.);  // <== Just in case these were modified for the real data!!
		MCTaskList2->AddToList(&averageEnergies);
	 }
	// Size cuts:
	MCTaskList2->AddToList(filterSize);

	// User cuts:
	if (usercuts)
		MCTaskList2->AddToList(usercuts);


	 if(collarea->IsMcRingAnalysis()){
		  // Now the hadronness, theta2 and offset cuts:
		  	    mfhadtheta2cuts->SetCutType(MFHadAlpha::kHadTheta2Offset);
	  }
	  else{
		  // Now the hadronness and theta2 cuts:
		    mfhadtheta2cuts->SetCutType(MFHadAlpha::kHadTheta2Stereo);
	  }
	  if (AeffAtmCorr)
		  // this ensures that the cuts are applied in reconstructed energy,
		  // and the effective area is that of reconstructed energy.
		  // The correction on effective area due to the shift from MEnergyEst to MEnergyEstAtmCorr is
		  // applied later on
		  mfhadtheta2cuts->SetNameEnergyEst("MEnergyEst");
	  else
		  mfhadtheta2cuts->SetNameEnergyEst(nameEnergyEst);
	  MCParList2->AddToList(&HadTheta2Cuts);
	  mfconthadtheta2->SetFilter(mfhadtheta2cuts);
	  MCTaskList2->AddToList(mfconthadtheta2);

	  // The same MMcCollectionAreaCalc task as before, now will fill MHMcCollectionArea::fHistSel
	  // and call MHMcCollectionArea::Calc in the PostProcess
	  MCTaskList2->AddToList(&collareaCalc);
	  MFillH migmatfill("MHMcEnergyMigration", "MMcEvt_1");
	  if(IsStandardAnalysis){
	  		  collareaCalcVsEest.SetUseEnergyEst(kTRUE);
	  		  MCTaskList2->AddToList(&collareaCalcVsEest);

	  		  // Now the energy migration matrix:
	  		  if (AeffAtmCorr)
	  			  // this ensures that the energy migration is calculated in terms of reconstructed energy.
	  			  // The correction on the energy migration due to the shift from MEnergyEst to MEnergyEstAtmCorr is
	  			  // applied later on
	  			  energyMigration->SetNameEnergyEst("MEnergyEst");
	  		  else
	  			  energyMigration->SetNameEnergyEst(nameEnergyEst);
	  		  energyMigration->SetHistCol(collarea->GetHist());
	  		  // => so that we have the same binning in the coll. area in Etrue (as needed by the unfolding)
	  		  MCParList2->AddToList(energyMigration);
	  		  // The binning below are needed by MHMcEnergyMigration:
	  		  MCParList2->AddToList(BinningDE);
	  		  MCParList2->AddToList(BinningImpact);
	  		  // task to fill the migration matrix:
	  		  migmatfill.SetBit(MFillH::kDoNotDisplay);
	  		  MCTaskList2->AddToList(&migmatfill);
	  	  }

	  MCParList2->AddToList(MCTaskList2);
	  MEvtLoop MCEvtLoop2;
	  MCEvtLoop2.SetParList(MCParList2);
	  MCEvtLoop2.SetDisplay(display);
	  if (!MCEvtLoop2.Eventloop())
		  return kFALSE;
	 MCTaskList2->PrintStatistics();

	  if (AeffAtmCorr)
	  {
		  if(IsStandardAnalysis){
		  TH3F* hde = (TH3F*) &((TH3F&)hdeltaetot_eest->GetHist());
		  collareaVsEest->Calc(assumedspectrum, &EffOnTimevsAzvsZd, hde);
		  }
		  TH3F* hdetrue = (TH3F*) &((TH3F&)hdeltaetot_etrue->GetHist());
		  collarea->Calc(assumedspectrum, &EffOnTimevsAzvsZd, hdetrue);
	  }
	  else
	  {
		  if(IsStandardAnalysis) collareaVsEest->Calc(assumedspectrum, &EffOnTimevsAzvsZd);
		  collarea->Calc(assumedspectrum, &EffOnTimevsAzvsZd);

	  }

	 delete MCTaskList;
	 delete MCTaskList2;
	 delete MinZdCut;
	 delete MaxZdCut;
	 return kTRUE;

	}
