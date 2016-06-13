// ROOT stuff
#include "TRandom.hh"
#include "TCanvas.h" 
#include "TFile.h" 
#include "TH1F.h" 
#include "TH2F.h" 
#include "TStyle.h" 
#include "TRandom3.hh" 
#include "TLegend.h" 
#include "TText.h" 
#include "TLine.h" 

// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include "Variable.hh" 
#include "PolynomialPdf.hh" 
#include "DalitzPlotPdf.hh" 
#include "DalitzVetoPdf.hh" 
#include "ResonancePdf.hh" 
#include "AddPdf.hh"
#include "ProdPdf.hh"
#include "GooPdf.hh" 
#include "FitManager.hh" 
#include "UnbinnedDataSet.hh"

using namespace std;

TCanvas* foo; 
TCanvas* foodal; 
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
UnbinnedDataSet* data = 0; 

Variable* m12 = 0;
Variable* m13 = 0;
Variable* eventNumber = 0; 
bool fitMasses = false; 
Variable* fixedRhoMass  = new Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable* fixedRhoWidth = new Variable("rho_width", 0.1503, 0.01, 0.1, 0.2); 

const fptype _mD0 = 1.86484; 
const fptype _mD02 = _mD0 *_mD0;
const fptype _mD02inv = 1./_mD02; 
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component. 
Variable* motherM = new Variable("motherM", _mD0);
Variable* chargeM = new Variable("chargeM", piPlusMass);
Variable* neutrlM = new Variable("neutrlM", piZeroMass);
Variable* massSum = new Variable("massSum", _mD0*_mD0 + 2*piPlusMass*piPlusMass + piZeroMass*piZeroMass); // = 3.53481 
Variable* constantOne = new Variable("constantOne", 1); 
Variable* constantZero = new Variable("constantZero", 0); 

GooPdf* kzero_veto = 0; 

fptype cpuGetM23 (fptype massPZ, fptype massPM) {
  return (_mD02 + piZeroMass*piZeroMass + piPlusMass*piPlusMass + piPlusMass*piPlusMass - massPZ - massPM); 
}

void getToyData (std::string toyFileName) {
  TH2F dalitzplot("dalitzplot", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit); 
  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  data = new UnbinnedDataSet(vars); 

  int len = 2048;
  char tmp[len];

  std::ifstream reader;
  reader.open(toyFileName.c_str()); 
  std::string buffer;
  while (!reader.eof()) {
    reader >> buffer;
    if (buffer == "====") break; 
    std::cout << buffer; 
  }

  double dummy = 0; 
  while (!reader.eof()) {
    reader.getline (tmp, len, '\n');
    /*
    reader >> dummy;
    reader >> dummy;      // m23, m(pi+ pi-), called m12 in processToyRoot convention. 
    reader >> m12->value; // Already swapped according to D* charge. m12 = m(pi+pi0)
    reader >> m13->value;

    // Errors on Dalitz variables
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    reader >> dummy; // Decay time
    reader >> dummy; // sigma_t

    reader >> dummy; // Md0
    reader >> dummy; // deltaM
    reader >> dummy; // ProbSig
    reader >> dummy; // Dst charge
    reader >> dummy; // Run
    reader >> dummy; // Event
    reader >> dummy; // Signal and four bkg fractions. 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
    // by throwing out events with a probability linear in m12. 
    // NB! This exercise continues below. 

    // EXERCISE 2: Instead of the above efficiency, impose a 
    // K0 veto, by throwing out events with 0.475 < m23 < 0.505. 

    // EXERCISE 3: Use both the above. 

    //eventNumber->value = data->getNumEvents(); 
    //data->addEvent(); 
    */

    //sscanf the buffer, 20 elements
    //sscanf (tmp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy, &dummy, &m12->value, &m13->value, &dummy, &dummy,
    //	&dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy);
    sscanf (tmp, "%lf %lf %lf %lf", &dummy, &dummy, &m12->value, &m13->value);

    std::vector <fptype> list;
    list.push_back (m12->value);
    list.push_back (m13->value);
    list.push_back (data->getNumEvents ());
    data->insertEventVector(list);

    dalitzplot.Fill(m12->value, m13->value); 
  }

#if 0
  dalitzplot.SetStats(false); 
  dalitzplot.Draw("colz");
  foodal->SaveAs("dalitzplot.png"); 
#endif
}

GooPdf* makeKzeroVeto () {
  if (kzero_veto) return kzero_veto; 

  VetoInfo* kVetoInfo = new VetoInfo();
  kVetoInfo->cyclic_index = PAIR_23; 
  kVetoInfo->minimum = new Variable("veto_min", 0.475*0.475);
  kVetoInfo->maximum = new Variable("veto_max", 0.505*0.505);
  vector<VetoInfo*> vetos; vetos.push_back(kVetoInfo); 
  kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos); 
  return kzero_veto;
}

DalitzPlotPdf* makeSignalPdf (GooPdf* eff = 0) {
  DecayInfo* dtop0pp = new DecayInfo();
  dtop0pp->motherMass  = _mD0; 
  dtop0pp->daug1Mass  = piZeroMass;
  dtop0pp->daug2Mass  = piPlusMass;
  dtop0pp->daug3Mass  = piPlusMass;
  dtop0pp->meson_radius  = 1.5; 
 
  ResonancePdf* rhop  = new ResonancePdf("rhop",
							     new Variable("rhop_amp_real", 1),
							     new Variable("rhop_amp_imag", 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_12);


  bool fixAmps = false;

  ResonancePdf* rhom  = new ResonancePdf("rhom", 
							     fixAmps ? new Variable("rhom_amp_real", 0.714) : 
							     new Variable("rhom_amp_real",  0.714, 0.001, 0, 0),
							     fixAmps ? new Variable("rhom_amp_imag", -0.025) :
							     new Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_13);

  ResonancePdf* rho0  = new ResonancePdf("rho0", 
							     fixAmps ? new Variable("rho0_amp_real", 0.565) : 
							     new Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
							     fixAmps ? new Variable("rho0_amp_imag", 0.164) :
							     new Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_23);

  Variable* sharedMass = new Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
  Variable* shareWidth = new Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0); 

  ResonancePdf* rhop_1450  = new ResonancePdf("rhop_1450", 
								  fixAmps ? new Variable("rhop_1450_amp_real", -0.174) : 
								  new Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1450_amp_imag", -0.117) :
								  new Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);

  ResonancePdf* rho0_1450  = new ResonancePdf("rho0_1450", 
								  fixAmps ? new Variable("rho0_1450_amp_real", 0.325) : 
								  new Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1450_amp_imag", 0.057) : 
								  new Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);

  ResonancePdf* rhom_1450  = new ResonancePdf("rhom_1450", 
								  fixAmps ? new Variable("rhom_1450_amp_real", 0.788) : 
								  new Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1450_amp_imag", 0.226) : 
								  new Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);

  sharedMass = new Variable("rhop_1700_mass",  1.720, 0.01, 1.6, 1.9);
  shareWidth = new Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0); 

  
  ResonancePdf* rhop_1700  = new ResonancePdf("rhop_1700", 
								  fixAmps ? new Variable("rhop_1700_amp_real", 2.151) : 
								  new Variable("rhop_1700_amp_real",  2.151, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1700_amp_imag", -0.658) : 
								  new Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);
  
  ResonancePdf* rho0_1700  = new ResonancePdf("rho0_1700", 
								  fixAmps ? new Variable("rho0_1700_amp_real",  2.400) : 
								  new Variable("rho0_1700_amp_real",  2.400, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1700_amp_imag", -0.734) : 
								  new Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);
  
  ResonancePdf* rhom_1700  = new ResonancePdf("rhom_1700", 
								  fixAmps ? new Variable("rhom_1700_amp_real",  1.286) : 
								  new Variable("rhom_1700_amp_real",  1.286, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1700_amp_imag", -1.532) : 
								  new Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);
  
  ResonancePdf* f0_980  = new ResonancePdf("f0_980", 
							       fixAmps ? new Variable("f0_980_amp_real",  0.008 * (-_mD02)) : 
							       new Variable("f0_980_amp_real",  0.008 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_980_amp_imag", -0.013 * (-_mD02)) : 
							       new Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_980_mass",     0.980, 0.01, 0.8, 1.2),
							       new Variable("f0_980_width",    0.044, 0.001, 0.001, 0.08),
							       0,
							       PAIR_23);
  
  ResonancePdf* f0_1370  = new ResonancePdf("f0_1370", 
								fixAmps ? new Variable("f0_1370_amp_real", -0.058 * (-_mD02)) : 
								new Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1370_amp_imag",  0.026 * (-_mD02)) : 
								new Variable("f0_1370_amp_imag",  0.026 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1370_mass",     1.434, 0.01, 1.2, 1.6),
								new Variable("f0_1370_width",    0.173, 0.01, 0.01, 0.4),
								0,
								PAIR_23);
  
  ResonancePdf* f0_1500  = new ResonancePdf("f0_1500", 
								fixAmps ? new Variable("f0_1500_amp_real", 0.057 * (-_mD02)) : 
								new Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1500_amp_imag", 0.012 * (-_mD02)) : 
								new Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1500_mass",     1.507, 0.01, 1.3, 1.7),
								new Variable("f0_1500_width",    0.109, 0.01, 0.01, 0.3),
								0,
								PAIR_23);
  
  ResonancePdf* f0_1710  = new ResonancePdf("f0_1710", 
								fixAmps ? new Variable("f0_1710_amp_real", 0.070 * (-_mD02)) : 
								new Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1710_amp_imag", 0.087 * (-_mD02)) : 
								new Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1710_mass",     1.714, 0.01, 1.5, 2.9), 
								new Variable("f0_1710_width",    0.140, 0.01, 0.01, 0.5),
								0,
								PAIR_23);
  
  ResonancePdf* f2_1270  = new ResonancePdf("f2_1270", 
								fixAmps ? new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
								fixAmps ? new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),  
								new Variable("f2_1270_mass",     1.2754, 0.01, 1.0, 1.5),
								new Variable("f2_1270_width",    0.1851, 0.01, 0.01, 0.4),
								2,
								PAIR_23);
  
  ResonancePdf* f0_600  = new ResonancePdf("f0_600", 
							       fixAmps ? new Variable("f0_600_amp_real", 0.068 * (-_mD02)) : 
							       new Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_600_amp_imag", 0.010 * (-_mD02)) : 
							       new Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_600_mass",     0.500, 0.01, 0.3, 0.7),
							       new Variable("f0_600_width",    0.400, 0.01, 0.2, 0.6), 
							       0,
							       PAIR_23);
  
  ResonancePdf* nonr  = new ResonancePdf("nonr",
							     fixAmps ? new Variable("nonr_amp_real", 0.5595 * (-1)) : 
							     new Variable("nonr_amp_real", 0.5595 * (-1),   0.001, 0, 0),
							     fixAmps ? new Variable("nonr_amp_imag", -0.108761 * (-1)) : 
							     new Variable("nonr_amp_imag", -0.108761* (-1), 0.1, 0, 0)); 

  dtop0pp->resonances.push_back(nonr); 
  dtop0pp->resonances.push_back(rhop);
  dtop0pp->resonances.push_back(rho0); 
  dtop0pp->resonances.push_back(rhom); 
  dtop0pp->resonances.push_back(rhop_1450); 
  dtop0pp->resonances.push_back(rho0_1450); 
  dtop0pp->resonances.push_back(rhom_1450); 
  dtop0pp->resonances.push_back(rhop_1700); 
  dtop0pp->resonances.push_back(rho0_1700); 
  dtop0pp->resonances.push_back(rhom_1700); 
  dtop0pp->resonances.push_back(f0_980); 
  dtop0pp->resonances.push_back(f0_1370); 
  dtop0pp->resonances.push_back(f0_1500); 
  dtop0pp->resonances.push_back(f0_1710); 
  dtop0pp->resonances.push_back(f2_1270); 
  dtop0pp->resonances.push_back(f0_600); 

  if (!fitMasses) {
    for (vector<ResonancePdf*>::iterator res = dtop0pp->resonances.begin(); res != dtop0pp->resonances.end(); ++res) {
      (*res)->setParameterConstantness(true); 
    }
  }

  if (!eff) {
    // By default create a constant efficiency. 
    vector<Variable*> offsets;
    vector<Variable*> observables;
    vector<Variable*> coefficients; 

    observables.push_back(m12);
    observables.push_back(m13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne); 
    eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  }

  return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

void runToyFit (std::string toyFileName) {
  m12 = new Variable("m12", 0, 3);
  m13 = new Variable("m13", 0, 3); 
  m12->numbins = 240;
  m13->numbins = 240;
  eventNumber = new CountVariable("eventNumber", 0, INT_MAX);
  getToyData(toyFileName);

  // EXERCISE 1 (real part): Create a PolynomialPdf which models
  // the efficiency you imposed in the preliminary, and use it in constructing
  // the signal PDF. 

  // EXERCISE 2: Create a K0 veto function and use it as the efficiency. 

  // EXERCISE 3: Make the efficiency a product of the two functions
  // from the previous exercises.

  DalitzPlotPdf* signal = makeSignalPdf(); 
  signal->setData(data); 
  signal->setDataSize(data->getNumEvents()); 
  FitManager datapdf(signal); 
  
  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);

  datapdf.setMaxCalls (10);
  datapdf.fit();
  datapdf.getMinuitValues();
  std::vector<Variable*> modParams;
  signal->getParameters(modParams);

  std::vector<double> expected;  //GooFit values
  expected.push_back(-5.59500e-01);
  expected.push_back(1.08761e-01);
  expected.push_back(1);
  expected.push_back(0);
  expected.push_back(5.65000e-01);
  expected.push_back(1.64000e-01);
  expected.push_back(7.14000e-01);
  expected.push_back(-2.50000e-02);;
  expected.push_back(-1.74000e-01);
  expected.push_back(-1.17000e-01);
  expected.push_back(3.25000e-01);
  expected.push_back(5.70000e-02);
  expected.push_back(7.88000e-01);
  expected.push_back(2.26000e-01);
  expected.push_back(2.15100);
  expected.push_back(-6.58000e-01);
  expected.push_back(2.40000);
  expected.push_back(-7.34000e-01);
  expected.push_back(1.28600);
  expected.push_back(-1.53200);
  expected.push_back(-2.78210e-02);
  expected.push_back(4.52092e-02);
  expected.push_back(2.01702e-01);
  expected.push_back(-9.04183e-02);
  expected.push_back(-1.98225e-01);
  expected.push_back(-4.17315e-02);
  expected.push_back(-2.43434e-01);
  expected.push_back(-3.02554e-01);
  expected.push_back(2.95316e-01);
  expected.push_back(4.65835e-02);
  expected.push_back(-2.36479e-01);
  expected.push_back(-3.47763e-02);
  expected.push_back(7.75800e-01);
  expected.push_back(1.50300e-01);
  expected.push_back(1.46500);
  expected.push_back(4.00000e-01);
  expected.push_back(1.72000);
  expected.push_back(2.50000e-01);
  expected.push_back(9.80000e-01);
  expected.push_back(4.40000e-02);
  expected.push_back(1.43400);
  expected.push_back(1.73000e-01);
  expected.push_back(1.50700);
  expected.push_back(1.09000e-01);
  expected.push_back(1.71400);
  expected.push_back(1.40000e-01);
  expected.push_back(1.27540);
  expected.push_back(1.85100e-01);
  expected.push_back(5.00000e-01);
  expected.push_back(4.00000e-01);
  expected.push_back(0);
  expected.push_back(1);  
  
  double variation;
  int count = 0;
  for (int i = 0; i < modParams.size(); i++) { 
    variation = fabs(expected[i] - modParams[i]->value); //expected - actual
    //check the variance of the generated parameter from its actual value and compare it to our epsilon of 0.001
    if (variation > 0.001) {
      std::cout << "\n" << modParams[i]->name << " value not in epsilon." << endl;
      std::cout << "Expected: " << expected[i] << endl;
      std::cout << "Actual: " << modParams[i]->value << endl;
      std::cout << "Variation: " << variation << endl;
      count++;
    }
  }
  
  std::cout << "\nTotal variances: " << count << endl;  
  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
}

int main (int argc, char** argv) {
#ifdef TARGET_MPI
  MPI_Init(&argc, &argv);

  //we have MPI, so lets do something slightly different here:
  int myId, numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

#ifndef TARGET_OMP
  //set the processes to gpus here
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  
  //No way to figure out how many processes per node, so we read the environment variable
  int nodes = atoi (getenv ("PBS_NUM_NODES"));
  if (nodes == 0)
    nodes = 1;
  int procsPerNode = numProcs/nodes;
  int localRank = myId % procsPerNode;

  if (deviceCount == 1 && localRank > 1)
  {
    printf ("Multi-process to one GPU!\n");
    cudaSetDevice (0);
  }
  else if (procsPerNode > 1 && deviceCount > 1)
  {
     if (localRank <= deviceCount)
     {
       printf ("setting multiple processes to multiple GPU's\n");
       cudaSetDevice (localRank);
     }
     else
     {
       printf ("More multi-processes than multi-gpu's!\n");
       cudaSetDevice (localRank % deviceCount);
     }
  }
  else
  {
    printf ("Multi-GPU's, using one process! %i, [%i,%i]\n", deviceCount, localRank, procsPerNode);
    cudaSetDevice (0);
  }
#endif
#endif

 //check to see that the file exists
 ifstream ifile(argv[1]);
 if (ifile) {

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1, 0);
  foo = new TCanvas(); 
  foodal = new TCanvas(); 
  foodal->Size(10, 10);

  runToyFit(argv[1]);

  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

#ifdef TARGET_MPI
  MPI_Finalize();
#endif

  delete m12;
  delete m13;
  delete eventNumber;
  delete foo;
  delete foodal;
  delete constantOne;
  delete constantZero;
  delete fixedRhoMass;
  delete fixedRhoWidth;
  delete motherM;
  delete chargeM;
  delete neutrlM;
  delete massSum;
  delete data;
 
 } else {
    std::cout << "**ERROR: File not found**" << endl;
 }
  return 0; 
}
