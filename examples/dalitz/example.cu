#include "Variable.hh" 
#include "../../PDFs/GaussianPdf.hh" 
#include "../../FitManager.hh" 
#include "../../UnbinnedDataSet.hh" 
#include "../../PDFs/GooPdf.hh"
#include "../../PdfBase.hh"

#include "TRandom.hh" 
#include "TH1F.h"
#include "TCanvas.h" 

#include <sys/time.h>
#include <sys/times.h>
#include <stdio.h>	

using namespace std; 

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


  Variable* xvar = new Variable("xvar", -5, 5); 

  int numbins = 10000;

  // Generate data
  TRandom donram(42); 
  UnbinnedDataSet data(xvar); // Stores events
  for (int i = 0; i < numbins; ++i) {
    fptype val = donram.Gaus(0.2, 1.1); // these are the values for mean and sigma
    if (fabs(val) > 5) {
	--i; 
	continue;
    } 
    data.addEvent(val); 
  }

  // Create the PDF
  Variable* mean = new Variable("mean", 0, 1, -10, 10);
  Variable* sigma = new Variable("sigma", 1, 0.5, 1.5); 
  GaussianPdf* gauss = new GaussianPdf("gauss", xvar, mean, sigma); 

  timeval startTime, stopTime, totalTime;

  // Run a fit
  gauss->setData(&data);
  FitManager datapdf(gauss); 
  gettimeofday(&startTime, NULL);
  datapdf.fit(); 
  
  // Redirect the output to 'output.txt'
  freopen("output.txt","w",stdout);

  fptype* host_output = new fptype[xvar->numbins];
  gauss->transformGrid(host_output);
  std::cout << "TRANSFORMGRID: \n" << std::endl;
  for (int i = 0; i < xvar->numbins; i++) {
    std::cout << host_output[i] << std::endl;
  }
    
  vector<fptype> res;
  gauss->evaluateAtPoints(xvar, res);
  std::cout << "\nEVALUATEATPOINTS: \n" << std::endl;
  for (int i = 0; i < res.size(); i++) {
      std::cout <<  res[i] << std::endl;
  }
     
  vector<vector<fptype> > values;
  gauss->getCompProbsAtDataPoints(values);
  std::cout << "\nGETCOMPPROBSATDATAPOINTS: \n" << std::endl;
  for (int i = 0; i < values.size(); i++) {
    for (int j = 0; j < values.size(); j++) {
      std::cout << values[i][j] << std::endl;
    }
  }
  
  gettimeofday(&stopTime, NULL);

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "\nWallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;

  return 0;
}
