#include "PdfBase.hh"
#include <mpi.h>
#include <typeinfo>


// This is code that belongs to the PdfBase class, that is, 
// it is common across all implementations. But it calls on device-side
// functions, and due to the nvcc translation-unit limitations, it cannot
// sit in its own object file; it must go in the CUDAglob.cu. So it's
// off on its own in this inline-cuda file, which GooPdf.cu 
// should include. 

#define CHECK_MPI(error) \
        if (error != MPI_SUCCESS) { \
					int length; \
					char message[MPI_MAX_ERROR_STRING]; \
          MPI_Error_string(error, message, &length); \
          printf("\n%.*s\n", length, message); \
          MPI_Abort(MPI_COMM_WORLD, 1);}

#ifdef CUDAPRINT
__host__ void PdfBase::copyParams (const std::vector<double>& pars) const {
  if (host_callnumber < 1) {
    std::cout << "Copying parameters: " << (long long) cudaArray << " ";
  }
  for (unsigned int i = 0; i < pars.size(); ++i) {
    host_params[i] = pars[i]; 
    
    if (host_callnumber < 1) {
      std::cout << pars[i] << " ";
    }
    
    if (isnan(host_params[i])) {
      std::cout << " agh, NaN, die " << i << std::endl;
      abortWithCudaPrintFlush(__FILE__, __LINE__, "NaN in parameter"); 
    }
  }
  
  if (host_callnumber < 1) {
    std::cout << std::endl; 
  }
  MEMCPY_TO_SYMBOL(cudaArray, host_params, pars.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
}
#else 
__host__ void PdfBase::copyParams (const std::vector<double>& pars) const {
  // copyParams method performs eponymous action! 

  for (unsigned int i = 0; i < pars.size(); ++i) {
    host_params[i] = pars[i]; 
    
    if (isnan(host_params[i])) {
      std::cout << " agh, parameter is NaN, die " << i << std::endl;
      abortWithCudaPrintFlush(__FILE__, __LINE__, "NaN in parameter"); 
    }
  }

  MEMCPY_TO_SYMBOL(cudaArray, host_params, pars.size()*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
}
#endif

__host__ void PdfBase::copyParams () {
  // Copies values of Variable objects
  parCont pars; 
  getParameters(pars); 
  std::vector<double> values; 
  for (parIter v = pars.begin(); v != pars.end(); ++v) {
    int index = (*v)->getIndex(); 
    if (index >= (int) values.size()) values.resize(index + 1);
    values[index] = (*v)->value;
  }
  copyParams(values); 
}

__host__ void PdfBase::copyNormFactors () const {
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  SYNCH(); // Ensure normalisation integrals are finished
}

__host__ void PdfBase::initialiseIndices (std::vector<unsigned int> pindices) {
  // Structure of the individual index array: Number of parameters, then the indices
  // requested by the subclass (which will be interpreted by the subclass kernel), 
  // then the number of observables, then the observable indices. Notice that the
  // observable indices are not set until 'setIndices' is called, usually from setData;
  // here we only reserve space for them by setting totalParams. 
  // This is to allow index sharing between PDFs - all the PDFs must be constructed 
  // before we know what observables exist. 

  if (totalParams + pindices.size() >= maxParams) {
    std::cout << "Major problem with pindices size: " << totalParams << " + " << pindices.size() << " >= " << maxParams << std::endl; 
  }

  assert(totalParams + pindices.size() < maxParams); 
  host_indices[totalParams] = pindices.size(); 
  for (int i = 1; i <= host_indices[totalParams]; ++i) {
    host_indices[totalParams+i] = pindices[i-1]; 
  }
  host_indices[totalParams + pindices.size() + 1] = observables.size(); 
  
  parameters = totalParams;
  totalParams += (2 + pindices.size() + observables.size()); 

  /* 
  std::cout << "host_indices after " << getName() << " initialisation : ";
  for (int i = 0; i < totalParams; ++i) {
    std::cout << host_indices[i] << " ";
  }
  
  std::cout << " | " 
	    << parameters << " " 
	    << totalParams << " " 
	    << cudaArray << " " 
	    << paramIndices << " "
	    << std::endl; 
  */
  MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice); 
}

/*
__host__ void PdfBase::setData (std::vector<std::map<Variable*, fptype> >& data) {
  // Old method retained for backwards compatibility 

  if (dev_event_array) {
    gooFree(dev_event_array);
    dev_event_array = 0; 
  }

  setIndices();
  int dimensions = observables.size();
  numEntries = data.size();
  numEvents = numEntries; 
  
  fptype* host_array = new fptype[data.size()*dimensions];
  for (unsigned int i = 0; i < data.size(); ++i) {
    for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
      assert(data[i].find(*v) != data[i].end()); 
      host_array[i*dimensions + (*v)->index] = data[i][*v]; 
    }
  }

  gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype)); 
  MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array; 
}
*/

__host__ void PdfBase::recursiveSetIndices () {
  for (unsigned int i = 0; i < components.size(); ++i) {
    components[i]->recursiveSetIndices(); 
  }

  int numParams = host_indices[parameters]; 
  int counter = 0; 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    host_indices[parameters + 2 + numParams + counter] = (*v)->index; 
    //std::cout << getName() << " set index of " << (*v)->name << " to " << (*v)->index << " " << (parameters + 2 + numParams + counter) << std::endl; 
    counter++; 
  }  
  generateNormRange(); 
}

__host__ void PdfBase::setIndices () {
  int counter = 0; 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    (*v)->index = counter++; 
  }
  recursiveSetIndices(); 
  MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams*sizeof(unsigned int), 0, cudaMemcpyHostToDevice); 

  //std::cout << "host_indices after " << getName() << " observable setIndices : ";
  //for (int i = 0; i < totalParams; ++i) {
  //std::cout << host_indices[i] << " ";
  //}
  //std::cout << std::endl; 
}

__host__ void PdfBase::setData (UnbinnedDataSet* data)
{
  //std::cout << "PdfBase::setData" << std::endl;
  if (dev_event_array) {
    gooFree(dev_event_array);
    SYNCH();
    dev_event_array = 0; 

    m_iEventsPerTask = 0;
  }

  setIndices();
  int dimensions = observables.size();
  numEntries = data->getNumEvents(); 
  numEvents = numEntries; 

#ifdef TARGET_MPI

  int world_size, world_rank;
	// Get the number of processes
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
	// Get the rank of the processes
  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank); 

  int num_events_per_proc = numEntries/world_size;

  if (world_rank == 0) {
    numEvents = data->getNumEvents();
    dimensions = observables.size();
  }

  MPI_Bcast(&numEvents, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dimensions, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);

  int *counts = new int[world_size];
  int *displacements = new int[world_size];

  //indexing for copying events over!
  for (int i = 0; i < world_size - 1; i++)
    counts[i] = num_events_per_proc;
  counts[world_size - 1] = numEntries - num_events_per_proc*(world_size - 1);
  
  //displacements into the array for indexing!
  displacements[0] = 0;
  for (int i = 1; i < world_size; i++)
    displacements[i] = displacements[i - 1] + counts[i - 1];

#endif

  fptype* host_array;

#ifdef TARGET_MPI

  if (world_rank == 0) {
    host_array = new fptype[numEntries*dimensions];
  }

  //This is an array to track if we need to redo indexing
  int fixme[observables.size ()];
  memset(fixme, 0, sizeof (int)*observables.size ());

  //printf ("Checking observables for Counts!\n");
  for (int i = 0; i < observables.size (); i++)
  {
    //printf ("%i - %s\n", i, observables[i]->name.c_str ());
    //cast this variable to see if its one we need to correct for
    CountVariable *c = dynamic_cast <CountVariable*> (observables[i]);
    //if it cast, mark it
    if (c)
    {
      fixme[i] = 1;
      //printf ("%i of %i - %s\n", i, observables.size (), c->name.c_str ());
    }
  }

	if (world_rank == 0) {
    //populate this array with our stuff
    for (int i = 0; i < numEntries; ++i)
    {
      for (obsIter v = obsBegin(); v != obsEnd(); ++v)
      {
        fptype currVal = data->getValue((*v), i);
        host_array[i*dimensions + (*v)->index] = currVal;
      }
    }

    printf("\nHost array: ");
    for (int i = 0; i < world_size; i++) {
      for (int j = 0; j < 6; j++) {	
    	  printf("%f ", host_array[displacements[i]+j]);
  	  }
    }
  }
	
#else

  host_array = new fptype[numEntries*dimensions];

  for (int i = 0; i < numEntries; ++i) 
  {
    for (obsIter v = obsBegin(); v != obsEnd(); ++v)
    {
      fptype currVal = data->getValue((*v), i);
      host_array[i*dimensions + (*v)->index] = currVal;
    }
  }

	printf("\nHost array: ");
	for (int i = 0; i < 6; i++) {
		printf("%f ", host_array[i]);
	}

#endif

#ifdef TARGET_MPI

  int mystart = displacements[world_rank];
  int myend = mystart + counts[world_rank];
  int mycount = myend - mystart;
	
	fptype *recv_buf = new fptype[mycount]; // The buffer to receive the scattered elements
	
  CHECK_MPI(MPI_Scatterv(host_array, counts, displacements, MPI_DOUBLE, recv_buf, mycount, MPI_DOUBLE, 0, MPI_COMM_WORLD));

	// Print the numbers scattered to each processor
  printf("\nProcessor rank %i out of %i processors: ", world_rank, world_size);

  for (int i = 0; i < 6; i++) {
    printf("%f ", recv_buf[i]);
  }

	// we need to fix our observables indexing to reflect having multiple cards
  for (int i = 1; i < world_size; i++)
  {
    for (int j = 0; j < counts[i]; j++)
    {
      //assumption is that the last observable is the index!
      for (int k = 0; k < dimensions; k++)
      {
        //Its counting, fix the indexing here
        if (fixme[k] > 0)
          recv_buf[j * dimensions + k] = float (j);
      }
    }
  }

#endif 

#ifdef TARGET_MPI 
	gooMalloc((void**) &dev_event_array, dimensions*mycount*sizeof(fptype));
  MEMCPY(dev_event_array, recv_buf, dimensions*mycount*sizeof(fptype), cudaMemcpyHostToDevice);
  MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
  if (world_rank == (world_size-1))  
    delete[] host_array;
  
  delete[] recv_buf; 
  printf("\ndev_event_array: ");

  for (int i = 0; i < 6; i++) {
    printf("%f ", dev_event_array[i]);
  }

  if (world_rank == (world_size - 1)) 
    printf("\n\n");

  //update everybody
  setNumPerTask(this, mycount);

  delete [] counts;
  delete [] displacements;
#else
  gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype)); 
  MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
  MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array;
#endif
}

__host__ void PdfBase::setData (BinnedDataSet* data)
{ 
  if (dev_event_array) { 
    gooFree(dev_event_array);
    dev_event_array = 0; 

    m_iEventsPerTask = 0;
  }

  setIndices();
  numEvents = 0; 
  numEntries = data->getNumBins(); 
  int dimensions = 2 + observables.size(); // Bin center (x,y, ...), bin value, and bin volume. 
  if (!fitControl->binnedFit()) setFitControl(new BinnedNllFit()); 

  fptype* host_array = new fptype[numEntries*dimensions]; 

#ifdef TARGET_MPI

	int world_size, world_rank;
  // Get the number of processes
  MPI_Comm_size (MPI_COMM_WORLD, &world_size);
  // Get the rank of the processes
  MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);

  //This is an array to track if we need to redo indexing
  int fixme[dimensions];
  memset(fixme, 0, sizeof (int)*dimensions);

  for (int i = 0; i < observables.size (); i++)
  {
    //cast this variable to see if its one we need to correct for
    CountVariable *c = dynamic_cast <CountVariable*> (observables[i]);
    //if it cast, mark it
    if (c)
      fixme[i] = 1;
  }

	if (world_rank == 0) { 
    // populate the array
    for (unsigned int i = 0; i < numEntries; ++i) {
      for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
        host_array[i*dimensions + (*v)->index] = data->getBinCenter((*v), i);
      }

      host_array[i*dimensions + observables.size() + 0] = data->getBinContent(i);
      host_array[i*dimensions + observables.size() + 1] = fitControl->binErrors() ? data->getBinError(i) : data->getBinVolume(i);
      numEvents += data->getBinContent(i);
    }

		printf("\nHost array: ");
  	for (int i = 0; i < 25; i++) {
    	printf("%f ", host_array[i]);
  	}
  }

#else

	// populate the array
  for (unsigned int i = 0; i < numEntries; ++i) {
  	  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    	  host_array[i*dimensions + (*v)->index] = data->getBinCenter((*v), i); 
    	}
		
			host_array[i*dimensions + observables.size() + 0] = data->getBinContent(i);
			host_array[i*dimensions + observables.size() + 1] = fitControl->binErrors() ? data->getBinError(i) : data->getBinVolume(i); 
			numEvents += data->getBinContent(i);	
	}

	printf("\nHost array: ");
  for (int i = 0; i < 25; i++) {
    printf("%f ", host_array[i]);
  }

#endif

#if TARGET_MPI

  int num_events_per_proc = numEvents/world_size;

  int *counts = new int[world_size];
  int *displacements = new int[world_size];

  //indexing for copying events over!
  for (int i = 0; i < world_size - 1; i++)
    counts[i] = num_events_per_proc;
  counts[world_size - 1] = numEvents - num_events_per_proc*(world_size - 1);

  //displacements into the array for indexing!
  displacements[0] = 0;
  for (int i = 1; i < world_size; i++)
    displacements[i] = displacements[i - 1] + counts[i - 1];

  int mystart = displacements[world_rank];
  int myend = mystart + counts[world_rank];
  int mycount = myend - mystart;
	
	fptype *recv_buf = new fptype[mycount];

	// check data type of fptype
	if (typeid(fptype) == typeid(float)) {
		printf("\nfloat");
	  CHECK_MPI(MPI_Scatterv(host_array, counts, displacements, MPI_FLOAT, recv_buf, mycount, MPI_FLOAT, 0, MPI_COMM_WORLD));
	} else {
		// fptype is of type double
		printf("\ndouble");
		CHECK_MPI(MPI_Scatterv(host_array, counts, displacements, MPI_DOUBLE, recv_buf, mycount, MPI_DOUBLE, 0, MPI_COMM_WORLD));
	}

  // Print the numbers scattered to each processor
  printf("\nProcessor rank %d"
           " out of %d processors: ", world_rank, world_size);

	// we need to fix our observables indexing to reflect having multiple cards
  for (int i = 1; i < world_size; i++)
  {
    for (int j = 0; j < counts[i]; j++)
    {
      //assumption is that the last observable is the index!
      for (int k = 0; k < dimensions; k++)
      {
        //Its counting, fix the indexing here
        if (fixme[k] > 0)
          recv_buf[(j + displacements[i])*dimensions + dimensions - k] = float (j);
      }
    }
  }

	for (int i = 0; i < 25; i++) {
    printf("%f ", recv_buf[i]);
  }

#endif

#ifdef TARGET_MPI
  gooMalloc((void**) &dev_event_array, dimensions*mycount*sizeof(fptype)); 
  MEMCPY(dev_event_array, recv_buf, dimensions*mycount*sizeof(fptype), cudaMemcpyHostToDevice); 
  MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  delete[] host_array;
	delete[] recv_buf;

  //update our displacements:
  for (int i = 0; i < world_size; i++)
    displacements[i] = 0;

  //update everybody
  setNumPerTask(this, mycount);

  delete [] counts;
  delete [] displacements;
#else
  gooMalloc((void**) &dev_event_array, dimensions*numEntries*sizeof(fptype)); 
  MEMCPY(dev_event_array, host_array, dimensions*numEntries*sizeof(fptype), cudaMemcpyHostToDevice);
#endif
}

__host__ void PdfBase::generateNormRange () {
  if (normRanges) gooFree(normRanges);
  gooMalloc((void**) &normRanges, 3*observables.size()*sizeof(fptype));
  
  fptype* host_norms = new fptype[3*observables.size()];
  int counter = 0; // Don't use index in this case to allow for, eg, 
  // a single observable whose index is 1; or two observables with indices
  // 0 and 2. Make one array per functor, as opposed to variable, to make
  // it easy to pass MetricTaker a range without worrying about which parts
  // to use. 
  for (obsIter v = obsBegin(); v != obsEnd(); ++v) {
    host_norms[3*counter+0] = (*v)->lowerlimit;
    host_norms[3*counter+1] = (*v)->upperlimit;
    host_norms[3*counter+2] = integrationBins > 0 ? integrationBins : (*v)->numbins;
    counter++; 
  }

  MEMCPY(normRanges, host_norms, 3*observables.size()*sizeof(fptype), cudaMemcpyHostToDevice);
  delete[] host_norms; 
}

void PdfBase::clearCurrentFit () {
  totalParams = 0; 
  gooFree(dev_event_array);
  dev_event_array = 0; 
}

__host__ void PdfBase::printProfileInfo (bool topLevel) {
#ifdef PROFILING
  if (topLevel) {
    cudaError_t err = MEMCPY_FROM_SYMBOL(host_timeHist, timeHistogram, 10000*sizeof(fptype), 0);
    if (cudaSuccess != err) {
      std::cout << "Error on copying timeHistogram: " << cudaGetErrorString(err) << std::endl;
      return;
    }
    
    std::cout << getName() << " : " << getFunctionIndex() << " " << host_timeHist[100*getFunctionIndex() + getParameterIndex()] << std::endl; 
    for (unsigned int i = 0; i < components.size(); ++i) {
      components[i]->printProfileInfo(false); 
    }
  }
#endif
}



gooError gooMalloc (void** target, size_t bytes) {
// Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_BACKEND==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP
  target[0] = malloc(bytes);
  if (target[0]) return gooSuccess;
  else return gooErrorMemoryAllocation; 
#else
  return (gooError) cudaMalloc(target, bytes); 
#endif
}

gooError gooFree (void* ptr) {
// Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_BACKEND==THRUST_DEVICE_BACKEND_OMP || THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP
  free(ptr);
  return gooSuccess;
#else
  return (gooError) cudaFree(ptr); 
#endif
}
