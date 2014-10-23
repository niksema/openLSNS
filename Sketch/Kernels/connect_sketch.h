//////////////////////////////////////////////////////////////////////////
//    Stage 2. Connectivity.
//    The connectivity matrix between elements within the network is specified
//    as the look-up-tables for weighted connections between targets and
//    sources. The sequences of sources and weights which specified in arrays
//    Srcs & Weights should be the most similar and continuous for sequential
//    targets. The metric of connections specified as look-up-table in Metric
//    and Connect arrays, where:
//        1. Metric is an array that specifies numbers of sources and offset for
//           initial element in Connect array for each target.
//        2. Connect is a look-up-table which specifies the the set of pairs
//           of weights and sourses in correspondent arrays for each target
//------------------------------------------------------------------------
//    The optimization for connectivity matrix must be done in the first
//    place (highest priority). Another possible optimization could be the
//    calculating of the partial sum for each target on single multiprocessor.

#ifndef __CONNECT_SKETCH_H
#define __CONNECT_SKETCH_H
//--- global parameters (READ-ONLY access, could be stored in the constant memory)
uint Ntrg;      // total number of targets (each one receives >=1 connections)
uint Nsrc;      // total number of sources (each one provides >= 1 connections);
uint Nc;        // total number of connections; Nsrc <= Nc
uint Nw;        // total number of weights; Nw <= Nc
//--- local parameters (READ-ONLY access, stored in global memory and could be
//    cached in the shared memory)
uint4 *Metric;  // L=Ntrg; specify the connectivity matrix for the network:
                // x - target index
                // y - number of sources for the target
                // z - index of initial element in the look-up-table for connections
                // w - reserved
uint2 *Connect; // L=Nc: look-up-table for all connections between targets and sources:
                // x - an index in sources array,
                // y - an index in weight array.
                // duplicated series of connections may/or may not be excluded
                // depending on optimization strategy
float *Weights; // L=Nw: list of weights/connectios (none zero)
                // duplicated weights may/or may not be excluded
                // depending on optimization strategy
//--- shared variables (READ/WRITE access, stored in global memory)
float *Srcs;    // L=Nsrc:  list of sources; the array 'PreSyn' (specified at stage I)
                // is a subset of array Srcs (Nsrc >= Nsout)
float *Trgs;    // L=Ntrg:  a list of targets (synapses)
                // results of synaptic summation
                // uses as arrays M and H for stage III

//--- <kernel 2.1>. synaptic summation
// The result is a weighted sum between correspondent sources and weights stored in
// array of targets.
void convolution( void )
{
    for( uint i = 0; i < Ntrg; ++i ){ // loop for all targets
        // load to shared memory (?)
        uint ti = Metric[i].x;    // make sure that ti < Ntrg
        uint n_src = Metric[i].y; // make sure that n_src+si < Nc
        uint si = Metric[i].z;    // make sure that n_src+si < Nc
        // math: to calculate the weighted sum for ti-th targer
        float sum = 0.;             // initialize the weighted sum for ti-th targer
        for( uint j = si; j < n_src+si; ++j ){// loop to calculate weighted sum for ti-th target
            // load to shared memory (?)
            uint xi = Connect[j].x;           // index of specific source; make sure that xi < Nsrc
            uint wi = Connect[j].y;           // index of specific weight; make sure that wi < Nw
            sum += Srcs[xi]*Weights[wi];      // calculate the partial sum for i-th target
        }
        // store weighted sum to global memory
        Trgs[ti] = sum;
    }
}

#endif // __CONNECT_SKETCH_H
