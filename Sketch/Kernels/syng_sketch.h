//////////////////////////////////////////////////////////////////////////
// Stage 4. Calculating the synaptic conductance
//
//------------------------------------------------------------------------
//

#ifndef __SYNG_SKETCH_H
#define __SYNG_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global variables
//------------------------------------------------------------------------
// constant memory 64 Kb
// for 1.x devices shared memory is 16 Kb per multiprocessor (warp of 32 threads),
// all devices >= 2.0 have 48 KB shared memory
//--- global parameters (READ-ONLY access, could be stored in the constant memory)


//////////////////////////////////////////////////////////////////////////
// <kernel 4.1> calculating the synaptic conductance
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and
// could be cached in the shared memory)
extern uint Nsyn;       // total number of synapses
extern uint4 *SynLUT;   // L=Nsyn: look-up-tables to calculate a synaptic conductance
                        // x - look-up-table of components for synaptic transmission
                        // y - look-up-table of components for presynaptic inhibition
                        // z - look-up-table of components for synaptic plasticity
                        // w - look-up-table of components for synaptic conductance
//--- shared variables (READ) (stored in global memory)
float *SynGates;        // L<=3*Nsyn: synaptic gate variables
                        // (synaptic transmission; presynaptic inhibition; synaptic plasticity)
extern float *Gmax;     // L=Nsyn: array of maximal conductance of synapses
//--- shared variables (WRITE) (stored in global memory)
extern float *Gsyn;     // L=Nsyn: array of synaptic conductance

//--- kernel 4.1 calculating the synaptic conductance
// The result is stored in array Gsyn (synaptic conductances).
// The variables in arrays SynGates & Gmax must be specified before.
void kernel_syn_g( uint id )
{
    // load to shared memory (?)
    uint mi = SynLUT[id].x;       // id-th component in look-up-table for synaptic transmission
    uint hi = SynLUT[id].y;       // id-th component in look-up-table for presynaptic inhibition
    uint zi = SynLUT[id].z;       // id-th component in look-up-table for synaptic plasticity
    uint si = SynLUT[id].w;       // synaptic conductance of id-th synapse
    // load from global memory
    float m = SynGates[mi];       // synaptic transmission
    float h = SynGates[hi];       // presynaptic inhibition
    float z = SynGates[zi];       // synaptic plasticity
    float g = Gmax[si];           // maximal conductance
    // math & store to global memory
    Gsyn[si] = g*m*z/(1.+h);      // synaptic conductance
    // ge = ???
}

//--- call <kernel 4.1> for all synapses
void syn_g( void )
{
    for( uint i = 0; i < Nsyn; ++i ){
        kernel_syn_g( i );
    }
}

#endif // __SYNG_SKETCH_H
