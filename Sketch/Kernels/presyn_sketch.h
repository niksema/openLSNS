//////////////////////////////////////////////////////////////////////////
// Stage 1. Common presynaptic operations (spikes discrimination
//          and transformation of membrane potentials of neurons) for signal
//          transmission between neurons and/or nerves
//------------------------------------------------------------------------
// CUDA specific: the size of constant memory 64 Kb
// for 1.x devices shared memory is 16 Kb per multiprocessor (warp of 32
// threads), all devices >= 2.0 have 48 KB shared memory
//------------------------------------------------------------------------

#ifndef __PRESYN_SKETCH_H
#define __PRESYN_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global parameters and variables
//--- global parameters (READ-ONLY access, could be stored in the constant
//    memory)
extern float Theshold;  // the threshold for spike discrimination
extern float Step;      // the integration step
extern uint Nv;         // total number of all compartments
extern uint Nout;       // total number of all synapses (either type) in the network
extern uint NStType;    // number of different types of "standars" synapses (AMPA/GABA)
extern float2 *StType;  // L=NStdSyn; parameters of "standard" synapses:
                        // x - exp(-Step/T) where T is a time constant
                        // y - amplitude jump
extern uint NSiType;    // number of different types of sigma synapses
extern float2 *SiType;  // L=NSiType; parameters of sigma synapses:
                        // x - half-voltage
                        // y - slope
extern uint NAType;     // number of different types of alpha synapses (AMPA/GABA)
extern float2 *AType1;  // L=NAType; constant parameters of alpha synapses:
                        // x - a0
                        // y - a
extern float4 *AType2;  // L=NAType; parameters of alpha synapses:
                        // x - alpha
                        // y - half-voltage
                        // z - slope
                        // w - time constant
extern uint NNmType;    // number of different types of nmda synapses
extern float2 *NmdaType1;// L=NNmdSyn; parameters of nmda synapses:
                        // x - alpha
                        // y - T
extern float4 *NmdaType2;// L=NNmdSyn; parameters of nmda synapses: (see stage 3)
                        // x - k
                        // y - half-voltage
                        // z - slope
                        // w - reserved
//--- shared variables (READ/WRITE access, stored in global memory)
extern uint *Spike;     // L=Nv; 1 if spike is generated, otherwise 0
extern uint *Spike_;    // L=Nv; 1 if V >= Threshold, otherwise 0
extern float *V;        // L=Nv; membrane potential of all compartments
extern float *PreSyn;   // L=Nout; presynaptic outputs

//////////////////////////////////////////////////////////////////////////
// <kernel 1.1> spike discrimination
//------------------------------------------------------------------------
// Optimization strategy: optimize the look-up-table in SpikeLUT array to
// provide coalescing memory accesses for V, Spike, Spike_ arrays.
//------------------------------------------------------------------------
// kernel specific parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could
//    be cached in the shared memory)
uint Nspikes;   // total number of neurons (axon compartments) involved
                // into spikes discrimination
uint *AxonLUT;  // L=Nspikes; look-up-table on axon/soma compartments
                // involved into spikes discrimination

//--- <kernel 1.1> if a spike generated then 1 otherwise 0.
// Spike = ( V[current] <= Threshold && V[previous] > Threshold )? 1: 0;
// the result is stored in the array 'Spike'
void kernel_spike_presyn( uint id )
{
    // load to shared memory (?)
    uint axon = AxonLUT[id];    // make sure axon < Nv
    // load from global memory
    uint lt_tr = uint( V[axon] < Theshold );
    // math & store to global memory
    Spike[axon] = lt_tr*Spike_[axon];
    Spike_[axon] = !lt_tr;
}

//--- call <kernel 1.1> for all compartments which involved in spike
// discrimination.
//------------------------------------------------------------------------
// Experiments: check if spike discrimination for all compartments faster
// then for soma(axon) compartments only.
//------------------------------------------------------------------------
void spike_presyn( void )
{
    for( uint i = 0; i < Nspikes; ++i ){
        kernel_spike_presyn( i );
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.2> presynaptic transformation of membrane potentials for
// sigma transformation of membrane potential
//------------------------------------------------------------------------
// Optimization strategy: optimize the look-up-table in SiSynLUT array to
// provide coalescing memory accesses for V & PreSyn arrays.
//------------------------------------------------------------------------
// kernel specific parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could
//    be cached in the shared memory)
uint Nsisyn;    // total number of sigma synapses
uint4 *SiSynLUT;// L=Nsisyn; look-up-tables for sigma synapses:
                // x - look-up-table of types of sigma synapses
                // y - look-up-table of axons provided sigma synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.2> presynaptic output is calculating according to <publication(s)>
// Out = 1./( 1.+exp(-( v-Hv )/Slp ))
// the result is stored in array PreSyn.
void kernel_sigma_presyn( uint id )
{
    // load to shared memory (?)
    uint type = SiSynLUT[id].x; // make sure type < NSiType
    uint axon = SiSynLUT[id].y; // make sure axon < Nv
    uint out = SiSynLUT[id].z;  // make sure out < Nout
    // load from constant memory (?)
    float hv = SiType[type].x;
    float slp = SiType[type].y;
    // load from global memory
    float v = V[axon];
    // math & store to global memory
    PreSyn[out] = 1./( 1.+exp(-( v-hv )/slp ));
}

//--- call <kernel 1.3> for all compartments which involved in presynaptic
// processing for sigma synapses
void sigma_presyn( void )
{
    for( uint i = 0; i < Nsisyn; ++i ){
        kernel_sigma_presyn( i );
    }
}


/////////////////////////////////////////--------------------------------
// not sure. TODO synaptic summation needs to be clarify


//////////////////////////////////////////////////////////////////////////
// <kernel 1.3> presynaptic transformation of membrane potentials for
// standard synapses
//------------------------------------------------------------------------
// Optimization strategy: optimize the look-up-table in SSynLUT array to
// provide coalescing memory accesses for Spike & PreSyn arrays.
//------------------------------------------------------------------------
// kernel specific parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could
//    be cached in the shared memory)
uint Nssyn;     // total number of "standard" synapses
uint4 *SSynLUT; // L=Nssyn; look-up-tables for "standard" synapses:
                // x - look-up-table of types of "standard" synapses
                // y - look-up-table of axons provided synaptic inputs to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.3> presynaptic output is calculating according to <publication(s)>
// Out = Out*exp(-step/T )+Spike*Jump;
// the result is stored in array PreSyn.
void kernel_std_presyn( uint id )
{
    // load to shared memory (?)
    uint type = SSynLUT[id].x;  // make sure type < NStType
    uint axon = SSynLUT[id].y;  // make sure axon < Nv
    uint out = SSynLUT[id].z;   // make sure out < Nout
    // load from constant memory (?)
    float expt = StType[type].x;
    float a = StType[type].y;
    // load from global memory
    float presyn = PreSyn[out];
    uint spike = Spike[axon];
    // math & store to global memory
    PreSyn[out] = presyn*expt+a*spike;
}

//--- call <kernel 1.2> for all compartments which involved in presynaptic
// processing for "standard" synapses
void std_presyn( void )
{
    for( uint i = 0; i < Nssyn; ++i ){
        kernel_std_presyn( i );
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.4> presynaptic transformation of membrane potentials for
// alpha synapses (AMPA/GABA)
//------------------------------------------------------------------------
// Optimization strategy: optimize the look-up-table in SSynLUT array to
// provide coalescing memory accesses for V & PreSyn arrays.
//------------------------------------------------------------------------
// kernel specific parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nasyn;     // total number of alpha synapses
uint4 *ASynLUT; // L=Nasyn; look-up-tables for alpha synapses:
                // x - look-up-table of types of alpha synapses
                // y - look-up-table of axons provided alpha synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.4> presynaptic output is calculating according to
// Wang X-J, Tegner J, Constantinidis C, Goldman-Rakic PS (2004)
// Division of labor among distinct subtypes of inhibitory neurons in a cortical microcircuit of working memory.
// Proc. Natl. Acad. Sci. (USA), 101, 1368-1373
// (see, http://www.cns.nyu.edu/wanglab/publications/pdf/wang04.pnas.pdf)
// T*dOut/dt = out0-Out*(out0+1), where out0 = alpha*T*f, f = a0+a/(1+exp(-( v-Hv )/Slp ))
// the result is stored in array PreSyn.
void kernel_alpha_presyn( uint id )
{
    // load to shared memory (?)
    uint type = ASynLUT[id].x;  // make sure type < NAType
    uint axon = ASynLUT[id].y;  // make sure axon < Nv
    uint out = ASynLUT[id].z;   // make sure out < Nout
    // load from constant memory (?)
    float a0 = AType1[type].x;
    float a = AType1[type].y;
    float alpha = AType2[type].x;
    float hv = AType2[type].y;
    float slp = AType2[type].z;
    float t = AType2[type].w;

    // load from global memory
    float presyn = PreSyn[out];
    float v = V[axon];
    // math
    float k = alpha*( a0+a/( 1.+exp(-( v-hv )/slp )))*t;
    float out0 = k/( k+1. );
    t = t/( k+1. );
    presyn = EXP_EULER( presyn, out0, step, t );
    // store to global memory
    PreSyn[out] = presyn;
}

void alpha_presyn( void )
{
    for( uint i = 0; i < Nasyn; ++i ){
        kernel_alpha_presyn( i );
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.5> presynaptic transformation of membrane potentials for
// nmda synapses
//------------------------------------------------------------------------
// kernel specific parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nnsyn;     // total number of nmda synapses
uint4 *NSynLUT; // L=Nnsyn; look-up-tables for nmda synapses:
                // x - look-up-table of types of nmda synapses
                // y - look-up-table of axons provided nmda synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - look-up-table of a-types of nmda synapses
float *PreSyn_; // L=Nnsyn; temporal variables needed for calculation

//--- <kernel 1.5> presynaptic output is calculating according to
// Wang X-J, Tegner J, Constantinidis C, Goldman-Rakic PS (2004)
// Division of labor among distinct subtypes of inhibitory neurons in a cortical microcircuit of working memory.
// Proc. Natl. Acad. Sci. (USA), 101, 1368-1373
// (see, http://www.cns.nyu.edu/wanglab/publications/pdf/wang04.pnas.pdf)
// T*dOut/dt = out0-Out*(out0+1), where out0 = alpha*T*f, f is the "output" function of
// alpha synapse. The result is stored to array PreSyn, the intermediate result stored
// to array PreSyn_.
void kernel_nmda_presyn( uint id )
{
    // load to shared memory (?)
    uint ntype = NSynLUT[id].x; // make sure type < NNmType
    uint axon = NSynLUT[id].y;  // make sure axon < Nv
    uint out = NSynLUT[id].z;   // make sure out < Nout
    uint atype = NSynLUT[id].w; // make sure type < NAType
    // load from constant memory (?)
    float a0 = AType1[atype].x;
    float a = AType1[atype].y;
    float alpha = AType2[atype].x;
    float hv = AType2[atype].y;
    float slp = AType2[atype].z;
    float t = AType2[atype].w;
    float nmda_alpha = NmdaType1[ntype].x;
    float nmda_t = NmdaType1[ntype].y;
    // load from global memory
    float v = V[axon];
    float presyn_ = PreSyn_[id];
    float presyn = PreSyn[out];
    // math for alpha synapse
    float k = alpha*( a0+a/( 1.+exp(-( v-hv )/slp )))*t;
    float out0 = k/( k+1. );
    t = t/( k+1. );
    presyn_ = EXP_EULER( presyn_, out0, step, t );
    // math for nmda synapse
    k = nmda_alpha*nmda_t*presyn_;
    out0 = k/( k+1. );
    nmda_t = nmda_t/( k+1. );
    presyn = EXP_EULER( presyn, out0, step, nmda_t );
    // store to global memory
    PreSyn_[id] = presyn_;
    PreSyn[out] = presyn;
}

//--- call <kernel 1.5> for all compartments which involved in presynaptic
// processing for nmda synapses
void kernel_nmda_presyn( void )
{
    for( uint i = 0; i < Nnsyn; ++i ){
        kernel_nmda_presyn( i );
    }
}

#endif // __PRESYN_SKETCH_H
