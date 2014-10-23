//////////////////////////////////////////////////////////////////////////
// Stage 1. Presynaptic operations (spikes discrimination
//          and transformation of membrane potentials) which
//          are common for signals transmission between neurons and nerves
//------------------------------------------------------------------------
// optimization strategy: optimize the look-up-table in SpikeLUT array to
// provide coalescing memory accesses for V, Spike, Spike_ arrays.
//------------------------------------------------------------------------
// experiments: check if spike discrimination for all compartments faster then
// for soma(axon) compartments.

#ifndef __PRESYN_SKETCH_H
#define __PRESYN_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global variables
//------------------------------------------------------------------------
// constant memory 64 Kb
// for 1.x devices shared memory is 16 Kb per multiprocessor (warp of 32 threads),
// all devices >= 2.0 have 48 KB shared memory
//--- global parameters (READ-ONLY access, could be stored in the constant memory)
float Theshold; // the threshold for spike discrimination
uint Nv;        // total number of all compartments
uint Nsout;     // total number of all synapses in the network
uint NStdSyn;   // number of different types of "standars" synapses
float2 *StdSyn; // L=NStdSyn; parameters of "standard" synapses:
                // x - time constant
                // y - amplitude jump
uint NSigSyn;   // number of different types of sigma synapses
float2 *SigSyn; // L=NSigSyn; parameters of sigma synapses:
                // x - half-voltage
                // y - slope
uint NAlpSyn;   // number of different types of alpha synapses

float2 *ASynCnt;// constant parameters of sigma synapses:
                // x - a0
                // y - a
float4 *ASyn;   // parameters of alpha synapses:
                // x - alpha
                // y - half-voltage
                // z - slope
                // w - time constant

uint NNmdSyn;   // number of different types of nmda synapses
float2 *NmdaSyn;// L=NNmdSyn; parameters of nmda synapses:
                // x - alpha
                // y - T

//--- shared variables (READ/WRITE access, stored in global memory)
uint *Spike;    // L=Nv; 1 if spike is generated, otherwise 0
uint *Spike_;   // L=Nv; 1 if V >= Threshold, otherwise 0
float *V;       // L=Nv; membrane potential of all compartments
float *PreSyn;  // L=Nsout; presynaptic outputs

//////////////////////////////////////////////////////////////////////////
// <kernel 1.1> spike discrimination
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nspikes;   // total number of neurons (axon compartments) involved
                // into spikes discrimination
uint *SpikeLUT; // L=Nspikes; look-up-table of soma compartments (array *V)
                // involved into spikes discrimination

//--- <kernel 1.1> if a spike generated then 1 otherwise 0.
// Spike = ( V[current] <= Threshold && V[previous] > Threshold )? 1: 0;
// the result is stored in the array 'Spike'
void spike_presyn( void )
{
    for( uint i = 0; i < Nspikes; ++i ){
        // load to shared memory (?)
        uint axon = SpikeLUT[i];

        // load from global memory
        uint lt_tr = uint( V[axon] < Theshold );
        // math & store to global memory
        Spike[axon] = lt_tr*Spike_[axon];
        Spike_[axon] = !lt_tr;
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.2> presynaptic transformation of membrane potentials for
// standard synapses
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nssyn;     // total number of "standard" synapses
uint4 *SSynLUT; // L=Nssyn; look-up-tables for "standard" synapses:
                // x - look-up-table of types of "standard" synapses
                // y - look-up-table of axons provided synaptic inputs to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.2> presynaptic output is calculating according to <publication(s)>
// Out = Out*exp(-step/T )+Spike*Jump;
// the result is stored in array PreSyn.
void std_presyn( void )
{
    for( uint syn = 0; syn < Nssyn; ++syn ){
        // load to shared memory (?)
        uint type = SSynLUT[i].x;
        uint axon = SSynLUT[i].y;
        uint out = SSynLUT[i].z;
        // load from constant memory (?)
        float t = StdSyn[type].x;
        float a = StdSyn[type].y;
        // load from global memory
        uint spike = Spike[axon];
        float presyn = PreSyn[out];
        // math & store to global memory
        PreSyn[out] = presyn*exp(-step/t )+a*spike;
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.3> presynaptic transformation of membrane potentials for
// sigma synapses
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nsisyn;    // total number of sigma synapses
uint4 *SiSynLUT;// L=Nsisyn; look-up-tables for sigma synapses:
                // x - look-up-table of types of sigma synapses
                // y - look-up-table of axons provided sigma synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.3> presynaptic output is calculating according to <publication(s)>
// Out = 1./( 1.+exp(-( v-Hv )/Slp ))
// the result is stored in array PreSyn.
void sigma_presyn( void )
{
    for( uint syn = 0; syn < Nsisyn; ++syn ){
        // load to shared memory (?)
        uint type = SiSynLUT[i].x;
        uint axon = SiSynLUT[i].y;
        uint out = SiSynLUT[i].z;
        // load from constant memory (?)
        float hv = SigSyn[type].x;
        float slp = SigSyn[type].y;
        // load from global memory
        float v = V[axon];
        // math & store to global memory
        PreSyn[out] = 1./( 1.+exp(-( v-hv )/slp ));
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.4> presynaptic transformation of membrane potentials for
// alpha synapses (AMPA/GABA)
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nasyn;     // total number of alpha synapses
uint4 *ASynLUT; // L=Nasyn; look-up-tables for alpha synapses:
                // x - look-up-table of types of alpha synapses
                // y - look-up-table of axons provided alpha synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - reserved

//--- <kernel 1.4> presynaptic output is calculating according to <publication(s)>
// T*dOut/dt = out0-Out*(out0+1), where out0 = alpha*T*f, f = a0+a/(1+exp(-( v-Hv )/Slp ))
// the result is stored in array PreSyn.
void alpha_presyn( void )
{
    for( uint syn = 0; syn < Nasyn; ++syn ){
        // load to shared memory (?)
        uint type = ASynLUT[i].x;
        uint axon = ASynLUT[i].y;
        uint out = ASynLUT[i].z;
        // load from constant memory (?)
        float a0 = ASynCnt[type].x;
        float a = ASynCnt[type].y;
        float alpha = ASyn[type].x;
        float hv = ASyn[type].y;
        float slp = Asyn[type].z;
        float t = Asyn[type].w;

        // load from global memory
        float v = V[axon];
        float k = alpha*( a0+a/( 1.+exp(-( v-hv )/slp )))*t;
        float out0 = k/( k+1. );
        float presyn = PreSyn[out];
        t = t/( k+1. );
        // math
        presyn = EXP_EULER( presyn, out0, step, t );
        // store to global memory
        PreSyn[out] = presyn;
    }
}

//////////////////////////////////////////////////////////////////////////
// <kernel 1.5> presynaptic transformation of membrane potentials for
// nmda synapses
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Nnsyn;     // total number of nmda synapses
uint4 *NSynLUT; // L=Nnsyn; look-up-tables for nmda synapses:
                // x - look-up-table of types of nmda synapses
                // y - look-up-table of axons provided nmda synapses to other neurons
                // z - look-up-table of presynaptic outputs
                // w - look-up-table of a-types of nmda synapses
float *PreOut;  // L=Nnsyn; temporal variables needed for calculation

//--- <kernel 1.5> presynaptic output is calculating according to <publication(s)>
// T*dOut/dt = out0-Out*(out0+1), where out0 = alpha*T*f, f is the "output" function of alpha synapse
// the result is stored to array PreSyn, the intermediate result stored to array NmdaOut.
void nmda_presyn( void )
{
    for( uint syn = 0; syn < Nnsyn; ++syn ){
        // load to shared memory (?)
        uint type = NSynLUT[i].x;
        uint axon = NSynLUT[i].y;
        uint out = NSynLUT[i].z;
        uint atype = NSynLUT[i].w;
        // load from constant memory (?)
        float a0 = ASynCnt[atype].x;
        float a = ASynCnt[atype].y;
        float alpha = ASyn[atype].x;
        float hv = ASyn[atype].y;
        float slp = Asyn[atype].z;
        float t = Asyn[atype].w;
        float nmda_alpha = NmdaSyn[type].x;
        float nmda_t = NmdaSyn[type].y;
        // load from global memory
        float v = V[axon];
        float pre_out = PreOut[syn];
        float nmda_out = PreSyn[out];
        // math
        float k = alpha*( a0+a/( 1.+exp(-( v-hv )/slp )))*t;
        float out0 = k/( k+1. );
        t = t/( k+1. );
        pre_out = EXP_EULER( pre_out, out0, step, t );
        k = nmda_alpha*pre_out*nmda_t;
        out0 = k/( k+1. );
        nmda_t = nmda_t/( k+1. );
        nmda_out = EXP_EULER( nmda_out, out0, step, nmda_t );
        // store to global memory
        PreOut[syn] = pre_out;
        PreSyn[out] = nmda_out;
    }
}

#endif // __PRESYN_SKETCH_H
