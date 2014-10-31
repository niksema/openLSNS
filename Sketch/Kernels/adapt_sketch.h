//////////////////////////////////////////////////////////////////////////
//    Stage 3. The synaptic fatigue (or depression), synaptic enhancement,
//    synaptic plastisity that involes the NMDA and/or AMPA glutamate
//    receptors should be implemented here later.

#ifndef __ADAPT_SKETCH_H
#define __ADAPT_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global parameters and variables
//--- global parameters (READ-ONLY access, could be stored in the constant
//    memory)
extern uint Nv;         // total number of all compartments
extern uint Nnmdasyn;   // total number of synapses with adaptation
extern uint NNmType;    // number of different types of nmda synapses
extern float4 *NmdaType2;// L=NNmdSyn; parameters of nmda synapses:
                        // x - k
                        // y - half-voltage
                        // z - slope
                        // w - reserved
//--- shared variables (READ/WRITE access, stored in global memory)
extern float *V;        // L=Nv; membrane potential of all compartments
extern float *MgIn;     // L=Nin; concentration of Mg ions inside cell


//////////////////////////////////////////////////////////////////////////
// <kernel 3.1>
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint4 *NmdaLUT; // L=Nnmdasyn; look-up-tables for
                // x - look-up-table of types of synaptic adaptation
                // y - look-up-table of dendrites which receive synaptic input
                // z - concentration of Mg+ inside the cell
                // w - look-up-table of outputs
float *PreSyn;  // L=Nadaptsyn; outputs of synaptic plasticity

//--- <kernel 3.1>
void adapt_nmda( uint id )
{
    // load from shared memory (?)
    uint type = NmdaLUT[id].x;  // make sure type < NNmType
    uint dendr = NmdaLUT[id].y; // make sure dendr < Nv
    uint ion = NmdaLUT[id].z;   // make sure dendr < Nin
    uint out = NmdaLUT[id].w;   // make sure out < Nnmdasyn
    // load from constant memory
    float km = NmdaType2[type].x;
    float hv = NmdaType2[type].y;
    float slp = NmdaType2[type].z;
    // load from global memory
    float mg = MgIn[ion];
    float v = V[dendr];
    // math & store to global memory
    PreSyn[out] = 1./( 1.+mg/km*exp(-( v-hv )/slp );
}

//--- call <kernel 3.1> for all nmda synapses
void adapt_nmda( void )
{
    for( uint i = 0; i < Nnmdasyn; ++i ){
        kernel_adapt_nmda( i );
    }
}

#endif // __ADAPT_SKETCH_H
