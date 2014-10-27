//////////////////////////////////////////////////////////////////////////
//    Stage 3. The synaptic fatigue (or depression), synaptic enhancement,
//    synaptic plastisity that involes the NMDA and/or AMPA glutamate
//    receptors should be implemented here later.

#ifndef __ADAPT_SKETCH_H
#define __ADAPT_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global variables
//--- global parameters (READ-ONLY access, could be stored in the constant memory)


//////////////////////////////////////////////////////////////////////////
// <kernel 3.1>
//------------------------------------------------------------------------
// kernel specific variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint Npsyn;     // total number of
uint4 *PSynLUT; // L=Npsyn; look-up-tables for
                // x - look-up-table of types of synaptic adaptation
                // y - look-up-table of dendrites which receive synaptic input
                // z - look-up-table of outputs
                // w - reserved

//--- <kernel 3.1>
void adapt_presyn( void )
{
    for( uint syn = 0; syn < Npsyn; ++syn ){
        // load to shared memory (?)
        uint type = PSynLUT[i].x;
        uint dendrite = PSynLUT[i].y;
        uint out = PSynLUT[i].z;
        //...
    }

}
#endif // __ADAPT_SKETCH_H
