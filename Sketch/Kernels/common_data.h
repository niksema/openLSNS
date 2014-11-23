//////////////////////////////////////////////////////////////////////////
// global parameters and variables for simulation kernels
#ifndef __COMMON_DATA_H
#define __COMMON_DATA_H
//--- global parameters (READ-ONLY access, could be stored in the constant
//    memory)

//--- general parameters
extern float Theshold;  // the threshold for spike discrimination
extern float Step;      // the integration step
extern uint Nv;         // total number of all compartments
extern uint Nout;       // total number of all synapses (either type) in the network

//--- parameters for calculating the normalize concentration of neuro-transmitter
//    (as sigmoid function) based on membrane potential of pre-synaptic neurons
extern uint NSiType;    // number of different models for calculating of normalize concentration
extern float2 *SiType;  // L=NSiType; parameters for calculation of transmitter concentration
                        // x - half-voltage
                        // y - slope

//--- parameters of pulse synapses
//    [##ref]
extern uint NStType;    // number of different types of pulse synapses
extern float2 *StType;  // L=NStSyn; parameters of "standard" synapses:
                        // x - exp(-Step/T) where T is a time constant
                        // y - alpha (maximal concentration of transmitter)

//--- parameters of AMPA/GABAa/NMDA synapses
//    http://senselab.med.yale.edu/modeldb/showmodel.asp?model=18500&file=\Channels_NEW\ampa.mod
//    http://senselab.med.yale.edu/modeldb/showmodel.asp?model=18500&file=\Channels_NEW\gabaa.mod
//    http://senselab.med.yale.edu/modeldb/showmodel.asp?model=18500&file=\Channels_NEW\nmda.mod
//    http://cns.iaf.cnrs-gif.fr/files/Dlgn96.pdf
//    http://www.ane.pl/pdf/6520.pdf
//    Kinetic of transmitter release for AMPA/GABA(a)/NMDA synapses:
//    T*dM/dt = -M+Alpha*F(Vpre)*(1-M): where F(Vpre) is normalize concentration of neurotransmetter
//    Synaptic current for AMPA/GABA(a):
//    Isyn = Gmax*M*(V-Esyn)
//    Synaptic current for NMDA:
//    Isyn = Gmax*M*Z(V-Esyn)
//    where Z = 1/(1+[Mg]*exp(-0.062V)/3.57) is magnezium block of NMDA synapse
extern uint NAType;     // number of different types of AMPA/GABAa synapses
extern float2 *AType;   // L=NAType; parameters of AMPA/GABAa synapses:
                        // x - Alpha (maximal concentration of transmitter)
                        // y - T (time constant)


//--- parameters of GABAb synapses
//    http://senselab.med.yale.edu/modeldb/showmodel.asp?model=18500&file=\Channels_NEW\gabab.mod
//    http://cns.iaf.cnrs-gif.fr/files/Dlgn96.pdf


//--- parameters of NMDA synapses
//    Kinetic of transmitter release for NMDA synapse:
//    T*dM/dt = -M+Alpha*F(Vpre)*(1-M): where F(Vpre) is normalize concentration of neurotransmetter
//    Synaptic current:
//    Isyn = Gmax*M*Z*(V-Esyn)
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

#endif // __COMMON_DATA_H
