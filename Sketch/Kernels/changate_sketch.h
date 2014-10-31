//////////////////////////////////////////////////////////////////////////
// Stage 5. Calculating the gate variables (activation and inactivation)
// for ion currents
//------------------------------------------------------------------------
//

#ifndef __CHANGATE_SKETCH_H
#define __CHANGATE_SKETCH_H
//////////////////////////////////////////////////////////////////////////
// global parameters and variables
//--- global parameters (READ-ONLY access, could be stored in the constant
//    memory)
extern uint Nv;         // total number of all compartments
extern uint Ngates;     // total number of all gate variables
//--- shared variables (READ/WRITE access, stored in global memory)
extern float *V;        // L=Nv: the array of membrane potentials
extern float *Gates;    // L=Ngates: array of gate variables for specific channels


//////////////////////////////////////////////////////////////////////////
// local parameters and variables
//--- local parameters (READ-ONLY access, stored in global memory and could be cached in the shared memory)
uint *V_LUT;    // L=N: look-up-table for compartments associated with the specific gate variable



        //--- <kernel ???>. pseudocode for calculating of gate variables
        // The result is a gate variable (activation ar inactivation)
        // stored in array Gates. Prop is a description of correspondent
        // type of gate variable (voltage dependent, z-type etc)
        void calc_gate( void )
        {
            for( uint i = 0; i < N; ++i ){  // loop for all specific gate variables
                uint vi = V_LUT[i];         // get the index of compartment associated with the gate variable to be processed
                Gates[i] = f( V[vi] );      // f specifyes operation to calculate a gate variable
                GatesP[i] = pow( Gates[i], Pow );
            }
        }

        //--- generic description of gate variable for voltage dependent channel
        // GENERIC_G = 1./( 1.+exp(-( v-Hv )/Slp ))
        // GENERIC_T = Tmax/( cosh(-( v-Hvt )/Slpt ))
        float f( float v )
        {
            if( instant_gate ){
                return GENERIC_G( 1., v, Hv, Slp );
            }
            else{
                float g0 = GENERIC_G( 1., v, Hv, Slp );
                float t = T0+GENERIC_T( Tmax, v, Hvt, Slpt );
                return EXP_EULER( Gate, gate0, step, t );
            }
        }

        //--- modified generic description of gate variable for voltage dependent channel
        // GENERIC_G = 1./( 1.+exp(( v-Hv )/Slp ))
        // GENERIC_T_MOD = 2.*Tmax/( exp(-( v-Hvt )/Slpt )+exp(( v-Hvmt )/Slpt1 ))
        float f( float v )
        {
            if( instant_gate ){
                return GENERIC_G( 1., v, Hv, Slp );
            }
            else{
                float g0 = GENERIC_G( 1., v, Hv, Slp );
                float t = T0+GENERIC_T_MOD( Tmax, v, Hvt, Slpt, Slpt1 );
                return EXP_EULER( Gate, gate0, step, t );
            }
        }

        //--- alpha/beta description of gate variable for voltage dependent channel
        // ALPHA_BETA = A*( B*v-Hv )/( C+exp(-( v-Hv )/Slp ))
        float f( float v )
        {
            if( instant_gate ){
                float alpha = ALPHA_BETA( v, Hva, Slpa, Aa, Ba, Ca );
                float beta = ALPHA_BETA( v, Hvb, Slpb, Ab, Bb, Cb );
                return alpha/( alpha+beta );
            }
            else{
                float alpha = ALPHA_BETA( v, Hva, Slpa, Aa, Ba, Ca );
                float beta = ALPHA_BETA( v, Hvb, Slpb, Ab, Bb, Cb );
                float t = 1./( alpha+beta );
                float gate0 = alpha*t;
                t = T0+t*Tmax;
                return EXP_EULER( Gate, gate0, step, t );
            }
        }

        //--- generic description of gate variable for ca-activated channel (z-channel)
        float f( float v )
        {
            if( instant_gate ){
                float gexp1 = exp( -( v-Hv )/Slp );
                float gexp2 = 1./gexp1;
                float in = In_ca*K1*gexp1;
                return in/( in+K2*gexp2 );
            }
            else{
                float gexp1 = exp( -( v-Hv )/Slp );
                float gexp2 = 1./gexp1;
                float in = In_ca*K1*gexp1;
                float gate0 = in/( in+K2*gexp2 );
                float t = T0+Tmax/( in+Gamma*gexp2 );
                return EXP_EULER( Gate, gate0, step, t );
            }
        }

        //--- alpha/beta description of gate variable for ca-activated channel (z-channel)
        float f( float v )
        {
            if( instant_gate ){
                float k = Alpha*pow( In_ca*Beta, Lymbda );
                return k/( k+1.);
            }
            else{
                float k = Alpha*pow( In_ca*Beta, Lymbda );
                float gate0 = k/( k+1.);
                float t = T0+Tmax/( Gamma*k+1. );
                return EXP_EULER( Gate, gate0, step, t );
            }
        }


#endif // __CHANGATE_SKETCH_H
