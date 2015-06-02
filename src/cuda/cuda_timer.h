#ifndef __CUDA_TIMER__
#define __CUDA_TIMER__

#include <cuda.h>

class cuda_timer{
	public:
		cuda_timer( void ) : IsStarted( false ){
			cudaEventCreate( &Start );
			cudaEventCreate( &Stop );
		}
		~cuda_timer( void ){
			cudaEventDestroy( Start );
			cudaEventDestroy( Stop );
		}
	public:
		void start( cudaStream_t s = 0 ){
			if( !IsStarted ){
				cudaEventRecord( Start, s ); 
				IsStarted = true;
			}
		};
		void stop( cudaStream_t s = 0 ){
			if( IsStarted ){
				cudaEventRecord( Stop, s ); 
				IsStarted = false;
			}
		};
		float elapsed( void ){
			if( IsStarted ){
				return -1.0; 
			}
			cudaEventSynchronize( Stop );
			float elapsed = 0;
			cudaEventElapsedTime( &elapsed, Start, Stop );
			return elapsed;
		};
	private:
		bool IsStarted;
		cudaEvent_t Start;
		cudaEvent_t Stop;
};

#endif /*__CUDA_TIMER__*/