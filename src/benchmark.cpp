// includes, system
#include <argp.h>

#include "info.h"
#include"kernel.h"

#define NUM_CALLS 100000
const char *argp_program_version = "benchmark version 1";
const char *argp_program_bug_address = "<your@email.address>";
static char doc[] = "Your program description.";
static char args_doc[] = "[FILENAME]...";
static struct argp_option options[] = { 
    { "info", 'i', 0, 0, "Display information about the target system"},
    { "kernel", 'k', 0, 0, "CUDA kernel benchmark"},
    { "word", 'w', 0, 0, "Compare words instead of characters."},
    { 0 } 
};

struct arguments{
    enum { CHARACTER_MODE, WORD_MODE, LINE_MODE } mode;
    bool isCaseInsensitive;
};

static error_t parse_opt( int key, char *arg, struct argp_state *state )
{
    float time = 0.;
    struct arguments *arguments = ( struct arguments *)state->input;
    switch( key ){
	  case 'i': 
	    arguments->isCaseInsensitive = true;
	    info( "benchmark" );
	    break;
	  case 'k': 
	    arguments->mode = arguments::LINE_MODE;
	    printf("Start the CUDA kernel benchmark. Please wait ...\n" );
	    time = call_kernel( NUM_CALLS, 16, 16, 192 );
	    printf("Kernel performance: total time %f msec, time per call %f msec\n", time, time/NUM_CALLS );
	    break;
	  case 'w': 
	    arguments->mode = arguments::WORD_MODE; 
	    printf("word\n" );
	    break;
	  case ARGP_KEY_ARG:
	    printf("ARGP_KEY_ARG\n" );
	    return 0;
	  default: 
	    return ARGP_ERR_UNKNOWN;
    }   
   
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

////////////////////////////////////////////////////////////////////////////////
// Program main
int main( int argc, char **argv )
{
    struct arguments arguments;
    arguments.mode = arguments::CHARACTER_MODE;
    arguments.isCaseInsensitive = false;
    return argp_parse( &argp, argc, argv, 0, 0, &arguments );
}
