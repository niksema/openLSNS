// convertor.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "load_nsm.h"

int _tmain(int argc, _TCHAR* argv[])
{
	nsm_model model;
	if( model.load_model("peripheralChemo.nsm")){
		model.save_model("new_peripheralChemo.nns");
	}
	getchar();
	return 0;
}

