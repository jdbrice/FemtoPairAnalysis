

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairDstAna/PairDstAna.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

	Logger::setGlobalLogLevel( "none" );

	TaskFactory::registerTaskRunner<PairDstAna>( "PairDstAna" );


	TaskEngine engine( argc, argv, "PairDstAna" );
	
	return 0;
}
