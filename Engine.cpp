

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairDstAna/PairDstAna.h"
#include "DeltaTofPdfs/DeltaTofPdfAna.h"
#include "DeltaTofPdfs/GeneratePdf.h"
#include "JPsi/JPsiHistogram.h"
#include "JPsi/JPsiMvaHistogram.h"
#include "JPsi/JPsiFitter.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

	Logger::setGlobalLogLevel( "none" );

	TaskFactory::registerTaskRunner<PairDstAna>( "PairDstAna" );
	TaskFactory::registerTaskRunner<DeltaTofPdfAna>( "DeltaTofPdfAna" );
	TaskFactory::registerTaskRunner<GeneratePdf>( "GeneratePdf" );
	TaskFactory::registerTaskRunner<JPsiHistogram>( "JPsiHistogram" );
	TaskFactory::registerTaskRunner<JPsiMvaHistogram>( "JPsiMvaHistogram" );
	TaskFactory::registerTaskRunner<JPsiFitter>( "JPsiFitter" );
	

	TaskEngine engine( argc, argv, "PairDstAna" );


	return 0;
}
