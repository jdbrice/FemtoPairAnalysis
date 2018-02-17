

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairDstAna/PairDstAna.h"
#include "PairDstAna/LikeSignRatioAna.h"
#include "DeltaTofPdfs/DeltaTofPdfAna.h"
#include "DeltaTofPdfs/GeneratePdf.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	// loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

	Logger::setGlobalLogLevel( "none" );

	TaskFactory::registerTaskRunner<PairDstAna>( "PairDstAna" );
	TaskFactory::registerTaskRunner<LikeSignRatioAna>( "LikeSignRatioAna" );
	TaskFactory::registerTaskRunner<DeltaTofPdfAna>( "DeltaTofPdfAna" );
	TaskFactory::registerTaskRunner<GeneratePdf>( "GeneratePdf" );
	

	TaskEngine engine( argc, argv, "PairDstAna" );


	return 0;
}
