

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "PairDstAna/PairDstAna.h"
#include "PairDstAna/LikeSignRatioAna.h"
#include "PairDstAna/MixedLSR.h"
#include "DeltaTofPdfs/DeltaTofPdfAna.h"
#include "DeltaTofPdfs/GeneratePdf.h"
#include "MixedEventAna.h"
#include "DcaPdfAna.h"
#include "SignalMaker.h"


#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	// loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);

	Logger::setGlobalLogLevel( "none" );

	TaskFactory::registerTaskRunner<PairDstAna>( "PairDstAna" );
	TaskFactory::registerTaskRunner<LikeSignRatioAna>( "LikeSignRatioAna" );
	TaskFactory::registerTaskRunner<MixedLSR>( "MixedLSR" );
	TaskFactory::registerTaskRunner<DeltaTofPdfAna>( "DeltaTofPdfAna" );
	TaskFactory::registerTaskRunner<GeneratePdf>( "GeneratePdf" );
	TaskFactory::registerTaskRunner<SignalMaker>( "SignalMaker" );
	TaskFactory::registerTaskRunner<MixedEventAna>( "MixedEventAna" );
	TaskFactory::registerTaskRunner<DcaPdfAna>( "DcaPdfAna" );
	

	TaskEngine engine( argc, argv, "PairDstAna" );


	return 0;
}
