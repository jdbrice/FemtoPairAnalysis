

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
#include "MixedPidAna.h"
#include "DcaPdfAna.h"
#include "SignalMaker.h"
#include "BackgroundMaker.h"
#include "JPsi/JPsiHistogram.h"
#include "JPsi/JPsiMvaHistogram.h"
#include "JPsi/JPsiFitter.h"
#include "ComparePid.h"

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
	TaskFactory::registerTaskRunner<BackgroundMaker>( "BackgroundMaker" );
	TaskFactory::registerTaskRunner<MixedEventAna>( "MixedEventAna" );
	TaskFactory::registerTaskRunner<MixedPidAna>( "MixedPidAna" );
	TaskFactory::registerTaskRunner<DcaPdfAna>( "DcaPdfAna" );
	TaskFactory::registerTaskRunner<JPsiHistogram>( "JPsiHistogram" );
	TaskFactory::registerTaskRunner<JPsiMvaHistogram>( "JPsiMvaHistogram" );
	TaskFactory::registerTaskRunner<JPsiFitter>( "JPsiFitter" );

    TaskFactory::registerTaskRunner<ComparePid>( "ComparePid" );
	

	TaskEngine engine( argc, argv, "PairDstAna" );


	return 0;
}
