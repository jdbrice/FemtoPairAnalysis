#ifndef COMPARE_PID_H
#define COMPARE_PID_H

#include "TreeAnalyzer.h"
#include "XmlRange.h"
#include "RooPlotLib.h"
#include "Reporter.h"
#include "FitConfidence.h"
#include "HistoBins.h"

#include "vendor/loguru.h"

#include "TNamed.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TLatex.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"

#include <algorithm>


class ComparePid : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float nSigmaDeltaY   = 3.0;
	float nSigmaDeltaZ   = 3.0;
	float nSigmaDeltaTof = 3.0;
	float deltaTofMin    = -0.75;
	float deltaTofMax    = 0.75;

	float pairPidMin = 1.30;

public:

	const int DEBUG = 1;
	ComparePid() {}
	~ComparePid() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		nSigmaDeltaY   = config.getFloat( "p.nSigmaDeltaY", 3.0 );
		nSigmaDeltaZ   = config.getFloat( "p.nSigmaDeltaZ", 3.0 );
		nSigmaDeltaTof = config.getFloat( "p.nSigmaDeltaTof", 3.0 );
		deltaTofMin    = config.getFloat( "p.DeltaTof:min", -1.0 );
		deltaTofMax    = config.getFloat( "p.DeltaTof:max", 1.0 );

		pairPidMin     = config.get<float>( "p.pairPidMin", 1.3 );

		LOG_F( INFO, "nSigmaDeltaY = %f", nSigmaDeltaY );
		LOG_F( INFO, "nSigmaDeltaZ = %f", nSigmaDeltaZ );
		// LOG_F( INFO, "nSigmaDeltaTof = %f", nSigmaDeltaTof );
		LOG_F( INFO, "%0.2f < deltaTof < %0.2f", deltaTofMin, deltaTofMax );
	}
protected:

	bool passDeltaY( float pt, float dy ){
		double sigy = -17.6867 + 18.4528*exp(0.637142/pt);
		double marg = 0.;
		if(pt>3.) marg = 0.5;

		if( fabs( dy ) <= ( nSigmaDeltaY + marg ) * sigy ) 
			return true;
		
		return false;
	}

	bool passDeltaZ( float pt, float dz ){
		double sigz = -32.6793 + 32.6034 * exp( 0.444217 / pt );
		double marg = 0.;
		if(pt>3.) 
			marg = 0.5;
		if( fabs( dz ) <= ( nSigmaDeltaZ + marg ) * sigz ) 
			return true;
		return false;
	}

	bool passDeltaTof( float pt, float dt ){
		double sigt = 0.1009 + 0.009116 * exp( 3.474 / pt );
		if( fabs( dt ) <= ( nSigmaDeltaTof ) * sigt ) 
			return true;
		return false;
	}

	virtual void analyzeEvent(){
		FemtoPair * pair = this->_fpr.get();

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// TRIGGER FLAG
		if ( pair->d1_mTriggerFlag <= 0 )
			return;
		if ( pair->d2_mTriggerFlag <= 0 )
			return;

		if ( pair->d1_mPt < 1.0 || pair->d2_mPt < 1.0 )
			return;

		if ( pair->d1_mDCA > 1.0 || pair->d2_mDCA > 1.0 )
			return;

		if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 )
			return;
		
		if ( fabs(lv.Rapidity()) > 0.5 )
			return;

		float pairPid = sqrt( pow( pair->d1_mPid, 2 ) + pow( pair->d2_mPid, 2 ) );
		if ( 0 == pair->mChargeSum && pairPid > pairPidMin ){
			book->fill( "mva", lv.M() );
		} // 0 == pair->mChargeSum


		if ( false == passDeltaY( pair->d1_mPt, pair->d1_mDeltaY ) || false == passDeltaY( pair->d2_mPt, pair->d2_mDeltaY ) )
			return;
		if ( false == passDeltaZ( pair->d1_mPt, pair->d1_mDeltaZ ) || false == passDeltaZ( pair->d2_mPt, pair->d2_mDeltaZ ) )
			return;

		if ( pair->d1_mMatchFlag != 1 || pair->d2_mMatchFlag != 1 )
			return;

		// if ( false == passDeltaTof( pair->d1_mPt, pair->d1_mDeltaTimeOfFlight )  )
		// 	return;
		// if ( false == passDeltaTof( pair->d2_mPt, pair->d2_mDeltaTimeOfFlight )  )
		// 	return;

		if ( pair->d1_mDeltaTimeOfFlight > deltaTofMax || pair->d1_mDeltaTimeOfFlight < deltaTofMin )
			return;
		if ( pair->d2_mDeltaTimeOfFlight > deltaTofMax || pair->d2_mDeltaTimeOfFlight < deltaTofMin )
			return;

		if ( pair->d1_mNSigmaPion < -1 || pair->d1_mNSigmaPion > 3 )
			return;
		if ( pair->d2_mNSigmaPion < -1 || pair->d2_mNSigmaPion > 3 )
			return;

		if ( 0 == pair->mChargeSum ){
			book->fill( "tra", lv.M() );
		} // 0 == pair->mChargeSum


	} //analyzeEvent


	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		TH1 * hmva = book->get( "mva" );
		TH1 * htra = book->get( "tra" );

		hmva->Sumw2();
		htra->Sumw2();

		hmva->Scale( 1.0, "width");
		htra->Scale( 1.0, "width");

		TAxis *x = hmva->GetXaxis();

		int x1 = x->FindBin( 1.1 );
		int x2 = x->FindBin( 2.5 );

		hmva->Scale( 1.0 / hmva->Integral( x1, x2 ) );
		htra->Scale( 1.0 / htra->Integral( x1, x2 ) );


		

	}




};


#endif

