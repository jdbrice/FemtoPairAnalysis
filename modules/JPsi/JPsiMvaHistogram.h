#ifndef JPSI_MVA_HISTOGRAM_H
#define JPSI_MVA_HISTOGRAM_H


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


class JPsiMvaHistogram : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float nSigmaDeltaY = 3.0;
	float nSigmaDeltaZ = 3.0;

	HistoBins ptBins;

	float pid = 1.0;

public:

	const int DEBUG = 1;
	JPsiMvaHistogram() {}
	~JPsiMvaHistogram() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		nSigmaDeltaY = config.getFloat( "p.nSigmaDeltaY", 3.0 );
		nSigmaDeltaZ = config.getFloat( "p.nSigmaDeltaZ", 3.0 );

		ptBins.load( config, "bins.pt" );

		pid = config.get<float>( "p.pid", 1.0 );
	}
protected:

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

		if ( pair->d1_mPt < 1.5 && pair->d2_mPt < 1.5 )
			return;

		if ( pair->d1_mPt < 1.3 || pair->d2_mPt < 1.3 )
			return;

		if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 )
			return;
		if ( fabs(lv.Rapidity()) > 0.5 )
			return;
		

		if ( pair->d1_mPid < pid || pair->d2_mPid < pid )
			return;

		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			book->fill( "uls", lv.M(), lv.Pt() );
		} // 0 == pair->mChargeSum


	} //analyzeEvent


	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		

	}




};


#endif

