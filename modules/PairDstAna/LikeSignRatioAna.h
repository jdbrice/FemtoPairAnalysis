#ifndef LIKE_SIGN_RATIO_ANA_H
#define LIKE_SIGN_RATIO_ANA_H

#include "TreeAnalyzer.h"


#include "TLorentzVector.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"

#include "vendor/loguru.h"

#include <algorithm>

class LikeSignRatioAna : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

public:

	const int DEBUG = 1;
	LikeSignRatioAna() {}
	~LikeSignRatioAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

	}
protected:

	void like_sign( float pidMin, float pidMax, TLorentzVector &lv1, TLorentzVector &lv2, TLorentzVector &lv ){
		book->fill( "ls", lv.M(), pidMin );

		mixed_pair( "ls_mixed_", 0.1, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.2, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.3, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.4, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.5, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.6, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.7, pidMin, pidMax, lv );
		mixed_pair( "ls_mixed_", 0.8, pidMin, pidMax, lv );
	}

	void unlike_sign( float pidMin, float pidMax, TLorentzVector &lv1, TLorentzVector &lv2, TLorentzVector &lv ){
		book->fill( "uls", lv.M(), pidMin );

		mixed_pair( "uls_mixed_", 0.1, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.2, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.3, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.4, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.5, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.6, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.7, pidMin, pidMax, lv );
		mixed_pair( "uls_mixed_", 0.8, pidMin, pidMax, lv );
		
	}

	void mixed_pair( string prefix, float X, float pidMin, float pidMax, TLorentzVector &lv ){
		// LOG_F( INFO, "filling X=%f --> %s", X, (prefix + dtes(X, "p") ).c_str() );
		if ( pidMin < X && pidMax >= X )
			book->fill( prefix + dtes( X, "p" ), lv.M() );
	}


	virtual void analyzeEvent(){
		FemtoPair * pair = this->_fpr.get();

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// PAIR kinematics
		if ( fabs( lv.Rapidity() ) > 0.5 )
			return;
		// if ( lv.Pt() < minPairPt )
		// 	return;

		// DAUGHTER kinematics
		if ( fabs( lv1.PseudoRapidity() ) > 0.5 )
			return;
		if ( lv1.Pt() < 1.1 )
			return;
		if ( fabs( lv2.PseudoRapidity() ) > 0.5 )
			return;
		if ( lv2.Pt() < 1.1 )
			return;

		// TRIGGER FLAG
		if ( pair->d1_mTriggerFlag <= 0 )
			return;
		if ( pair->d2_mTriggerFlag <= 0 )
			return;

		float pidMin = std::min( pair->d1_mPid, pair->d2_mPid );
		float pidMax = std::max( pair->d1_mPid, pair->d2_mPid );


		if ( 0 != pair->mChargeSum ){
			like_sign( pidMin, pidMax, lv1, lv2, lv );
		} else {
			unlike_sign( pidMin, pidMax, lv1, lv2, lv );
		}
	}



	virtual void postMake(){


		TH1* hlsr = book->get("lsr");

		for ( float i :  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8} ){
			LOG_F( INFO, "Ratio @ %f", i );
			int bin = hlsr->GetXaxis()->FindBin( i );

			TH1 * huls = book->get( "uls_mixed_" + dtes( i ) );
			TH1 * hls = book->get( "ls_mixed_" + dtes( i ) );

			int binM1 = huls->GetXaxis()->FindBin(0.0);
			int binM2 = huls->GetXaxis()->FindBin(2.6);

			double Euls = 0, Els = 0;
			float Iuls = huls->IntegralAndError( binM1, binM2, Euls );
			LOG_F( INFO, "Iuls = %f +/- %f", Iuls, Euls );
			float Ils = hls->IntegralAndError( binM1, binM2, Els );
			LOG_F( INFO, "Ils = %f +/- %f", Ils, Els );

			LOG_F(INFO, "");
			float v = Iuls / Ils;
			hlsr->SetBinContent( bin, v );
			hlsr->SetBinError( bin, v * (( Euls / Iuls ) + ( Els / Ils )) );


		}

	}


};

#endif