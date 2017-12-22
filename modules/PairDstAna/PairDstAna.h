#ifndef PAIR_DST_ANA_H
#define PAIR_DST_ANA_H


#include "TreeAnalyzer.h"
#include "XmlRange.h"
#include "RooPlotLib.h"

#include "vendor/loguru.h"

#include "TNamed.h"
#include "TTree.h"
#include "TLorentzVector.h"


#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"

#include <algorithm>



template <>
TString XmlConfig::get<TString>( string path ) const {
	TString r( getString( path ) );
	return r;
}

template <>
TString XmlConfig::get<TString>( string path, TString dv ) const {
	if ( !exists( path ) )
		return dv;
	TString r( getString( path ) );
	return r;
}

class PairDstAna : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;

	float pidCut = -10;
	float bgMax = 0.0;
	float bgMin = 0.0;
	float minPairPt = 0.0;
	float bgScale = 1.0;
public:

	const int DEBUG = 1;
	PairDstAna() {}
	~PairDstAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		pidCut = config.getFloat( "p.PidCut", 0.8 );
		bgMax = config.getFloat( "p.BgMax", 0.8 );
		bgMin = config.getFloat( "p.BgMin", 0.2 );
		minPairPt = config.getFloat( "p.MinPairPt", 0.0 );
		bgScale = config.getFloat( "p.BgScale", 1.0 );

	}
protected:


	virtual void analyzeEvent(){
		FemtoPair * pair = this->_fpr.get();

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// PAIR kinematics
		if ( fabs( lv.Rapidity() ) > 0.5 )
			return;
		if ( lv.Pt() < minPairPt )
			return;

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


		////////////////////////////////////////////////////////////////////////
		/// Like-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 != pair->mChargeSum ){
			
			if ( true ){
				book->fill( "shapecorr_D", lv.M() );
				book->fill( "shapecorr_D_rb", lv.M() );
				if ( pidMax > pidCut ){
					book->fill( "shapecorr_N", lv.M() );
					book->fill( "shapecorr_N_rb", lv.M() );
				}
			}

			TH2 * hdMdPidMin = (TH2*)book->get("lsdMdPidMin");
			int pidBin = hdMdPidMin->GetYaxis()->FindBin( pidMin );
			int mBin = hdMdPidMin->GetXaxis()->FindBin( lv.M() );

			for ( int iPidBin = 1; iPidBin <= pidBin; iPidBin++ ){
				float iPid = hdMdPidMin->GetYaxis()->GetBinCenter( iPidBin );
				hdMdPidMin->Fill( lv.M(), iPid );
			}

			if (pair->d1_mPid > pidCut && pair->d2_mPid > pidCut) {
				book->fill("lscand", lv.M() );
				book->fill("lscand_rb", lv.M() );
			}


		} // 0 != pair->mChargeSum


		////////////////////////////////////////////////////////////////////////
		/// Opposite-Sign
		////////////////////////////////////////////////////////////////////////
		if ( 0 == pair->mChargeSum ){
			book->fill( "uls", lv.M() );
			book->fill( "uls_rb", lv.M() );

			if (pair->d1_mPid > pidCut || pair->d2_mPid > pidCut) {
				book->fill("scand", lv.M() );
				book->fill("scand_rb", lv.M() );
			}

			if (pair->d1_mPid > pidCut && pair->d2_mPid > pidCut) {
				book->fill("pcand", lv.M() );
				book->fill("pcand_rb", lv.M() );
			}

			


			if ( 	pidMax > bgMin && pidMin < bgMax
					// && (pair->d1_mDCA > 0.5 || pair->d2_mDCA > 0.5 )
					// && (fabs(pair->d1_mDeltaY) > 20 || fabs(pair->d2_mDeltaY) > 20 )
					// && ( fabs(pair->d1_mDeltaTimeOfFlight) > 1 || fabs( pair->d2_mDeltaTimeOfFlight ) > 1 )
				) {
				
				book->fill("bgshape", lv.M() );
				book->fill("bgshape_rb", lv.M() );
			}


			{
				TH2 * hdMdPidMin = (TH2*)book->get("ulsdMdPidMin");
				float pidMin = std::min( pair->d1_mPid, pair->d2_mPid );
				int pidBin = hdMdPidMin->GetYaxis()->FindBin( pidMin );
				int mBin = hdMdPidMin->GetXaxis()->FindBin( lv.M() );

				for ( int iPidBin = 1; iPidBin <= pidBin; iPidBin++ ){
					float iPid = hdMdPidMin->GetYaxis()->GetBinCenter( iPidBin );
					hdMdPidMin->Fill( lv.M(), iPid );
				}
			}

			if ( pidMax > pidCut && pidMin < bgMax 
				// && ( fabs(pair->d1_mDeltaTimeOfFlight) > 0.75 || fabs( pair->d2_mDeltaTimeOfFlight ) > 0.75 )
				// && (pair->d1_mDCA > 0.5 || pair->d2_mDCA > 0.5 )
				){
				TH2 * hdMdPidMin = (TH2*)book->get("mulsdMdPidMin");

				int pidBin = hdMdPidMin->GetYaxis()->FindBin( pidMin );
				int mBin = hdMdPidMin->GetXaxis()->FindBin( lv.M() );

				for ( int iPidBin = 1; iPidBin <= pidBin; iPidBin++ ){
					float iPid = hdMdPidMin->GetYaxis()->GetBinCenter( iPidBin );
					hdMdPidMin->Fill( lv.M(), iPid );
				}
			}

		} // 0 == pair->mChargeSum


	} //analyzeEvent

	float histMeanY1D( TH1 * h ){
		float avg = 0;
		int navg = 0;
		for ( int bi = 1; bi <= h->GetXaxis()->GetNbins(); bi++ ){
			avg += h->GetBinContent( bi );
			navg ++;
		}
		if ( 0 == navg ) return 0;
		return avg / navg;
	}

	TH1 * derivativeHist( TH1 * h ){
		TH1 * hd = (TH1*)h->Clone( TString::Format( "%s_dx", h->GetName() ) );
		hd->Reset();
		for ( int i = 2; i <= h->GetXaxis()->GetNbins(); i++ ){
			float y1 = h->GetBinContent( i - 1 );
			float y2 = h->GetBinContent( i );

			float x1 = h->GetXaxis()->GetBinCenter( i - 1 );
			float x2 = h->GetXaxis()->GetBinCenter( i );

			hd->SetBinContent( i, (y2-y1) / (x2-x1) );
		}
		return hd;
	}

	TH2 * derivativeHist2D( TH2 * h ){
		TH2 * hd = (TH2*)h->Clone( TString::Format( "%s_dy", h->GetName() ) );
		hd->Reset();

		for ( int x = 1; x <= h->GetXaxis()->GetNbins(); x++ ){
			float t = h->GetBinContent( h->GetBin( x, 1 ) );
			for ( int i = 2; i <= h->GetYaxis()->GetNbins(); i++ ){
				int b1 = h->GetBin( x, i - 1 );
				int b2 = h->GetBin( x, i );
				float y1 = h->GetBinContent( b1 );
				float y2 = h->GetBinContent( b2 );

				float x1 = h->GetXaxis()->GetBinCenter( b1 );
				float x2 = h->GetXaxis()->GetBinCenter( b2 );
				if ( t <= 0 ) continue;
				hd->SetBinContent( b2, (-(y2-y1) / (x2-x1)) / t );
			}	
		}
		return hd;
	}


	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		vector<string> hnames = { "uls_rb", "scand_rb", "pcand_rb", "bgshape_rb", "lscand_rb" };

		for ( string hn : hnames ){
			book->get( hn )->Scale( 1.0, "width" );
		}

		book->get( "bgshape" )->SetLineColor(kRed);
		book->get( "bgshape_rb" )->SetLineColor(kRed);


		book->clone( "pcand_rb", "sig_rb" );
		TH1 * hsig_rb = book->get( "sig_rb" );

		book->clone( "pcand", "sig" );
		TH1 * hsig = book->get( "sig" );

		book->clone( "pcand_rb", "sig_wosc_rb" );
		TH1 * hsig_wosc_rb = book->get( "sig_wosc_rb" );

		book->clone( "pcand", "sig_wosc" );
		TH1 * hsig_wosc = book->get( "sig_wosc" );

		book->clone( "bgshape_rb", "bgs_wosc_rb");
		TH1 * hbgs_wosc_rb = book->get( "bgs_wosc_rb" );
		book->clone( "bgshape", "bgs_wosc");
		TH1 * hbgs_wosc = book->get( "bgs_wosc" );

		book->clone( "bgshape_rb", "bgs_rb");
		TH1 * hbgs_rb = book->get( "bgs_rb" );

		book->clone( "bgshape", "bgs");
		TH1 * hbgs = book->get( "bgs" );

		LOG_F( INFO, "Background Scaler = %f", bgScale );
		
		hbgs_wosc_rb->Scale( bgScale / hbgs->Integral() );
		hbgs_wosc->Scale( bgScale / hbgs->Integral() );
		
		
		


		hsig_wosc_rb->Add( hbgs_rb, -1 );
		hsig_wosc->Add( hbgs, -1 );

		book->clone( "shapecorr_N", "shapecorr" );
		book->clone( "shapecorr_N_rb", "shapecorr_rb" );
		book->get("shapecorr")->Divide( book->get( "shapecorr_D" ) );
		TH1 * hsc_rb = book->get("shapecorr_rb");
		hsc_rb->Divide( book->get( "shapecorr_D_rb" ) );

		// apply the shape correction
		float avgsc = histMeanY1D( hsc_rb );
		// avgsc = avgsc / (float)navgsc;
		LOG_F(INFO, "avg shape corr= %f", avgsc);
		hsc_rb->Scale( 1.0 / avgsc );


		TH2 * hmdMdPidMin = (TH2*)book->get("mulsdMdPidMin");
		TH2 * hmdMdPidMin_dy = derivativeHist2D( hmdMdPidMin );
		TH1 * hmixedcorr = book->get( "mixedcorr" );
		for ( int ib = 1; ib <= hmdMdPidMin_dy->GetXaxis()->GetNbins(); ib++ ){
			TH1 *hp = hmdMdPidMin_dy->ProjectionY( "temp", ib, ib );
			if ( 0 >= hp->Integral() ) continue;
			hmixedcorr->SetBinContent( ib, hp->Integral( 1, 15 ) );
		}


		TH1 * hdMdPidMin = book->get("dMdPidMin");
		TH2 * hulsdMdPidMin = (TH2*)book->get("lsdMdPidMin");
		int pidBin = hulsdMdPidMin->GetXaxis()->FindBin( pidCut );
		for ( int ib = 1; ib <= hdMdPidMin->GetXaxis()->GetNbins(); ib++ ){
			TH1 *hp = hulsdMdPidMin->ProjectionY( "temp", ib, ib );
			if ( 0 >= hp->Integral() ) continue;
			hdMdPidMin->SetBinContent( ib, hp->Integral( pidBin, -1 ) / hp->Integral() );
		}

		derivativeHist( hulsdMdPidMin->ProjectionY( "d_lmr", 1, 10 ) );
		derivativeHist( hulsdMdPidMin->ProjectionY( "d_imr", 13, 16 ) );
		derivativeHist( hulsdMdPidMin->ProjectionY( "d_hmr", 16, -1 ) );
		derivativeHist2D( hulsdMdPidMin );
		

		float avgrc = histMeanY1D( hdMdPidMin );
		hdMdPidMin->Scale( 1.0 / avgrc );

		avgrc = histMeanY1D( hmixedcorr );
		hmixedcorr->Scale( 1.0 / avgrc );


		TF1 * fcorr = new TF1( "fcorr", "pol3", 0, 5 );
		// fcorr->SetParameters( 0.430253, 1.14899, -0.650197, 0.125849 );
		fcorr->SetParameters( 0.639019, 0.827669, -0.507094, 0.0993058 );
		for ( int bi = 1; bi <= hbgs_rb->GetXaxis()->GetNbins(); bi++ ){
			float m = hsig_rb->GetBinCenter( bi );
		
			float corr = 1.0;
			// if ( m < 0.6 ) corr = 0.9;
			// if ( m < 0.5 ) corr = 0.85;
			// if ( m < 0.4 ) corr = 0.8;
			// if ( m > 0.8 ) corr = 0.95;
			// if ( m > 0.9 ) corr = 1.05;
			// if ( m > 1.5 && m < 3.0) corr = 1.12;
			// if ( m > 2.5 && m <= 3.6 ) m = 2.5;
			
			if ( m < 0.75 )
				corr = fcorr->Eval( m );

			if ( m > 1.5 && m < 3.0) corr = 1.08;
			if ( m > 3.5 ) corr = 0.70;
			LOG_F( INFO, "M=%f -> x%f", m, corr );
			hbgs_rb->SetBinContent( bi, hbgs_rb->GetBinContent(bi) * (1.0*corr) );
		}

		for ( int bi = 1; bi <= hbgs->GetXaxis()->GetNbins(); bi++ ){
			float m = hsig->GetBinCenter( bi );
		
			float corr = 1.0;
			// if ( m < 0.6 ) corr = 0.9;
			// if ( m < 0.5 ) corr = 0.85;
			// if ( m < 0.4 ) corr = 0.8;
			// if ( m > 0.8 ) corr = 0.95;
			// if ( m > 0.9 ) corr = 1.05;
			// if ( m > 1.5 && m < 3.0) corr = 1.12;
			// if ( m > 2.5 && m <= 3.6 ) m = 2.5;
			
			if ( m < 0.75 )
				corr = fcorr->Eval( m );

			if ( m > 1.5 && m < 3.0) corr = 1.08;
			// if ( m > 3.5 ) corr = 0.70;
			LOG_F( INFO, "M=%f -> x%f", m, corr );
			hbgs->SetBinContent( bi, hbgs->GetBinContent(bi) * (1.0*corr) );
		}

		// TH1 * hCorr = hsc_rb;
		// for ( int bi = 1; bi <= hbgs_rb->GetXaxis()->GetNbins(); bi++ ){
		// 	float m = hsig_rb->GetBinCenter( bi );
		// 	int cbin = hCorr->GetXaxis()->FindBin( m );
		// 	float corr = hCorr->GetBinContent( cbin );
			
		// 	if ( 0 >= corr || m > 50.5 )
		// 		continue;
		// 	LOG_F( INFO, "M=%f -> x%f", m, corr );
		// 	hbgs_rb->SetBinContent( bi, hbgs_rb->GetBinContent(bi) * (1.0*corr) );
		// }

		TH1 * lscand = book->get("lscand");
		lscand->SetLineColor( kBlack );
		float bgNormMin = config.getFloat( "p.BgNormMin", 2.0 );
		LOG_F( INFO, "Normalizing background in M > %f", bgNormMin );
		int lsBin1 = lscand->GetXaxis()->FindBin( bgNormMin );
		bgScale = bgScale * lscand->Integral( lsBin1, -1 );

		hbgs_rb->Scale( bgScale / hbgs->Integral( lsBin1, -1 ) );
		hbgs->Sumw2();
		hbgs->Scale( bgScale / hbgs->Integral( lsBin1, -1 ) );


		hsig_rb->Add( hbgs_rb, -1 );
		hsig->Add( hbgs, -1 );



		rpl.style( "lscand" ).set( config, "style.lscand" );
		rpl.style( "pcand" ).set( config, "style.pcand" );
		rpl.style( "bgs" ).set( config, "style.bgs" );



	}

};


#endif

