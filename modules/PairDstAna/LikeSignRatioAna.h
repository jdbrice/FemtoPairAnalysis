#ifndef LIKE_SIGN_RATIO_METHOD_H
#define LIKE_SIGN_RATIO_METHOD_H

#include "TPaveStats.h"
#include "TRandom3.h"

#include "TLatex.h"

class LikeSignRatioAna : public TreeAnalyzer
{
protected:

	BranchReader<FemtoPair> _fpr;
	float pidCut = 1.3;
	float minPairPt = 0;
	float maxPairPt = 100;

	bool use_cum = false;

	vector<double> pidVals;


	vector<TLorentzVector> buffer_pos_lv, buffer_neg_lv;
	vector<float> buffer_pos_pid, buffer_neg_pid;

public:
	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		book->cd();

		pidCut = config.get<float>( "p.pidCut", 1.3 );
		minPairPt = config.get<float>("p.minPairPt", 2.2 );
		maxPairPt = config.get<float>("p.maxPairPt", 100 );

		use_cum = config.getBool( "p.cum", use_cum );
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
		if ( lv.Pt() < minPairPt || lv.Pt() > maxPairPt )
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

		float pidMin  = std::min( pair->d1_mPid, pair->d2_mPid );
		float pidMax  = std::max( pair->d1_mPid, pair->d2_mPid );
		float pairPid = sqrt( pow(pair->d1_mPid, 2) + pow(pair->d2_mPid, 2) );

		// pidVals.push_back( pairPid);


		////////////////////////////////////////////////////////////////////////
		/// Like-Sign
		////////////////////////////////////////////////////////////////////////
		// if ( pairPid < 1.0 )
		// 	return;

		buffer_pos_lv.push_back( lv1 );
		buffer_pos_pid.push_back( pair->d1_mPid );
		buffer_neg_lv.push_back( lv2 );
		buffer_neg_pid.push_back( pair->d2_mPid );

		// if ( pidMin < 0.2 )
		// 	return;

		if ( 0 != pair->mChargeSum ){
			book->fill( "all_ls", lv.M() );
		} else if ( 0 == pair->mChargeSum ){
			book->fill( "all_uls", lv.M() );
		}


		if ( pairPid > pidCut ){
			if ( 0 != pair->mChargeSum ){
				book->fill( "ls", lv.M() );
			} else if ( 0 == pair->mChargeSum ){
				book->fill( "uls", lv.M() );
			}
		}  


		if ( pairPid < pidCut || false == use_cum ){
			if ( 0 != pair->mChargeSum ){
				book->fill( "ls_pid_mass", lv.M(), pairPid );
			} else if ( 0 == pair->mChargeSum ){
				book->fill( "uls_pid_mass", lv.M(), pairPid );
			}
		}
	}


	void spliceInX( TH2 * h2, TH1 * h, int ix ){
		TAxis * x = h2->GetXaxis();
		TAxis * y = h2->GetYaxis();

		TAxis * hx = h->GetXaxis();

		if ( y->GetNbins() != hx->GetNbins() ){
			LOG_F( ERROR, "Number of bins mismatch y=%d != hx=%d", y->GetNbins(), hx->GetNbins() );
		}

		for ( size_t iy = 1; iy <= y->GetNbins(); iy++ ){
			h2->SetBinContent( ix, iy, h->GetBinContent( iy ) );
		}
	}


	void draw_signal( Reporter &rp ){

		RooPlotLib rpl;
		rpl.link( book );

		TH1 * huls = book->get( "uls" );
		TH1 * hls = book->get( "ls" );
		TH1 * hbg = book->get( "bg" );

		if ( hbg->GetXaxis()->GetBinWidth(4) != hbg->GetXaxis()->GetBinWidth(12) ){
			// huls->Scale( 1.0, "width" );
			// hls->Scale( 1.0, "width" );
			// hbg->Scale( 1.0, "width" );
		}

		rpl.style( "uls" ).set( config, "style.uls" ).set("draw", "hp" ).set( "logy", 1 ).set( "min", 1e1 ).draw();
		rpl.style( "ls" ).set( config, "style.ls" ).set( "draw", "same hp" ).draw();
		rpl.style( "bg" ).set( config, "style.bg" ).set( "draw", "same hp" ).draw();

		book->clone( "uls", "rsig" );

		book->get( "rsig" )->Add( hbg, -1 );

		rpl.style( "rsig" ).set( "lc", "#222" ).set( "lw", 2 ).set( "draw", "same hp" ).draw();

		rp.savePage();
	}

	TH1 * rolling_ratio( string name, TH1 * hnum, TH1 * hden, short k1 = 0, short k2 = 0 ){
		TH1 *hratio = (TH1*)hnum->Clone( name.c_str() );

		if ( hnum->GetXaxis()->GetNbins() != hden->GetXaxis()->GetNbins() ){
			LOG_F( ERROR, "bin mismatch " );
			return hratio;
		}

		int n_bins = hratio->GetXaxis()->GetNbins();
		for ( int ix = 1; ix <= n_bins; ix++ ){

			int ix1 = ix - k1;
			int ix2 = ix + k2;

			if ( ix1 < 1 )
				ix1 = 1;
			if ( ix2 > n_bins )
				ix2 = n_bins;

			float num = hnum->GetBinContent( ix );
			float den = hden->GetBinContent( ix );
			
			// for ( int ixx = ix1; ixx <= ix2; ixx++ ){
			// 	num += hnum->GetBinContent( ixx );
			// 	den += hden->GetBinContent( ixx );
			// }

			den = hden->Integral();
			if ( 0 == den ){
				hratio->SetBinContent( ix, 0 );
				continue;
			}

			
			hratio->SetBinContent( ix, num / den );
			hratio->SetBinError( ix, 0 );
		}


		return hratio;
	}


	void processSlices(  ){

		RooPlotLib rpl;
		rpl.link( book );
		gStyle->SetOptStat( 0 );
		gStyle->SetOptFit( 111 );


		Reporter rp( "efits.pdf", 300, 300 );


		TH2 * hls2  = (TH2*)book->get( "ls_pid_mass" );
		TH2 * huls2 = (TH2*)book->get( "uls_pid_mass" );

		book->clone( "ls", "bg" );
		
		TH1 *huls = book->get( "uls" );
		TH1 *hls = book->get( "ls" );
		TH1 *hbg = book->get( "bg" );

		TH1 * hf = (TH1*)book->get( "cfactor" )->Clone( "tmpfactor" );
		hf->Reset();

		TAxis * x = hls2->GetXaxis();
		TAxis * y = hls2->GetYaxis();

		for ( size_t ix = 1; ix <= hf->GetXaxis()->GetNbins(); ix++ ){
			size_t a = ix;
			size_t b = ix;
			// if ( b >= x->GetNbins() ){
			// 	a = ix - 4;
			// 	b = ix;
			// }
			TH1 * hpuls = huls2->ProjectionY( TString::Format("huls_%lu", ix), a, b );
			TH1 * hpls = hls2->ProjectionY( TString::Format("hls_%lu", ix) , a, b );

			float mlow = x->GetBinLowEdge(ix);
			float mhigh = x->GetBinUpEdge(ix);

			
			if ( use_cum ){
				hpuls = hpuls->GetCumulative(true);
				hpls = hpls->GetCumulative(true);
			}
			

			hpuls->Sumw2();
			hpls->Sumw2();

			TH1 * hr = (TH1*)hpuls->Clone( TString::Format("hratio_%lu", ix) );
			// hr->SetDirectory(0);

			hr->Divide( hpls );
			// 
			// TH1 *hr = rolling_ratio( "tmpratio", hpuls, hpls );

			hr->SetTitle( TString::Format( "%f<m<%f", mlow, mhigh ) );
			hpuls->SetTitle( TString::Format( "%f<m<%f", mlow, mhigh ) );
			hpls->SetTitle( TString::Format( "%f<m<%f", mlow, mhigh ) );

			rpl.style( hpuls ).set(config, "style.hpuls" );
			rpl.style( hpls ).set(config, "style.hpls" );

			hpuls->Draw( "hpe" );
			hpuls->SetMinimum(1);
			hpls->Draw("same hpe");
			gPad->SetLogy(1);

			rp.savePage();


			rpl.style( hr ).set( config, "style.ratio" ).draw();
			gPad->SetLogy(0);

			TF1 * fpol0 = new TF1( "fpol0", config.getString( "p.fit:formula", "pol0" ).c_str() );
			hr->Fit( fpol0, "R", "", config.get<float>("p.fit:min", 1.0), config.get<float>("p.fit:max", 1.2) );

			// hr->SetMinimum( fpol0->Eval( 1.4 ) - 1.0 );
			hr->SetMinimum( 0.0 );
			hr->SetMaximum( fpol0->Eval( 1.4 ) + 3.0 );

			hbg->SetBinContent( ix, hbg->GetBinContent( ix ) * fpol0->Eval(1.4) );
			hf->SetBinContent( ix, fpol0->Eval( 1.4 ) );

			
			TLatex tl;
			tl.SetTextSize( 12.0 / 360.0 );
			tl.DrawLatexNDC( 0.4, 0.40, TString::Format("r(pairPID=1.3) = %0.4f ", fpol0->Eval( 1.3 ) ) );
			tl.DrawLatexNDC( 0.4, 0.35, TString::Format("r(pairPID=1.4) = %0.4f ", fpol0->Eval( 1.4 ) ) );

			// fpol0->Draw("")

			gPad->Update();
			TPaveStats *st = (TPaveStats*)hr->FindObject("stats");
			if (st ){
				st->SetX1NDC( 0.25 ); st->SetX2NDC( 0.5 );
				st->SetY1NDC( 0.8 ); st->SetY2NDC( 0.9 );
			}
			
			gPad->Modified(); 

			rp.savePage();
		}

		hf->Draw();
		rp.savePage();

		draw_signal( rp );
	}

	virtual void postMake(){
		TreeAnalyzer::postMake();

		RooPlotLib rpl;
		rpl.link( book );

		processSlices();
		// mix_like_sign();
		return;

		book->clone( "ls", "bg" );
		book->clone( "uls_pid_mass", "ratio_pid_mass" );
		TH2 * h2 = (TH2*)book->get( "ratio_pid_mass" );
		h2->Divide( (TH2*)book->get( "ls_pid_mass" ) );


		TH2 * hls  = (TH2*)book->get( "ls_pid_mass" );
		TH2 * huls = (TH2*)book->get( "uls_pid_mass" );
		TH2 * hrC  = (TH2*)book->get( "ratio_cum" );
		TH2 * hr   = (TH2*)book->get( "ratio" );

		TAxis * x = h2->GetXaxis();
		TAxis * y = h2->GetYaxis();
		for ( size_t ix = 1; ix <= x->GetNbins(); ix++ ){
			TH1 * hpuls = huls->ProjectionY( "tmpuls", ix, ix );
			TH1 * hpulsC = hpuls->GetCumulative( true );
			hpuls->Sumw2();
			spliceInX( (TH2*)book->get( "uls_cum" ), hpulsC, ix );


			TH1 * hpls = hls->ProjectionY( "tmpls", ix, ix );
			TH1 * hplsC = hpls->GetCumulative( true );
			hpls->Sumw2();
			spliceInX( (TH2*)book->get( "ls_cum" ), hplsC, ix );

			hpulsC->Divide( hplsC );
			hpuls->Divide( hpls );

			spliceInX( hrC, hpulsC, ix );
			spliceInX( hr, hpuls, ix );
		}



		Reporter rp( "fits.pdf", 300, 300 );
		
		TH1 * hbg = book->get( "bg" );

		for ( size_t ix = 1; ix <= x->GetNbins(); ix++ ){
			TH1 * hp = hr->ProjectionY( "py", ix, ix );
			hp->Draw( "hp" );
			TF1 * fpol0 = new TF1( "fpol0", "pol0" );
			hp->Fit( fpol0, "R", "", config.get<float>("p.fmin", 1.0), config.get<float>("p.fmax", 1.2) );

			hbg->SetBinContent( ix, hbg->GetBinContent( ix ) * fpol0->GetParameter(0) );
			gPad->Update();
			TPaveStats *st = (TPaveStats*)hp->FindObject("stats");
			if (st ){
				st->SetX1NDC( 0.25 ); st->SetX2NDC( 0.5 );
				st->SetY1NDC( 0.8 ); st->SetY2NDC( 0.9 );
			}
			
			gPad->Modified(); 

			rp.savePage();
		}

		book->get( "uls" )->Sumw2();
		book->get( "ls" )->Sumw2();
		book->get( "bg" )->Sumw2();

		if ( hbg->GetXaxis()->GetBinWidth(4) != hbg->GetXaxis()->GetBinWidth(12) ){
			book->get("uls")->Scale( 1.0, "width" );
			book->get("ls")->Scale( 1.0, "width" );
			book->get("bg")->Scale( 1.0, "width" );
		}

		rpl.style( "uls" ).set( config, "style.uls" ).set("draw", "hp" ).set( "logy", 1 ).set( "min", 10 ).draw();
		rpl.style( "ls" ).set( config, "style.ls" ).set( "draw", "same hp" ).draw();
		rpl.style( "bg" ).set( config, "style.bg" ).set( "draw", "same hp" ).draw();

		book->clone( "uls", "rsig" );

		book->get( "rsig" )->Add( hbg, -1 );

		rpl.style( "rsig" ).set( "draw", "same hp" ).draw();



		rp.savePage();



		vector<double> bin_edges = HistoBins::makeQuantileBins( pidVals, 100, -1, -1 );
		cout << vts( bin_edges ) << endl;

	}


	virtual void mix_like_sign(){
		return;
		size_t n_pos = buffer_pos_pid.size();
		size_t n_neg = buffer_neg_pid.size();
		LOG_F( INFO, "buffer_neg_pid.size() = %lu", n_pos );
		LOG_F( INFO, "buffer_pos_pid.size() = %lu", n_neg );

		TRandom3 r3;
		r3.SetSeed( 100  );

		TH1 * h = book->get( "mixed_ls" );

		for ( size_t i = 0; i < 1000000; i++ ){

			size_t i_pos = r3.Integer( n_pos );
			size_t i_neg = r3.Integer( n_neg );

			float p1 = buffer_pos_pid[ i_pos ];
			float p2 = buffer_neg_pid[ i_neg ];
			TLorentzVector lv = buffer_pos_lv[ i_pos ] + buffer_neg_lv[ i_neg ];

			float pairPid = sqrt( pow( p1,2 ) + pow( p2, 2 ) );
			h->Fill( lv.M(), pairPid );
		}
	}

};


#endif


