#ifndef MIXED_PID_ANA_H
#define MIXED_PID_ANA_H

#include "TreeAnalyzer.h"

#include "vendor/loguru.h"

#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"
#include "TriggerPatchMapper.h"


#include "TRandom3.h"
#include "TLatex.h"

#include "MixedEventAna.h"


class MixedPidAna : public TreeAnalyzer
{
protected:
	BranchReader<FemtoPair> _fpr;
	size_t eventIndex = 0;
	size_t mixN;
	size_t nMixed;
	size_t nMixedWeighted;



	vector<TrackObj> buffer_pos;
	vector<TrackObj> buffer_neg;

	TRandom3 r3;

	TH1 *huls_mix, *huls_same, *hls_mix, *hls_same;
	TH1 *huls_mix_dpid, *huls_same_dpid, *hls_mix_dpid, *hls_same_dpid;
	TH2 * hmap_hit;
	double Ihmap_hit;
	// int chargeSum;
	vector<string> named_weights;

	string rpName;


	// params
	float minPt, maxPt;
	float MIN_MASS, MAX_MASS;
	map<string, float> norm_const;
	map<string, TH1*> hWeight;

	float Iuls_same, Ils_same;
	float Iuls_mix, Ils_mix;

public:
	MixedPidAna() {}
	~MixedPidAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		r3.SetSeed(0);


		// chargeSum = config.get<int>( "params.chargeSum", 0 );
		
		minPt = config.get<float>( "params.minPt", 0 );
		maxPt = config.get<float>( "params.maxPt", 100 );

		mixN = config.get<unsigned int>( "params.mixN", 5 );

		MIN_MASS = config.get<float>( "params.minMass", 0 );
		MAX_MASS = config.get<float>( "params.maxMass", 0.1 );
		LOG_F( INFO, "MASS = ( %f, %f)", MIN_MASS, MAX_MASS );
		LOG_F( INFO, "mixN=%lu", mixN );

		rpName = string(TString::Format( "rpPidMixed_mass%f_%f.pdf", MIN_MASS, MAX_MASS ).Data());

	}

protected:

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

		// if ( 0 == iEventLoop  ){
		// 	hsame = book->get( "uls_same" );
		// 	hmix = book->get( "uls_mix" );
		// 	hsame = book->get( "ls_same" );
		// 	hmix = book->get( "ls_mix" );
		// }


		if ( iEventLoop <= 1 )
			nMixed = 0;
		if ( iEventLoop <= 2 )
			nMixedWeighted = 0;
	}

	virtual void addToBuffer( TrackObj &to ){

		if ( to.mNHitsFit > 0 )
			buffer_pos.push_back( to );
		else 
			buffer_neg.push_back( to );
	}

	string cprefix( int chargeSum ){
		if ( 0 == chargeSum )
			return "uls_";
		return "ls_";
	}

	virtual void SameEventLoop(  ){

		FemtoPair * pair = this->_fpr.get();

		if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 ) return;
		

		TLorentzVector lv1, lv2, lv;
		lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
		lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
		lv = lv1 + lv2;

		// PAIR kinematics
		if ( fabs( lv.Rapidity() ) > 0.5 )
			return;
		if ( lv.Pt() < minPt || lv.Pt() > maxPt )
			return;

		// DAUGHTER kinematics
		if ( fabs( lv1.PseudoRapidity() ) > 0.5 )
			return;
		if ( lv1.Pt() < 0.1 )
			return;
		if ( fabs( lv2.PseudoRapidity() ) > 0.5 )
			return;
		if ( lv2.Pt() < 0.1 )
			return;

		// TRIGGER FLAG
		if ( pair->d1_mTriggerFlag <= 0 )
			return;
		if ( pair->d2_mTriggerFlag <= 0 )
			return;

		if ( lv.M() > MAX_MASS || lv.M() < MIN_MASS  ) return;

		if ( pair->d1_mPid > 1.0 || pair->d2_mPid > 1.0 ) return;
		if ( pair->d1_mPid < -0.2 || pair->d2_mPid < -0.2 ) return;

		float pairPid = sqrt( pow( pair->d1_mPid, 2 ) + pow( pair->d2_mPid, 2 ) );

		TrackObj to1( eventIndex, lv1, pair, 1 );
		TrackObj to2( eventIndex, lv2, pair, 2 );

		addToBuffer( to1 );
		addToBuffer( to2 );

		book->fill( cprefix( pair->mChargeSum ) + "same", pairPid );
		book->fill( cprefix( pair->mChargeSum ) + "same_deltaPid", abs( pair->d1_mPid - pair->d2_mPid ) );
		book->fill( cprefix( pair->mChargeSum ) + "same_deltaPid_pid", TMath::Min( pair->d1_mPid, pair->d2_mPid ), abs( pair->d1_mPid - pair->d2_mPid ) );

		eventIndex ++;
	}

	virtual void MixedEventLoop( int chargeSum ){
		TLorentzVector lv1, lv2, lv;
			
		vector<TrackObj> *col1 = &buffer_pos;
		vector<TrackObj> *col2 = &buffer_pos;

		if ( chargeSum == 0 ){
			col1 = &buffer_pos;
			col2 = &buffer_neg;
		} else if ( chargeSum == 2 ){
			col1 = &buffer_pos;
			col2 = &buffer_pos;
		} else if ( chargeSum == -2 ){
			col1 = &buffer_neg;
			col2 = &buffer_neg;
		}

		string prefix = cprefix( chargeSum );
		
		for ( size_t i = 0; i < mixN; i++ ){
			size_t it1 = r3.Integer( col1->size() );
			TrackObj &to1 = (*col1)[it1];
			for ( size_t j = 0; j < mixN; j++ ){
				size_t it2 = r3.Integer( col2->size() );
				TrackObj &to2 = (*col2)[it2];
				
				if ( to1.eventIndex == to2.eventIndex )
					continue;

				lv = to1.lv + to2.lv;
				if ( lv.Pt() < minPt || lv.Pt() > maxPt ) continue;
				if ( to1.tp == to2.tp ) continue;

				float pairPid = sqrt( pow( to1.mPid, 2 ) + pow( to2.mPid, 2 ) );
				book->fill( prefix + "mix", pairPid );
				book->fill( prefix + "mix_deltaPid", abs( to1.mPid - to2.mPid ) );
				book->fill( prefix + "mix_deltaPid_pid", TMath::Min( to1.mPid, to2.mPid ), abs( to1.mPid - to2.mPid ) );
				nMixed++;
			} // loop j
		} // loop i
	}

	virtual void WeightedMixedEventLoop( int chargeSum  ){
		TLorentzVector lv1, lv2, lv;
			
		vector<TrackObj> *col1 = &buffer_pos;
		vector<TrackObj> *col2 = &buffer_pos;

		float Isame = Iuls_same;
		float Imix = Iuls_mix;
		string prefix = cprefix( chargeSum );
		TH2 * hwsame = (TH2*)book->get( prefix + "same_deltaPid_pid" );
		TH2 * hwmix = (TH2*)book->get( prefix + "mix_deltaPid_pid" );

		if ( chargeSum == 0 ){
			col1 = &buffer_pos;
			col2 = &buffer_neg;
		} else if ( chargeSum == 2 ){
			col1 = &buffer_pos;
			col2 = &buffer_pos;
			Isame = Ils_same;
			Imix  = Ils_mix;
		} else if ( chargeSum == -2 ){
			col1 = &buffer_neg;
			col2 = &buffer_neg;
			Isame = Ils_same;
			Imix  = Ils_mix;
		}

		for ( size_t i = 0; i < mixN; i++ ){
			size_t it1 = r3.Integer( col1->size() );
			TrackObj &to1 = (*col1)[it1];
			for ( size_t j = 0; j < mixN; j++ ){
				size_t it2 = r3.Integer( col2->size() );
				TrackObj &to2 = (*col2)[it2];
				
				if ( to1.eventIndex == to2.eventIndex )
					continue;

				lv = to1.lv + to2.lv;
				if ( lv.Pt() < minPt || lv.Pt() > maxPt ) continue;
				if ( to1.tp == to2.tp ) continue;

				float pairPid = sqrt( pow( to1.mPid, 2 ) + pow( to2.mPid, 2 ) );
				float dPid = fabs( to1.mPid - to2.mPid );
				float pidmin = TMath::Min( to1.mPid, to2.mPid );

				float wsame = fweight( hwsame, Isame, pidmin, dPid );
				float wmixed = fweight( hwmix, Imix, pidmin, dPid );

				if ( wmixed > 0 )
					wmixed = 1.0 / wmixed;
				else 
					continue;

				book->fill( prefix + "mixw", pairPid, wmixed * wsame );
				book->fill( prefix + "mixw_deltaPid", dPid, wmixed * wsame );
				book->fill( prefix + "mixw_deltaPid_pid", TMath::Min( to1.mPid, to2.mPid ), abs( to1.mPid - to2.mPid ), wmixed * wsame );

				nMixedWeighted++;
			} // loop j
		} // loop i
	}

	virtual void analyzeEvent(){

		if ( 0 == iEventLoop){
			SameEventLoop();
		}

		if ( 1 == iEventLoop ){
			MixedEventLoop(0);
			MixedEventLoop(2);
			MixedEventLoop(-2);
		} // iEventLoop == 1

		if ( 2 == iEventLoop ){
			WeightedMixedEventLoop(0);
			WeightedMixedEventLoop(2);
			WeightedMixedEventLoop(-2);
		} // iEventLoop == 1

	} // analyzeEvent


	float fweight( TH1 * h1, float I, float v ){
		if ( nullptr == h1 )
			return 0;

		int iBin1 = h1->GetXaxis()->FindBin( v );
		float c = h1->GetBinContent( iBin1 );
		if ( c <= 0 )
			return 0;
		return c / I;
	}

	float fweight( TH1 * h1, float I, float hit1, float hit2 ){
		TH2 * hh = dynamic_cast<TH2*>(h1);
		if ( nullptr == hh )
			return 0;

		int iBin1 = hh->GetXaxis()->FindBin( hit1 );
		int iBin2 = hh->GetYaxis()->FindBin( hit2 );
		float v = hh->GetBinContent( iBin1, iBin2 );
		if ( v <= 0 )
			return 0;//5.06842372022301050e-05;
		return v / I;
	}

	// float fweight( TH1 * h1, float I, float hit1, float hit2 ){
	// 	TH2 * hh = dynamic_cast<TH2*>(h1);
	// 	if ( nullptr == hh )
	// 		return 0;

	// 	int iBin1 = hh->GetXaxis()->FindBin( hit1 );
	// 	int iBin2 = hh->GetYaxis()->FindBin( hit2 );
	// 	float v = hh->GetBinContent( iBin1, iBin2 );
	// 	if ( v <= 0 )
	// 		return 0;//5.06842372022301050e-05;
	// 	return v / I;
	// }
	// float fweight( TH1 * h1, float I, float hit1, float hit2, float hit3 ){

	// 	TH3 * hh = dynamic_cast<TH3*>(h1);
	// 	if ( nullptr == hh )
	// 		return 0;

	// 	int iBin1 = hh->GetXaxis()->FindBin( hit1 );
	// 	int iBin2 = hh->GetYaxis()->FindBin( hit2 );
	// 	int iBin3 = hh->GetZaxis()->FindBin( hit3 );
	// 	int iBin = hh->GetBin( iBin1, iBin2, iBin3 );
	// 	// LOG_F( INFO, "bin=(%d, %d, %d)", iBin1, iBin2, iBin3 );
	// 	float v = hh->GetBinContent( iBin );
	// 	if ( v <= 0 )
	// 		return 0;//5.06842372022301050e-05;
	// 	return v / I;
	// }
	// void compareSlice( float pt1, float pt2, TCanvas * can ){

		
	// 	int ptBin1 = hsame->GetYaxis()->FindBin( pt1 );
	// 	int ptBin2 = hsame->GetYaxis()->FindBin( pt2 );


	// 	TH1 * hs = hsame->ProjectionX( TString::Format("hs_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );
	// 	TH1 * hm = hmix->ProjectionX( TString::Format("hm_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );
	// 	TH1 * hmw = hmixw->ProjectionX( TString::Format("hmw_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );


	// 	int hsbin1 = hs->GetXaxis()->FindBin( 2.0 );
	// 	int hsbin2 = hs->GetXaxis()->FindBin( 2.4 );

	// 	int hmbin1 = hm->GetXaxis()->FindBin( 2.0 );
	// 	int hmbin2 = hm->GetXaxis()->FindBin( 2.4 );

	// 	float hmbw = hm->GetBinWidth( 10 );
	// 	float hsbw = hs->GetBinWidth( 10 );

	// 	float Ihs = hs->Integral( hsbin1, hsbin2 );
	// 	float Ihm = hm->Integral( hmbin1, hmbin2 );
		

	// 	hm->Sumw2();
	// 	hmw->Sumw2();


	// 	hm->Scale( Ihs / Ihm );
	// 	hm->Scale( hsbw / hmbw  );


	// 	hmbin1 = hm->GetXaxis()->FindBin( 0.3 );
	// 	hmbin2 = hm->GetXaxis()->FindBin( 4.0 );

	// 	hsbin1 = hs->GetXaxis()->FindBin( 0.3 );
	// 	hsbin2 = hs->GetXaxis()->FindBin( 4.0 );

	// 	Ihs = hs->Integral( hsbin1, hsbin2 );
	// 	float Ihmw = hmw->Integral( hmbin1, hmbin2 );
		
	// 	hmw->Scale( Ihs / Ihmw );
	// 	hmw->Scale( hsbw / hmbw  );
		

	// 	hs->SetTitle( TString::Format( "%0.2f < p_{T}^{#mu#mu} < %0.2f (GeV/c); M_{#mu#mu} (GeV/c^{2}); dN/dM_{#mu#mu} (GeV/c^{2})^{-1}", pt1, pt2 ) );
		
	// 	hs->Draw("pe");
	// 	hm->Draw("same hpe");
	// 	hmw->Draw("same hpe");
		
	// 	gPad->SetLogy(0);
	// 	can->Print( rpName.c_str() );


	// 	HistoBins ratio_bins( config, "bins.mass_ratio" );
	// 	// Compute the ratios
	// 	TH1 * hsratio = (TH1*) hs->Rebin( ratio_bins.nBins(), TString::Format("hsr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() );

	// 	TH1 * hmr = (TH1*)hm->Rebin( ratio_bins.nBins(), TString::Format("hmr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() ); 

	// 	for ( int i = 1; i < hmr->GetXaxis()->GetNbins()+1; i++ ){
	// 		float v1 = hmr->GetBinContent( i );
	// 		float v2 = hsratio->GetBinContent( i );
	// 		float e2 = hsratio->GetBinError( i );
	// 		if ( 0 == v2 ) continue;
	// 		hmr->SetBinContent( i, (v1/v2) * (hmbw / hsbw) );
	// 		hmr->SetBinError( i, (e2 / v2) * (hmbw / hsbw) );
	// 	}

	// 	TH1 * hmwr = (TH1*)hmw->Rebin( ratio_bins.nBins(), TString::Format("hmwr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() ); 
	// 	// TH1 * hmwr = (TH1*)hmw->Clone( TString::Format("hmwr_%d_%d", ptBin1, ptBin2 ) ); 
		
	// 	for ( int i = 1; i < hmwr->GetXaxis()->GetNbins()+1; i++ ){
	// 		float v1 = hmwr->GetBinContent( i );
	// 		float v2 = hsratio->GetBinContent( i );
	// 		float e2 = hsratio->GetBinError( i );
	// 		if ( 0 == v2 || 0 == v1  ) continue;
	// 		hmwr->SetBinContent( i, (v1 / v2) * (hmbw / hsbw) );
	// 		hmwr->SetBinError( i, (e2 / v2) * (hmbw / hsbw) );
	// 	}

	// 	hmr->Draw( "" );
		
	// 	TF1 * fp0 = new TF1( "fp0", "pol0" );
	// 	hmr->Fit( fp0, "R", "", 2.0, 5.0 );
	// 	hmr->SetMinimum( 0.5 );
	// 	hmr->SetMaximum( 1.5 );
	// 	gPad->SetLogy(0);
		
	// 	TLatex tl;
	// 	tl.SetTextSize( 12.0 / 360.0 );
	// 	tl.DrawLatexNDC( 0.2, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", fp0->GetChisquare(), fp0->GetNDF(), fp0->GetChisquare() / (float)fp0->GetNDF() ) );
	// 	tl.DrawLatexNDC( 0.2, 0.8, TString::Format("decorr = %0.3f", fp0->GetParameter( 0 )) );
		
		
	// 	TF1 * fp1 = new TF1( "fp1", "pol0" );
	// 	fp1->SetRange( 0.2, 3.0 );
	// 	hmwr->Fit( fp1, "NR", "", 0.2, 3.0 );
	// 	hmwr->Draw("same");
	// 	fp1->SetLineColor(kBlack);
	// 	fp1->SetRange( 0, 4.0 );
	// 	fp1->Draw("same");
	// 	tl.DrawLatexNDC( 0.2, 0.75, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", fp1->GetChisquare(), fp1->GetNDF(), fp1->GetChisquare() / (float)fp1->GetNDF() ) );
	// 	tl.DrawLatexNDC( 0.2, 0.7, TString::Format("corr = %0.3f", fp1->GetParameter( 0 )) );
		
	// 	can->Print( rpName.c_str() );
	// }
	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		LOG_F( INFO, "# of Mixed = %lu", nMixed );
		LOG_F( INFO, "iEventLoop=%d", iEventLoop );

		if ( iEventLoop == 0 ){
			LOG_F( INFO, "buffer_pos.size()=%lu", buffer_pos.size() );
			LOG_F( INFO, "buffer_neg.size()=%lu", buffer_neg.size() );

			Iuls_same = book->get( "uls_same_deltaPid_pid" )->Integral();
			Ils_same = book->get( "ls_same_deltaPid_pid" )->Integral();
		}

		if ( iEventLoop == 1 ){
			Iuls_mix = book->get( "uls_mix_deltaPid_pid" )->Integral();
			Ils_mix = book->get( "ls_mix_deltaPid_pid" )->Integral();
		}
	}

};


#endif