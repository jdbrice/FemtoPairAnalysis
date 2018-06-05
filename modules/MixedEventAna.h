#ifndef MIXED_EVENT_ANA_H
#define MIXED_EVENT_ANA_H

#include "TreeAnalyzer.h"

#include "vendor/loguru.h"

#include "FemtoDstFormat/BranchReader.h"
#include "FemtoPairFormat/FemtoPair.h"
#include "TriggerPatchMapper.h"


#include "TRandom3.h"
#include "TLatex.h"

class TrackObj{
public:
	TLorentzVector lv;
	size_t eventIndex;
	float pid;
	int hit;
	int tp;

	float mDeltaY;
	float mDeltaZ;
	float mDeltaTimeOfFlight;
	int mMatchFlag;
	int mTriggerFlag;
	int mCell;
	int mModule;
	int mBackleg;
	float mPt;
	float mEta;
	float mPhi;
	int mId;
	int mNHitsFit;
	int mNHitsMax;
	int mNHitsDedx;
	float mNSigmaPion;
	float mDCA;
	float mPid;

	TrackObj( size_t _eventIndex, TLorentzVector &_lv, FemtoPair * fp, int index ){
		lv = _lv;
		this->eventIndex = _eventIndex;

		if ( 1 == index ){
			mDeltaY = fp->d1_mDeltaY;
			mDeltaZ = fp->d1_mDeltaZ;
			mDeltaTimeOfFlight = fp->d1_mDeltaTimeOfFlight;
			mMatchFlag = fp->d1_mMatchFlag;
			mTriggerFlag = fp->d1_mTriggerFlag;
			mCell = fp->d1_mCell;
			mModule = fp->d1_mModule;
			mBackleg = fp->d1_mBackleg;
			mPt = fp->d1_mPt;
			mEta = fp->d1_mEta;
			mPhi = fp->d1_mPhi;
			mId = fp->d1_mId;
			mNHitsFit = fp->d1_mNHitsFit;
			mNHitsMax = fp->d1_mNHitsMax;
			mNHitsDedx = fp->d1_mNHitsDedx;
			mNSigmaPion = fp->d1_mNSigmaPion;
			mDCA = fp->d1_mDCA;
			mPid = fp->d1_mPid;
		} else {
			mDeltaY = fp->d2_mDeltaY;
			mDeltaZ = fp->d2_mDeltaZ;
			mDeltaTimeOfFlight = fp->d2_mDeltaTimeOfFlight;
			mMatchFlag = fp->d2_mMatchFlag;
			mTriggerFlag = fp->d2_mTriggerFlag;
			mCell = fp->d2_mCell;
			mModule = fp->d2_mModule;
			mBackleg = fp->d2_mBackleg;
			mPt = fp->d2_mPt;
			mEta = fp->d2_mEta;
			mPhi = fp->d2_mPhi;
			mId = fp->d2_mId;
			mNHitsFit = fp->d2_mNHitsFit;
			mNHitsMax = fp->d2_mNHitsMax;
			mNHitsDedx = fp->d2_mNHitsDedx;
			mNSigmaPion = fp->d2_mNSigmaPion;
			mDCA = fp->d2_mDCA;
			mPid = fp->d2_mPid;
		}

		this->hit = mBackleg * 60 + mModule * 12 + mCell;
		this->tp = TriggerPatchMapper::findTriggerPatch( hit );
	}
};

class MixedEventAna : public TreeAnalyzer
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

	TH2 *hmix, *hsame, *hmixw;
	TH2 * hmap_hit;
	double Ihmap_hit;
	int chargeSum;
	vector<string> named_weights;

	string rpName;


	// params
	float minPt, maxPt;
	float MIN_PID, MAX_PID, CUT_PID;
	bool usePairPID = true;
	map<string, float> norm_const;
	map<string, TH1*> hWeight;

public:
	MixedEventAna() {}
	~MixedEventAna() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		this->_fpr.setup( this->chain, "Pairs" );
		
		book->cd();

		r3.SetSeed(0);


		chargeSum = config.get<int>( "params.chargeSum", 0 );
		
		minPt = config.get<float>( "params.minPt", 0 );
		maxPt = config.get<float>( "params.maxPt", 100 );

		named_weights = config.getStringVector( "params.weights" );
		mixN = config.get<unsigned int>( "params.mixN", 5 );

		MIN_PID = config.get<float>( "params.minPID", 0 );
		MAX_PID = config.get<float>( "params.maxPID", 0.1 );
		CUT_PID = config.get<float>( "params.cutPID", 0.1 );
		usePairPID = config.getBool( "params.usePairPID", true );
		LOG_F( INFO, "PID = ( %f, %f)", MIN_PID, MAX_PID );
		LOG_F( INFO, "mixN=%lu", mixN );

		rpName = string(TString::Format( "rpMixed_cs%d_pid%f_%f.pdf", chargeSum, MIN_PID, MAX_PID ).Data());

	}

protected:

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

		if ( 0 == iEventLoop  ){
			hsame = (TH2*)book->get( "hsame" );
			hmix = (TH2*)book->get( "hmix" );
			hmixw = (TH2*)book->get( "hmixw" );

			hmap_hit = (TH2*)book->get( "map_hit" );

			for ( string n : named_weights ){
				if ( book->exists( n ) )
					hWeight[n] = (TH2*)book->get( n );
			}
		}


		if ( iEventLoop <= 1 )
			nMixed = 0;
		if ( iEventLoop <= 2 )
			nMixedWeighted = 0;
	}

	int deltaBL( int bl1, int bl2 ){
		if ( bl1 > bl2 )
			return bl1 - bl2;
		return bl2 - bl1;
	}

	virtual void addToBuffer( TrackObj &to ){

		if ( to.mNHitsFit > 0 )
			buffer_pos.push_back( to );
		else 
			buffer_neg.push_back( to );
	}

	double meanContent( TH1 * h1 ){

		TH2 * hh = dynamic_cast<TH2*>(h1);
		if ( nullptr == hh )
			return 0;
		double total = 0.0;
		double n_bins = 0;
		for ( int i = 1; i < hh->GetXaxis()->GetNbins()+1; i++ ){
			for ( int j = 1; j < hh->GetYaxis()->GetNbins()+1; j++ ){
				double v = hh->GetBinContent( i, j );
				if ( v <= 1 ) continue;
				total += v;
				n_bins += 1;
			}
		}
		if ( n_bins < 0.99 ){
			return 0;
		}
		return total / n_bins;
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
		

		float pairPid = sqrt( pow( pair->d1_mPid, 2 ) + pow( pair->d2_mPid, 2 ) );
		float minPid = TMath::Min( pair->d1_mPid, pair->d2_mPid );
		float maxPid = TMath::Max( pair->d1_mPid, pair->d2_mPid );

		if ( usePairPID ){
			if ( pairPid > MAX_PID || pairPid < MIN_PID  ) return;
		} else {
			if ( minPid > MAX_PID || maxPid < MIN_PID ) return;
			if ( minPid < CUT_PID ) return;
		}
		

		if ( abs(pair->mChargeSum) == chargeSum ){
			TrackObj to1( eventIndex, lv1, pair, 1 );
			TrackObj to2( eventIndex, lv2, pair, 2 );

			addToBuffer( to1 );
			addToBuffer( to2 );
			
			if ( lv.Pt() > minPt && lv.Pt() < maxPt && lv.M() > 0.0){
				book->fill( "hsame", lv.M(), lv.Pt() );
				

				float dEta = fabs( lv1.PseudoRapidity() - lv2.PseudoRapidity() );
				float dPhi = fabs(lv1.DeltaPhi( lv2 ));
				float angle = to1.lv.Angle( to2.lv.Vect() );
				((TH2*)hWeight[ "same_pt_angle" ])->Fill( angle, lv.Pt() );
				((TH2*)hWeight[ "same_dEta_dPhi" ])->Fill( dPhi, dEta );
			}
			
		} // if chargeSum

		eventIndex ++;
	}

	virtual void MixedEventLoop(){
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
				
				float dPhi = fabs(to1.lv.DeltaPhi( to2.lv ));
				float dEta = fabs( to1.lv.PseudoRapidity() - to2.lv.PseudoRapidity() );
				float angle = to1.lv.Angle( to2.lv.Vect() );

				((TH2*)hWeight[ "mix_pt_angle" ])->Fill( angle, lv.Pt() );
				((TH2*)hWeight[ "mix_dEta_dPhi" ])->Fill( dPhi, dEta );

				hmix->Fill( lv.M(), lv.Pt() );
				nMixed++;
			} // loop j
		} // loop i
	}

	virtual void WeightedMixedEventLoop(){
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
				
				float dPhi = fabs(to1.lv.DeltaPhi( to2.lv ));
				float dEta = fabs( to1.lv.PseudoRapidity() - to2.lv.PseudoRapidity() );
				float angle = to1.lv.Angle( to2.lv.Vect() );

				float wmixed = fweight( hWeight[ "mix_pt_angle" ], norm_const["mix_pt_angle"], angle, lv.Pt() );
				float wsame = fweight( hWeight[ "same_pt_angle" ], norm_const["same_pt_angle"], angle, lv.Pt() );

				if ( wmixed > 0 )
					wmixed = 1.0 / wmixed;
				else 
					continue;

				((TH2*)hWeight[ "mixwa_pt_angle" ])->Fill( angle, lv.Pt(), wmixed );
				((TH2*)hWeight[ "mixwb_pt_angle" ])->Fill( angle, lv.Pt(), wmixed * wsame );

				((TH2*)hWeight[ "mixw_dEta_dPhi" ])->Fill( dPhi, dEta, wmixed );

				hmixw->Fill( lv.M(), lv.Pt(), wmixed * wsame );
				nMixedWeighted++;
			} // loop j
		} // loop i
	}

	virtual void analyzeEvent(){

		if ( 0 == iEventLoop){
			SameEventLoop();
		}

		if ( 1 == iEventLoop ){
			MixedEventLoop();
		} // iEventLoop == 1

		if ( 2 == iEventLoop ){
			WeightedMixedEventLoop();
		} // iEventLoop == 2

	} // analyzeEvent

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
	float fweight( TH1 * h1, float I, float hit1, float hit2, float hit3 ){

		TH3 * hh = dynamic_cast<TH3*>(h1);
		if ( nullptr == hh )
			return 0;

		int iBin1 = hh->GetXaxis()->FindBin( hit1 );
		int iBin2 = hh->GetYaxis()->FindBin( hit2 );
		int iBin3 = hh->GetZaxis()->FindBin( hit3 );
		int iBin = hh->GetBin( iBin1, iBin2, iBin3 );
		// LOG_F( INFO, "bin=(%d, %d, %d)", iBin1, iBin2, iBin3 );
		float v = hh->GetBinContent( iBin );
		if ( v <= 0 )
			return 0;//5.06842372022301050e-05;
		return v / I;
	}
	void compareSlice( float pt1, float pt2, TCanvas * can ){

		
		int ptBin1 = hsame->GetYaxis()->FindBin( pt1 );
		int ptBin2 = hsame->GetYaxis()->FindBin( pt2 );


		TH1 * hs = hsame->ProjectionX( TString::Format("hs_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );
		TH1 * hm = hmix->ProjectionX( TString::Format("hm_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );
		TH1 * hmw = hmixw->ProjectionX( TString::Format("hmw_%d_%d", ptBin1, ptBin2 ), ptBin1, ptBin2 );


		int hsbin1 = hs->GetXaxis()->FindBin( config.get<float>( "params.norm:min", 2.0 ) );
		int hsbin2 = hs->GetXaxis()->FindBin( config.get<float>( "params.norm:max", 5.0 ) );

		int hmbin1 = hm->GetXaxis()->FindBin( config.get<float>( "params.norm:min", 2.0 ) );
		int hmbin2 = hm->GetXaxis()->FindBin( config.get<float>( "params.norm:max", 5.0 ) );

		float hmbw = hm->GetBinWidth( 10 );
		float hsbw = hs->GetBinWidth( 10 );

		float Ihs = hs->Integral( hsbin1, hsbin2 );
		float Ihm = hm->Integral( hmbin1, hmbin2 );
		

		hm->Sumw2();
		hmw->Sumw2();


		hm->Scale( Ihs / Ihm );
		hm->Scale( hsbw / hmbw  );


		hmbin1 = hm->GetXaxis()->FindBin( config.get<float>( "params.weightedNorm:min", 0.2 ) );
		hmbin2 = hm->GetXaxis()->FindBin( config.get<float>( "params.weightedNorm:max", 4.0 ) );

		hsbin1 = hs->GetXaxis()->FindBin( config.get<float>( "params.weightedNorm:min", 0.2 ) );
		hsbin2 = hs->GetXaxis()->FindBin( config.get<float>( "params.weightedNorm:max", 4.0 ) );

		Ihs = hs->Integral( hsbin1, hsbin2 );
		float Ihmw = hmw->Integral( hmbin1, hmbin2 );
		
		hmw->Scale( Ihs / Ihmw );
		hmw->Scale( hsbw / hmbw  );
		

		hs->SetTitle( TString::Format( "%0.2f < p_{T}^{#mu#mu} < %0.2f (GeV/c); M_{#mu#mu} (GeV/c^{2}); dN/dM_{#mu#mu} (GeV/c^{2})^{-1}", pt1, pt2 ) );
		
		hs->Draw("pe");
		hm->Draw("same hpe");
		hmw->Draw("same hpe");
		
		gPad->SetLogy(1);
		can->Print( rpName.c_str() );


		HistoBins ratio_bins( config, "bins.mass_ratio" );
		// Compute the ratios
		TH1 * hsratio = (TH1*) hs->Rebin( ratio_bins.nBins(), TString::Format("hsr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() );

		TH1 * hmr = (TH1*)hm->Rebin( ratio_bins.nBins(), TString::Format("hmr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() ); 

		for ( int i = 1; i < hmr->GetXaxis()->GetNbins()+1; i++ ){
			float v1 = hmr->GetBinContent( i );
			float v2 = hsratio->GetBinContent( i );
			float e2 = hsratio->GetBinError( i );
			if ( 0 == v2 ) continue;
			hmr->SetBinContent( i, (v1/v2) * (hmbw / hsbw) );
			hmr->SetBinError( i, (e2 / v2) * (hmbw / hsbw) );
		}

		TH1 * hmwr = (TH1*)hmw->Rebin( ratio_bins.nBins(), TString::Format("hmwr_%d_%d", ptBin1, ptBin2 ), ratio_bins.bins.data() ); 
		// TH1 * hmwr = (TH1*)hmw->Clone( TString::Format("hmwr_%d_%d", ptBin1, ptBin2 ) ); 
		
		for ( int i = 1; i < hmwr->GetXaxis()->GetNbins()+1; i++ ){
			float v1 = hmwr->GetBinContent( i );
			float v2 = hsratio->GetBinContent( i );
			float e2 = hsratio->GetBinError( i );
			if ( 0 == v2 || 0 == v1  ) continue;
			hmwr->SetBinContent( i, (v1 / v2) * (hmbw / hsbw) );
			hmwr->SetBinError( i, (e2 / v2) * (hmbw / hsbw) );
		}

		hmr->Draw( "" );
		
		TF1 * fp0 = new TF1( "fp0", "pol0" );
		hmr->Fit( fp0, "R", "", 2.0, 5.0 );
		hmr->SetMinimum( 0.5 );
		hmr->SetMaximum( 1.5 );
		gPad->SetLogy(0);
		
		TLatex tl;
		tl.SetTextSize( 12.0 / 360.0 );
		tl.DrawLatexNDC( 0.2, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", fp0->GetChisquare(), fp0->GetNDF(), fp0->GetChisquare() / (float)fp0->GetNDF() ) );
		tl.DrawLatexNDC( 0.2, 0.8, TString::Format("decorr = %0.3f", fp0->GetParameter( 0 )) );
		
		
		TF1 * fp1 = new TF1( "fp1", "pol0" );
		fp1->SetRange( 0.2, 3.0 );
		hmwr->Fit( fp1, "NR", "", 0.2, 3.0 );
		hmwr->Draw("same");
		fp1->SetLineColor(kBlack);
		fp1->SetRange( 0, 4.0 );
		fp1->Draw("same");
		tl.DrawLatexNDC( 0.2, 0.75, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", fp1->GetChisquare(), fp1->GetNDF(), fp1->GetChisquare() / (float)fp1->GetNDF() ) );
		tl.DrawLatexNDC( 0.2, 0.7, TString::Format("corr = %0.3f", fp1->GetParameter( 0 )) );
		
		can->Print( rpName.c_str() );
	}
	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		LOG_F( INFO, "# of Mixed = %lu", nMixed );

		LOG_F( INFO, "iEventLoop=%d", iEventLoop );
		if ( iEventLoop <= 1 ){
			LOG_F( INFO, "buffer_pos.size()=%lu", buffer_pos.size() );
			LOG_F( INFO, "buffer_neg.size()=%lu", buffer_neg.size() );
			


			for ( string n : named_weights ){
				// norm_const[ n ] = meanContent( hWeight[ n ] );
				norm_const[ n ] = hWeight[ n ]->Integral();
				LOG_F( INFO, "Norm[%s] = %f", n.c_str(), norm_const[ n ] );
			}

			hsame->SetLineColor(kRed);
			hmixw->SetLineColor(kBlack);
		} else if ( 2 == iEventLoop ){



			TCanvas * can = new TCanvas( "can", "can", 2000, 1000 );
			gStyle->SetOptStat(0);
			can->Print( (rpName+"[").c_str() );
			can->SetTopMargin( 0.01 );
			can->SetLeftMargin( 0.1 );
			can->SetRightMargin( 0.05 );

			book->get( "same_pt_angle" )->Draw("colz");
			can->Print( rpName.c_str() );
			book->get( "mix_pt_angle" )->Draw("colz");
			can->Print( rpName.c_str() );
			book->get( "mixwa_pt_angle" )->Draw("colz");
			can->Print( rpName.c_str() );
			book->get( "mixwb_pt_angle" )->Draw("colz");
			can->Print( rpName.c_str() );

			compareSlice( 0, 30, can );
			float step = config.get<float>( "params.step", 0.25 );
			for ( float i = 0; i < 10; i+=step ){
				compareSlice( i, i+step, can );
			}
			
			TH2 * hs = ((TH2*)book->get( "same_dEta_dPhi" ));
			TH2 * hm = ((TH2*)book->get( "mix_dEta_dPhi" ));
			hm->SetLineColor( kRed);

			int idPhi1 = hs->GetXaxis(  )->FindBin( 1.0 );
			int idPhi2 = hs->GetXaxis(  )->FindBin( 1.5 );

			int idEta1 = hs->GetYaxis(  )->FindBin( 0.0 );
			int idEta2 = hs->GetYaxis(  )->FindBin( 1.0 );

			double n = hs->Integral( idPhi1, idPhi2, idEta1, idEta2 );
			double d = hm->Integral(idPhi1, idPhi2, idEta1, idEta2);
			if ( d > 0 && n > 0 ){
				LOG_F ( INFO, "n=%f, d=%f", n, d );
				hm->Scale( n / d );
			}
			else {
				LOG_F( ERROR, "Cannot scale hm" );
			}


			TH1 *hsx = hs->ProjectionX( "hsx" );
			TH1 *hsy = hs->ProjectionY( "hsy" );

			TH1 *hmx = hm->ProjectionX( "hmx" );
			TH1 *hmy = hm->ProjectionY( "hmy" );

			TH2 * hc = (TH2*)hs->Clone( "corr_dEta_dPhi" );
			hc->Add( hm, -1 );




			can->Print( (rpName + "]").c_str() );

		}
	}

};


#endif