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


	vector<TrackObj> buffer_pos;
	vector<TrackObj> buffer_neg;

	TRandom3 r3;

	TH2 *hmix, *hsame, *hmixw;
	TH2 * hmap_hit;
	double Ihmap_hit;
	int chargeSum;
	vector<string> named_weights;



	// params
	float minPt, maxPt;

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

	virtual void analyzeEvent(){

		if ( 0 == iEventLoop){

			FemtoPair * pair = this->_fpr.get();

			if ( fabs(pair->d1_mEta) > 0.5 || fabs(pair->d2_mEta) > 0.5 ) return;
			

			TLorentzVector lv1, lv2, lv;
			lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.105 );
			lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.105 );
			lv = lv1 + lv2;
			// if ( lv.M() < 2.0 ) return;

			float pairPid = sqrt( pow( pair->d1_mPid, 2 ) + pow( pair->d2_mPid, 2 ) );

			if ( pairPid > 0.2  ) return;

			if ( abs(pair->mChargeSum) == chargeSum ){
				TrackObj to1( eventIndex, lv1, pair, 1 );
				TrackObj to2( eventIndex, lv2, pair, 2 );

				addToBuffer( to1 );
				addToBuffer( to2 );

				
				if ( lv.Pt() > minPt && lv.Pt() < maxPt ){
					book->fill( "hsame", lv.M(), lv.Pt() );
					hWeight[ "map_hit"   ]->Fill( to1.hit, to2.hit );

					float dEta = fabs( lv1.PseudoRapidity() - lv2.PseudoRapidity() );
					float dPhi = fabs(lv1.DeltaPhi( lv2 ));
					((TH2*)hWeight[ "dEta_dPhi" ])->Fill( dPhi, dEta );
					((TH2*)hWeight[ "pt_dEta" ])->Fill( dEta, lv.Pt() );
					((TH2*)hWeight[ "pt_dPhi" ])->Fill( dPhi, lv.Pt() );
					((TH3*)hWeight[ "pt_dEta_dPhi" ])->Fill( dPhi, dEta, lv.Pt() );
					
					((TH2*)hWeight[ "pt_dpt" ])->Fill( pair->d1_mPt, lv.Pt() );
					((TH2*)hWeight[ "pt_dpt" ])->Fill( pair->d2_mPt, lv.Pt() );

					((TH2*)hWeight[ "pt_dBL" ])->Fill( deltaBL(to1.mBackleg, to2.mBackleg), lv.Pt() );

					book->fill( "same_dEta_dPhi", dPhi, dEta );

					// book->fill( "same_pt", lv.Pt() );
					// book->fill( "same_eta", lv.PseudoRapidity() );
					// book->fill( "same_phi", lv.Phi() );
					// book->fill( "same_dCell", abs(to1.mCell - to2.mCell), lv.M() );
					// book->fill( "same_dMod", abs(to1.mModule - to2.mModule), lv.M() );
					// book->fill( "same_dBL", deltaBL(to1.mBackleg, to2.mBackleg), lv.Pt() );

					book->fill( "same_pt_angle", to1.lv.Angle( to2.lv.Vect() ), lv.Pt() );
				}
				
			} // if chargeSum

			eventIndex ++;
		}


		if ( 1 == iEventLoop ){
			
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
					// if ( deltaBL(to1.mBackleg, to2.mBackleg ) > 15  ) continue;
					// LOG_F( INFO, "dBL=%d", deltaBL(to1.mBackleg, to2.mBackleg ) );
					// if ( lv.M() < 2.0 ) return;
					// hmix->Fill( lv.M(), lv.Pt() );

					float dPhi = fabs(to1.lv.DeltaPhi( to2.lv ));
					float dEta = fabs( to1.lv.PseudoRapidity() - to2.lv.PseudoRapidity() );
					float angle = to1.lv.Angle( to2.lv.Vect() );
					// float wdd = fweight( hWeight[ "dEta_dPhi" ], norm_const["dEta_dPhi"], to1.lv.DeltaPhi( to2.lv ), fabs( to1.lv.PseudoRapidity() - to2.lv.PseudoRapidity() ) );
					// float whit = fweight( hmap_hit, norm_const[ "map_hit" ], to1.hit, to2.hit );
					// float wdPhi = fweight( (TH2*)hWeight[ "pt_dPhi" ], norm_const[ "pt_dPhi" ], dPhi, lv.Pt() );
					// float wdEta = fweight( (TH2*)hWeight[ "pt_dEta" ], norm_const[ "pt_dEta" ], dEta, lv.Pt() );
					// float wptdd = fweight( (TH3*)hWeight[ "pt_dEta_dPhi" ], norm_const[ "pt_dEta_dPhi" ], dPhi, dEta, lv.Pt() );
					// int hit1 = to1.hit;
					// int hit2 = to2.hit;
					// if ( to2.hit < to1.hit ){
					// 	hit1 = to2.hit;
					// 	hit2 = to1.hit;
					// }
					// float wpthit = fweight( (TH3*)hWeight[ "pt_hit_map" ], norm_const[ "pt_hit_map" ], hit1, hit2, lv.Pt() );
					// LOG_F( INFO, "weight dPhi ( pt=%f, dPhi=%f, dEta=%f ) = %f", lv.Pt(), dPhi, dEta, wptdd );

					// float wptdpt1 = fweight( (TH2*)hWeight[ "pt_dpt" ], norm_const[ "pt_dpt" ], to1.lv.Pt(), lv.Pt() );
					// float wptdpt2 = fweight( (TH2*)hWeight[ "pt_dpt" ], norm_const[ "pt_dpt" ], to2.lv.Pt(), lv.Pt() );

					// float wdbl = fweight( (TH2*)hWeight[ "pt_dBL" ], norm_const[ "pt_dBL" ], deltaBL(to1.mBackleg, to2.mBackleg), lv.Pt() );

					float wpta = fweight( hWeight[ "same_pt_angle" ], norm_const[ "same_pt_angle" ], angle, lv.Pt() );

					// book->fill( "mix_pt", lv.Pt() );
					// book->fill( "mix_eta", lv.PseudoRapidity() );
					// book->fill( "mix_phi", lv.Phi() );
					// book->fill( "mix_dCell", abs(to1.mCell - to2.mCell), lv.M() );
					// book->fill( "mix_dMod", abs(to1.mModule - to2.mModule), lv.M() );
					// book->fill( "mix_dBL", deltaBL(to1.mBackleg, to2.mBackleg), lv.Pt() );

					book->fill( "mix_pt_angle", angle, lv.Pt() );
					book->fill( "mixw_pt_angle", angle, lv.Pt(), wpta );

					book->fill( "mix_dEta_dPhi", dPhi, dEta );
					book->fill( "mixw_dEta_dPhi", dPhi, dEta,  wpta );
					// if ( wdPhi != wdPhi )
					// 	LOG_F( INFO, "wdEta=%f, wdPhi=%f", wdEta, wpta );
					

					hmix->Fill( lv.M(), lv.Pt() );
					if ( to1.hit != to2.hit )
						hmixw->Fill( lv.M(), lv.Pt(),  wpta );
				}
			}
		} // iEventLoop == 1
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

		int bin1 = hs->GetXaxis()->FindBin( 2.0 );
		int bin2 = hs->GetXaxis()->FindBin( 4.0 );

		float Ihs = hs->Integral( bin1, bin2 );
		float Ihm = hm->Integral( bin1, bin2 );
		float Ihmw = hmw->Integral( bin1, bin2 );


		hm->Scale( Ihs / Ihm );


		bin1 = hs->GetXaxis()->FindBin( 0.2 );
		bin2 = hs->GetXaxis()->FindBin( 1.5 );

		Ihs = hs->Integral( bin1, bin2 );
		Ihmw = hmw->Integral( bin1, bin2 );
		hmw->Scale( Ihs / Ihmw );
		
		// TH1 * hsum = (TH1*)hmw->Clone( TString::Format("hsum_%d_%d", ptBin1, ptBin2 ) );
		// hsum->Add( hm );

		hs->SetTitle( TString::Format( "%0.2f < p_{T}^{#mu#mu} < %0.2f (GeV/c); M_{#mu#mu} (GeV/c^{2})", pt1, pt2 ) );
		
		hs->Draw("pe");
		hm->Draw("same");
		if ( (pt1 + pt2)/2.0 > 2.01 )
			hmw->Draw("same");
		// hsum->SetLineColor( kGreen );
		// hsum->Draw("same");
		gPad->SetLogy(1);
		can->Print( "rpMixed.pdf" );

		TH1 * hmr = (TH1*)hm->Clone( TString::Format("hmr_%d_%d", ptBin1, ptBin2 ) ); 
		for ( int i = 1; i < hmr->GetXaxis()->GetNbins()+1; i++ ){
			float v1 = hmr->GetBinContent( i );
			float v2 = hs->GetBinContent( i );
			float e2 = hs->GetBinError( i );
			if ( 0 == v2 ) continue;
			hmr->SetBinContent( i, v1/v2 );
			hmr->SetBinError( i, e2 / v2 );
		}

		TH1 * hmwr = (TH1*)hmw->Clone( TString::Format("hmwr_%d_%d", ptBin1, ptBin2 ) ); 
		if ( (pt1 + pt2)/2.0 > 2.01 ){
			for ( int i = 1; i < hmwr->GetXaxis()->GetNbins()+1; i++ ){
				float v1 = hmwr->GetBinContent( i );
				float v2 = hs->GetBinContent( i );
				float e2 = hs->GetBinError( i );
				if ( 0 == v2 || 0 == v1  ) continue;
				hmwr->SetBinContent( i, v1 / v2 );
				hmwr->SetBinError( i, e2 / v2 );
			}
		}
		// hmr->Divide( hs );
		// hmw->Divide( hs );

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
		
		
		if ( (pt1 + pt2)/2.0 > 2.01 ){
			TF1 * fp1 = new TF1( "fp1", "pol0" );
			hmwr->Fit( fp1, "NR", "", 0.0, 1.2 );
			hmwr->Draw("same");
			fp1->Draw("same");
			tl.DrawLatexNDC( 0.2, 0.75, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", fp1->GetChisquare(), fp1->GetNDF(), fp1->GetChisquare() / (float)fp1->GetNDF() ) );
			tl.DrawLatexNDC( 0.2, 0.7, TString::Format("corr = %0.3f", fp1->GetParameter( 0 )) );
		}
		can->Print( "rpMixed.pdf" );
	}
	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		LOG_F( INFO, "iEventLoop=%d", iEventLoop );
		if ( 0 == iEventLoop ){
			LOG_F( INFO, "buffer_pos.size()=%lu", buffer_pos.size() );
			LOG_F( INFO, "buffer_neg.size()=%lu", buffer_neg.size() );
			


			for ( string n : named_weights ){
				// norm_const[ n ] = meanContent( hWeight[ n ] );
				norm_const[ n ] = hWeight[ n ]->Integral();
				LOG_F( INFO, "Norm[%s] = %f", n.c_str(), norm_const[ n ] );
			}

			hsame->SetLineColor(kRed);
			hmixw->SetLineColor(kBlack);
		} else if ( 1 == iEventLoop ){

			TCanvas * can = new TCanvas( "can", "can", 2000, 1000 );
			gStyle->SetOptStat(0);
			can->Print( "rpMixed.pdf[" );
			can->SetTopMargin( 0.01 );
			can->SetLeftMargin( 0.1 );
			can->SetRightMargin( 0.05 );

			for ( float i = 0; i < 10; i+=0.25 ){
				compareSlice( i, i+0.25, can );
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




			can->Print( "rpMixed.pdf]" );

		}
	}

};


#endif