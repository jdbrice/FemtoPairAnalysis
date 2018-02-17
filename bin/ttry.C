



void ttry( double scale = 0.98, double escale = 100.0 ){
	TFile *f = new TFile( "res.root" );
	TFile *fct = new TFile( "MTD_Scaled_STAR.root" );


	TCanvas * c1 = new TCanvas( "c1" );

	TH1 * rs = (TH1*)f->Get("pcand_rb")->Clone("raw_signal");
	TH1 * bg = (TH1*)f->Get("bgs_rb")->Clone( "background" );

	TH1 * ss = (TH1*)rs->Clone( "corr_signal" );
	// ss->Scale( 1.0, "width" );
	
	gStyle->SetOptStat(0);

	rs->Sumw2();
	ss->Sumw2();
	bg->Sumw2();

	rs->Draw("");
	rs->Draw("hsame");
	bg->Draw("same");
	bg->Draw("hsame");

	bg->Scale( scale / bg->Integral() );
	cout << "bg I = " << bg->Integral() <<endl;

	ss->Add( bg, -1 );
	ss->Scale( 1/(3.54e11) * escale );
	

	TCanvas * c2 = new TCanvas( "c2", "", 1200, 900 );
	ss->Draw();
	ss->Draw("hsame");

	ss->GetXaxis()->SetRangeUser( 0, 3.4 );
	ss->SetMaximum( 1e-6 );
	ss->SetMinimum( 1e-12 );

	TH1 * Scaled_jpsi_mumu = (TH1*)fct->Get( "Scaled_jpsi_mumu" );
	Scaled_jpsi_mumu->SetLineColor(kRed);
	Scaled_jpsi_mumu->Draw("samec");

	TH1 * hSum = (TH1*)Scaled_jpsi_mumu->Clone( "hsum" );

	TH1 * Scaled_omega_mumu = (TH1*)fct->Get( "Scaled_omega_mumu" );
	Scaled_omega_mumu->SetLineColor(kRed);
	Scaled_omega_mumu->Draw("samec");
	hSum->Add( Scaled_omega_mumu );

	TH1 * Scaled_phi_mumu = (TH1*)fct->Get( "Scaled_phi_mumu" );
	Scaled_phi_mumu->SetLineColor(kRed);
	Scaled_phi_mumu->Draw("samec");
	hSum->Add( Scaled_phi_mumu );

	TH1 * Scaled_ccbar_mumu = (TH1*)fct->Get( "Scaled_ccbar_mumu" );
	TH1 * Smooth_ccbar_mumu = (TH1*)Scaled_ccbar_mumu->Clone( "Smooth_ccbar_mumu" );
	Smooth_ccbar_mumu->Reset();
	TRandom3 r3;
	r3.SetSeed(0);
	for ( int i = 0; i <1000000; i++ ){
		float v = (Scaled_ccbar_mumu->GetRandom() + r3.Gaus( 0, 0.1 ));
		if ( v < 0.210 ) 
			continue;
		Smooth_ccbar_mumu->Fill( v );
	}
	Smooth_ccbar_mumu->Scale( Scaled_ccbar_mumu->Integral() / Smooth_ccbar_mumu->Integral() );
	Smooth_ccbar_mumu->SetLineColor(kRed);
	Smooth_ccbar_mumu->Draw("samec");
	hSum->Add( Smooth_ccbar_mumu );

	TH1 * Scaled_rho_mumu = (TH1*)fct->Get( "Scaled_rho_mumu" );
	Scaled_rho_mumu->SetLineColor(kRed);
	Scaled_rho_mumu->Draw("samec");
	hSum->Add( Scaled_rho_mumu );

	TH1 * Scaled_eta_mumu = (TH1*)fct->Get( "Scaled_eta_mumu" );
	Scaled_eta_mumu->SetLineColor(kRed);
	Scaled_eta_mumu->Draw("samec");
	hSum->Add( Scaled_eta_mumu );

	TH1 * Scaled_eta_gammamumu = (TH1*)fct->Get( "Scaled_eta_gammamumu" );
	Scaled_eta_gammamumu->SetLineColor(kRed);
	Scaled_eta_gammamumu->Draw("samec");
	hSum->Add( Scaled_eta_gammamumu );

	TH1 * Scaled_etaprime_gammamumu = (TH1*)fct->Get( "Scaled_etaprime_gammamumu" );
	Scaled_etaprime_gammamumu->SetLineColor(kRed);
	Scaled_etaprime_gammamumu->Draw("samec");
	hSum->Add( Scaled_etaprime_gammamumu );

	hSum->SetLineColor(kBlack);
	hSum->Draw("samec");


	gPad->SetLogy(1);

	TFile * fexport = new TFile( "signal_cocktail_export.root", "RECREATE" );
	fexport->cd();

	Scaled_omega_mumu->Write();
	Scaled_phi_mumu->Write();
	Smooth_ccbar_mumu->Write();
	Scaled_rho_mumu->Write();
	Scaled_eta_mumu->Write();
	Scaled_eta_gammamumu->Write();
	Scaled_etaprime_gammamumu->Write();
	hSum->Write();

	rs->Write();
	bg->Write();
	ss->Write();

	fexport->Write();

}