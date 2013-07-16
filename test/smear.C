{
  gSystem->Load("/afs/cern.ch/user/m/meridian/work/tests/KinFitter/lib/Linux26SL5_x86_64_gcc462/5.32.00/libTKinFitter.so");
  TH1F aSmear("aSmear","aSmear",200,-2,2);
  TH1F bSmear("bSmear","bSmear",200,-2,2);
  TH1F hSmear("hSmear","hSmear",200,-2,2);

  TH1F aMagSmear("aMagSmear","aMagSmear",600,0., 300);
  TH1F bMagSmear("bMagSmear","bMagSmear",600,0., 300);
  TH1F hMagSmear("hMagSmear","hMagSmear",600,0., 300);

  TH1F aMag("aMag","aMag",600,0., 300);
  TH1F bMag("bMag","bMag",600,0., 300);
  TH1F hMag("hMag","hMag",600,0., 300);

  TH1F aMagFit("aMagFit","aMagFit",200,-2,2);
  TH1F bMagFit("bMagFit","bMagFit",200,-2,2);
  TH1F hMagFit("hMagFit","hMagFit",200,-2,2);

  TH1F aPt("aPt","aPt",600,0., 300);
  TH1F bPt("bPt","bPt",600,0., 300);
  TH1F hPt("hPt","hPt",600,0., 300);

  TH1F thetaDir("theta","theta",300,0., 2*TMath::Pi());
  TH1F phiDir("phi","phi",300, -TMath::Pi()/2, TMath::Pi()/2);

  TH1F ptBoost("ptBoost","ptBoost",300,0., 300);
  TH1F etaBoost("etaBoost","etaBoost",100,-5., 5.);
  TH1F phiBoost("phiBoost","phiBoost",100,-TMath::Pi()/2., TMath::Pi()/2.);

  TH1F fitChi2("fitChi2","fitChi2",1000,0.,100.);
  TH1F fitPChi2("fitPChi2","fitPChi2",100,0.,1.);

  TH1F aDiff("aDiff","aDiff",200,-1.,1.);
  TH1F bDiff("bDiff","bDiff",200,-1.,1.);

  TH1F aPull("aPull","aPull",200,-5.,5.);
  TH1F bPull("bPull","bPull",200,-5.,5.);
  // Isotropic 2 body decay
  TF1 theta("theta","1",0,2*TMath::Pi());
  TF1 phi("phi","1",-TMath::Pi()/2.,TMath::Pi()/2.);

  //Resonance parameters, mass, pt, eta
  TF1 mass("mass","1/(TMath::Pi()*[1]* (1+TMath::Power((x-[0])/[1],2)))",10.,200.);
  mass.SetParameter(0,91.181);
  mass.SetParameter(1,2.4);

  TF1 pt_boost("pt_boost","1",0.,100.);
  TF1 eta_boost("eta_boost","1",-2.,2.);

  for (int i=0;i<10000;++i)
    {
      //      std::cout << "++++" << std::endl; 
      //      float my_mass=mass.GetRandom();
      float my_mass=91.181;

      float my_theta=theta.GetRandom();
      float my_phi=phi.GetRandom();
      thetaDir.Fill(my_theta);
      phiDir.Fill(my_phi);

      float my_pt_boost=pt_boost.GetRandom();
      float my_eta_boost=eta_boost.GetRandom();
      float my_phi_boost=phi.GetRandom();
      ptBoost.Fill(my_pt_boost);
      etaBoost.Fill(my_eta_boost);
      phiBoost.Fill(my_phi_boost);

      TVector3 a_v;
      TVector3 b_v;

      TLorentzVector h;
      h.SetPtEtaPhiM(my_pt_boost,my_eta_boost,my_phi_boost,my_mass);
      TVector3 boost=h.BoostVector();
      //      boost.Print();

      a_v.SetMagThetaPhi(my_mass/2,my_theta,my_phi);
      b_v.SetMagThetaPhi(my_mass/2,TMath::Pi()-my_theta,TMath::Pi()+my_phi);

//       a_v.Print();
//       b_v.Print();

      TLorentzVector a(a_v,a_v.Mag());
      TLorentzVector b(b_v,b_v.Mag());

      a.Boost(boost);
      b.Boost(boost);
      
      h=a+b;

      //      std::cout << h.M() << std::endl;
      
      float aS=0.1*gRandom->Gaus(0.,1.);
      float bS=0.1*gRandom->Gaus(0.,1.);

      TVector3 smearVect_a = a.Vect();
      smearVect_a.SetMag(a.Vect().Mag()*(1+aS));
      TLorentzVector a_smear(smearVect_a,smearVect_a.Mag());

      TVector3 smearVect_b = b.Vect();
      smearVect_b.SetMag(b.Vect().Mag()*(1+bS));
      TLorentzVector b_smear(smearVect_b,smearVect_b.Mag());

      TLorentzVector h_smear=a_smear+b_smear;
      // std::cout << h_smear.M() << std::endl;

      TMatrixD covJet_a(3,3);
      TMatrixD covJet_b(3,3);

      covJet_a(0,0)=TMath::Power(0.1,2);
      covJet_a(1,1)=1e-9;
      covJet_a(2,2)=1e-9;

      covJet_b(0,0)=TMath::Power(0.1,2);
      covJet_b(1,1)=1e-9;
      covJet_b(2,2)=1e-9;

      TFitParticleRelPtEtaPhi fit_a("fit_a","fit_a",&smearVect_a,0.,&covJet_a);
      TFitParticleRelPtEtaPhi fit_b("fit_b","fit_b",&smearVect_b,0.,&covJet_b);

      //BW mass constraint
      //      TFitConstraintMBW Zmass("Zmass","Zmass",0,0,91.181,2.4);
      //standard (delta) mass constraint
      TFitConstraintM Zmass("Zmass","Zmass",0,0,91.181);
      Zmass.addParticles1(&fit_a,&fit_b);

      TKinFitter fit;
      fit.addMeasParticle(&fit_a);
      fit.addMeasParticle(&fit_b);
      fit.addConstraint(&Zmass);
      fit.setMaxNbIter(50);
      fit.setMaxDeltaS(5e-5);
      fit.setMaxF(1e-4);
      fit.setVerbosity(0);
      fit.fit();

      if ( fit.getStatus() == 0 ) {
	fitPChi2.Fill( TMath::Prob( fit.getS(), fit.getNDF() ) );
	fitChi2.Fill( fit.getS());
	const TMatrixD* parfit_a = fit_a.getParCurr();
	const TMatrixD* covMatrixFit_a = fit_a.getCovMatrixFit();
	const TMatrixD* parfit_b = fit_b.getParCurr();
	const TMatrixD* covMatrixFit_b = fit_b.getCovMatrixFit();
	aDiff.Fill( (*parfit_a)(0,0) - 1/(1+aS));
	aPull.Fill( ((*parfit_a)(0,0) - 1/(1+aS))/TMath::Sqrt((*covMatrixFit_a)(0,0)));
	bDiff.Fill( (*parfit_b)(0,0) - 1/(1+bS));
	bPull.Fill( ((*parfit_b)(0,0) - 1/(1+bS))/TMath::Sqrt((*covMatrixFit_b)(0,0)));
	aMagFit.Fill((*parfit_a)(0,0) -1);
	bMagFit.Fill((*parfit_b)(0,0) -1);
	TLorentzVector a_fit=*fit_a.getCurr4Vec();
	TLorentzVector b_fit=*fit_b.getCurr4Vec();
	hMagFit.Fill((a_fit+b_fit).M()/h.M()-1.);
      }

      aSmear.Fill(aS);
      bSmear.Fill(bS);
      hSmear.Fill(h_smear.M()/h.M()-1);

      aMagSmear.Fill(a_smear.Vect().Mag());
      bMagSmear.Fill(b_smear.Vect().Mag());
      hMagSmear.Fill(h_smear.M());

      aMag.Fill(a.Vect().Mag());
      bMag.Fill(b.Vect().Mag());
      hMag.Fill(h.M());

      aPt.Fill(a.Pt());
      bPt.Fill(b.Pt());
      hPt.Fill(h.Pt());
      
    }

  TFile *f=new TFile("smear.root","RECREATE");
  aSmear.Write();
  bSmear.Write();
  hSmear.Write();

  aMagSmear.Write();
  bMagSmear.Write();
  hMagSmear.Write();

  aMagFit.Write();
  bMagFit.Write();
  hMagFit.Write();

  aMag.Write();
  bMag.Write();
  hMag.Write();

  aPt.Write();
  bPt.Write();
  hPt.Write();

  thetaDir.Write();
  phiDir.Write();

  ptBoost.Write();
  etaBoost.Write();
  phiBoost.Write();

  fitChi2.Write();
  fitPChi2.Write();

  aDiff.Write();
  aPull.Write();
  bDiff.Write();
  bPull.Write();
  
  f->Write();
  f->Close();
}
