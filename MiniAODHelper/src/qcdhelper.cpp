/*
 * qcdhelper.cpp
 * 
 * Copyright 2017 Michael Wassmer
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * see https://gitlab.cern.ch/ttH/reference/blob/master/definitions/Moriond17.md#69-lj-channel-qcd-estimation
 */


#include "MiniAOD/MiniAODHelper/interface/qcdhelper.hpp"


QCDHelper::QCDHelper(TString path_to_sf_file_)
{
	LoadFile(path_to_sf_file_);
}

void QCDHelper::LoadFile(TString path_to_sf_file_)
{
	path_to_sf_file = path_to_sf_file_;
	if(path_to_sf_file!="")
	{
		scalefactor_file = TFile::Open(path_to_sf_file);
	}
	if(scalefactor_file)
	{
		El_SF = (TH2D*)scalefactor_file->Get("El_FakeSF");
		Mu_SF = (TH2D*)scalefactor_file->Get("Mu_FakeSF");
		initialized = true;
	}
}

void QCDHelper::Reset()
{
	scalefactor_file = 0;
	path_to_sf_file = "";
	Mu_SF = 0;
	El_SF = 0;
	initialized = false;
}

double QCDHelper::GetScaleFactor(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons)
{
	if(!initialized) return 0.;
	int bin = -1;
	double sf = 0.;
	if(n_isoinverted_electrons==1&&n_isoinverted_muons==0) 
	{
		bin = El_SF->GetBin(n_jets,n_btags);
		sf = El_SF->GetBinContent(bin);
	}
	else if(n_isoinverted_electrons==0&&n_isoinverted_muons==1)
	{
		bin = Mu_SF->GetBin(n_jets,n_btags);
		sf = Mu_SF->GetBinContent(bin);
	}
	else 
	{
		return 0.;
	}
	return sf<0. ? 0. : sf;
}

double QCDHelper::GetScaleFactorError(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons)
{
	if(!initialized) return 0.;
	int bin = -1;
	double sf_err = 0.;
	if(n_isoinverted_electrons==1&&n_isoinverted_muons==0) 
	{
		bin = El_SF->GetBin(n_jets,n_btags);
		sf_err = El_SF->GetBinError(bin);
	}
	else if(n_isoinverted_electrons==0&&n_isoinverted_muons==1)
	{
		bin = Mu_SF->GetBin(n_jets,n_btags);
		sf_err = Mu_SF->GetBinError(bin);
	}
	else 
	{
		return 0.;
	}
	return sf_err<0. ? 0. : sf_err;
}

double QCDHelper::GetScaleFactorErrorUp(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons)
{
	if(!initialized) return 0.;
	double sf = GetScaleFactor(n_jets,n_btags,n_isoinverted_electrons,n_isoinverted_muons);
	double sf_err = GetScaleFactorError(n_jets,n_btags,n_isoinverted_electrons,n_isoinverted_muons);
	return sf+sf_err <0. ? 0. : sf+sf_err;
}

double QCDHelper::GetScaleFactorErrorDown(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons)
{
	if(!initialized) return 0.;
	double sf = GetScaleFactor(n_jets,n_btags,n_isoinverted_electrons,n_isoinverted_muons);
	double sf_err = GetScaleFactorError(n_jets,n_btags,n_isoinverted_electrons,n_isoinverted_muons);
	return sf-sf_err <0. ? 0. : sf-sf_err;
}

QCDHelper::~QCDHelper()
{
	if(initialized) scalefactor_file->Close();
}


