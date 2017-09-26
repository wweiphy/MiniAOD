/*
 * qcdhelper.h
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


#ifndef QCDHELPER_H
#define QCDHELPER_H

#include <TFile.h>
#include <TString.h>
#include <TH2D.h>

class QCDHelper
{
	public:
		QCDHelper(TString path_to_sf_file_);
		~QCDHelper();
		double GetScaleFactor(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons);
		double GetScaleFactorError(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons);
		double GetScaleFactorErrorUp(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons);
		double GetScaleFactorErrorDown(int n_jets, int n_btags, int n_isoinverted_electrons, int n_isoinverted_muons);
		void Reset();
		void LoadFile(TString path_to_sf_file_);
			
	private:
		TFile* scalefactor_file = 0;
		TString path_to_sf_file = "";
		TH2D* Mu_SF = 0;
		TH2D* El_SF = 0;
		bool initialized = false;
};

#endif /* QCDHELPER_H */ 
