ivim_tools2
===========

This software package provides the necessary MATLAB scripts, phase 
distributions and data sets to reproduce the figures in the abstract 
"A Framework to calculate the IVIM signal for different diffusion gradient
profiles" by A. Wetscherek and F. B. Laun, submitted to the ISMRM 23rd Annual 
Meeting & Exhibition. 

The provided framework is based on the concept of using normalized phase distributions 
as described in an article recently accepted for publication in Magnetic Resonance in
Medicine, the scripts of which are available as ivim_tools:

A. Wetscherek,  B. Stieltjes and F. B. Laun: "Flow-Compensated Intravoxel Incoherent 
Motion Diffusion Imaging" (doi: 10.1002/mrm.25410).

ivim_tools2 includes phase distributions for a larger variety of gradient profiles,
such as sine and cosine gradients and the interface of the functions to calculate 
the IVIM signal attenuation is more user friendly and the code more readable. This 
currently comes at the price of speed, such that fitting of parameter maps might take
much longer, if the methods of ivim_tools2 are used.

---------------------------------------------------------------------------

Contents:

  0. Copyright

  1. Phase distributions
  
  2. Script part1: attenuation curves for different attenuation 

  3. Script part2: model fit to liver / pancreas ROI data

---------------------------------------------------------------------------

0. Copyright
------------

Included in this software distribution are the following files, which are 
covered by seperate BSD licenses: 
- jacobianest.m from the DERIVESTsuite (c) 2007 by John D'Errico

This software is covered by the following BSD license:
---------------------------------------------------------------------------

Copyright (c) 2014, Andreas Wetscherek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

---------------------------------------------------------------------------



1. Phase distributions
----------------------

The phase distributions are found in the file generic.dat, which contains the
following variables:

-	N:		Ratios Tv/l for which phase distributions are available (between 1 and 100).
			for 0<1<N linear interpolation is exact.

-	phis:	centers of the bins for the phase distributions (x-values).

-	bip: bipolar gradient profile (used in MRM manuscript).
-   fc0: simple flow compensated profile .
-	fc1: symmetric flow compensated profile (used in MRM manuscript).
-	sin_x: sine profile with x full periods during diffusion time T. (x = 1, 2, 4, 6, 10, 16)
-	cos0_x: simple (theoretical) cosine profile with x full periods during diffusion time T. (x = 1, 2, 4, 6, 10, 16)
-	cos1_x: symmetric cosine profile with x full periods during diffusion time T. (x = 2, 4, 6, 10, 16)
 
 

2.  Script part1: attenuation curves for different attenuation 
------------------------------------------------------

- 	script reproduces the figures from the ISMRM abstract and in addition produces an attenuation
	curve for venules (v = 5-10 mm/s, l = 2mm).

-	get_norm_phd(N, profile) fetches the normalized phase distribution for the specified profile
	(e.g. 'bip' or 'cos0_6') and ratio N = Tv/l and is called by get_IVIM_signal.

-	get_IVIM_signal(b, T, l, v, profile) calculates the IVIM signal for the specified parameters,
	where b, T, l and v must be matrices of equal size. The input must be shaped such that along the 
	1st dimension N = Tv/l is constant. This is chosen, because then the same phase distribution applies
	and it doesn't need to be fetched for every data point. The function returns an equally sized matrix
	containing the signal attenuation.
	
-	get_IVIM_laminar(b, T, l, v, profile) calculates the IVIM signal for the specified parameters, but
	assumes a parabolic velocity profile as in laminar flow, where v equals the average velocity. The
	shape of the matrices is not relevant here, since the function reshapes the data and passes it on
	to get_IVIM_signal
	
	
	
 3. Script part2: model fit to liver / pancreas ROI data
 ------------------------------------------------
 
-	script performs a fit of the IVIM model to the data set provided and returns plots of the best fit
	solution. Parameter values and fit uncertainties are quoted in the standard output.
 
-	ROI data is from the study in the manuscript and contained in the file 'data.mat'. The acquisition
	settings are specified in the equally sized variables b, T and betta (the latter for the gradient 
	profile). The variables pancreas / s_pancreas and liver / s_liver contain the averaged signal in the
	ROI and the standard deviation, respectively.
	
-	generic_model_laminar(x, b, T, betta) calculates the IVIM signal for data sets where either the 'bip'
	profile (betta == 1) or the 'fc1' profile (betta == 2) is used. b, T and betta must be equally shaped.
	The parameters x are the parameters of the IVIM model: 
	
		D = x(1) : tissue diffusion coefficient (10 ^-3 mm²/s)
		f = x(2) : perfusion fraction
		l = x(3) : characteristic length of vessel segments (mm)
		v = x(4) : characteristic velocity (mm/s)