 This repository contains files that I created to do bispectrum analysis of wave data from surface
 and adv files from the observations made at Fire Island, MVCO and Matanzas Inlet. 

 My workflow consisted of separating the data based on the skewness events (positive and negative).
 The spectra was bifurcated into two parts based on the sign of skewness. Then, I created an
 average spectra file. Once that was created, I used the Stingray python package to get the 3rd
 order cummulant and bispectrum magnitude. 

Matlab codes are used to read in the skewness and spectrum data and write out ".mat files that can 
be read in python notebooks to perform the bispectrum analysis. 

 The idea was that i can see some patterns similar to Crawford and Hay's 2001 paper. I did not have full
 understanding of taking the wave data and doing a transformation based on wave theory as the paper explains.
 So, i just take the data and try to perform bispectrum analysis straight away. 

 This folder has 3 subfolders.

 1. adv_data -> where i used ADV near bottom data to bifurcate the spectrum. 
 2. workhorse_data-> used the workhorse surface spectra
 3. individual_data-> used the spectra data for individual events such as the 14 Feb storm in Fire Island. 


  I tried to quantify skewness by taking the values under the diagonal of a bispectrum magnitude and 3rd order
  cummulant. Both did not help me quantify the skewness neither were the patterns very clear. 

  I think the transformation that i don't know could be a reason. 

  The workhorse data had the *best* analysis. 

 The presentation contains some notes on what I learnt. 
