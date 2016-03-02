GEM
===

GxExM project

Collaboration with Bonn University funded by the NZ Royal Society

Aim is to quantify uncertanty in climate change simulations from the choice of hybrid and sowing date.

The raw data (.out) is at: C:\Dropbox\2014-GEM_paper\rawData or Ed's hard drive

NIWA's ERA-40 data across NZ (5 x 5 km grid)
Filtered for arable lands (~2100 pixels)
5 sowing dates
5 hybrids


23 December 2015: Moving to new machine. Many supporting data at C:\GitHubRepos\GEM_files

Analysis steps:

1) Met files (APSIM) created with netCDFToMet_GEM based on F:\\ERA-40_Reanalysis\\ from NIWA filtered 
2) Obs: Arable filter with CreateFilters.R at C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\gisLayers\\luFilter\\LandUseFilter_Allcats.txt LandUse = 3 (arable lands only)
3)



