function [kmlStr] = roms2kml_isobaths(romsGrdName,kmlFileName,isobaths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [kmlStr] = roms2kml_isobaths(romsGrdName,kmlFileName,isobaths)
%
%   E.g.: roms2kml_isobaths('roms_grd.nc','demo.kml',[200 500 1000])
%
% This function exports the isobaths of the ROMS configuration to KML or KMZ (zipped) format 
% to use with Google Earth.
% The output kml or kmz file can be dropped directly to Google Earth.
% 
% Input arguments:
%	romsGrdName:  provide the directory if the file is not in the same working path
%	kmlFileName: provide the name of the output *.kml or *.kmz file
%	isobaths: vector containing the isobaths to plot
%
% Created by Pablo Otero (Oct 2009) based on the matlab tools available at:
% http://code.google.com/p/googleearthtoolbox
%
% Some remarks: 'mat2gray.m' and 'gra2ind.m' should be added to the matlab path. The
% ge_imagesc.m were also modified to properly plot the figure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc=netcdf(romsGrdName);
lon=nc{'lon_rho'}(:); lat=nc{'lat_rho'}(:);
mask=nc{'mask_rho'}(:);
h=nc{'h'}(:);
close(nc)

[c,kk]=contour(lon,lat,h,isobaths);
kmlStr=ge_contour(lon,lat,h,'numLevels',isobaths);
ge_output(kmlFileName,kmlStr,'name','ROMS isobaths');