function [kmlStr] = roms2kml_currents(romsFileName,romsGrdName,kmlFileName,timestep,level,spares)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [kmlStr] = roms2kml_currents(romsFileName,romsGrdName,kmlFileName,timestep,level,spares)
%
%   E.g.: roms2kml_currents('roms_his.nc','roms_grd.nc','demo.kml',1,-10,3)
%
% This function exports a ROMS current vector map to KML or KMZ (zipped) format to use with Google Earth
% The output kml or kmz file can be dropped directly to Google Earth.
%
% The scaling of the arrow is dependent on the distance between adjacent points. The minimum distance between
% rho-points of the model is assigned to correspond with the maximum speed. Additionally, this script also
% writes to the kml file scale arrows at the corners of the domain. The scale arrow is predefined to be 0.2m/s,
% which should accomplish main requirements on regional ocean modeling.
% 
% Input arguments:
% 	romsFileName: provide the directory if the file is not in the same working path
%	romsGrdName:  provide the directory if the file is not in the same working path
%	kmlFileName: provide the name of the output *.kml or *.kmz file
%	timestep: ROMS time step
%	level: positive value indicates the s-layer and negative values depth in meters
%	spares: spatial resolution of the output, defined as the gap (N-1) between cells of the model.
%		1-> |x|x|x|x|x| 2-> |x| |x| |x| 3-> |x| | |X| | 4-> |x| | | |x| and so on...
%
% Some parameters can be modified into the code. Please, visit the help page:
%       showdemo ge_quiver
% to gain insight about other options
%
% Created by Pablo Otero (Oct 2009) based on the matlab tools available at:
% http://code.google.com/p/googleearthtoolbox
%
% Some remarks: 'mat2gray.m' and 'gra2ind.m' should be added to the matlab path. Note that some
% functions like ge_output.m, ge_quiver.m or authoptions.m have been modified from the original source.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<5)
 disp(['Please, check the number of arguments']);
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Extract data from ROMS file %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 nc=netcdf(romsGrdName);
 lon=nc{'lon_rho'}(:); lat=nc{'lat_rho'}(:);  mask=nc{'mask_rho'}(:);
 lonu=nc{'lon_u'}(:);  latu=nc{'lat_u'}(:);   masku=nc{'mask_u'}(:);
 lonv=nc{'lon_v'}(:);  latv=nc{'lat_v'}(:);   maskv=nc{'mask_v'}(:);
 close(nc)

 nc=netcdf(romsFileName);
 time=nc{'scrum_time'}(timestep)./86400;
 close(nc);

 % Get the horizontal slice
 datau=get_hslice(romsFileName,romsGrdName,'u',timestep,level,'u');
 isee=find(masku~=0);
 datau=griddata(lonu(isee),latu(isee),datau(isee),lon,lat,'nearest');

 datav=get_hslice(romsFileName,romsGrdName,'v',timestep,level,'v');
 isee=find(maskv~=0);
 datav=griddata(lonv(isee),latv(isee),datav(isee),lon,lat,'nearest');

 mask(mask==0)=nan;
 datau=datau.*mask; datav=datav.*mask;

 %Reduce the spatial resolution
 datau=squeeze(datau(1:spares:end,1:spares:end));
 datav=squeeze(datav(1:spares:end,1:spares:end));
 mask=squeeze(mask(1:spares:end,1:spares:end));
 lon=squeeze(lon(1:spares:end,1:spares:end));
 lat=squeeze(lat(1:spares:end,1:spares:end));

 %Scale depending on the maximum speed and the distance between points
 %The factor 'maxL' indicates how many cells are allowed to be occupied by 
 %the maximum vector speed. Note that this operation is performed after 
 %reducing the spatial resolution.
 % E.g. (1 cell)	(3 cell)	(5 cell)
 % |x|x|x|x|x|x|	|x|x|x|x|x|x|	|x|x|x|x|x|x
 % |>|			|---->|		|-------->|
 maxL=3;
 distlon=min(diff(lon(1,:)));
 distlat=min(diff(lat(:,1)));
 mindist=min(distlon,distlat);
 speed=sqrt(datau.^2+datav.^2);
 maxspeed=max(max(speed));
 scaleL=mindist*maxL/maxspeed;

 %%%%%%%%%%%%%%%%%%%
 % Create the plot %
 %%%%%%%%%%%%%%%%%%%
 figure
 isee=find(~isnan(mask));
 arrows(lon(isee),lat(isee),complex(datau(isee),datav(isee)).*mindist./maxspeed,1,'k');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Write the kml structure %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Create the description of the file
 kk=sprintf('%3.4g',time);
 if(level>0)
   snippet=['Currents at the ',num2str(level), ' s-level during model day',kk];
 else
   snippet=['Currents at ',num2str(level), ' m-depth during model day ',kk];
 end
 description=['Currents from the Regional Ocean Modelling System (ROMS) configuration ran by the Ocean Modeling Group of the Instituto Espanol de Oceanografia (IEO) at the Centro Oceanografico de A Coruna, Spain.'];

 % Write the structure of the kml file
 kmlStr = ge_quiver(lon(isee),lat(isee),datau(isee).*scaleL,datav(isee).*scaleL,'magnitudeScale',1,'snippet',snippet,'description',description);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Write the reference arrow at each corner of the domain %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % First, the arrow
 kmlStr1 = ge_quiver(max(max(lon)),max(max(lat)),0.2*scaleL,0,'name','0.2 m/s','description','Arrow scale','magnitudeScale',1,'lineColor','ffff0000','lineWidth',2);
 kmlStr2 = ge_quiver(max(max(lon)),min(min(lat)),0.2*scaleL,0,'name','0.2 m/s','description','Arrow scale','magnitudeScale',1,'lineColor','ffff0000','lineWidth',2);
 kmlStr3 = ge_quiver(min(min(lon)),max(max(lat)),0.2*scaleL,0,'name','0.2 m/s','description','Arrow scale','magnitudeScale',1,'lineColor','ffff0000','lineWidth',2);
 kmlStr4 = ge_quiver(min(min(lon)),min(min(lat)),0.2*scaleL,0,'name','0.2 m/s','description','Arrow scale ','magnitudeScale',1,'lineColor','ffff0000','lineWidth',2);

 % And the text of the scale. We only want the text and not associated icon of the placemark. Reducing
 % the size of the icon (with the scale option) overlaps the text over the arrow scale. Humm...one possible
 % solution is to make transparent the icon. Making the icon fully transparent also displaces the text. Hum...
 % Then, we can select a very low opacity per cent (e.g. 2%). To know the hexadecimal value, we can use the
 % Calculator application of our operating system. Insert 2/256 in decimal mode and after permute to hex mode.

 iconURLStr='http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png';
 iconColorStr='02ffffff';

 kmlStr5 = ge_point(max(max(lon)),max(max(lat)),10,'name','0.2 m/s','iconURL',iconURLStr,'iconColor',iconColorStr);
 kmlStr6 = ge_point(max(max(lon)),min(min(lat)),10,'name','0.2 m/s','iconURL',iconURLStr,'iconColor',iconColorStr);
 kmlStr7 = ge_point(min(min(lon)),max(max(lat)),10,'name','0.2 m/s','iconURL',iconURLStr,'iconColor',iconColorStr);
 kmlStr8 = ge_point(min(min(lon)),min(min(lat)),10,'name','0.2 m/s','iconURL',iconURLStr,'iconColor',iconColorStr);


 %%%%%%%%%%%%%%%%%%%%%%
 % Write the kml file %
 %%%%%%%%%%%%%%%%%%%%%%
 ge_output(kmlFileName,[kmlStr kmlStr1 kmlStr2 kmlStr3 kmlStr4, kmlStr5, kmlStr6, kmlStr7, kmlStr8],'name','Model currents','description',description,'snippet',snippet);


