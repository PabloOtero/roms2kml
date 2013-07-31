function [kmlStr] = roms2kml_slice(romsFileName,romsGrdName,kmlFileName,timestep,level,vname,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [kmlStr] = roms2kml_slice(romsFileName,romsGrdName,kmlFileName,timestep,level,vname,varargin)
%
%   E.g.: roms2kml_slice('roms_his.nc','roms_grd.nc','demo.kml',1,-10,'temp')
%
% This function exports a ROMS slice to KML or KMZ (zipped) format to use with Google Earth
% Do not forget to export also the output PNG image of the ROMS slice.
% The output kml or kmz file can be dropped directly to Google Earth.
% 
% Input arguments:
% 	romsFileName: provide the directory if the file is not in the same working path
%	romsGrdName:  provide the directory if the file is not in the same working path
%	kmlFileName: provide the name of the output *.kml or *.kmz file
%	timestep: ROMS time step
%	level: positive value indicates the s-layer and negative values depth in meters
%	vname: choose between (2D variables) 'zeta','ubar','vbar'...choose level=0
%                           (3D variables) 'temp','salt','u','v','w','omega','speed','AKs','AKt'
%
% Optional input parameters:
%	Optional parameterse are the same than those provided to 'ge_imagesc'. Please,
%	visit the help page:
%       showdemo ge_imagesc
%
% Created by Pablo Otero (Oct 2009) based on the matlab tools available at:
% http://code.google.com/p/googleearthtoolbox
%
% Some remarks: 'mat2gray.m' and 'gra2ind.m' should be added to the matlab path. The
% ge_imagesc.m should be also modified to properly plot the figure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<6)
 disp(['Please, check the number of arguments']);
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Extract data from ROMS file %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc=netcdf(romsGrdName);
lon=nc{'lon_rho'}(:); lat=nc{'lat_rho'}(:);
mask=nc{'mask_rho'}(:);
if(strcmp(vname,'u') | strcmp(vname,'ubar'))
 lonu=nc{'lon_u'}(:); latu=nc{'lat_u'}(:); masku=nc{'mask_u'}(:);
elseif(strcmp(vname,'v') | strcmp(vname,'vbar'))
 lonv=nc{'lon_v'}(:); latv=nc{'lat_v'}(:); maskv=nc{'mask_v'}(:);
elseif(strcmp(vname,'speed'))
 lonu=nc{'lon_u'}(:); latu=nc{'lat_u'}(:); masku=nc{'mask_u'}(:);
 lonv=nc{'lon_v'}(:); latv=nc{'lat_v'}(:); maskv=nc{'mask_v'}(:);
end
close(nc)

% Get the horizontal slice
if(strcmp(vname,'temp') | strcmp(vname,'salt') | strcmp(vname,'zeta'))
 data=get_hslice(romsFileName,romsGrdName,vname,timestep,level,'r');
elseif (strcmp(vname,'w') | strcmp(vname,'omega') | strcmp(vname,'AKt') | strcmp(vname,'Aks') )
 data=get_hslice(romsFileName,romsGrdName,vname,timestep,level,'w');
elseif (strcmp(vname,'u') | strcmp(vname,'ubar'))
 data=get_hslice(romsFileName,romsGrdName,vname,timestep,level,'u');
 isee=find(masku~=0);
 data=griddata(lonu(isee),latu(isee),data(isee),lon,lat,'nearest');
elseif (strcmp(vname,'v') | strcmp(vname,'vbar'))
 data=get_hslice(romsFileName,romsGrdName,vname,timestep,level,'v');
 isee=find(maskv~=0);
 data=griddata(lonv(isee),latv(isee),data(isee),lon,lat,'nearest');
elseif (strcmp(vname,'speed'))
 datau=get_hslice(romsFileName,romsGrdName,'u',timestep,level,'u');
 isee=find(masku~=0);
 datau=griddata(lonu(isee),latu(isee),datau(isee),lon,lat,'nearest');
 datav=get_hslice(romsFileName,romsGrdName,'v',timestep,level,'v');
 isee=find(maskv~=0);
 datav=griddata(lonv(isee),latv(isee),datav(isee),lon,lat,'nearest');
 data=sqrt(datau.^2+datav.^2);
 clear datau, datav;
end
 masknan=mask; masknan(mask==0)=nan;
 data=data.*masknan;
 nc=netcdf(romsFileName);
 time=nc{'scrum_time'}(timestep)./86400;
 close(nc);

 % The colormap is expanded to the complete range of data
 isee=find(mask==1);
 cLimLow = min(min(data(isee)));
 cLimHigh = max(max(data(isee)));

 %%%%%%%%%%%%%%%%%%%%%%%
 % Optional parameters %
 %%%%%%%%%%%%%%%%%%%%%%%
 AuthorizedOptions = authoptions( 'ge_imagesc' );

 % Predefined options
    altitudeMode = 'clampToGround';
    cBarFormatStr = '%-0.2g';
    visibility = 1;
    pngFileName = strcat(kmlFileName(1:end-3),'png');
    altitude = 10;

 % New options introduced as input arguments
 if ~isempty(varargin)
  for k = 1:2:length(varargin(:))
    if ~strcmp(varargin{k}, AuthorizedOptions)
        s=dbstack;
        error(['Unauthorized parameter name ' 39 varargin{k} 39 ' in ' 10,...
            'parameter/value passed to ' 39 s(end-1).name 39 '.']);
    end
    eval([varargin{k},'=varargin{',num2str(k+1),'};'])
  end
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Write the kml structure %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Create the description of the file
 if(strcmp(vname,'zeta'))
  vnameStr='Sea Surface elevation (m)';
 elseif(strcmp(vname,'ubar'))
  vnameStr='West-East batropic speed (m/s)';
 elseif(strcmp(vname,'vbar'))
  vnameStr='South-North batropic speed (m/s)';
 elseif(strcmp(vname,'temp'))
  vnameStr='Temperature (ºC)';
 elseif(strcmp(vname,'salt'))
  vnameStr='Salinity';
 elseif(strcmp(vname,'u'))
  vnameStr='West-East speed (m/s)';
 elseif(strcmp(vname,'v'))
  vnameStr='South-North speed (m/s)';
 elseif(strcmp(vname,'omega'))
  vnameStr='Vertical velocity (m/s)';
 elseif(strcmp(vname,'speed'))
  vnameStr='Speed (m/s)';
 elseif(strcmp(vname,'AKs'))
  vnameStr='Vertical salinity diffusion coeffcient (AKs)';
 elseif(strcmp(vname,'AKt'))
  vnameStr='Vertical thermal diffusion coeffcient (AKt)';
 end

 timeStr=sprintf('%3.4g',time);

 if(level>0)
   snippet=[vnameStr,' at the ',num2str(level), ' s-level during model day ',timeStr];
 elseif(level<0)
   snippet=[vnameStr,' at ',num2str(level), ' m-depth during model day ',timeStr];
 else
   snippet=[vnameStr, ' during model day ',timeStr];
 end
 description=['Data from the Regional Ocean Modelling System (ROMS) configuration run by the Ocean Modeling Group of the Instituto Espanol de Oceanografia (IEO) at the Centro Oceanografico de A Coruna, Spain.'];

 kmlStr = ge_imagesc(lon(1,:),lat(:,1),data,'imgURL',pngFileName,'alphaMatrix',mask,'cLimLow',cLimLow,'cLimHigh',cLimHigh,'altitudeMode',altitudeMode,'altitude',altitude);

 % Create the colorbar
 output = ge_colorbar(lon(end),lat(1),data,...
                            'numUnits',20,...
                             'cLimLow',cLimLow,...
                            'cLimHigh',cLimHigh,...
                       'cBarFormatStr',cBarFormatStr,...
                       'visibility',visibility);

 %%%%%%%%%%%%%%%%%%%%%%
 % Write the kml file %
 %%%%%%%%%%%%%%%%%%%%%%
 ge_output(kmlFileName,[output kmlStr],'name',kmlFileName);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Create the plot following pcolor and overwrite %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % The function 'ge_imagesc' creates a PNG image using 'imwrite', 
 % but we have experienced some problems related to the scaling
 % when exported to Google Earth (GE). GE only works with images
 % that have the same relation "number of pixels/geographical degree"
 % in both directions. In the experienced problem, we had one output
 % of the ROMS model with similar number of cells in both directions,
 % but representing a very different geographical extension. The result
 % was a latitudinal distortion when loaded in GE. We have tried to 
 % change the input resolution in 'imwrite' without success. The solution
 % taken here is to remove the image created by 'imwrite' and create 
 % other new.
 
 % The previous image is deleted
 delete(pngFileName);
 
 % A new figure is created with the same proportion in both directions
 fig=figure;
 pcolor(lon,lat,data.*mask); shading flat; axis image; axis off;

 % Examine the proportion of the box used in the plotting and close.
 posbox=get(gca,'PlotBoxAspectRatio')
 close(fig)

 % Another figure is created keeping the same relative size. We
 % depart from a 600x600 resolution image.
 figure('Position',[0 0 600*posbox(1)/posbox(2) 600]);
 pcolor(lon,lat,data.*mask); shading flat; axis image; axis off;

 % Expand axes to fill figure window
 set(gca,'Position',[0 0 1 1])

 % Set the color in figure and axes to none to keep a transparent background
 set(gca,'Color','none');
 set(gcf,'Color','none');
 export_fig(pngFileName,'-png')
