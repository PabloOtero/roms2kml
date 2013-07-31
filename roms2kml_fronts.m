function [kmlStr] = roms2kml_fronts(romsFileName,romsGrdName,kmlFileName,time_step,variable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [kmlStr] = roms2kml_fronts(romsFileName,romsGrdName,kmlFileName,time_step,variable)
%
%   E.g.: roms2kml_fronts('roms_his.nc','roms_grd.nc','demo.kml',1,'MLD')
%
% This function exports the isobaths of the ROMS configuration to KML or KMZ (zipped) format 
% to use with Google Earth.
% The output kml or kmz file can be dropped directly to Google Earth.
% 
% Input arguments:
%	romsGrdName:  provide the directory if the file is not in the same working path
%	kmlFileName: provide the name of the output *.kml or *.kmz file
%	time_step: vector containing the isobaths to plot
%
% Created by Pablo Otero (Jun 2012) based on the matlab tools available at:
% http://code.google.com/p/googleearthtoolbox
%
% Some remarks: 'mat2gray.m' and 'gra2ind.m' should be added to the matlab path. The
% ge_imagesc.m were also modified to properly plot the figure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all
%romsGrdName = '/data/Roms_simula/SimulaRaia/op_conf/Raia_grd3_masked.nc.1';
%romsFileName = '/data/Roms_simula/SimulaRaia/mar/avg/Raia_20120416_avg.nc.1';
%time_step=1;
%variable='SSS'
%kmlFileName='kk_SSS.kml';

if(nargin<5)
 disp(['Please, check the number of arguments']);
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Extract data from ROMS file %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 nc=netcdf(romsGrdName);
 lon=nc{'lon_rho'}(:); lat=nc{'lat_rho'}(:);  mask=nc{'mask_rho'}(:); mask(mask==0)=nan;
 close(nc)

 nc=netcdf(romsFileName);
 time=nc{'scrum_time'}(time_step)./86400;
 t = nc{'temp'}(time_step,:,:,:);
 s = nc{'salt'}(time_step,:,:,:);
 close(nc);


if(strcmp(variable,'MLD'))

  %%% Depth calculation %%%
  z = zatr(romsFileName);
  if(length(z)==3)
     z = squeeze(z(time_step,:,:,:));
  end
  dzw = diff(z);
  zw = zatw(romsFileName);
  if(length(z)==3)
    zw = squeeze(zw(time_step,:,:,:));
  end
  dz = diff(zw);

  %%% RHO VARIABLE %%%
  r = potden80(s,t,z); %From Ocean's pack

  %%% VERTICAL GRADIENT OF RHO
  r = reshape(r,size(t));

  % Voy a prescindir del valor absoluto para evitar momento de alta salinidad superficial que provocan alta estratificacion.
  %gradrho = abs(diff(r)./dzw);
  gradrho = diff(r)./dzw;

  %%% Compute MLD based on:
  %%%
  %%% - a temperature change from the ocean surface of 0.5 degree Celsius,
  %%% - a density change from the ocean surface of 0.125 (sigma units), 
  %%% - and a variable density change from the ocean surface corresponding to a temperature change of 0.5 degree Celsius.
  %%%
  %%% The MLD based on the variable density criterion is designed to account for the large variability of the coefficient of thermal expansion that characterizes seawater. 
  %%%
  %%% See http://www.nodc.noaa.gov/OC5/WOA94/mix.html
  for i=1:size(gradrho,2)
   for j=1:size(gradrho,3)
    %level = find(gradrho(:,i,j)==min(gradrho(:,i,j))); 			% Previous criterion ICES JMS
    %level = find( abs( t(:,i,j)-t(end,i,j) ) >= 0.5 );  			%Criterion a)
    level = find( r(:,i,j) >= squeeze(r(end,i,j))+0.125 ); 			%Criterion b)
    %level = find( r(:,i,j) >= squeeze(potden80(s(end,i,j),t(end,i,j)-0.5,0)) ); 	%Criterion c)
    dz(1:max(level),i,j) = 0;
   end
  end

  %%% MIX LAYER THICKNESS
  h1 = squeeze(sum(dz)).*mask;

elseif(strcmp(variable,'SST'))

  h1 = squeeze(t(end,:,:)).*mask;

elseif(strcmp(variable,'SSS'))

  h1 = squeeze(s(end,:,:)).*mask;

else

  disp(['Bad name of the variable']);

end


  %%% LOOKING FOR THE EDGE OF THE FRONT - Sobel edge filter
  Sobel = [1 2 1; 0 0 0; -1 -2 -1];
  H = conv2(h1,Sobel);
  V = conv2(h1,Sobel');
  edge = sqrt(H.^2 + V.^2);

  % THRESHOLDING TO THE EDGE
  alfa=0.1;
  I_max=max(max(edge));
  I_min=min(min(edge));
  level=alfa*(I_max-I_min)+I_min;
  Ibw=max(edge,level.*ones(size(edge)));
  Ibw=Ibw-level;
  Ibw=squeeze(Ibw(2:end-1,2:end-1));
  Ibw(1,:)=0; Ibw(end,:)=0; Ibw(:,1)=0; Ibw(:,end)=0;
  Ibw=Ibw.*mask;

  [n,m]= size(h1);
  xx= 2:m-1;
  yy= 2:n-1;

  filt= [1 2 1;0 0 0; -1 -2 -1]/8;  tv= 2;
  imo= conv2(h1, rot90(filt), 'same').^2 + conv2(h1, filt, 'same').^2;
  thresh= sqrt( tv* nanmean(nanmean( imo(yy,xx) ))  ); % Estimate threshold

  % The filters are defined for sqrt(imo), but since we calculated imo, compare
  %  to thresh ^2
  imout= ( imo >= thresh^2 );
  %imout= ( imo >= (thresh^2/20) );

  % Thin the wide edges
  xpeak= imo(yy,xx-1) <= imo(yy,xx) & imo(yy,xx) > imo(yy,xx+1) ;
  ypeak= imo(yy-1,xx) <= imo(yy,xx) & imo(yy,xx) > imo(yy+1,xx) ;
  imout(yy,xx)= imout(yy,xx) & ( xpeak | ypeak );
  imout(1:2,:)=0; imout(n-1:n,:)=0; imout(:,1:2)=0; imout(:,m-1:m)=0;

  % Link edge pixels together into lists of sequential edge points, one
  % list for each edge contour.  Discard contours less than 5 pixels long.
  [edgelist, labelededgeim] = pixellink(imout, 5, complex(lat,lon));

  % Display the labeled edge image with separate colours for each
  % distinct edge (choose your favorite colourmap!)
  labelededgeim(labelededgeim~=0)=1;
  %labelededgeim(labelededgeim==0)=nan;

  figure,    pcolor(lon,lat,labelededgeim); colormap(vga); shading flat; axis image; go_coast;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%% KML OPTIONS  %%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %AuthorizedOptions = authoptions(mfilename);


if(strcmp(variable,'MLD'))
    name = 'MLD fronts';
    lineColor = 'FFdaa520'; %Mostaza
elseif(strcmp(variable,'SST'))
    name = 'SST fronts';
    lineColor = 'FFb22222'; %Rojo
elseif(strcmp(variable,'SSS'))
    name = 'SSS fronts';
    lineColor = 'FF66cdaa'; %Verde	
end
           id = 'contour';
        idTag = 'id';
    timeStamp = ' ';
timeSpanStart = ' ';
 timeSpanStop = ' ';
  description = '';
   visibility = 1;
    
    lineWidth = 2.0;
      snippet = ' ';
     altitude = 1.0;
      extrude = 0;
   tessellate = 1;
 altitudeMode = 'relativeToGround';
  msgToScreen = false;
  forceAsLine = true;
     numLevels = 1;
   lineValues = [];
    lineAlpha = 'FF';
       cMap   = 'jet';
       region = ' ';

% parsepairs %script that pars


 output = '';

 for count=1:length(edgelist)
        output = [output ge_plot(edgelist{count}(:,2),edgelist{count}(:,1),...
                'name',name, ...
                'snippet', snippet,...
                'description',description, ...
                'region', region, ...
                'timeStamp',timeStamp, ...
                'timeSpanStart',timeSpanStart, ...
                'timeSpanStop',timeSpanStop, ...
                'visibility',visibility, ...
                'lineColor',lineColor,...
                'lineWidth',lineWidth,...
                'altitude',altitude,...
                'altitudeMode',altitudeMode,...
                'extrude',extrude,...
                'tessellate',tessellate) ];

 end


if(strcmp(variable,'MLD'))
 ge_output(kmlFileName,output,'name','ROMS - MLD fronts');
elseif(strcmp(variable,'SST'))
 ge_output(kmlFileName,output,'name','ROMS - SST fronts');
elseif(strcmp(variable,'SSS'))
 ge_output(kmlFileName,output,'name','ROMS - SSS fronts');
end
 


















