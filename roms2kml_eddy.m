function [kmlStr] = roms2kml_eddy(romsFileName,romsGrdName,kmlFileName,vlevel,year_output,month_output,day_output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [kmlStr] = roms2kml_eddy(romsFileName,romsGrdName,kmlFileName,vlevel,year_output,month_output,day_output)
%
%   E.g.: roms2kml_eddy('/data/Roms_simula/SimulaRaia/op_out/avg/Raia_','/data/Roms_simula/SimulaRaia/op_conf/Raia_grd3_masked.nc','demo.kml',40,2012,4,16)
%
% This function locate eddies in ROMS outputs, calsify them as cyclonic or anticyclonic, tracks (if possible)
% during previous dates and export all the visualization to KML or KMZ (zipped) format to use with Google Earth.
% The output kml or kmz file can be dropped directly to Google Earth.
% 
% Input arguments:
%       romsFileName: This script is adpated to work with the operational model of the ROMS-IEO
%                     It is aim to work with daily averaged files, one time_step per file. These files have
%                     the typical name: Raia_yyyymmdd_avg.nc. Thus, you must provide the complete route 
%                     of the average file and the prefix, for example: '/data/Raia_'
%       romsGrdName:  provide the directory if the file is not in the same working path
%       kmlFileName:  provide the name of the output *.kml or *.kmz file
%       vlevel:       vertical level of the ROMS file (positive means s-level, negative depth in meters)
%       year_output:  yyyy
%       month_output: mm
%       day_output:   dd
%
% As defined by Henson and Thomas (2007), an eddy consists of a region of high vorticity (the core), 
% surrounded by a circulation cell (the ring), which experiences high rates of strain. Such regions can be 
% detected with the Okubo-Weiss parameter (Okubo, 1970; Weiss 1991), commonly abbreviated W. 
% The Okubo-Weiss parameter estimates vorticity and the normal and sheer components of strain and combines 
% them in a way that allows easy identification of regions that are either strain dominated (W > 0) 
% or vorticity-dominated (W < 0). For a detailed description, please see Henson and Thomas (2007).
%
% Flag clusters of cells where W is significantly negative. These cells are vorticity-dominated 
% and therefore likely to be eddy cores. We distinguish two main methods:
% 1) The threshold value of W to be specified as the number of standard deviations 
% from the mean (following Henson and Thomas (2007) and Isern-Fontanet et al. (2003, 2004, 2006), 
% which used the value -0.2). 
% 2) or as an absolute value (following Chelton et al. (2007), which used -2.0e-022)
%                   W_method (select 1 or 2 on the user defined options)
%
% Created by Pablo Otero (Jun 2012) based on the matlab tools available at:
% http://code.google.com/p/googleearthtoolbox
% and the trackin tool of Jean-Yves Tinevez (simple tracker)
% http://www.mathworks.com/matlabcentral/fileexchange/34040
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Lines to test the script
%
%clear all
%close all
%romsGrdName = '/data/Roms_simula/SimulaRaia/op_conf/Raia_grd3_masked.nc';
%romsFileName = '/data/Roms_simula/SimulaRaia/op_out/avg/Raia_';
%year_output = 2012;
%month_output = 4;
%day_output = 16;
%kmlFileName = 'kk_eddy.kml';
%vlevel=40;
%coef=1;



%------------------------------ PREDEFINED OPTIONS ------------------------------------------ %

%-- ROMS file

malla = '.nc';
time_step=1;		% If history file, then select the time_step to process.

%-- PROCESSING OPTIONS

coef=1;                 % Coefficient to multiply the Okubo_Weiss and the vorticity fields 
track=1; 		% To track the pathway of the eddy. If "0", the kml only includes the position
			% of the current eddies without their pathways
radius_eddy_core=15; 	% Ignore eddies with radius (km) below this threshold
W_method = 1;		% Select 1 or 2
daysbefore = 6;		% Number of days previous to the selected date to track the eddy

%--KML

%           id = 'polygon';
%        idTag = 'id';
    timeStamp = ' ';
timeSpanStart = ' ';
 timeSpanStop = ' ';
  description = '';
   visibility = 1;

    lineWidth = 0.5;
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
       region = ' ';

         name = 'Eddy core';
description_cyclonic = 'Cyclonic eddy core detected from the surface velocity field of the ROMS model';
description_anticyclonic = 'Anti-cyclonic eddy core detected from the surface velocity field of the ROMS model';  
    lineColor = '7fbdbdbd';
polyColor_red = '7fb43104';
polyColor_blue = '7f5858fa'; 

     name_track = 'Eddy track';
lineColor_track = 'FFdba901';
lineWidth_track = 2.0;
 altitude_track = 100.0; % Trackes are higher than polygons to facilitate the visualization
description_track = 'track during previous time steps';

%--TRACKING TOOL

max_linking_distance = 0.2;
max_gap_closing = 1;
debug = true;

%---------------------------------- END OF USER DEFINITIONS ------------------------------------

end_date=julian(year_output,month_output,day_output);
days=[end_date-daysbefore:end_date];

for i_frame=1:length(days)

  kk=gregorian(days(i_frame))
  year=num2str(kk(1));
  if(kk(2)<10)
    month=strcat('0',num2str(kk(2)));
  else
    month=num2str(kk(2));
  end
  if(kk(3)<10)
    day=strcat('0',num2str(kk(3)));
  else
    day=num2str(kk(3));
  end

  daystr=strcat(year,month,day)

  filein=strcat(romsFileName, daystr, '_avg',malla);

  [lat,lon,mask,W]=get_okubo(filein,romsGrdName,time_step,vlevel,coef);
  W=W.*mask;
  W2=W;

  [lat,lon,mask,xi]=get_vort(filein,romsGrdName,time_step,vlevel,coef);

  if (W_method==1)
        isee=find(W<0);
        Wmean=mean(W(isee));
        Wsd=std(W(isee));
        limit_core=Wmean-0.2*Wsd; %Aunque Henson and Thomas dicen literalmente "limit_core=-0.2*Wsd"
        W(W>=limit_core)=0;
  elseif (W_method==2)
        W(W>=-2.0e-022)=0;
  else
        disp(['Error: Select 1 or 2 in the W_method']);
  end
  
  % Identify the core eddies
  c=contourc(lon(1,:),lat(:,1),W,[0 0]);
  isee=find(c(1,:)==0);
  c(1,isee)=nan; c(2,isee)=nan;

  % Ignore the small ones
  area_eddy_core=pi*(radius_eddy_core/111)^2;

  neddy=0;
  output = '';
  isee=find(isnan(c(1,:)));
  for count=1:length(isee)
   if(count~=length(isee))
     [area,cx,cy]=polycenter(c(1,isee(count)+1:isee(count+1)-1),c(2,isee(count)+1:isee(count+1)-1)); 
   else
     [area,cx,cy]=polycenter(c(1,isee(count)+1:end),c(2,isee(count)+1:end));
   end
     
   if(area>area_eddy_core)
      icenter=find( abs(lon(1,:)-cx) == min(min(abs(lon-cx)))   );
      jcenter=find( abs(lat(:,1)-cy) == min(min(abs(lat-cy)))   );

      neddy=neddy+1;
      xcenter(neddy)=cx; ycenter(neddy)=cy;
    
    if( i_frame==length(days) )  
      if(xi(jcenter,icenter)<0)
        output = [output ge_poly(c(1,isee(count)+1:isee(count+1)-1),c(2,isee(count)+1:isee(count+1)-1),...
                'name',name, ...
                'polyColor',polyColor_red,...
                'snippet', snippet,...
                'description',description_anticyclonic, ...
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
      else
        output = [output ge_poly(c(1,isee(count)+1:isee(count+1)-1),c(2,isee(count)+1:isee(count+1)-1),...
                'name',name, ...
                'polyColor',polyColor_blue,...
                'snippet', snippet,...
                'description',description_cyclonic, ...
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
      end % End painting red or blue in cyclonic or anticyclonic eddies
     end % End constructing kml features of the last day
    end % End selecting the bigger eddies
   end % End looping in the total perimeter of the contoured eddies

  %Save the centroid location
  points{i_frame} = [xcenter' ycenter']
  clear xcenter ycenter
end


%Track points and make a plot
figure(1)
clf
hold on

for i_frame = 1 : length(days)  
    str = num2str(i_frame);
    for j_point = 1 : size(points{i_frame}, 1)
        pos = points{i_frame}(j_point, :);
        plot(pos(1), pos(2), 'x')
        text('Position', pos, 'String', str)
    end   
end

[ tracks adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

n_tracks = numel(tracks);
colors = hsv(n_tracks);

all_points = vertcat(points{:});



 

for i_track = 1 : n_tracks 
    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.  
   
  if( ~isnan(tracks{i_track}(end)) )
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    
    plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :))   

    output = [output ge_plot(track_points(:,1),track_points(:, 2),...
                'name',name_track, ...
                'snippet', snippet,...
                'description',description_track, ...
                'region', region, ...
                'timeStamp',timeStamp, ...
                'timeSpanStart',timeSpanStart, ...
                'timeSpanStop',timeSpanStop, ...
                'visibility',visibility, ...
                'lineColor',lineColor_track,...
                'lineWidth',lineWidth_track,...
                'altitude',altitude_track,...
                'altitudeMode',altitudeMode,...
                'extrude',extrude,...
                'tessellate',tessellate) ];
  end
end

pcolor(lon,lat,xi); shading flat; hold on;
contour(lon,lat,W,[0 0],'k');

% Me queda guardar los tracks en el kml
ge_output(kmlFileName,output,'name','ROMS - Eddy cores and tracks');




