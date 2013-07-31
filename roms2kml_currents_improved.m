function [kmlStr] = roms2kml_currents_improved(romsFileName,romsGrdName,kmlFileName,timestep,level,spares)
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


%TEST
%romsFileName='/data/Roms_simula/SimulaRaia/op_out/his/Raia_20121118_his.nc';
%romsGrdName='/data/Roms_simula/SimulaRaia/op_conf/Raia_grd3_masked.nc';
%kmlFileName='/data/Roms_simula/SimulaRaia/op_out/kml/CurrentsSurf_20121118_test3.kml';
%timestep=1;
%level=-1;
%spares=4;


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

 altitude = 1.0;

 %%%%%%%%%%%%%%%%%%%
 % Create the plot %
 %%%%%%%%%%%%%%%%%%%
 figure
 isee=find(~isnan(mask));
 arrows(lon(isee),lat(isee),complex(datau(isee),datav(isee)).*mindist./maxspeed,1,'k');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Write the kml structure %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(kmlFileName,'w');

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<kml xmlns="http://earth.google.com/kml/2.0">\n');
fprintf(fid,'<Document>\n');

part1 = sprintf('<Placemark>\n<open>0</open>\n<Style>\n<LineStyle>\n<color>ff000000</color>\n<width>1</width>\n</LineStyle>\n</Style>\n<LineString>\n<coordinates>');
part2 = sprintf('</coordinates>\n</LineString>\n</Placemark>\n');


 % Pinto las referencias en las esquinas
 DX=0.2*scaleL;
 DY=0*scaleL;
 direction = 90;
 magnitude = DX;
 alpharad = deg2rad( [direction, (direction + 12), (direction - 12)]  );           %[rad]
 divPer= [magnitude, (3*magnitude/4), (3*magnitude/4) ]';

 %how to do this step once in parrallel?
        x(1) = sin( alpharad(1) ) * divPer(1);
        y(1) = cos( alpharad(1) ) * divPer(1);
        x(2) = sin( alpharad(2) ) * divPer(2);
        y(2) = cos( alpharad(2) ) * divPer(2);
        x(3) = sin( alpharad(3) ) * divPer(3);
        y(3) = cos( alpharad(3) ) * divPer(3);


 lon_ref=(max(max(lon))+min(min(lon)))/2;
 lat_ref=max(max(lat))+0.1;

 % bug fixes & addtions from Brett Grant
 x = x/(cos(lat_ref*(pi/180)));

 arrow_coords = [ num2str(lon_ref,'%11.7f'),',',num2str(lat_ref,'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon_ref+x(1),'%11.7f'),',',num2str(lat_ref+y(1),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon_ref+x(2),'%11.7f'),',',num2str(lat_ref+y(2),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon_ref+x(3),'%11.7f'),',',num2str(lat_ref+y(3),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon_ref+x(1),'%11.7f'),',',num2str(lat_ref+y(1),'%11.7f'),',',num2str(altitude),' '];

 fprintf(fid,'<Placemark>\n');
 fprintf(fid,'<name>0.2 m/s</name>\n');
 fprintf(fid,'<description>Flecha de escala</description>\n');
 fprintf(fid,'<Style>\n');
 fprintf(fid,'<LineStyle>\n');
 fprintf(fid,'<color>ff8A0886</color>\n');
 fprintf(fid,'<width>2</width>\n');
 fprintf(fid,'</LineStyle>\n');
 fprintf(fid,'</Style>\n');
 fprintf(fid,' <visibility>1</visibility>\n');
 fprintf(fid,'<LineString>\n');
 fprintf(fid,'<extrude>0</extrude>\n');
 fprintf(fid,'<altitudeMode>clampToGround</altitudeMode>\n');
 fprintf(fid,'<tessellate>1</tessellate>\n');
 fprintf(fid,'<coordinates>\n');
 fprintf(fid,'%s',arrow_coords);  
 fprintf(fid,'</coordinates>\n');
 fprintf(fid,'</LineString>\n');
 fprintf(fid,'</Placemark>\n');

%Decidimos quitarle la etiqueta "punto"
if(0)
 fprintf(fid,'<Placemark>\n');
 fprintf(fid,'<name>0.2 m/s</name>\n');
 fprintf(fid,'<description>Flecha de escala</description>\n');
 fprintf(fid,'<Style>\n');
 fprintf(fid,'<IconStyle>\n');
 fprintf(fid,'<color>02ffffff</color>\n');
 fprintf(fid,'<scale>1</scale>\n');
 fprintf(fid,'<Icon>\n');
 fprintf(fid,'<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n');
 fprintf(fid,'</Icon>\n');
 fprintf(fid,'</IconStyle>\n');
 fprintf(fid,'</Style>\n');
 fprintf(fid,'<Point id="point_point">\n');
 fprintf(fid,'<altitudeMode>relativeToGround</altitudeMode>\n');
 fprintf(fid,'<tessellate>1</tessellate>\n');
 fprintf(fid,'<extrude>0</extrude>\n');
 fprintf(fid,'<coordinates>\n');
 fprintf(fid,'%s',[num2str(lon_ref,'%11.7f'),',',num2str(lat_ref,'%11.7f'),',',num2str(10)]);  
 fprintf(fid,'</coordinates>\n');
 fprintf(fid,'</Point>\n');
 fprintf(fid,'</Placemark>\n');
end


% Añado la parte de ge_quiver
[row_count, col_count] = size(lon);
output = cell((row_count * col_count), 1);

lon=lon;
lat=lat;
DX=datau.*scaleL;
DY=datav.*scaleL;


for row = 1:row_count
    for col = 1:col_count

        if (DX(row,col) == 0) && (DY(row,col) == 0)

            direction = 0;

        elseif DX(row,col) == 0

            if DY(row,col) > 0
                direction = 0;
            else
                direction = 180;
            end

        elseif DY(row,col) == 0

            if DX(row,col) > 0
                direction = 90;
            else
                direction = 270;
            end

        else
            direction = rad2deg( atan2(  DX(row,col) , DY(row,col)  ) );

        end

        magnitude = sqrt( DX(row,col) .^ 2 + DY(row,col) .^ 2 );

       %continued ugliness
        alpharad = deg2rad( [direction, (direction + 12), (direction - 12)]  );           %[rad]
        divPer= [magnitude, (3*magnitude/4), (3*magnitude/4) ]';

        %how to do this step once in parrallel?
        %     x = sin( alpharad ) * divPer;
        %     y = cos( alpharad ) * divPer;
        x(1) = sin( alpharad(1) ) * divPer(1);
        y(1) = cos( alpharad(1) ) * divPer(1);
        x(2) = sin( alpharad(2) ) * divPer(2);
        y(2) = cos( alpharad(2) ) * divPer(2);
        x(3) = sin( alpharad(3) ) * divPer(3);
        y(3) = cos( alpharad(3) ) * divPer(3);


        % bug fixes & addtions from Brett Grant
% % %   %added for aspect ratio
        x = x/(cos(lat(row, col)*(pi/180)));
% % % modified output accuracy
        arrow_coords = [ num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col)+x(1),'%11.7f'),',',num2str(lat(row, col)+y(1),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col)+x(2),'%11.7f'),',',num2str(lat(row, col)+y(2),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col)+x(3),'%11.7f'),',',num2str(lat(row, col)+y(3),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col)+x(1),'%11.7f'),',',num2str(lat(row, col)+y(1),'%11.7f'),',',num2str(altitude),' '];

        lolo = strfind(arrow_coords,'NaN');
        if lolo
            arrow_coords = [ num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' ',...
                        num2str(lon(row, col),'%11.7f'),',',num2str(lat(row, col),'%11.7f'),',',num2str(altitude),' '];
        end
        clear lolo

   fprintf(fid,'%s',part1);
   fprintf(fid,'%s',arrow_coords); 
   fprintf(fid,'%s',part2); 

 end 
end

fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');
fclose(fid);


