%% Flowline Analysis
%Date created:      
%29 July 2018, updated 21 December 2020

%Copyright:         
%Annika Horlings, University of Washington
%with some code included from Chade Greene, University of Texas at Austin 
%(available on Mathworks) and Nicholas Holschuh, Amherst University.

%Purpose:
%To extract coordinates of ice flowines, speed, and horizontal divergence
%along ice flow lines in Antarctica or Greenland.

%% Read in the data
%Interactively read in MEASURES velocity data
cd('/Users/annikahorlings/Documents/Projects/Layer_thinning/Paper/Final_RunsMARCH2020/Data/Velocities_and_strain_rates');
[xg, yg, x_vel] = grdread('antarctic_ice_vel_phase_map_v01.nc');%x velocities 1996 - 2018
[xg, yg, y_vel] = grdread('antarctic_ice_vel_phase_map_v01.nc');%y velocities 1996 - 2018
x = double(xg);
y = double(yg);
vx=double(x_vel); 
vy=double(y_vel);

%Load REMA, for plotting
cd('/Users/annikahorlings/Documents/Projects/Layer_thinning/Paper/Final_RunsMARCH2020/Horizontal_strain_rates_maps/rema');
[xr,yr,Zr] = rema_data('xy');
xr; % xim: the x axis of the underlying image
yr; % yim: the y axis of the underlying image
Zr; %get(gcf, 'Colormap'); % im: the n x m x 1 values for the ArcticDEM hillshade
minc = 0;% minc/maxc: the limits for black and white in your hillshade
maxc = 2000;

%Use Nick's code to use multiple colormaps by converting underlying
%hillshade to a "truecolor" image instead of a "falsecolor" image (ie it
%needs to be n x m x 3, containing r,g,b values for each pixel, not n x ,x
%1, mapped to a colorbar). The script colorlock_mat.m will take in an image
%and a colormap, to produce a truecolor image.
%No hillshade
im_rgb = colorlock_mat(Zr,gray,[minc maxc]);

%% Convert x and y to lat and long 
[x_lat, y_lon_f] = ps2ll(x, y(1:length(x))); 
[x_lat_f, y_lon] = ps2ll([x' ones(1, (length(y) - length(x)))]', y);

%% Calculate speed
sp = sqrt(x_vel.^2+y_vel.^2); %calculate speed
speed = double(sp);     %convert to double precision

%Smooth the data for strain rate calculations. (Too noisy if you don't smooth.)
vx_s = smooth_ndh(vx, 20, 0); %smooth the data
vy_s = smooth_ndh(vy, 20, 0);

%% Calculate the total horizontal strain rate
e_horiz = (gradient(vx_s, 450) + gradient(vy_s, 450));  %the total lhorizontal strain rate 
                            %is the sum of the horizontal strain rates in 
                            %the x direction and y directions 

%% Make the grid
[xgrd, ygrd] = meshgrid(x,y); %this creates the grid

DomainSize = 450; %spatial resolutino of measures data, in km

%% Define starting coordinates
%%%%%Option 1:%%%%%%%

%If you don't know what the field looks like...
% figure(1);
% streamslice(x, y, vx, vy);  %display a slice of what it looks like
%                             %note that streamslice only takes double
%                             %precision numbers
% %in polar sterographic coordinates
% xstart = -1.6*10^6;    %in meters
% ystart = -2*10^5;       %in meters


%%%%Option 2:%%%%%%%

%Alternatively, interactively get points from figure ...
% message = sprintf('Double click on the following map to select starting coordinates.');
% m = msgbox(message);
% uiwait(m)
% 
% f = figure(1);
% imagesc(x, y, speed);
% xlabel('Easting (km)');
% ylabel('Northing (km)');
% c = colorbar;
% caxis([0 500]);
% ylabel(c, 'Speed m yr^{-1})');
% [xstart,ystart] = getpts;
% close(f);


%%%%Option 3:%%%%%%%

%Define the starting coordinates:
xstart = [-1.2908e06 -1.2908e06 -1.2908e06 -1.2908e06 -1.2908e06 -1.2908e06 -1.58e06 -1.6e06 -1.7e06 -1.6992e06 -1.5456e06 -1.3297e06 -1.2538e06]; 
ystart = [-5.4838e05 -4.5e05 -4e05 -3.5e05 -3e05 -2.5e05 -0.5e05 -0.5e05 -0.5e05 -1.0736e05 -1.3418*10^5 -5.8712e05 -4.4409e05];


%% Calculate streamlines
%This portion includes Chad Greene's Code 

%Compute streamlines from vector input
XY = stream2(x,y,vx,vy,xstart,ystart,[0.1 floor(DomainSize/0.045)]);    %stream2 only takes
                                                                        %double precision numbers
%If you want to view a preliminary plot of the streamline, then:
% figure(2);
% streamline(x,y,vx,vy,xstart,ystart,[0.1 floor(DomainSize/0.045)]);  %streamline only takes double
                                                                     %double precision numbers                                                                        
% Preallocate outputs: 
lat = cell(size(xstart)); 
lon = cell(size(xstart)); 

% Convert stream lines to georeferenced coordinates: 
for k = 1:length(XY)
    try % Some seed locations may result in empty outputs 
        [lat{k},lon{k}] = ps2ll(XY{k}(:,1),XY{k}(:,2)); 
    end
end

%% Calculate associated speed, strain rate, distance vectors along flowline
%This portion includes Chad Greene's code.
%Uses cell arrays, which is advantageous when more than one flowline is
%being calculated.

% Preallocate outputs: 
d = cell(size(xstart)); 
v = cell(size(xstart)); 
t = cell(size(xstart)); 
e_h = cell(size(xstart));
a = cell(size(xstart));
st = cell(size(xstart));

for k = 1:length(XY)
            d{k} = pathdist(lat{k},lon{k},'units','meters'); % distance, units of meters
            v{k} = interp2(x,y,speed,XY{k}(:, 1),XY{k}(:, 2),'linear');  % speed, units of meters/year
                            %for one flowline; figure out how to do with
                            %multiple flowlines
            t{k} = cumsum([0;diff(d{k})]./v{k}); % time, units of years
            e_h{k} = interp2(x, y, e_horiz, XY{k}(:, 1),XY{k}(:, 2),'linear'); %total horizontal strain
                            %rate in units of yr^-1 along the flowline.
end
    
% Convert cells to matrices;
d = cell2mat(d); 
v = cell2mat(v);
t = cell2mat(t); 
e_h = cell2mat(e_h);
XY = cell2mat(XY);

%% Plotting
%If you want to plot the outputs:

%Plot the flowlines
figure;
cmap = intense; %I made this colormap
colormap(cmap) % or whatever colormap you want.
n = size(cmap,1); %# size of colormap
dmap=(0.1-0)/n; %# color step

f = imagesc(xr/1000,yr/1000, im_rgb);
f.AlphaData = ~isnan(Zr); % makes NaNs (ocean) transparent
set(gca, 'Ydir', 'Normal');
set(gca,'visible','off');
hold all; 

im_handle = imagesc(x/1000, y/1000, speed);
hold on;
for i = 2:2:length(XY(1, :))
    plot((XY(:, i-1))/1000, (XY(:, i))/1000, '-', 'Color', 'k', 'LineWidth', 2); %plot the flowline
    hold on;
end

for i = 2:2:length(XY(1, :))
    plot((XY(1, i-1))/1000, (XY(1, i))/1000, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1 1 1]); %plot the start of the flowline
    hold on;
end
c = colorbar;
caxis([0 800]);
ylabel(c, 'Speed (m yr^{-1})');
xlabel('Easting (km)');
ylabel('Northing (km)');
set(gca,'FontSize',14);

%Plot the speed and divergence rate as a function of distance along the flowline

figure;
subplot(2, 1, 1); %Speed
for i = 1:length(v(1, :))
    plot(d(:, i)/1000, v(:, i), 'LineWidth', 2);
    hold on;
end
ylabel('Speed (m yr{-1})');
xlabel('Distance along flowline (km)');

subplot(2, 1, 2); %Total horizontal divergence rate
for i = 1:length(e_h(1, :))
    plot(d(:, i)/1000, e_h(:, i), 'LineWidth', 2);
    hold on;
end
ylabel('Total Horizontal Divergence Rate (yr^{-1})');
xlabel('Distance along flowline (km)');


