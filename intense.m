function cm_data=intense(m)
%sea colormap
%Annika Horlings
%24 April 2020
%This colormap uses the customcolormap function downloaded from mathworks


%% Make sure to use the custom colormap function
cm = customcolormap([0 0.1 0.2 0.4 0.6 0.8 1], ([102 51 0; 153 76 0; 153 153 0; 0 153 0; 0 153 153; 0 0 153; 255 255 255]./255));


%[0 0.20 0.4 0.6 0.7 0.8 0.9 1], [255/255 99/255 71/255; 255/255 20/255 147/255; 0 0 255/255; 106/255 90/255 205/255; 30/255 144/255 255/255; 72/255 209/255 204/255; 175/255 238/255 238/255; 1 1 1]);

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);
  
end
end