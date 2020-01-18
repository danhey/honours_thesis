function avg_array = bwdmean(center_array, dir)
%% Input Parameters
% center_array: 2D array of values defined at cell centers
% w: 'x' or 'y', direction in which average is taken

%% Out Parameter
% avg_array: 2D array of averaged values

center_shifted = circshift(center_array, -dir);
avg_array = (center_shifted + center_array) / 2;