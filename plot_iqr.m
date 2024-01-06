function p1 = plot_iqr(X, data, centerline, cmap, five_to_ninetyfive, Alpha)
%plot_iqr creates a lineplot of the interquartile range (IQR) given data
%   Usage: plot_iqr(X, data, centerline, cmap, five_to_ninetyfive)
% 
%   Inputs:
%       X: A data domain against which to plot the measure of interest
%           (defaults to 1:M)
%       data: An N x M matrix, with N being individuals in a population,
%       and M being the measure of interest.
%       centerline: Whether the center of the line should be median or mean
%           (defaults to median)
%       cmap: The color for the shaded region denoting the IQR
%       five_to_ninetyfive: true (1) or false (0). If true, plots the 5 to
%           95th percentile as well as the IQR. (defaults to false)
%       Alpha: the base transparency 
% 
%   Output: a line plot of the IQR.

    if ~exist('centerline','var') || isempty(centerline) || strcmp(centerline, 'median')
        Y = median(data,1);
    elseif strcmp(centerline, 'mean')
        Y = mean(data,1);
    end

    if ~exist('cmap','var') || isempty(cmap)
        cmap = [0.5, 0.5, 0.5];
    end

    if ~exist('five_to_ninetyfive','var') || isempty(five_to_ninetyfive)
        five_to_ninetyfive = false;
    end

    if ~exist('Alpha','var') || isempty(Alpha)
        Alpha = 0.25;
    end

    if ~exist('X','var') || isempty(X)
        X = 1:size(Y,2);
    end

    err75 = prctile(data,75);
    err25 = prctile(data,25);
    hold on;
    p2 = plot(X, err25, 'w-');
    p3 = plot(X, err75, 'w-');
    p4 = fill([X'; flipud(X')], [err25'; flipud(err75')], cmap, 'EdgeColor', 'none', 'FaceAlpha', Alpha);

    if five_to_ninetyfive
        err95 = prctile(data,95);
        err05 = prctile(data,5);
        p5 = plot(X, err05, 'w-');
        p6 = plot(X, err95, 'w-');
        p7 = fill([X'; flipud(X')], [err05'; flipud(err95')], cmap, 'EdgeColor', 'none', 'FaceAlpha', Alpha/2); 
    end

    
    p1 = plot(X,Y, 'k-');
    p1.LineWidth = 1;
    hold off;


end