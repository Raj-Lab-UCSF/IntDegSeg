function [ev_changepts, ev_sfit, ev_1deriv, ev_2deriv] = find_ev_IntDegSeg(ev, method, knots, order)
%find_ev_IntDegSeg Finds the two indices bounding Int Deg Seg regimes
%   This function fits a spline to the eigenvalues and then finds change
%   points in the eigenvalue spectrum.
% 
%   Usage: ev_changepts = find_ev_IntDegSeg(ev, method, knots, order)
% 
%   Inputs:
%       ev: The eigenvalues of the Laplacian
%       method:
%           '1deriv' : evaluates the analytical first derivative then uses
%           the findchangepts function to find two points in the first
%           derivative that partition the data into the three segments by
%           the root mean square 'rms' statistic.
%           '2deriv' : evaluates the analytical second derivative then
%           finds its first maxima and last minima as the regime bounds.
%       knots: the number of knots in the spline (default = 3)
%       order: The polynomial order for the spline function (default = 16)
% 
%   Outputs:
%       ev_changepts: the identified change points in the eigen-spectrum.
%       ev_sfit: the spline fit to the eigenvalues
%       ev_1deriv: the analytical first derivative
%       ev_2deriv: the analytical second derivative
% 
%   Written by Benjamin Sipes, June 2023

if ~exist('knots','var') || isempty(knots)
    knots = 3;
end

if ~exist('order','var') || isempty(order)
    order = 10;
end

nroi = length(ev);

%%%% Determine Int Deg Seg Regimes from Template:
splinefit = spap2(knots, order, (1:nroi)', ev);
ev_sfit = fnval(splinefit, (1:nroi)');
ev_fit_deriv = fnder(splinefit);
ev_1deriv = fnval(ev_fit_deriv, (1:nroi)');

switch method
    case '1deriv'
        ev_changepts = findchangepts(ev_1deriv, "MaxNumChanges",2,"Statistic",'rms');

    case '2deriv'
        ev_fit_2deriv = fnder(ev_fit_deriv);
        ev_2deriv = fnval(ev_fit_2deriv, (1:nroi)');
        ev_changepts = zeros(2,1);
        ev_changepts(1) = min(find(islocalmax(ev_2deriv)));
        ev_changepts(2) = max(find(islocalmin(ev_2deriv)));
end

end