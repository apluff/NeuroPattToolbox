function index = windingNumberAngles(vx, vy, loc, radius, type)
% Estimate the winding number (Poincare index) around a circle of radius
% RADIUS centred on the the row and column coordinates in LOC in the vector
% field defined by VX and VY, by computing the angular differences between
% all consecutive vectors. TYPE is a string indicating where the centre of
% the path should be: "centre" (default value) indicates the centre of the
% cell containing LOC, "point" indicates the nearest grid point, "both"
% indicates that both should be computed and returned.

if nargin < 5
    type = "centre";
end

% First check that calculation is not off the edge of the array
if any(floor(loc(1:2))-radius < 0) || floor(loc(1))+radius>size(vx,1) ||...
        floor(loc(2))+radius>size(vx,2)
    index = nan;
    return
end

% Convert vectors to complex values
vf = vx + 1i*vy;

% Compute the indices of the counterclockwise circuit
circs = cell(1,2);
% Compute path based on the centre of the cell containing LOC
if strcmp(type, "centre") || strcmp(type, "both")
    circRowCent = floor(loc(1)) + [-radius+1:radius, radius:-1:-radius+1];
    circColCent = floor(loc(2)) + [0:-1:-radius+1, -radius+1:radius, radius:-1:1];
    circs{1} = sub2ind(size(vf), circRowCent, circColCent);
end
% Compute paths based on the nearest grid point
if strcmp(type, "point") || strcmp(type, "both")
    circRowPoint = round(loc(1)) + [-radius:radius-1, radius:-1:-radius+1];
    circColPoint = round(loc(2)) + [0:-1:-radius+1, -radius:radius-1, radius:-1:1];
    circs{2} = sub2ind(size(vf), circRowPoint, circColPoint);
end

% Calculate the total angle difference
index = nan(1, length(circs));
for icirc = 1:length(circs)
    if ~isempty(circs{icirc})
        ivf = vf(circs{icirc});
        totAngDiff = sum(anglesubtract(angle([ivf(2:end) ivf(1)]), angle(ivf)));
        index(icirc) = round(totAngDiff/(2*pi));
    end
end

% Output the index
index = index(~isnan(index));
if isempty(index)
    index = nan;
end

