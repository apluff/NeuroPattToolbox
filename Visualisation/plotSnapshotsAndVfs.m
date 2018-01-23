function plotSnapshotsAndVfs(snaps, vfx, vfy, allPatternLocs, isPhase, plotSeparate, nsnaps)
% Function to plot consecutive snapshots with their corresponding velocity
% fields. VFS can either be a complex matrix where real gives x-coords and 
% imag gives y-coords, or can be input as two varibles. Default is to overlay
% vectors on snapshots, but setting options="separate" will plot them 
% separately.

% Parse inputs and set defaults
if ~exist('isPhase', 'var') || isempty(isPhase)
    isPhase = false;
end
if ~exist('plotSeparate', 'var')
    plotSeparate = false;
end

% Choose colour maps
if isPhase
    clrmap = pmkmp_new;
    clims = [-pi, pi];
else
    clrmap = parula;
    clims = [min(snaps(:)), max(snaps(:))];
end

if ~exist('nsnaps', 'var')
    nsnaps = size(snaps, 3);
end
[nx, ny, nvfs] = size(vfx);

% Alternately plot snapshots and vector fields
for isnap = 1:nsnaps
    if plotSeparate
        subplot(2, nsnaps, isnap)
    else
        subplot(1, nsnaps, isnap)
    end
    
    imagesc(snaps(:,:,isnap), clims)
    colormap(clrmap)
    axis xy
    
    if isnap <= nvfs
        if plotSeparate
            if nvfs < nsnaps
                subplot(2, 2*nsnaps, 2*(nsnaps+isnap) + [0 1])
            else
                subplot(2, nsnaps, nsnaps+isnap)
            end
        else
            hold on
        end
        
        quiver(vfx(:,:,isnap), vfy(:,:,isnap))
        if exist('allPatternLocs', 'var') && ~isempty(allPatternLocs)
            hold on
            % Add critical points
            thisTime = allPatternLocs(:,3) == isnap;
            thisLoc = allPatternLocs(thisTime, 1:2);
            scatter(thisLoc(:,2), thisLoc(:,1), 'filled')
        end
        
        xlim([-0.5, ny+1.5])
        ylim([-0.5, nx+1.5])
        xticks([])
        yticks([])
        
        hold off
    end
end