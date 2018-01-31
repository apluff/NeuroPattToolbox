function [hFigsSorted] = sortfighandlearray(hFigs,varargin)
%SORTFIGHANDLEARRAY Due to changes in Matlab 2014b graphics system, figure
%handles are no longer doubles, but rather an graphics object. Hence, these
%objects can not be sorted directly. This function accepts an array of
%figure handles, and returns an array of similar length, with figure
%handles sorted with respect to their numeric property 'Number'.
%
% *** Inputs:
%   hFigs - Figure handles array. 
%           An array who holds the handles to a series of figures.
%
%   varargin{1} - String.
%           Should match 'ascend' or 'descend', depending on  how the
%           figure handles should be sorted. By default, figure handles are
%           sorted in an ascending order.
%
% *** Outputs:
%   hFigSorted- Figure handles array.
%           Figure handles are sorted with respect to their property
%           'Number'.
% ------------------------------------------------------------------------
%
% By Andreas Berg, from MATLAB Central

narginchk(1,2)

if nargin == 1
   sortStr = 'ascend';
else
   if strcmpi(varargin{1},'ascend') || strcmpi(varargin{1},'descend')
      sortStr = varargin{1};
   else
      error('Bad input. 2nd input argument should match string "ascend" or "descend"')
   end
end

% Allocate output
nFigs = length(hFigs);
hFigsSorted = gobjects(nFigs,1);

% Capture an array of figure numbers
figNumbersArray = [hFigs.Number];

% Sort and find indices of the sorted figures
[~,I] = sort(figNumbersArray,sortStr);

% Assemble output array
hFigsSorted(1:nFigs) = hFigs(I);

end

