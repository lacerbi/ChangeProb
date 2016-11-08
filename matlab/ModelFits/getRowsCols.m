function [R,C] = getRowsCols(N,longDim)
% [R,C] = GETROWSCOLS(N) is a function to determine the number of rows and
% columns for creating a figure with subplots that best fits a square 
% plotting area. 'N' is the number of plots, 'R' is the number of rows, and
% 'C' is the number of columns. The optional argument 'longDim' selects the
% dimension with the greater number of plots if plots do not perfectly fit
% into a square (enter 'R' for rows or 'C' for columns). By default, this 
% is columns so width > height.
%
% Created by SML Aug 2016

% Defaults:
if nargin < 2
    longDim = 'C';
end

% Checks:
assert((N>0)&(mod(N,1)==0), 'N needs to be a postive integer.')
assert((longDim=='C')|(longDim=='R'), 'longDim argument must be either R or C.')

% Get dimensions of plot:
d1 = ceil(sqrt(N));
d2 = ceil(N/d1);

% Assign to rows/columns:
switch longDim
    case 'C' % cols > rows
        R = min(d1,d2);
        C = max(d1,d2);
    case 'R' % rows > cols
        R = max(d1,d2);
        C = min(d1,d2);
end

end