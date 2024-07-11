function [A,b,x,theta,p,R,dw,sd] = ...
    ct_fanlinear(image,N,theta,p,R,dw,sd)
%FANLINEARTOMO  Creates 2D fan-beam linear-detector tomography test problem
%
%   [A,b,x,theta,p,R,d] = fanlineartomo(N)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p,R)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p,R,dw)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p,R,dw,sd)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p,R,dw,sd,isDisp)
%   [A,b,x,theta,p,R,d] = fanlineartomo(N,theta,p,R,dw,sd,isDisp,isMatrix)
%
% This function uses the "line model" to create a 2D X-ray tomography test 
% problem with an N-times-N pixel domain, using p rays in fan formation for 
% each angle in the vector theta. Unlike fancurvedtomo, in which the detector 
% is curved, in fanlineartomo the detector is linear, corresponding to the 
% central slice of a flat-panel detector.
%
% Input:
%   N         Scalar denoting the number of discretization intervals in 
%             each dimesion, such that the domain consists of N^2 cells.
%   theta     Vector containing the angles in degrees. Default: theta =
%             0:2:358.
%   p         Number of rays for each angle. Default: p = round(sqrt(2)*N).
%   R         The distance from the source to the center of the domain
%             is R*N. Default: R = 2.
%   dw        Total width of linear detector. Same units as R.
%   sd        Source to detector distance. Same units as R.
%   isDisp    If isDisp is non-zero it specifies the time in seconds 
%             to pause in the display of the rays. If zero (the default), 
%             no display is shown.
%   isMatrix  If non-zero, a sparse matrix is returned in A (default).
%             If zero, instead a function handle is returned.
%
% Output:
%   A         If input isMatrix is 1 (default): coefficient matrix with
%             N^2 columns and length(theta)*p rows.
%             If input isMatrix is 0: A function handle representing a
%             matrix-free version of A in which the forward and backward
%             operations are computed as A(x,'notransp') and A(y,'transp'),
%             respectively, for column vectors x and y of appropriate size.
%             The size of A can be retrieved using A([],'size'). The matrix
%             is never formed explicitly, thus saving memory.
%   b         Vector containing the rhs of the test problem.
%   x         Vector containing the exact solution, with elements
%             between 0 and 1.
%   theta     Vector containing the used angles in degrees.
%   p         The number of used rays for each angle.
%   R         The radius in side lengths. 
%   d         The span of the rays.
%
% See also: paralleltomo, fancurvedtomo, seismictomo, seismicwavetomo,
%           sphericaltomo.

% Reference: A. C. Kak and M. Slaney, Principles of Computerized
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Default value of sd.
if nargin < 7 || isempty(sd)
    sd = 3;
end
sd = sd*N;

% Default value of dw.
if nargin < 6 || isempty(dw)
    dw = 2.5;
end
dw = dw*N;

% Default value of R.
if nargin < 5 || isempty(R)
    R = 2;
end
R = R*N;

% Default value of the number of rays p.
if nargin < 4 || isempty(p)
    p = round(sqrt(2)*N);
end

% Default value of the angles theta.
if nargin < 3 || isempty(theta)
    theta = 0:2:358;
end

% Make sure theta is double precison to prevent round-off issues caused by
% single input.
theta = double(theta);

% Input check. The source must lie outside the domain.
if R < sqrt(2)/2*N
    error('R must be greater than half squareroot 2')
end

% Construct either matrix or function handle

% Anonymous function rotation matrix
Omega_x = @(omega_par) [cosd(omega_par) -sind(omega_par)];
Omega_y = @(omega_par) [sind(omega_par)  cosd(omega_par)];

% nA denotes the number of angles.
nA = length(theta);

% The starting values both the x and the y coordinates.
x0 = 0;
y0 = R;
xy0 = [x0;y0];

% Width of detector element
dew = dw/p;

% Set angles to match linear detector
if mod(p,2)
    depos = (1:(p-1)/2)*dew;
    deposall = [-depos(end:-1:1),0,depos];
    omega = 90-atand(sd./depos);
    omega = [-omega(end:-1:1),0,omega];
else
    depos = ((1:p/2)-0.5)*dew;
    deposall = [-depos(end:-1:1),depos];
    omega = 90-atand(sd./depos);
    omega = [-omega(end:-1:1),omega];
end

% The intersection lines.
x = (-N/2:N/2)';
y = x;
    
% Initialize vectors that contains the row numbers, the column numbers
% and the values for creating the matrix A efficiently.
rows = zeros(2*N*nA*p,1);
cols = rows;
vals = rows;
idxend = 0;

II = 1:nA;
JJ = 1:p;

% Loop over the chosen angles of the source.
for i = II
    
    % The starting points for the current angle theta.
    x0theta = Omega_x(theta(i))*xy0;
    y0theta = Omega_y(theta(i))*xy0;
    
    % The starting (center) direction vector (opposite xytheta) and
    % normalized.
    xytheta = [x0theta; y0theta];
    abtheta = -xytheta/R;
    
    % Loop over the rays.
    for j = JJ
        
        % The direction vector for the current ray.
        a = Omega_x(omega(j))*abtheta;
        b = Omega_y(omega(j))*abtheta;
        
        xdtheta = Omega_x(theta(i))*[deposall(j); R-sd];
        ydtheta = Omega_y(theta(i))*[deposall(j); R-sd];
        
        % Use the parametrisation of line to get the y-coordinates of
        % intersections with x = k, i.e., x constant.
        tx = (x - x0theta)/a;
        yx = b*tx + y0theta;
        
        % Use the parametrisation of line to get the x-coordinates of
        % intersections with y = k, i.e., y constant.
        ty = (y - y0theta)/b;
        xy = a*ty + x0theta;
        
        % Collect the intersection times and coordinates.
        t = [tx; ty];
        xxy = [x; xy];
        yxy = [yx; y];
        
        % Sort the coordinates according to intersection time.
        [~,I] = sort(t);
        xxy = xxy(I);
        yxy = yxy(I);
        
        % Skip the points outside the box.
        I = (xxy >= -N/2 & xxy <= N/2 & yxy >= -N/2 & yxy <= N/2);
        xxy = xxy(I);
        yxy = yxy(I);
        
        % Skip double points.
        I = (abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
        xxy(I) = [];
        yxy(I) = [];
        
        % Calculate the length within cell and determines the number of
        % cells which is hit.
        aval = sqrt(diff(xxy).^2 + diff(yxy).^2);
        col = [];
        
        % Store the values inside the box.
        if numel(aval) > 0
            
            % Calculates the midpoints of the line within the cells.
            xm = 0.5*(xxy(1:end-1)+xxy(2:end)) + N/2;
            ym = 0.5*(yxy(1:end-1)+yxy(2:end)) + N/2;
            
            % Translate the midpoint coordinates to index.
            col = floor(xm)*N + (N - floor(ym));

        end
        
        if ~isempty(col)
            % Create the indices to store the values to vector for
            % later creation of A matrix.
            idxstart = idxend + 1;
            idxend = idxstart + numel(col) - 1;
            idx = idxstart:idxend;
            
            % Store row numbers, column numbers and values.
            rows(idx) = (i-1)*p + j;
            cols(idx) = col;
            vals(idx) = aval;
        end
        
    end % end j
    
end % end i

% Truncate excess zeros.
rows = rows(1:idxend);
cols = cols(1:idxend);
vals = vals(1:idxend);

% Create sparse matrix A from the stored values.
A = sparse(rows,cols,vals,p*nA,N^2);

if nargout > 2
    % Create phantom head as a reshaped vector.
    x = image;
    x = x(:);
    % Create rhs.
    b = A*x;
end

if nargout > 5
    R = R/N;
end