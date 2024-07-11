clc;
clear;

addpath('AIRToolsII/')
AIRToolsII_setup('temporary')

N_pixels = 225;
theta_inc = 1;
theta = 0:theta_inc:360;

seed = 1;
eta = 0.008;

imdata = imread('brain.jpeg');
% imdata = imread('AIRToolsII/testprobs/head.png');
im = rgb2gray(imdata);
imgray = im2double(im);

%%

[A,b,x] = ct_image(imgray, N_pixels,theta);

x_init = x;

% Build RHS and print matrix parameters

M = size(A,1)
N = size(A,2)

% Error

rng(seed);
sigma = sqrt(eta*eta*norm(b)*norm(b)/M);
error = normrnd(0,sigma,M,1);
b_error = b + error;

% Indices of columns to be deleted

is_pixel_zero = zeros(N,1);
for i=1:M
    if b_error(i) == 0
        for j=1:N
            if full(A(i,j)) ~= 0
                is_pixel_zero(j) = 1;
            end
        end
    end
end
idx_non_zero_pixels = [];
for j=1:N
    if is_pixel_zero(j) == 0
        idx_non_zero_pixels(end+1) = j-1;
    end
end

% Delete cols

idx = any(is_pixel_zero ~= 1, 2);
A = A(:,idx);
x = x(idx);
N = size(A,2);

% Delete rows with b values set to zero

idx = any(vecnorm(A,2,2) ~= 0, 2);
b = b(idx,:);
b_error = b_error(idx,:);
A = A(idx,:);
M = size(A,1);

% Delete rows of zeros

idx = any(vecnorm(A,2,2) ~= 0, 2);
b = b(idx,:);
A = A(idx,:);

M = size(A,1)
N = size(A,2)

x_ls = lsqr(A,b_error,1E-18,10000);
fprintf("CGLS solution compleated!\n");

x_ls_full = zeros(1,N_pixels*N_pixels);
for i=1:N
    x_ls_full(idx_non_zero_pixels(i)+1) = x_ls(i);
end

max_pixel = max(max(x_init), max(x_ls));
min_pixel = min(min(x_init), min(x_ls));

figure_init = figure(1);
imagesc(reshape(x_init,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_pixel max_pixel])
filename_fig = "ct_real/x_" + int2str(M) + "_" + int2str(N) + ".png";
saveas(figure_init,filename_fig);

figure_ls = figure(2);
imagesc(reshape(x_ls_full,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_pixel max_pixel])
filename_fig = "ct_real/x_ls_" + int2str(M) + "_" + int2str(N) + ".png";
saveas(figure_ls,filename_fig);

% CSR

NNZ = size(nonzeros(A),1);

row_idx = zeros(M+1,1);
row_idx(1) = 0;
cols = zeros(NNZ,1);
values = zeros(NNZ,1);

counter = 0;
element_per_row = 0;
for i=1:M
    for j=1:N
        if A(i,j) ~= 0
            element_per_row = element_per_row + 1;
            counter = counter + 1;
            cols(counter,1) = j-1;
            values(counter,1) = A(i,j);
        end
    end
    row_idx(i+1) = row_idx(i) + element_per_row;
    element_per_row = 0;
end

filename_idx_non_zero_pixels = "ct_real/idx_non_zero_pixels_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_idx_non_zero_pixels = fopen(filename_idx_non_zero_pixels,'w');
fwrite(file_idx_non_zero_pixels, N_pixels*N_pixels, 'int');
fwrite(file_idx_non_zero_pixels, idx_non_zero_pixels, 'int');
fclose(file_idx_non_zero_pixels);

filename_is_pixel_zero = "ct_real/is_pixel_zero_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_is_pixel_zero = fopen(filename_is_pixel_zero,'w');
fwrite(file_is_pixel_zero, N_pixels*N_pixels, 'int');
fwrite(file_is_pixel_zero, is_pixel_zero, 'int');
fclose(file_is_pixel_zero);

filename_error = "ct_real/error_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_error = fopen(filename_error,'w');
fwrite(file_error, error, 'double');
fclose(file_error);

filename_x = "ct_real/x_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_x = fopen(filename_x,'w');
fwrite(file_x, x, 'double');
fclose(file_x);

filename_x_ls = "ct_real/x_ls_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_x_ls = fopen(filename_x_ls,'w');
fwrite(file_x_ls, x_ls, 'double');
fclose(file_x_ls);

filename_b = "ct_real/b_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_b = fopen(filename_b,'w');
fwrite(file_b, b, 'double');
fclose(file_b);

filename_b_error = "ct_real/b_error_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_b_error = fopen(filename_b_error,'w');
fwrite(file_b_error, b_error, 'double');
fclose(file_b_error);

filename_row_idx = "ct_real/row_idx_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_row_idx = fopen(filename_row_idx,'w');
fwrite(file_row_idx, row_idx, 'int');
fclose(file_row_idx);

filename_cols = "ct_real/cols_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_cols = fopen(filename_cols,'w');
fwrite(file_cols, NNZ, 'int');
fwrite(file_cols, cols, 'int');
fclose(file_cols);

filename_values = "ct_real/values_" + int2str(M) + "_" + int2str(N) + "_" + int2str(seed) + ".bin";
file_values = fopen(filename_values,'w');
fwrite(file_values, NNZ, 'int');
fwrite(file_values, values, 'double');
fclose(file_values);

function [A,b,x,theta,p,d] = ct_image(image,N,theta,p,d)
    %PARALLELTOMO  Creates a 2D parallel-beam tomography test problem
    %
    %   [A,b,x,theta,p,d] = paralleltomo(N)
    %   [A,b,x,theta,p,d] = paralleltomo(N,theta)
    %   [A,b,x,theta,p,d] = paralleltomo(N,theta,p)
    %   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d)
    %   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp)
    %   [A,b,x,theta,p,d] = paralleltomo(N,theta,p,d,isDisp,isMatrix)
    %
    % This function uses the "line model" to create a 2D X-ray tomography test
    % problem with an N-times-N pixel domain, using p parallel rays for each
    % angle in the vector theta.
    %
    % Input:
    %   N         Scalar denoting the number of discretization intervals in
    %             each dimesion, such that the domain consists of N^2 cells.
    %   theta     Vector containing the projetion angles in degrees.
    %             Default: theta = 0:1:179.
    %   p         Number of rays for each angle. Default: p = round(sqrt(2)*N).
    %   d         Scalar denoting the distance from the first ray to the last.
    %             Default: d = p-1.
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
    %             elements in A.
    %   b         Vector containing the rhs of the test problem.
    %   x         Vector containing the exact solution, with elements
    %             between 0 and 1.
    %   theta     Vector containing the used angles in degrees.
    %   p         The number of used rays for each angle.
    %   d         The distance between the first and the last ray.
    %
    % See also: fancurvedtomo, fanlineartomo, seismictomo, seismicwavetomo,
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
    
    % Default value of the projection angles theta.
    if nargin < 3 || isempty(theta)
        theta = 0:179;
    end
    
    % Make sure theta is double precison to prevent round-off issues caused by
    % single input.
    theta = double(theta);
    
    % Default value of the number of rays.
    if nargin < 4 || isempty(p)
        p = round(sqrt(2)*N);
    end
    
    % Default value of d.
    if nargin < 5 || isempty(d)
        d = p-1;
    end
    
    % Define the number of angles.
    nA = length(theta);
    
    % Radii for the circles.
    radii  = linspace(0,2,p+1);
    radii  = radii(2:end);
    
    % Image coordinates.
    centerImg = ceil(N/2);
    
    % Determine the quarature parameters.
    dx = sqrt(2)/N;
    nPhi = ceil((4*pi/dx)*radii);
    dPhi = 2*pi./nPhi;
        
    % Initialize vectors that contains the row numbers, the column numbers
    % and the values for creating the matrix A effiecently.
    rows = zeros(2*N*nA*p,1);
    cols = rows;
    vals = rows;
    idxend = 0;
    
    II = 1:nA;
    JJ = 1:p;
    
    % Loop over theta.
    for m = II
        
        % Angular position of source.
        xix = cosd(theta(m));
        xiy = sind(theta(m));
        
        % Loop over the circles.
        for n = JJ
            
            % (x,y) coordinates of circle.
            k = (0:nPhi(n))*dPhi(n);
            xx = (xix + radii(n)*cos(k))/dx + centerImg;
            yy = (xiy + radii(n)*sin(k))/dx + centerImg;
            
            % Round to get pixel index.
            col = round( xx );
            row = round( yy );
            
            % Discard if outside box domain.
            IInew = find(col>0 & col<=N & row>0 & row<=N);
            row = row(IInew);
            col = col(IInew);
            J = (N-row+1) + (col-1)*N;
            
            % Convert to linear index and bin
            Ju = unique(J);
            w = histc(J,Ju);
            
            % Determine rows, columns and weights.
            i = (m-1)*p + n;
            ii = repmat(i,length(Ju),1);
            jj = Ju(:);
            aa = (2*pi*radii(n)/nPhi(n))*w(:);
            
            % Store the values, if any.
            if ~isempty(jj)
                % Create the indices to store the values to vector for
                % later creation of A matrix.
                idxstart = idxend + 1;
                idxend = idxstart + numel(jj) - 1;
                idx = idxstart:idxend;
    
                % Store row numbers, column numbers and values.
                rows(idx) = ii;
                cols(idx) = jj;
                vals(idx) = aa;
            end
        end
    end
    
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

end