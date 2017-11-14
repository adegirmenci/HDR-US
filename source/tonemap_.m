function RGBldr = tonemap_(RGBhdr, varargin)
%TONEMAP   Render high dynamic range image for viewing.
%   RGB = TONEMAP(HDR) performs tone mapping on the high dynamic range
%   image HDR to a lower dynamic range image RGB suitable for display.
%
%   RGB = TONEMAP(HDR, 'AdjustLightness', [LOW HIGH], ...) adjusts the
%   overall lightness of the rendered image by passing the luminance
%   values of the low dynamic range image to IMADJUST with the values LOW
%   and HIGH, which are in the range [0, 1].
%
%   RGB = TONEMAP(HDR, 'AdjustSaturation', SCALE, ...) adjusts the
%   saturation of colors in the rendered image.  When SCALE is greater
%   than 1, the colors are more saturated.  A SCALE value in the range
%   [0, 1) results in less saturated colors.
%
%   RGB = TONEMAP(HDR, 'NumberOfTiles', [ROWS COLS], ...) sets the number
%   of tiles used during the adaptive histogram equalization part of the
%   tone mapping operation.  ROWS and COLS specify the number of tile rows
%   and columns.  Both ROWS and COLS must be at least 2.  The total number
%   of image tiles is equal to ROWS * COLS.  A larger number of tiles
%   results in an image with greater local contrast.  The default for ROWS
%   and COLS is 4.
%
%   Class Support
%   -------------
%   The high dynamic range image HDR must be a m-by-n-by-3 single or
%   double array.  The output image RGB is an m-by-n-by-3 uint8 image.
%
%   Example
%   -------
%   Load a high dynamic range image, convert it to a low dynamic range
%   image while deepening shadows and increasing saturation, and display
%   the results.
%
%       hdr = hdrread('office.hdr');
%       imshow(hdr)
%       rgb = tonemap(hdr, 'AdjustLightness', [0.1 1], ...
%                     'AdjustSaturation', 1.5);
%       figure;
%       imshow(rgb);
%
%   See also ADAPTHISTEQ, HDRREAD, STRETCHLIM.

%   Copyright 2007-2010 The MathWorks, Inc.

% Parse and validate input arguments.
validateattributes(RGBhdr, {'single', 'double'}, {'real'}, mfilename, 'RGBhdr', 1);
if ((ndims(RGBhdr) ~= 3) || (size(RGBhdr, 3) ~= 3))
    error(message('images:tonemap:badImageDimensions'))
end
options = parseArgs(varargin{:});

% Transform the HDR image to a new HDR image in the range [0,1] by taking
% the base-2 logarithm and linearly scaling it.
[RGBlog2Scaled, hasNonzero] = lognormal(RGBhdr);

% Convert the image to a low dynamic range image by adaptive histogram
% equalization.
if (hasNonzero)
    RGBldr = toneOperator(RGBlog2Scaled, ...
                          options.LRemap, ...
                          options.saturation, ...
                          options.numtiles);
else
    % "HDR" image only has zeros.  Return another image of zeros.
    RGBldr = RGBlog2Scaled;
end

%RGBldr = uint8(RGBldr * 255);



function options = parseArgs(varargin)
% Get user-provided and default options.

% Create a structure with default values, and map actual param-value pair
% names to convenient names for internal use.
knownParams = {'AdjustLightness',  'LRemap',     [0 1], {'nonempty', 'vector', 'real', 'nonnan', 'positive'};
               'AdjustSaturation', 'saturation',     1, {'scalar', 'real', 'nonnan', 'positive'};
               'NumberOfTiles',    'numtiles',   [4 4], {'nonempty', 'vector', 'integer', 'real', 'finite', 'positive', 'nonzero'}};
options = cell2struct(knownParams(:,3), knownParams(:,2), 1);

if (rem(nargin, 2) ~= 0)
    error(message('images:tonemap:paramValuePairs'))
end

% Loop over the P-V pairs.
for p = 1:2:nargin
    % Get the parameter name.
    paramName = varargin{p};
    if (~ischar(paramName))
        error(message('images:tonemap:badParamName'))
    end
    
    % Look for the parameter amongst the possible values.
    idx = strmatch(lower(paramName), lower(knownParams(:,1)));
    if (isempty(idx))
        error(message('images:tonemap:unknownParamName', paramName));
    elseif (numel(idx) > 1)
        error(message('images:tonemap:ambiguousParamName', paramName));
    end

    % Validate the value.
    options.(knownParams{idx, 2}) = varargin{p+1};
    validateattributes(varargin{p+1}, ...
                  {'double'}, ...
                  knownParams{idx,4}, ...
                  mfilename, ...
                  knownParams{idx,1}, ...
                  p+2);  % p+2 = Param name + 1 + offset to first arg
end



function [RGBlog2Scaled, hasNonzero] = lognormal(RGBhdr)
% Take the base-2 logarithm of an HDR image and return another HDR in [0,1].

% Remove 0's from each channel.  This can change color quality, but it's
% unlikely to have a big impact and prevents log(0) --> -inf.  That's worse.
minNonzero = min(RGBhdr(RGBhdr ~= 0));

if (isempty(minNonzero))
    
    RGBlog2Scaled = zeros(size(RGBhdr), class(RGBhdr));
    hasNonzero = false;
    
else
    
    RGBhdr(RGBhdr == 0) = minNonzero;

    % Ward's method equalizes the log-luminance histogram.
    RGBlog2 = log2(RGBhdr);
    RGBlog2Scaled = convertToImage(RGBlog2); % Normalize to [0,1]
    hasNonzero = true;
    
end



function RGBldr = toneOperator(RGBlog2Scaled, LRemap, saturation, numtiles)
% Convert the image from HDR to LDR.

% Colorspaces for HDR imagery is tricky.  For simplicity, assign the
% log-luminance image to be in sRGB.
Lab = sRGB2Lab(RGBlog2Scaled);

% Tone map the L* values from the RGB HDR to preserve overall color as much
% as possible.  This decreases global saturation, which can be reintroduced
% by scaling the a* and b* channels.
Lab(:,:,1) = Lab(:,:,1) ./ 100;
Lab(:,:,1) = adapthisteq(Lab(:,:,1), 'NumTiles', numtiles);
Lab(:,:,1) = imadjust(Lab(:,:,1), LRemap, [0 1]) * 100;
Lab(:,:,2) = Lab(:,:,2) * saturation;
Lab(:,:,3) = Lab(:,:,3) * saturation;

% Convert the image back to sRGB.
RGBldr = Lab2sRGB(Lab);



function y = convertToImage(x)
% Rescale an image with arbitrary range to [0,1].

xMin = min(x(:));
xMax = max(x(:));

if (xMin == xMax)
    
    % Avoid dividing by zero.
    if (xMin == 0)
        y = x;
    else
        y = x ./ xMin;
    end
    
else
    % Linearly scale the values to [0,1].
    y = (x - xMin) ./ (xMax - xMin);
end



function Lab = sRGB2Lab(rgb)
% Convert sRGB values in the range [0,1] to Lab via XYZ.

dims = size(rgb);
rgb = reshape(permute(rgb, [3 1 2]), [3, dims(1) * dims(2)]);
Lab = XYZ2Lab(sRGB2XYZ(rgb));
Lab = permute(reshape(Lab, [3, dims(1), dims(2)]), [2 3 1]);



function rgb = Lab2sRGB(lab)
% Convert Lab values to sRGB values in the range [0,1] via XYZ.

dims = size(lab);
lab = reshape(permute(lab, [3 1 2]), [3, dims(1) * dims(2)]);
rgb = XYZ2sRGB(Lab2XYZ(lab));
rgb = permute(reshape(rgb, [3, dims(1), dims(2)]), [2 3 1]);


function xyz = sRGB2XYZ(rgb)
% Convert sRGB values to XYZ assuming a D65 whitepoint.

% See <http://brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html> for
% implementation details.
M = [0.412424    0.212656    0.0193324  
     0.357579    0.715158    0.119193
     0.180464    0.0721856   0.950444];

mask = (rgb > 0.04045);

rgb(mask) = ((rgb(mask) + 0.055) ./ 1.055) .^ 2.4;
rgb(~mask) = rgb(~mask) ./ 12.92;

xyz = M' * rgb;


function lab = XYZ2Lab(xyz)
% Convert XYZ to Lab using a D65 whitepoint.

% See <http://brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html> for
% implementation details and explanation of E and K.
E = 216/24389;
K = 24389/27;

% Whitepoint adjustment.
xyz(1,:) = xyz(1,:) ./ 0.95047;
xyz(3,:) = xyz(3,:) ./ 1.08883;

mask = (xyz > E);

f = xyz;
f(mask) = f(mask) .^ (1/3);
f(~mask) = (K * f(~mask) + 16) ./ 116;

lab = [116 * f(2,:) - 16;
       500 * (f(1,:) - f(2,:));
       200 * (f(2,:) - f(3,:))];



function xyz = Lab2XYZ(lab)
% Convert Lab to XYZ using a D65 whitepoint.

% See <http://brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html> for
% implementation details and explanation of E and K.
E = 216/24389;
K = 24389/27;

f = lab;
xyz = f;

mask = lab(1,:) > K*E;
xyz(2,mask) = ((lab(1,mask) + 16) ./ 116) .^ 3;
xyz(2,~mask) = lab(1,~mask) ./ K;

mask = (xyz(2,:) > E);
f(2,mask) = (lab(1,mask) + 16) ./ 116;
f(2,~mask) = (K * xyz(2,~mask) + 16) ./ 116;

f(1,:) = lab(2,:) ./ 500 + f(2,:);
f(3,:) = f(2,:) - lab(3,:) ./ 200;

tmp = f(1,:) .^ 3;
mask = (tmp > E);
xyz(1,mask) = tmp(mask);
xyz(1,~mask) = (116 * f(1,~mask) - 16) ./ K;

tmp = f(3,:) .^ 3;
mask = (tmp > E);
xyz(3,mask) = tmp(mask);
xyz(3,~mask) = (116 * f(3,~mask) - 16) ./ K;

% Whitepoint adjustment.
xyz(1,:) = xyz(1,:) * 0.95047;
xyz(3,:) = xyz(3,:) * 1.08883;



function rgb = XYZ2sRGB(xyz)
% Convert XYZ to sRGB using a D65 whitepoint.

% See <http://brucelindbloom.com/index.html?Eqn_XYZ_to_RGB.html> for
% implementation details.
M = [0.412424    0.212656    0.0193324  
     0.357579    0.715158    0.119193
     0.180464    0.0721856   0.950444];

rgb = M' \ xyz;

mask = (rgb > 0.0031308);
rgb(mask) = ((1.055 * rgb(mask)) .^ (1 / 2.4)) - 0.055;
rgb(~mask) = 12.92 * rgb(~mask);
