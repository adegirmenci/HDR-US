function plotHistogram(image, title_, saveResults, workingDir)
%PLOTHISTOGRAM Plot the histogram
%
% Example:
%   plotHistogram(stack(:,:,1,1), 'Power: -26.7 dB', saveResults);
%
% Copyright (C) 2016-2017, by Alperen Degirmenci
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Histogram edges
edges = linspace(0,1,256);
edges = edges - edges(1)/2;
% yMax = round(numel(image)/256);
yMax = 1400;
% yMax = floor(log(numel(image)));

figure('units','normalized','position',[.25 .25 .45 .25])
hold on
ax = histogram(image, edges,...
    'Normalization','count',...
    'FaceColor', 'k',...
    'FaceAlpha', 1);
xlim([-0.01,1.01])
ylim([0,yMax])
axis off
legend(title_)

vals = ax.Values;
% Convert to logarithmic scale
%vals = log(ax.Values);
% Set max value to yMax
vals(vals > yMax) = yMax;
% Initialize histogram image
histImg = ones(128, 256);
[height, ~] = size(histImg);
% Fill out histogram
for i = 1:numel(ax.Values)
    val = vals(i);
    pctg = val/yMax;
    nPix = round(height*pctg);
    histImg(1:(height-nPix),i) = 0;
end
figure
imshow(histImg)

if(saveResults)
    fileName = [workingDir,filesep,'logHist_',title_,'.png'];
    fprintf('Saving histogram to:\n\t%s\n',fileName)
    imwrite(histImg, fileName)
end

end