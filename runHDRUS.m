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

close all; clear all; clc

saveResults = false;

%% Step 0: Choose Working Directory
cd ..
root = pwd; % current directory
cd source
workingDir = uigetdir([root,filesep,'data']); % This opens a dialog
% workingDir = [pwd,filesep,'Data',filesep,'exVivo']; % Alternatively, set working dir by hand
relativePath = workingDir(length(root)+2:end);
study = strsplit(relativePath,filesep);
study = study{end};

%% Step 1: Load Image Data

[FileName,PowerdB] = importPowerList([workingDir,filesep,study,'_power.txt']);

numImg = length(FileName);
numPower = length(unique(PowerdB));

fprintf('Loaded %d entries from %s_power.txt\n', numImg, study);

images(numImg) = struct('image',0.0,...
                        'path',' ',...
                        'power',0.0);

% load images
for i = 1:numImg
    images(i).path = ['.',filesep,relativePath,filesep,FileName{i}];
    images(i).power = PowerdB(i);
    images(i).image = imread(images(i).path);
end

exps_ = 10.0.^(([images(:).power])/20.0); % convert dB to power ratio

%%

%      Use HDR Toolbox to create an HDR radiance map for ultrasound:
%	   1) Read a stack of LDR images
%	   2) Read exposure values from the EXIF
%	   3) Estimate the Camera Response Function (CRF)
%	   4) Build the radiance map using the stack and stack_exposure
%	   5) Save the radiance map in .hdr format (optional)
%	   6) Show the tone mapped version of the radiance map

%% 1
format = 'jpg';

disp('1) Read a stack of LDR images');
[stack, norm_value] = ReadLDRStack(workingDir, format, 1);

%% 2
disp('2) Read exposure values');
%stack_exposure = ReadLDRStackInfo(name_folder, format); % from exif
stack_exposure = exps_;

%%
disp('3) Estimate the Response Curve');
%[lin_fun, ~] = DebevecCRF(stack, stack_exposure);
[lin_fun, ~, lEGray, zGray] = DebevecCRF(stack, stack_exposure);

h = figure(1);
set(h, 'Name', 'Response Curve');
plot(x,0:255);
semilogx(lin_fun(:,1), 0:255, 'r', lin_fun(:,2), 0:255, 'g', lin_fun(:,3), 0:255, 'b');

gGray = log(lin_fun(:,1));
zGray = squeeze(zGray(:,:,1));
Xposure = log(exp(lEGray(:,1))*stack_exposure);

figure('units','normalized','position',[.25 .25 .45 .3])
plot(Xposure(:),zGray(:),'.','Color',ones(1,3)*.6);
hold on
plot(gGray,0:255,'k','LineWidth',3);
xticks(linspace(-10,10,5))
yticks(round(linspace(0,255,5)))
box on; grid on;
xlim([-5,5])
ylim([-5,260])
xlabel('Log Pressure')
ylabel('Pixel Intensity')
legend('Data','Recovered','Location','northwest')
set(gca,'fontsize',16)

% % export_fig responseFunc -pdf -depsc

%%
disp('4) Build the radiance map using the stack and stack_exposure');
imgHDR = BuildHDR(stack, stack_exposure, 'LUT', lin_fun, 'Robertson', 'log', 1);

% disp('5) Save the radiance map in the .hdr format');
% hdrimwrite(imgHDR, 'stack_hdr_image.exr');

%% Show HDR luminance

figure('units','pixels','position',[400 200 size(imgHDR,2)+75 size(imgHDR,1)])
hdr1 = squeeze(imgHDR(:,:,1));

logHDR = log10(hdr1);
c = logspace(min(logHDR(:)),max(logHDR(:)),7);

imagesc(logHDR);
hold on
axis equal; colormap jet; % hsv, pink
caxis(log10([c(1) c(end)]));
colorbar('FontSize',11,'YTick',log10(c),'YTickLabel',round(c,2),'TickLabelInterpreter','latex');
% colorbar('FontSize',11,'YTick',c,'YTickLabel',c);
axis([0 size(imgHDR,2), 0 size(imgHDR,1)])
set(gca,'XTick',[])
set(gca,'YTick',[])

% saveas(gcf,['.',filesep,relativePath,filesep,'luminanceMap_wColorbar'],'png')
% export_fig(['.',filesep,relativePath,filesep,'luminanceMap_wColorbar'],'-png');

%%

disp('6) Show the tone mapped version of the radiance map with gamma encoding');
% h = figure;
% set(h, 'Name', 'Tone mapped version of the built HDR image');

imgOutReinhardGlobal = ReinhardTMO(double(imgHDR), 0.095, 0, 'global');
imgOutReinhardLocal = ReinhardTMO(double(imgHDR), 0.095, 0, 'local');
imgOutReinhardBilateral = ReinhardTMO(double(imgHDR), 0.095, 0, 'bilateral');
imgOutDurand = DurandTMO(double(imgHDR),24);
imgOutAdaptHist = tonemap_(imgHDR, 'AdjustLightness', [0.1 1.0], ...
                                  'AdjustSaturation', 0.6, ...
                                  'NumberOfTiles', [8,8]); % CLAHE
imgOutAdaptHist = ClampImg(imgOutAdaptHist, 0, 1);

% Gamma correction

imgOutReinhardGlobal = GammaTMO(imgOutReinhardGlobal, 2.2, 0.0, 0);
imgOutReinhardLocal = GammaTMO(imgOutReinhardLocal, 2.2, 0.0, 0);
imgOutReinhardBilateral = GammaTMO(imgOutReinhardBilateral, 2.2, 0.0, 0);
imgOutDurand = GammaTMO(imgOutDurand, 2.2, 0.0, 0);
imgOutAdaptHist = GammaTMO(imgOutAdaptHist, 1.0, 0.0, 0);

% Rescale

imgOutReinhardGlobal = rescale(imgOutReinhardGlobal);
imgOutReinhardLocal = rescale(imgOutReinhardLocal);
imgOutReinhardBilateral = rescale(imgOutReinhardBilateral);
imgOutDurand = rescale(imgOutDurand);
imgOutAdaptHist = rescale(imgOutAdaptHist);

% imgOutReinhardGlobal = (imgOutReinhardGlobal - min(imgOutReinhardGlobal(:))) / (max(imgOutReinhardGlobal(:)) - min(imgOutReinhardGlobal(:)));
% imgOutReinhardLocal = (imgOutReinhardLocal - min(imgOutReinhardLocal(:))) / (max(imgOutReinhardLocal(:)) - min(imgOutReinhardLocal(:)));
% imgOutReinhardBilateral = (imgOutReinhardBilateral - min(imgOutReinhardBilateral(:))) / (max(imgOutReinhardBilateral(:)) - min(imgOutReinhardBilateral(:)));
% imgOutDurand = (imgOutDurand - min(imgOutDurand(:))) / (max(imgOutDurand(:)) - min(imgOutDurand(:)));
% imgOutAdaptHist = (imgOutAdaptHist - min(imgOutAdaptHist(:))) / (max(imgOutAdaptHist(:)) - min(imgOutAdaptHist(:)));

% figure; imshow(imgOutReinhardGlobal);
% figure; imshow(imgOutReinhardLocal)
% figure; imshow(imgOutReinhardBilateral)
% figure; imshow(imgOutDurand)
% figure; imshow(imgOutAdaptHist)

%%
if(saveResults)
    imwrite(imgOutReinhardGlobal, [workingDir,filesep,'hdrReinhard_global','.png']);
    imwrite(imgOutReinhardLocal, [workingDir,filesep,'hdrReinhard_local','.png']);
    imwrite(imgOutReinhardBilateral, [workingDir,filesep,'hdrReinhard_bilateral','.png']);
    imwrite(imgOutDurand, [workingDir,filesep,'hdrDurand','.png']);
    imwrite(imgOutAdaptHist, [workingDir,filesep,'hdrAdaptHist.png']);
end
