% Copyright (C) 2016-2018, by Alperen Degirmenci
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

close all; clear all; clc

saveResults = false;

%% Step 0: Choose Working Directory
root = pwd; % current directory
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
% !! Make sure that the file order from ReadLDRStack matches expsure order
% For example *_1,jpg, *_2.jpg, ... *_11.jpg will be ordered 1,11,2, so
% renaming to 01,02,...11 would mitigate this.

%%
disp('3) Estimate the Response Curve');
%[lin_fun, ~] = DebevecCRF(stack, stack_exposure);
[lin_fun, ~, lEGray, zGray] = DebevecCRF(stack, stack_exposure);

gGray = log(lin_fun(:,1));
zGray = squeeze(zGray(:,:,1));
Xposure = log(exp(lEGray(:,1))*stack_exposure);

figure('units','normalized','position',[.25 .25 .45 .25])
plot(Xposure(:),zGray(:),'.','Color',ones(1,3)*.55);
hold on
plot(gGray,0:255,'k','LineWidth',3);
box on; grid on;
xticks(linspace(-10,10,5))
yticks(round(linspace(0,255,3)))
xlim([-5,5])
ylim([-5,260])
xlabel('Log Pressure')
ylabel('Pixel Intensity')
legend('Sampled Pixels','Recovered Curve','Location','northwest')
ax = gca;
ax.LineWidth = 1;
ax.FontName = 'Times New Roman';
set(ax,'FontSize',16)

if(saveResults)
    % Save as Variable
    save([workingDir,filesep,'responseFunction.mat'],'lin_fun')
    
    % Save as PDF
    set(gcf,'units','inch')
    ppos = get(gcf, 'Position');
    set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
    print([workingDir,filesep,'responseFunc.pdf'],'-dpdf')
    
    % % If you have the export_fig package, this is easier to use
    % export_fig([workingDir,filesep,'responseFunc'],'-pdf','-depsc');
end

% % To reuse a previously calculated response function:
% load([workingDir,filesep,'responseFunction.mat'])

%%
disp('4) Build the radiance map using the stack and stack_exposure');
imgHDR = BuildHDR(stack, stack_exposure, 'LUT', lin_fun, 'Robertson', 'log', 1);

if(saveResults)
    disp('5) Save the radiance map in the .hdr format');
    hdrimwrite(imgHDR, [workingDir,filesep,'stack_hdr_image.hdr']);
end

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
title('HDR Image')

% saveas(gcf,[workingDir,filesep,'luminanceMap_wColorbar'],'png')
% export_fig([workingDir,filesep,'luminanceMap_wColorbar'],'-png');

%% Tone Mapping

disp('6) Show the tone mapped version of the radiance map with gamma encoding');

nTMOs = 6;

hdrAdaptHist = tonemap_(imgHDR, 'AdjustLightness', [0.1 1.0], ...
                                  'AdjustSaturation', 0.6, ...
                                  'NumberOfTiles', [8,8]); % CLAHE
hdrAdaptHist = ClampImg(hdrAdaptHist, 0, 1);

hdrSimpleCompression = (imgHDR-min(imgHDR(:)))./(max(imgHDR(:)-min(imgHDR(:))));
hdrReinhardGlobal = ReinhardTMO(double(imgHDR), 0, 0, 'global');
hdrReinhardLocal = ReinhardTMO(double(imgHDR), 0, 0, 'local');
hdrReinhardBilateral = ReinhardTMO(double(imgHDR), 0, 0, 'bilateral');
hdrDurand = DurandTMO(double(imgHDR),22);

% Some other TMO options
% hdrWardHistAdj = WardHistAdjTMO(double(imgHDR));
% hdrTumblin = TumblinTMO(double(imgHDR));
% hdrRDevlin = ReinhardDevlinTMO(double(imgHDR));
% hdrRaman = RamanTMO([], stack);
% hdrMertens = MertensTMO([], stack);
% hdrBruce = BruceExpoBlendTMO([], stack, 29, 6);

% Gamma correction
hdrAdaptHist = GammaTMO(hdrAdaptHist, 1.0, 0.0, 0);
gamma = 2.2;
hdrSimpleCompression = GammaTMO(hdrSimpleCompression, gamma, 0.0, 0);
hdrReinhardGlobal = GammaTMO(hdrReinhardGlobal, gamma, 0.0, 0);
hdrReinhardLocal = GammaTMO(hdrReinhardLocal, gamma, 0.0, 0);
hdrReinhardBilateral = GammaTMO(hdrReinhardBilateral, gamma, 0.0, 0);
hdrDurand = GammaTMO(hdrDurand, gamma, 0.0, 0);

% Rescale
hdrAdaptHist = rescale(hdrAdaptHist);
hdrSimpleCompression = rescale(hdrSimpleCompression);
hdrReinhardGlobal = rescale(hdrReinhardGlobal);
hdrReinhardLocal = rescale(hdrReinhardLocal);
hdrReinhardBilateral = rescale(hdrReinhardBilateral);
hdrDurand = rescale(hdrDurand);

% Grayscale
hdrAdaptHist = rgb2gray(hdrAdaptHist);
hdrSimpleCompression = rgb2gray(hdrSimpleCompression);
hdrReinhardGlobal = rgb2gray(hdrReinhardGlobal);
hdrReinhardLocal = rgb2gray(hdrReinhardLocal);
hdrReinhardBilateral = rgb2gray(hdrReinhardBilateral);
hdrDurand = rgb2gray(hdrDurand);

%% Display Results

figure('units','pixels','position',[200 200 size(imgHDR,2)*6 size(imgHDR,1)])
subplot(1,6,1); imshow(hdrAdaptHist); title('Adaptive Hist Eq')
subplot(1,6,2); imshow(hdrSimpleCompression); title('Simple Compression')
subplot(1,6,3); imshow(hdrReinhardGlobal); title('Reinhard Global')
subplot(1,6,4); imshow(hdrReinhardLocal); title('Reinhard Local')
subplot(1,6,5); imshow(hdrReinhardBilateral); title('Reinhard Bilateral')
subplot(1,6,6); imshow(hdrDurand); title('Durand')

%% Save Results

if(saveResults)
    imwrite(hdrAdaptHist, [workingDir,filesep,'hdrAdaptHist.png']);
    imwrite(hdrSimpleCompression, [workingDir,filesep,'hdrSimpleCompression','.png']);
    imwrite(hdrReinhardGlobal, [workingDir,filesep,'hdrReinhard_global','.png']);
    imwrite(hdrReinhardLocal, [workingDir,filesep,'hdrReinhard_local','.png']);
    imwrite(hdrReinhardBilateral, [workingDir,filesep,'hdrReinhard_bilateral','.png']);
    imwrite(hdrDurand, [workingDir,filesep,'hdrDurand','.png']);
end

%% TMO Evaluation

for i = 1:numel(images)
    images(i).image = rgb2gray(im2double(images(i).image));
end

%% Set sliding window (patch) size

squareSize = 40; % patch size
stepSize = 10; % steps between patch locations

% patch locations
xList = 1 : stepSize : (size(imgHDR, 2)-squareSize);
yList = 1 : stepSize : (size(imgHDR, 1)-squareSize);

% Grid
[Xlist, Ylist] = ndgrid(xList,yList);
XY = [Xlist(:),Ylist(:)];
numXY = size(XY,1);

%% Prealloc

psnrScores = zeros(nTMOs, numel(images), numXY);

% for parfor
localPSNRscores = zeros(nTMOs, numel(images));

%% Evaluate TMOs

reverseStr = ''; % for status updates

for k = 1:numXY
    roi_ = [XY(k,:), squareSize, squareSize]; % current patch
    
    % Crop each TMO
    hdrAdaptHistCrop = imcrop(hdrAdaptHist, roi_);
    hdrSimpleCompressionCrop = imcrop(hdrSimpleCompression, roi_);
    hdrReinhardGlobalCrop = imcrop(hdrReinhardGlobal, roi_);
    hdrReinhardLocalCrop = imcrop(hdrReinhardLocal, roi_);
    hdrReinhardBilateralCrop = imcrop(hdrReinhardBilateral, roi_);
    hdrDurandCrop = imcrop(hdrDurand, roi_);
    
    % Execute in parallel
    parfor j = 1:numel(images)
        refImage = imcrop(images(j).image, roi_);
        
        localPSNRscores(:,j) = [psnr(hdrAdaptHistCrop,refImage);
                                psnr(hdrSimpleCompressionCrop, refImage);
                                psnr(hdrReinhardGlobalCrop, refImage);
                                psnr(hdrReinhardLocalCrop,refImage);
                                psnr(hdrReinhardBilateralCrop,refImage);
                                psnr(hdrDurandCrop,refImage)];
    end
    
    % add local parfor results to global
    psnrScores(:,:,k) = localPSNRscores;
    
    % Update status
    if(~mod(k-1,5))
        msg = sprintf('Processed %d/%d', k, numXY);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
end
fprintf('\n')

%% Process PSNR grid scores
% close all

powers = PowerdB';

hdrImages = zeros(size(imgHDR,1), size(imgHDR,2), 5);
hdrImages(:,:,1) =  hdrAdaptHist;
hdrImages(:,:,2) =  hdrSimpleCompression;
hdrImages(:,:,3) =  hdrReinhardGlobal;
hdrImages(:,:,4) =  hdrReinhardLocal;
hdrImages(:,:,5) =  hdrReinhardBilateral;
hdrImages(:,:,6) =  hdrDurand;

% Plot prep
figure('units','pixels','position',[200 200 size(imgHDR,2)*5 size(imgHDR,1)])
titles = {'Adaptive Hist Eq', 'Simple Compress', 'Reinhard Global', 'Reinhard Local', 'Reinhard Bilateral', 'Durand'};

% find arg max to get the reference image index that corresponds most with
% that window
maxMap = zeros(size(psnrScores,1),size(psnrScores,3));
argMaxMap = zeros(size(psnrScores,1),size(psnrScores,3));
for i = 1:size(psnrScores,1)
    localPSNR = squeeze(psnrScores(i,:,:)); % 15  X  52689
    [localMax, localArgMax] = max(localPSNR,[],1); % 1  X  52689

    maxMap(i,:) = localMax;
    argMaxMap(i,:) = localArgMax;

    argMaxImage = reshape(localArgMax,size(Xlist));
    
    argMaxDB = arrayfun(@(x) powers(x),argMaxImage);
    
    subplot(1,nTMOs,i);
    [C,h] = contourf(Xlist, Ylist, fliplr(argMaxDB),powers,'LineStyle','none');%,'ShowText','on')
    title(titles{i})

    colormap(jet)
    axis equal
    %axis off
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
    xlim([Xlist(1)-squareSize/2 Xlist(end)+squareSize/2])
    ylim([Ylist(1)-squareSize/2 Ylist(end)+squareSize/2])
end

if(saveResults)
    % Save as PDF
    set(gcf,'units','inch')
    ppos = get(gcf, 'Position');
    set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
    print([workingDir,filesep,'TMO_similarityMaps.pdf'],'-dpdf')
    
    % % If you have the export_fig package, this is easier to use
    % export_fig([workingDir,filesep,'TMO_similarityMaps'],'-pdf','-depsc');
end

%% Histogram
i_s = [1,5,9,15];
saveResults = false;

for i = 1:numel(i_s)
    plotHistogram(stack(:,:,1,i_s(i)),...
              sprintf('%.2fdB',PowerdB(i_s(i))),...
              saveResults, workingDir);
end

plotHistogram(hdrDurand,...
              'HDR_Durand',...
              saveResults, workingDir);
return
%% Study the effect of number of reference images on HDR quality

% Don't need tone mapping
% Compare HDR to HDR directly
% Take the 15 image version as ground truth

% set flags
normalize = true;
computeSSIM = true; % skipping this will make computations faster, MSE and SSIM have similar trends
computeMSE = true;
computeDynRanges = false; % not that informative
saveOutput = true; % save MAT files

% Generate combinations with 2 to 15 images
nImgs = 1:15;
imgIdxs = cell(numel(nImgs),1);
for i = 1:numel(imgIdxs)
    nImg = nImgs(i);
    if(i == 1)
        imgIdxs{i} = 8;
    else
        imgIdxs{i} = round(linspace(1,15,nImg));
    end
    % disp(imgIdxs{i})
end

% compute HDR
HDRimgs = zeros(size(imgHDR,1),size(imgHDR,2),size(imgHDR,3),numel(nImgs));
for i = 1:numel(imgIdxs)
    idxs = imgIdxs{i};
    subStack = stack(:,:,:,idxs);
    subExposures = stack_exposure(idxs);
    localHDR = BuildHDR(subStack, subExposures, 'LUT', lin_fun, 'Robertson', 'log', 1);
    HDRimgs(:,:,:,i) = localHDR;
end

% compute similarity
if(computeSSIM)
    ssimScores = zeros(numel(imgIdxs),1);
end
if(computeMSE)
    mseScores = zeros(numel(imgIdxs),1);
end
groundTruth = HDRimgs(:,:,1,end);
if(normalize)
    max_ = max(groundTruth(:));
    min_ = min(groundTruth(:));
    groundTruth = 255.*(groundTruth - min_) ./ (max_ - min_);
end
for i = 1:numel(imgIdxs)
    thisImage = HDRimgs(:,:,1,i);
    if(normalize)
        thisImage = 255.*(thisImage - min_) ./ (max_ - min_);
    end
    if(computeSSIM)
        ssimScores(i) = ssim(thisImage, groundTruth);
    end
    if(computeMSE)
        mseScores(i) = immse(thisImage, groundTruth);
    end
end

%% Plot
figure('units','normalized','position',[.25 .25 .3 .25])

nImgs = 1:15;
ax = gca;

yyaxis left
plot(nImgs, ssimScores, '-o', 'LineWidth', 2)
ylabel('SSIM')
xlabel('Number of images')
ylim([-0.1,1.1])
xlim([nImgs(1)-0.2,nImgs(end)+0.2])
xticks([1,5,10,15])
yticks(linspace(0,1,5))

yyaxis right
plot(nImgs, mseScores, '-.+', 'LineWidth', 2)
ylabel('MSE')
ax.YDir = 'reverse';
ylim([-0.1,1.1]*max(mseScores))
yticks(round(linspace(0,round(max(mseScores)),4)))

box on; grid on;

legend('SSIM','MSE','Location','southeast')
ax.LineWidth = 1;
ax.FontName = 'Times New Roman';
set(ax,'FontSize',16)

if(saveResults)
    % Save as PDF
    set(gcf,'units','inch')
    ppos = get(gcf, 'Position');
    set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
    print([workingDir,filesep,'nImagesStudy.pdf'],'-dpdf')
    
    % % If you have the export_fig package, this is easier to use
    % export_fig([workingDir,filesep,'nImagesStudy'],'-pdf','-depsc');
end

uniformSSIM = ssimScores;
uniformMSE = mseScores;
%% More complex: consider all combinations of images, not just linspace

% Generate combinations with 2 to 15 images
nImgs = 1:15;
imgIdxs = cell(numel(nImgs),1);
for i = 1:numel(imgIdxs)
    nImg = nImgs(i);
    imgIdxs{i} = combnk(nImgs, nImg); % n choose k
    % disp(imgIdxs{i})
end

% compute HDR
% too expensive to hold all in memory, compute ssim/mse and discard
groundTruth = imgHDR(:,:,1);
if(normalize)
    max_ = max(groundTruth(:));
    min_ = min(groundTruth(:));
    groundTruth = 255.*(groundTruth - min_) ./ (max_ - min_);
end

if(computeSSIM)
    ssimScores = cell(numel(nImgs),1);
end
if(computeMSE)
    mseScores = cell(numel(nImgs),1);
end
if(computeDynRanges)
    dynamicRanges = cell(numel(nImgs),2);
end
reverseStr = ''; % for status updates
for i = 1:numel(imgIdxs)
    if(computeSSIM)
        localssimScores = zeros(size(imgIdxs{i},1),1); % prealloc for use in parfor
    end
    if(computeMSE)
        localmseScores = zeros(size(imgIdxs{i},1),1); % prealloc for use in parfor
    end
    if(computeDynRanges)
        localDynRangesMin = zeros(size(imgIdxs{i},1),1); % prealloc for use in parfor
        localDynRangesMax = zeros(size(imgIdxs{i},1),1); % prealloc for use in parfor
    end
    subImgIdxs = imgIdxs{i};
    
    fprintf('Processing (%d choose %d) %d combinations', nImgs(end), nImgs(i), size(imgIdxs{i},1))
    parfor j = 1:size(imgIdxs{i},1)
        idxs = subImgIdxs(j,:);
        subStack = stack(:,:,:,idxs);
        subExposures = stack_exposure(idxs);
        thisHDR = BuildHDR(subStack, subExposures, 'LUT', lin_fun, 'Robertson', 'log', 1);
        thisHDR = thisHDR(:,:,1);
        
        % keep track of the dynamic range
        if(computeDynRanges)
            localDynRangesMin(j) = min(thisHDR(:));
            localDynRangesMax(j) = max(thisHDR(:));
        end
        
        % normalize dynamic range using ground truth's min/max
        if(normalize)
            thisHDR = 255.*(thisHDR - min_) ./ (max_ - min_);
        end
        
        % compute similarity
        if(computeSSIM)
            localssimScores(j) = ssim(thisHDR, groundTruth);
        end
        if(computeMSE)
            localmseScores(j) = immse(thisHDR, groundTruth);
        end

    end

    if(computeSSIM)
        ssimScores{i} = localssimScores;
    end
    if(computeMSE)
        mseScores{i} = localmseScores;
    end
    if(computeDynRanges)
        dynamicRanges{i,1} = localDynRangesMin;
        dynamicRanges{i,2} = localDynRangesMax;
    end
    fprintf('...done.\n')
    reverseStr = '';
end
fprintf('\n')

if(saveOutput)
    if(computeSSIM)
        save([workingDir,filesep,'ssimScores.mat'],'ssimScores')
    end
    if(computeMSE)
        save([workingDir,filesep,'mseScores.mat'],'mseScores')
    end
    if(computeDynRanges)
        save([workingDir,filesep,'dynamicRanges.mat'],'dynamicRanges')
    end
end

%%
if(computeDynRanges)
    dynRanges = cellfun(@(mins,maxs) maxs./mins, dynamicRanges(:,1), dynamicRanges(:,2), 'UniformOutput', false);
    [maxDynRanges, argmaxDynRanges] = cellfun(@max, dynRanges, 'UniformOutput', false);
    maxDynRanges = cell2mat(maxDynRanges);
end
%% Find the argmax of scores
if(computeSSIM)
    [maxSSIM, argmaxSSIM] = cellfun(@max, ssimScores, 'UniformOutput', false);
    [minSSIM, argminSSIM] = cellfun(@min, ssimScores, 'UniformOutput', false);
    maxSSIM = cell2mat(maxSSIM);
    minSSIM = cell2mat(minSSIM);
end
if(computeMSE)
    [maxMSE, argmaxMSE] = cellfun(@max, mseScores, 'UniformOutput', false);
    [minMSE, argminMSE] = cellfun(@min, mseScores, 'UniformOutput', false);
    maxMSE = cell2mat(maxMSE);
    minMSE = cell2mat(minMSE);
end

%% Plot max scores
if(computeMSE && computeSSIM)
    nImgs = 2:15;
    
    figure('units','normalized','position',[.25 .25 .3 .25])
    
    hold on
    
    ax = gca;
    
    yyaxis left
    semilogx(nImgs, maxSSIM(nImgs), '-o', 'LineWidth', 2)
    semilogx(nImgs, minSSIM(nImgs), '--x', 'LineWidth', 2)
    ylabel('SSIM')
    xlabel('Number of images')
    ylim([-0.1,1.1])
    xlim([nImgs(1)-0.1,nImgs(end)+1])
    xticks([1,5,10,15])
    yticks(linspace(0,1,5))
    
    yyaxis right
    semilogx(nImgs, minMSE(nImgs), '-.+', 'LineWidth', 2)
    semilogx(nImgs, maxMSE(nImgs), ':^', 'LineWidth', 2)
    ylabel('MSE')
    ax.YDir = 'reverse';
    ylim([-0.1*maxMSE(nImgs(1)),1.1*maxMSE(nImgs(1))])
    yticks(linspace(0,round(maxMSE(nImgs(1)),1),5))
    
    box on; grid on;
    
    legend('SSIM max','SSIM min','MSE min','MSE max','Location','southeast')
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    set(ax,'FontSize',16)
    
    if(saveResults)
        % Save as PDF
        set(gcf,'units','inch')
        ppos = get(gcf, 'Position');
        set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
        print([workingDir,filesep,'nImagesStudy.pdf'],'-dpdf')
        
        % % If you have the export_fig package, this is easier to use
        % export_fig([workingDir,filesep,'nImagesStudy'],'-pdf','-depsc');
    end
end
%% Only plot MSE
if(computeMSE)
    nImgs = 2:15;
    
    figure('units','normalized','position',[.25 .25 .3 .25])
    
    ax = gca;
    
    plot(nImgs, minMSE(nImgs), '-+', 'LineWidth', 2, 'MarkerSize', 9)
    hold on
    plot(nImgs, maxMSE(nImgs), '-.^', 'LineWidth', 2, 'MarkerSize', 9)
    plot(nImgs, uniformMSE(nImgs), ':ok', 'LineWidth', 2, 'MarkerSize', 9)
    ylabel('MSE')
    xlabel(sprintf('Number of images used out of %d', nImgs(end)))
    xlim([nImgs(1)-0.2, nImgs(end)+0.1])
    ylim([-0.1*maxMSE(nImgs(1)),1.1*maxMSE(nImgs(1))])
    xticks(round(linspace(nImgs(1),nImgs(end),4)))
%     yticks(round(linspace(0,round(maxMSE(nImgs(1))/5)*5,5)))
    yticks(linspace(0,100,5))
    yticks(linspace(0,180,5))
    
    box off
    
    legend('Best','Worst','Uniform','Location','northeast')
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    set(ax,'FontSize',16)
    
    if(saveResults)
        % Save as PDF
        set(gcf,'units','inch')
        ppos = get(gcf, 'Position');
        set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
        print([workingDir,filesep,'nImagesStudy_MSE.pdf'],'-dpdf')
        
        % % If you have the export_fig package, this is easier to use
        % export_fig([workingDir,filesep,'nImagesStudy_MSE'],'-pdf','-depsc');
    end
end
%% Only plot SSIM
if(computeSSIM)
    nImgs = 2:15;
    
    figure('units','normalized','position',[.25 .25 .3 .25])
    
    ax = gca;
    
    h1 = plot(nImgs, maxSSIM(nImgs), '-+', 'LineWidth', 2, 'MarkerSize', 9);
    hold on
    hline = refline([0,0.9]);
    hline.Color = 'k';
    hline.LineStyle = '--';
    h2 = plot(nImgs, minSSIM(nImgs), '-.^', 'LineWidth', 2, 'MarkerSize', 9);
    h3 = plot(nImgs, uniformSSIM(nImgs), ':ok', 'LineWidth', 2, 'MarkerSize', 9);
    
    xlabel(sprintf('Number of images used out of %d', nImgs(end)))
    xlim([nImgs(1)-0.2,nImgs(end)+0.2])
    ylabel('SSIM')
    ylim([-0.1, 1.1])
    xticks(round(linspace(nImgs(1),15,4)))
    yticks(linspace(0,1,5))
    
    box off
    %box on; grid on;
    
    legend([h1,h2,h3],{'Best','Worst','Uniform'},'Location','southeast')
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    set(ax,'FontSize',16)
    
    if(saveResults)
        % Save as PDF
        set(gcf,'units','inch')
        ppos = get(gcf, 'Position');
        set(gcf, 'PaperSize', [ppos(3) ppos(4)]);
        print([workingDir,filesep,'nImagesStudy_SSIM.pdf'],'-dpdf')
        
        % % If you have the export_fig package, this is easier to use
        % export_fig([workingDir,filesep,'nImagesStudy_SSIM'],'-pdf','-depsc');
    end
end
