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

saveFlag = false;

folderNums = 15;

hdrAdaptHist = zeros(373,293,numel(folderNums));
hdrReinhard_global = zeros(373,293,numel(folderNums));
hdrReinhard_local = zeros(373,293,numel(folderNums));
hdrReinhard_bilateral = zeros(373,293,numel(folderNums));
hdrDurand = zeros(373,293,numel(folderNums));


for i = 1:numel(folderNums)
    
    currFolder = folderNums(i);
    
    if(ispc)
        name_folder = ['.\Data\20170214_190148624\forHDRtoolbox\',num2str(currFolder)];
    elseif(ismac)
        name_folder = ['./Data/20170214_190148624/forHDRtoolbox/',num2str(currFolder)];
    end
    
    files = dir([name_folder,filesep,'hdr*.png']);
    
    assert(numel(files) == 5)
    
    tmp = im2double(imread([name_folder,filesep,files(1).name]));
    hdrAdaptHist(:,:,i) = tmp(:,:,1);
    tmp = im2double(imread([name_folder,filesep,files(2).name]));
    hdrDurand(:,:,i) = tmp(:,:,1);
    tmp = im2double(imread([name_folder,filesep,files(3).name]));
    hdrReinhard_bilateral(:,:,i) = tmp(:,:,1);
    tmp = im2double(imread([name_folder,filesep,files(4).name]));
    hdrReinhard_global(:,:,i) = tmp(:,:,1);
    tmp = im2double(imread([name_folder,filesep,files(5).name]));
    hdrReinhard_local(:,:,i) = tmp(:,:,1);
end

refImageList = dir([name_folder,filesep,'rgb*.jpg']);
refImages = zeros(373,293,numel(refImageList));
for i = 1:numel(refImageList)
    tmp = im2double(imread([name_folder,filesep,refImageList(i).name]));
    refImages(:,:,i) = tmp(:,:,1);
end

%% Set window size
%h = figure;
%imshow(squeeze(refImages(:,:,end)))
%roi_ = getrect;
%close(h);

% 20, 10 is pretty good for visualizing the catheter + tissue surface
squareSize = 40;
stepSize = 10;

xList = 1 : stepSize : (size(hdrAdaptHist, 2)-squareSize);
yList = 1 : stepSize : (size(hdrAdaptHist, 1)-squareSize);

[Xlist, Ylist] = ndgrid(xList,yList);
XY = [Xlist(:),Ylist(:)];
numXY = size(XY,1);

%%
ssimScores = zeros(numel(folderNums),5,numel(refImageList), size(XY,1));
psnrScores = zeros(numel(folderNums),5,numel(refImageList), size(XY,1));
mutualInfoScores = zeros(numel(folderNums),5,numel(refImageList), size(XY,1));

localSSIMscores1 = zeros(1, numel(refImageList));
localSSIMscores2 = zeros(1, numel(refImageList));
localSSIMscores3 = zeros(1, numel(refImageList));
localSSIMscores4 = zeros(1, numel(refImageList));
localSSIMscores5 = zeros(1, numel(refImageList));

localPSNRscores1 = zeros(1, numel(refImageList));
localPSNRscores2 = zeros(1, numel(refImageList));
localPSNRscores3 = zeros(1, numel(refImageList));
localPSNRscores4 = zeros(1, numel(refImageList));
localPSNRscores5 = zeros(1, numel(refImageList));

localMIscores1 = zeros(1, numel(refImageList));
localMIscores2 = zeros(1, numel(refImageList));
localMIscores3 = zeros(1, numel(refImageList));
localMIscores4 = zeros(1, numel(refImageList));
localMIscores5 = zeros(1, numel(refImageList));

exponents = [1,1,1];%[0.25,2,1];

reverseStr = '';
for i = 1:numel(folderNums)
    fprintf('\nIterating through %d folder(s):\n',numel(folderNums));
    for k = 1:numXY
        roi_ = [XY(k,:), squareSize, squareSize];
            
        hdrAdaptHistCrop = imcrop(hdrAdaptHist(:,:,i), roi_);
        hdrDurandCrop = imcrop(hdrDurand(:,:,i), roi_);
        hdrReinhard_globalCrop = imcrop(hdrReinhard_global(:,:,i), roi_);
        hdrReinhard_localCrop = imcrop(hdrReinhard_local(:,:,i), roi_);
        hdrReinhard_bilateralCrop = imcrop(hdrReinhard_bilateral(:,:,i), roi_);

        parfor j = 1:numel(refImageList)
            refImage = squeeze(refImages(:,:,j));
            refImage = imcrop(refImage, roi_);

            %currFolder = folderNums(i);
            
            localSSIMscores1(j) = ssim(hdrAdaptHistCrop,refImage, 'Exponents', exponents);
            localSSIMscores2(j) = ssim(hdrDurandCrop,refImage, 'Exponents', exponents);
            localSSIMscores3(j) = ssim(hdrReinhard_globalCrop,refImage, 'Exponents', exponents);
            localSSIMscores4(j) = ssim(hdrReinhard_localCrop,refImage, 'Exponents', exponents);
            localSSIMscores5(j) = ssim(hdrReinhard_bilateralCrop,refImage, 'Exponents', exponents);

            localPSNRscores1(j) = psnr(hdrAdaptHistCrop,refImage);
            localPSNRscores2(j) = psnr(hdrDurandCrop,refImage);
            localPSNRscores3(j) = psnr(hdrReinhard_globalCrop, refImage);
            localPSNRscores4(j) = psnr(hdrReinhard_localCrop,refImage);
            localPSNRscores5(j) = psnr(hdrReinhard_bilateralCrop,refImage);
            
            localMIscores1(j) = compMutualInfo(hdrAdaptHistCrop,refImage, 3);
            localMIscores2(j) = compMutualInfo(hdrDurandCrop,refImage, 3);
            localMIscores3(j) = compMutualInfo(hdrReinhard_globalCrop, refImage, 3);
            localMIscores4(j) = compMutualInfo(hdrReinhard_localCrop,refImage, 3);
            localMIscores5(j) = compMutualInfo(hdrReinhard_bilateralCrop,refImage, 3);
            
%             immseScores(i,1,j,k) = immse(hdrAdaptHistCrop,refImage);
%             immseScores(i,2,j,k) = immse(hdrDurandCrop,refImage);
%             immseScores(i,3,j,k) = immse(hdrReinhard_globalCrop, refImage);
%             immseScores(i,4,j,k) = immse(hdrReinhard_localCrop,refImage);
%             immseScores(i,5,j,k) = immse(hdrReinhard_bilateralCrop,refImage);
        end
        
        localSSIMscores = [localSSIMscores1;
                           localSSIMscores2;
                           localSSIMscores3;
                           localSSIMscores4;
                           localSSIMscores5];
                      
        localPSNRscores = [localPSNRscores1;
                           localPSNRscores2;
                           localPSNRscores3;
                           localPSNRscores4;
                           localPSNRscores5];
        
        localMIscores = [localMIscores1;
                           localMIscores2;
                           localMIscores3;
                           localMIscores4;
                           localMIscores5];
                       
        ssimScores(i,:,:,k) = localSSIMscores;
        psnrScores(i,:,:,k) = localPSNRscores;
        mutualInfoScores(i,:,:,k) = localMIscores;
        
        if(~mod(k-1,5))
            msg = sprintf('Processed %d/%d', k, numXY);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        
    end
    fprintf('.')
end
fprintf('\n')

% save('.\Data\20170214_190148624\ssimScores_window50_20171106', 'ssimScores')
% save('.\Data\20170214_190148624\psnrScores_window50_20171106', 'psnrScores')

%%
% 
% figure
% h = plot(1,1);
% colr = get(gca,'colororder');
% delete(h)
% hold on
% for j = 1:numel(refImageList)
%     for i = 1:5
%     plot3(folderNums, squeeze(ssimScores(:,i,j)), ones(size(folderNums))*j, '-', 'Color', colr(i,:))
%     end
% end
% xlabel('N Images Combined'); ylabel('Score'); zlabel('Ref Image Idx')
% title('Structural Similarity Index')
% legend('Adaptive Hist','Durand','Global','Local','Bilateral','Location','nw')
% xlim([2,15])
% box on; grid on;
% % export_fig(['.\Data\20170214_190148624\forHDRtoolbox\','ssimScores'],'-pdf','-depsc');
% 
% figure
% hold on
% for j = 1:numel(refImageList)
%     for i = 1:5
%     plot3(folderNums, squeeze(psnrScores(:,i,j)), ones(size(folderNums))*j, '-', 'Color', colr(i,:))
%     end
% end
% xlabel('N Images Combined'); ylabel('Score'); zlabel('Ref Image Idx')
% title('Peak Signal-to-Noise Ratio (PSNR)')
% legend('Adaptive Hist','Durand','Global','Local','Bilateral','Location','nw')
% xlim([2,15])
% box on; grid on;
% % export_fig(['.\Data\20170214_190148624\forHDRtoolbox\','psnrScores'],'-pdf','-depsc');
% 
% figure
% hold on
% for j = 1:numel(refImageList)
%     for i = 1:5
%     plot3(folderNums, squeeze(immseScores(:,i,j)), ones(size(folderNums))*j, '-', 'Color', colr(i,:))
%     end
% end
% xlabel('N Images Combined'); ylabel('Score'); zlabel('Ref Image Idx')
% title('Mean-Squared Error')
% legend('Adaptive Hist','Durand','Global','Local','Bilateral','Location','ne')
% xlim([2,15])
% box on; grid on;
% % export_fig(['.\Data\20170214_190148624\forHDRtoolbox\','immseScores'],'-pdf','-depsc');

%% Process SSIM grid scores

close all

powers = [-26.7,-24.6,-22.2,-20.1,-18,-15.9,-14.1,-12.0,-10.5,-8.4,-6.3,-4.2,-2.4,-0.9,0.0];

hdrImages = zeros(size(hdrAdaptHist,1), size(hdrAdaptHist,2), 5);
hdrImages(:,:,1) =  squeeze(hdrAdaptHist(:,:,1));
hdrImages(:,:,2) =  squeeze(hdrDurand(:,:,1));
hdrImages(:,:,3) =  squeeze(hdrReinhard_global(:,:,1));
hdrImages(:,:,4) =  squeeze(hdrReinhard_local(:,:,1));
hdrImages(:,:,5) =  squeeze(hdrReinhard_bilateral(:,:,1));

% find arg max to get the reference image index that corresponds most with
% that window
ssimMaxRef = squeeze(ssimScores); % 5  X  15  X  52689
maxMap = zeros(size(ssimMaxRef,1),size(ssimMaxRef,3));
argMaxMap = zeros(size(ssimMaxRef,1),size(ssimMaxRef,3));
for i = 1:size(ssimMaxRef,1)
    localSSIM = squeeze(ssimMaxRef(i,:,:)); % 15  X  52689
    [localMax, localArgMax] = max(localSSIM,[],1); % 1  X  52689
    maxMap(i,:) = localMax;
    argMaxMap(i,:) = localArgMax;
%     figure
    argMaxImage = reshape(localArgMax,size(Xlist));
%     %argMaxImage = flipud(reshape(localArgMax,size(Xlist))');
%     %imagesc(argMaxImage)
%     %imagesc((argMaxImage-1)/14.0)
%     colormap bone
%     axis equal
%     xlim([0 size(argMaxImage,2)])
%     ylim([0 size(argMaxImage,1)])
%     colorbar('FontSize',11);
%     
%     thisImage = squeeze(hdrImages(:,:,i));
%     
%     combinedImage = repmat(thisImage*.5,[1,1,3]);
%     combinedImage(51:end-50,51:end-50,1) = combinedImage(51:end-50,51:end-50,1) + (argMaxImage-6)/7.*.5;
%     combinedImage(51:end-50,51:end-50,2) = combinedImage(51:end-50,51:end-50,2) + (argMaxImage-7)/7.*.5;
%     imshow(combinedImage)

    figure
    %contourf(flipud(argMaxImage),1:15,'ShowText','on')
    %contourf(Xlist, Ylist, fliplr(argMaxImage),1:15,'ShowText','on')
    imagesc(repmat(squeeze(hdrImages(:,:,i)),[1,1,3]))
    hold on
    argMaxDB = arrayfun(@(x) powers(x),argMaxImage);
    [C,h] = contourf(Xlist+squareSize/2, Ylist+squareSize/2, argMaxDB, powers,'LineWidth',1.0,'Color','k');%,'ShowText','on')
    set(h,'LineColor','none')
    colormap(jet)
    axis equal
    %axis off
    %xlim([Xlist(1)-squareSize/2 Xlist(end)+squareSize/2])
    %ylim([Ylist(1)-squareSize/2 Ylist(end)+squareSize/2])
    export_fig(sprintf('.\\Data\\20170214_190148624\\forHDRtoolbox\\SSIMcontourfsNL_TMO_%d.pdf',i),'-pdf','-depsc');
    
end

%% Process PSNR grid scores
close all

powers = [-26.7,-24.6,-22.2,-20.1,-18,-15.9,-14.1,-12.0,-10.5,-8.4,-6.3,-4.2,-2.4,-0.9,0.0];

hdrImages = zeros(size(hdrAdaptHist,1), size(hdrAdaptHist,2), 5);
hdrImages(:,:,1) =  squeeze(hdrAdaptHist(:,:,1));
hdrImages(:,:,2) =  squeeze(hdrDurand(:,:,1));
hdrImages(:,:,3) =  squeeze(hdrReinhard_global(:,:,1));
hdrImages(:,:,4) =  squeeze(hdrReinhard_local(:,:,1));
hdrImages(:,:,5) =  squeeze(hdrReinhard_bilateral(:,:,1));

% find arg max to get the reference image index that corresponds most with
% that window
psnrMaxRef = squeeze(psnrScores); % 5  X  15  X  52689
maxMap = zeros(size(psnrMaxRef,1),size(psnrMaxRef,3));
argMaxMap = zeros(size(psnrMaxRef,1),size(psnrMaxRef,3));
for i = 1:size(psnrMaxRef,1)
    localPSNR = squeeze(psnrMaxRef(i,:,:)); % 15  X  52689
    [localMax, localArgMax] = max(localPSNR,[],1); % 1  X  52689
    maxMap(i,:) = localMax;
    argMaxMap(i,:) = localArgMax;
%     figure
    argMaxImage = reshape(localArgMax,size(Xlist));
%     %argMaxImage = flipud(reshape(localArgMax,size(Xlist))');
%     %imagesc(argMaxImage)
%     %imagesc((argMaxImage-1)/14.0)
%     colormap bone
%     axis equal
%     xlim([0 size(argMaxImage,2)])
%     ylim([0 size(argMaxImage,1)])
%     colorbar('FontSize',11);
%     
%     thisImage = squeeze(hdrImages(:,:,i));
%     
%     combinedImage = repmat(thisImage*.5,[1,1,3]);
%     combinedImage(51:end-50,51:end-50,1) = combinedImage(51:end-50,51:end-50,1) + (argMaxImage-6)/7.*.5;
%     combinedImage(51:end-50,51:end-50,2) = combinedImage(51:end-50,51:end-50,2) + (argMaxImage-7)/7.*.5;
%     imshow(combinedImage)
    
    argMaxDB = arrayfun(@(x) powers(x),argMaxImage);
    figure
    %contourf(flipud(argMaxImage),1:15,'ShowText','on')
    %imagesc(repmat(squeeze(hdrImages(:,:,i)),[1,1,3]))
    %hold on
    [C,h] = contourf(Xlist, Ylist, fliplr(argMaxDB),powers,'LineWidth',0.01,'Color','k');%,'ShowText','on')
    %[C,h] = contourf(Xlist+squareSize/2, Ylist+squareSize/2, argMaxDB,powers,'LineWidth',1.0,'Color','k');%,'ShowText','on')
    %set(h,'LineColor','none')
    %colormap(flipud(jet))
    colormap(jet)
    %hc = contourfcmap(Xlist, Ylist, fliplr(argMaxDB),powers+.001,jet(14), 'cbarloc', 'eastoutside', 'method', 'calccontour');
%     clabel(C,h,powers,'FontWeight','bold',...
%                        'FontSize',11,...
%                        'Color','k',...
%                        'Interpreter','latex') % 'FontSize',10, 'EdgeColor','k', 'BackgroundColor','w',...
    axis equal
    %axis off
    xlim([Xlist(1)-squareSize/2 Xlist(end)+squareSize/2])
    ylim([Ylist(1)-squareSize/2 Ylist(end)+squareSize/2])
%     export_fig(sprintf('.\\Data\\20170214_190148624\\forHDRtoolbox\\PSNRcontourfsNL_win40_step10_TMO_%d.pdf',i),'-pdf','-depsc');
end

%% Process MI grid scores
close all

powers = [-26.7,-24.6,-22.2,-20.1,-18,-15.9,-14.1,-12.0,-10.5,-8.4,-6.3,-4.2,-2.4,-0.9,0.0];

hdrImages = zeros(size(hdrAdaptHist,1), size(hdrAdaptHist,2), 5);
hdrImages(:,:,1) =  squeeze(hdrAdaptHist(:,:,1));
hdrImages(:,:,2) =  squeeze(hdrDurand(:,:,1));
hdrImages(:,:,3) =  squeeze(hdrReinhard_global(:,:,1));
hdrImages(:,:,4) =  squeeze(hdrReinhard_local(:,:,1));
hdrImages(:,:,5) =  squeeze(hdrReinhard_bilateral(:,:,1));

% find arg max to get the reference image index that corresponds most with
% that window
miMaxRef = squeeze(mutualInfoScores); % 5  X  15  X  52689
maxMap = zeros(size(miMaxRef,1),size(miMaxRef,3));
argMaxMap = zeros(size(miMaxRef,1),size(miMaxRef,3));
for i = 1:size(miMaxRef,1)
    localMI = squeeze(miMaxRef(i,:,:)); % 15  X  52689
    [localMax, localArgMax] = min(localMI,[],1); % 1  X  52689
    maxMap(i,:) = localMax;
    argMaxMap(i,:) = localArgMax;
%     figure
    argMaxImage = reshape(localArgMax,size(Xlist));
%     %argMaxImage = flipud(reshape(localArgMax,size(Xlist))');
%     %imagesc(argMaxImage)
%     %imagesc((argMaxImage-1)/14.0)
%     colormap bone
%     axis equal
%     xlim([0 size(argMaxImage,2)])
%     ylim([0 size(argMaxImage,1)])
%     colorbar('FontSize',11);
%     
%     thisImage = squeeze(hdrImages(:,:,i));
%     
%     combinedImage = repmat(thisImage*.5,[1,1,3]);
%     combinedImage(51:end-50,51:end-50,1) = combinedImage(51:end-50,51:end-50,1) + (argMaxImage-6)/7.*.5;
%     combinedImage(51:end-50,51:end-50,2) = combinedImage(51:end-50,51:end-50,2) + (argMaxImage-7)/7.*.5;
%     imshow(combinedImage)
    
    argMaxDB = arrayfun(@(x) powers(x),argMaxImage);
    figure
    %contourf(flipud(argMaxImage),1:15,'ShowText','on')
    %imagesc(repmat(squeeze(hdrImages(:,:,i)),[1,1,3]))
    %hold on
    [C,h] = contourf(Xlist, Ylist, fliplr(argMaxDB),powers,'LineWidth',0.01,'Color','k');%,'ShowText','on')
    %[C,h] = contourf(Xlist+squareSize/2, Ylist+squareSize/2, argMaxDB,powers,'LineWidth',1.0,'Color','k');%,'ShowText','on')
    %set(h,'LineColor','none')
    %colormap(flipud(jet))
    colormap(jet)
    %hc = contourfcmap(Xlist, Ylist, fliplr(argMaxDB),powers+.001,jet(14), 'cbarloc', 'eastoutside', 'method', 'calccontour');
%     clabel(C,h,powers,'FontWeight','bold',...
%                        'FontSize',11,...
%                        'Color','k',...
%                        'Interpreter','latex') % 'FontSize',10, 'EdgeColor','k', 'BackgroundColor','w',...
    axis equal
    %axis off
    xlim([Xlist(1)-squareSize/2 Xlist(end)+squareSize/2])
    ylim([Ylist(1)-squareSize/2 Ylist(end)+squareSize/2])
%     export_fig(sprintf('.\\Data\\20170214_190148624\\forHDRtoolbox\\PSNRcontourfsNL_win40_step10_TMO_%d.pdf',i),'-pdf','-depsc');
end