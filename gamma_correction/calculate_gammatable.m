%% my code
% This code is working very fine.
clear all
% 
load('D:\RF_fMRI\matlab_stim\gamma_correction\Luminance_projector1206.mat')
luminance=cell2mat(luminance)
INV_CLUT1=[];
for ii=1:3
    normLum=(luminance(ii,:)-luminance(ii,1))/(max(luminance(ii,:))-luminance(ii,1))
    aa=interp1([0:0.05:1],normLum,[0:255]'/255); % input  --> output
    bb=interp1(aa(:,1),[0:255]'/255,[0:255]'/255);%<-------- inverse function: output --> input
    figure;plot([aa bb])
    INV_CLUT1=[INV_CLUT1 bb];
end


%% fitting gamma function with curve-fitting.
% In this code, there is no inverse gamma function for the curve-fitting.
% The inverse gamma function was obtained from Joana.

clear all
% 
load('D:\RF_fMRI\matlab_stim\gamma_correction\Luminance_projector1206.mat')
luminance=cell2mat(luminance)
% 
Fits=[];
for ich=1:4
    indexValues=[0:0.05:1]*255;
    luminanceMeasurements = luminance(ich,:);


% 
%     figure; clf;
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'Position',[30 140 800 800]);
%     subplot(2,2,1);
%     plot(indexValues,luminanceMeasurements,'+');
%     hold on;
%     xlabel('Pixel Values');
%     ylabel('Luminance (cd/m2)');
%     strTitle{1}='Sampled Luminance Function';
%     strTitle{2}='Phase-1 Linear CLUT';
%     title(strTitle);
%     axis([0 256 0 max(luminanceMeasurements)]);
%     axis('square');
%     hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate and plot best-fit power function to sampled data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %zero-correct sampled luminance values
    lums=luminanceMeasurements-luminanceMeasurements(1);
    %normalize sampled luminance values
    normalizedLum=lums./max(lums);
    %trim zero level
    pixels=indexValues(2:end);
    normalizedLum=normalizedLum(2:end);

    %curve fit empirical luminance values 
    fitType = 7;  %curve-fitting
    outputx = [0:255];
    [Fit,~]=FitGamma(pixels',normalizedLum',outputx',fitType);
    Fits=[Fits Fit];

    %plot sampled luminance and curve fit results
    %figure(2);clf;hold on;
%     subplot(2,2,2); hold on;
%     plot(pixels,normalizedLum,'+'); %sampled luminance
%     plot(outputx,Fit,'r');  %curve fit results
%     axis([0 256 0 1]);
    
end


INV_CLUT2 =textread('D:\RF_fMRI\matlab_stim\gamma_correction\LUT_projector1206.txt') % <--- This inverse gamma function was calculated from Joana code.
figure;
for ii=1:3
    subplot(2,2,ii);
    plot([0:255],[Fits(:,ii)*255 INV_CLUT(:,ii) [0:255]']);
    axis equal
end

increments = 0:0.05:1;


%%




clear all
% 
% load('D:\RF_fMRI\matlab_stim\Luminance_projector.mat')
luminance=cell2mat(luminance)

tmp={};
for ich=1:4
    indexValues=[0:0.05:1]*255;
    luminanceMeasurements = luminance(ich,:);



    figure; clf;
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Position',[30 140 800 800]);
    subplot(2,2,1);
    plot(indexValues,luminanceMeasurements,'+');
    hold on;
    xlabel('Pixel Values');
    ylabel('Luminance (cd/m2)');
    strTitle{1}='Sampled Luminance Function';
    strTitle{2}='Phase-1 Linear CLUT';
    title(strTitle);
    axis([0 256 0 max(luminanceMeasurements)]);
    axis('square');
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate and plot best-fit power function to sampled data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %zero-correct sampled luminance values
    lums=luminanceMeasurements-luminanceMeasurements(1);
    %normalize sampled luminance values
    normalizedLum=lums./max(lums);
    %trim zero level
    pixels=indexValues(2:end);
    normalizedLum=normalizedLum(2:end);

    %curve fit empirical luminance values 
    fitType = 2;  %extended power function
    outputx = [0:255];
    [extendedFit,extendedX]=FitGamma(pixels',normalizedLum',outputx',fitType);

    %plot sampled luminance and curve fit results
    %figure(2);clf;hold on;
    subplot(2,2,2); hold on;
    plot(pixels,normalizedLum,'+'); %sampled luminance
    plot(outputx,extendedFit,'r');  %curve fit results
    axis([0 256 0 1]);
    xlabel('Pixel Values');
    ylabel('Normalized Luminance');
    strTitle{1}='Power Function Fit to Sampled Luminance Readings';
    strTitle{2}=['Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
    title(strTitle);
    hold off;
    
    maxLum = max(luminanceMeasurements);
    luminanceRamp=[0:1/255:1];
    pow=extendedX(1);
    offset=extendedX(2);
    invertedRamp=((maxLum-offset)*(luminanceRamp.^(1/pow)))+offset; %invert gamma w/o rounding
    %normalize inverse gamma table
    invertedRamp=invertedRamp./max(invertedRamp);
    %plot inverse gamma function
    %figure(4); clf; hold on;
    subplot(2,2,3); hold on;
    pels=[0:255];
    plot(pels,invertedRamp,'r');
    axis('square');
    axis([0 255 0 1]);
    xlabel('Pixel Values');
    ylabel('Inverse Gamma Table');
    strTitle{1}='Inverse Gamma Table Function';
    strTitle{2}=['for Exponent = ',num2str(extendedX(1)),'; Offset = ',num2str(extendedX(2))];
    title(strTitle);
    hold off;
    
    tmp{ich} = repmat(invertedRamp',1,3);
end








