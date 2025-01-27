 %% Opening the Macbeth chart
[file, path] = uigetfile('*.png');
filename =[path,file];
I = imread(filename);
figure(1);imshow(I);
%% Cropping the region of interest and displaying the histograms
J = imcrop(I);
figure(2);
subplot(1,3,1); imhist(J(:,:,1)); %for Red channel
subplot(1,3,2); imhist(J(:,:,2)); %for Green channel
subplot(1,3,3); imhist(J(:,:,3));  %for Blue channel
%% get cropped image BW
R = I(:,:,1); G = I(:,:,2); B = I(:,:,3);
%BW = (R >145) & (R<173) & (G>50) & (G<71) & (B>108) & (B<148);
BW = (R >177) & (R<255) & (G>171) & (G<253) & (B>98) & (B<182);
figure(3); imshow (BW);
%% get cropped image RGB
BW = ~BW;
BW = im2uint8(BW); 
BW = BW.*(1/255);
Icrop = BW.*I(:,:,:);
figure();imshow(Icrop);
%% PARAMETRIC SEGMENTATION
[file, path] = uigetfile('*.png');
filename =[path,file];
I = imread(filename);
figure(1);imshow(I);
%% Crop region of interest
J = imcrop(I);
%% Applying Gaussian distribution function per pixel
R = J(:,:,1); G = J(:,:,2); B = J(:,:,3);
mu_R = mean2(R); mu_G = mean2(G); mu_B = mean2(B);
mu = [mu_R,mu_G,mu_B];
R_1 = R(:); G_1 = G(:); B_1 = B(:);
RGB_1 = double([R_1 G_1 B_1]);
cov_RGB = cov(RGB_1);
fac = 1./(sqrt(2*pi.*det(cov_RGB)));
p = zeros(size(I,1),size(I,2));
for i = 1:size(I,1)
    for j = 1:size(I,2)
        x_r_mu = double([I(i,j,1) I(i,j,2) I(i,j,3)]) - mu;
        p(i,j) = fac.*exp(-0.5*x_r_mu*inv(cov_RGB)*x_r_mu'); 
    end
end
figure(); imagesc(p); colormap(gray);
% get segmented image in RGB
p = p./max(p);
%p = (p>0.1); %scaling
p = (p<0.1);
p = im2uint8(p); 
p = p./255;
Icrop = p.*I(:,:,:);
figure();imshow(Icrop);
%% Applying 1-D Gaussian
R = J(:,:,1); G = J(:,:,2); B = J(:,:,3);
mu_R = mean2(R); mu_G = mean2(G); mu_B = mean2(B);
R_1 = R(:); G_1 = G(:); B_1 = B(:);
R_1 = double(R_1); G_1 = double(G_1); B_1 = double(B_1); 
cov_R = cov(R_1); cov_G = cov(G_1); cov_B = cov(B_1);
fac_R = 1./(sqrt(2*pi.*det(cov_R)));
fac_G = 1./(sqrt(2*pi.*det(cov_G)));
fac_B = 1./(sqrt(2*pi.*det(cov_B)));
p_r = zeros(size(I,1),size(I,2));
p_g = zeros(size(I,1),size(I,2));
p_b = zeros(size(I,1),size(I,2));
for i = 1:size(I,1)
    for j = 1:size(I,2)
        x_r_mu = double([I(i,j,1)]) - mu_R;
        p_r(i,j) = fac_R.*exp(-0.5*x_r_mu*inv(cov_R)*x_r_mu'); 
        x_g_mu = double([I(i,j,2)]) - mu_G;
        p_g(i,j) = fac_G.*exp(-0.5*x_g_mu*inv(cov_G)*x_g_mu'); 
        x_b_mu = double([I(i,j,3)]) - mu_B;
        p_b(i,j) = fac_B.*exp(-0.5*x_b_mu*inv(cov_B)*x_b_mu'); 
    end
end
p = p_r(:,:).* p_g(:,:) .* p_b(:,:);
p = p./max(p);
figure(); imagesc(p); colormap(gray);
% get segmented image in RGB
p = (p<0.1); %scaling
p = im2uint8(p); 
p = p./255;
Icrop = p.*I(:,:,:);
figure();imshow(Icrop);
%% NON PARAMETRIC code from Maam Jing
clear all ; close all;
BINS = 32;
[filename,pathname] = uigetfile('*.jpg');
J = imread([pathname,filename]);
[I, rect] = imcrop(J);
%% Get the r g of the whole image
J = double(J);
R = J(:,:,1); G = J(:,:,2); B = J(:,:,3); % dissect the channels
I= R + G + B; %Equation 1
I(find(I==0))=100000; %to prevent NaNs
rJ = R./ I; gJ = G./I; %Equation 2
%% Crop the region of interest in the rg space
r = imcrop(rJ, rect); g = imcrop(gJ, rect);
rint = round(r*(BINS-1) + 1);
gint = round(g*(BINS-1) + 1);
colors = gint(:) + (rint(:)-1)*BINS;
%% Compute rg-histogram
% This is the 1-d version of a 2-d histogram
hist = zeros(BINS*BINS,1);
for row = 1:BINS
    for col = 1:(BINS-row+1)
        hist(col+(row-1)*BINS) = length(find(colors==(((col + (row-1)*BINS)))));
    end
end
%% Backproject histogram
rJint = round(rJ*(BINS-1)+1);
gJint = round(gJ*(BINS-1)+1);
colorsJ = gJint(:) + (rJint(:)-1)*BINS;
HB = hist(colorsJ);
HBImage = reshape(HB,size(J,[1,2]));
%% Making it colored
HBImage = im2uint8(HBImage); 
HBImage = HBImage./max(HBImage);
[filename,pathname] = uigetfile('*.jpg');
J = imread([pathname,filename]);
Icrop = HBImage.*J(:,:,:);
figure();imshow(Icrop);
%% Making it colored for foreground
HBImage = im2uint8(HBImage); 
HBImage = HBImage./max(HBImage);
BW = (HBImage<0.9); %logical operation
BW = im2uint8(BW); 
BW = BW.*(1/255);
[filename,pathname] = uigetfile('*.png');
J = imread([pathname,filename]);
Icrop = BW.*J(:,:,:);
figure();imshow(Icrop);