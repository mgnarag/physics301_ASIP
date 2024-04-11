%% Step 1a. Reading the light source and Munsell
P_lambda = readmatrix('4_PCA 2/D65_and_A.xls'); 
Lamda = 400:5:700; 
P_lambda = interp1(P_lambda(:,1), P_lambda(:,2), Lamda'); 
load('4_PCA 2/munsell400_700_5.mat');
R_lambda = munsell; %reflectance of each color chip
% Step 1b. Getting the color signal
C_lambda = P_lambda.*R_lambda(:,:);
%% Step 1'
P_lambda = readmatrix('LED_Array.xls');
range = P_lambda(:, 1) > 399 & P_lambda(:, 1) < 701;
P_lambda = P_lambda(range, :);
%P_lambda = readmatrix('Philipps_LED.xls');
Lamda = [400:5:700]; %400 to 700 nm in intervals of 5nm
Lamda = Lamda'; %transposing Lambda
P_lambda = interp1(P_lambda(:,1), P_lambda(:,2), Lamda); 
P_lambda = P_lambda.*(1/max(P_lambda));
load('4_PCA 2/munsell400_700_5.mat');
R_lambda = munsell; %reflectance of each color chip
% Step 1b. Getting the color signal
C_lambda = P_lambda.*R_lambda(:,:);
%% 2. Computing eigenspectra of ensemble color signals
[coeff, score, latent, tsquared, explained,mu] = pca(C_lambda'); 
eigen = coeff(:,1); %get the two PCs
%% 3. Computing transformation matrix T
%FOR CANON D5
[data] = textread('4_PCA 2/camera_8_spectra.txt','','delimiter', ' ');
Lambda_raw = data(:,1); %1st column of data is the wavelength
S_R_raw = data(:,2);
S_G_raw = data(:,3);
S_B_raw = data(:,4); 
Lamda = 400:5:700; 
S_R = interp1(Lambda_raw, S_R_raw, Lamda'); 
S_G = interp1(Lambda_raw, S_G_raw, Lamda');
S_B = interp1(Lambda_raw, S_B_raw, Lamda');
S_rgb = [S_R, S_G, S_B];
T = S_rgb'*eigen;
%% 3' Olympus
S_R = readmatrix('red camera sensitivity.xls'); S_R = S_R(:,2);
S_G = readmatrix('green camera sensitivity.xls'); S_G = S_G(:,2);
S_B = readmatrix('blue camera sensitivity.xls'); S_B = S_B(:,2);
S_rgb = [S_R, S_G, S_B];
T = S_rgb'*eigen;
%% 4. Let's get values of q first
R_lambda_macbeth = readmatrix('MacbethColorChecker.xls');
range = R_lambda_macbeth(:, 1) > 399 & R_lambda_macbeth(:, 1) < 701;
R_lambda_macbeth = R_lambda_macbeth(range, :);
q_r = []; q_g = []; q_b = [];
for i = 2:25 %index of Macbeth patches, starts at i=2 
    PRS_red = P_lambda.* R_lambda_macbeth(:,i).*S_R(:,1); 
    PRS_red_sum = sum(PRS_red);
    PRS_green = P_lambda.* R_lambda_macbeth(:,i).*S_G(:,1);
    PRS_green_sum = sum(PRS_green);
    PRS_blue = P_lambda.* R_lambda_macbeth(:,i).*S_B(:,1);
    PRS_blue_sum = sum(PRS_blue);
    q_r = [q_r, PRS_red_sum];
    q_g = [q_g, PRS_green_sum];
    q_b = [q_b, PRS_blue_sum];
end
%% 4' From actual data
R_lambda_macbeth = readmatrix('MacbethColorChecker.xls');
range = R_lambda_macbeth(:, 1) > 399 & R_lambda_macbeth(:, 1) < 701;
R_lambda_macbeth = R_lambda_macbeth(range, :);
q_RGB = readmatrix(['GIMP Macbeth patch.xls']); 
q_r = q_RGB(:,2); q_r = q_r';
q_g = q_RGB(:,3); q_g = q_g';
q_b = q_RGB(:,4); q_b = q_b';
%% Applying the Wiener estimation
patch = 1;
C1 = inv(T'*T);
C2 = C1*T';
q = vertcat(q_r(:,patch),q_g(:,patch),q_b(:,patch));
a = C2*q;
%% For all patches
t = tiledlayout(4,6);
for patch = 1:24
    C1 = inv(T'*T);
    C2 = C1*T'; 
    q = vertcat(q_r(:,patch),q_g(:,patch),q_b(:,patch));
    a = C2*q;
    %Color signal using Equation 1:
    color_spectra_macbeth= P_lambda.*R_lambda_macbeth(:,patch+1);
    %Reconstructed color signal using Equation 3:
    C3 = eigen*a;
    nexttile; plot(color_spectra_macbeth,'-r'); hold;plot(C3,'-b');
    title(patch);
end
%% For all patches'
t = tiledlayout(4,6);
for patch = 1:24
    C1 = inv(T'*T);
    C2 = C1*T'; 
    q = vertcat(q_r(:,patch),q_g(:,patch),q_b(:,patch));
    a = C2*q;
    %Color signal using Equation 1:
    color_spectra_macbeth= P_lambda.*R_lambda_macbeth(:,patch+1);
    %Reconstructed color signal using Equation 3:
    C3 = eigen*a;
    nexttile; plot((color_spectra_macbeth./max(color_spectra_macbeth)),'-r'); hold;plot((C3./max(C3)),'-b');
    title(patch);
end
