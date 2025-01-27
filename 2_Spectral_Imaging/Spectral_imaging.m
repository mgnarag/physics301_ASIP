%%reading the spectral sensitivity, S(lambda), of Canon D5 camera
[data] = textread('2_Spectral_Imaging/camera_8_spectra.txt','','delimiter', ' ');
Lambda_raw = data(:,1); %1st column of data is the wavelength
S_R_raw = data(:,2); %2nd column contains the sensitivity on R channel
S_G_raw = data(:,3); %3rd column contains the sensitivity on G channel
S_B_raw = data(:,4); %4th column contains the sensitivity on B channel

%majority of the spectral power distribution P(lambda) of light sources and
%reflectance of an object R(lambda) are in intervals d_lambda = 5nm.
%The S(lambada) that we have are in intervals of 4nm, thus, we need to
%interpolate the data:
Lamda = [400:5:700]; %400 to 700 nm in intervals of 5nm
Lamda = Lamda'; %transposing Lambda
%final spectral sensitivity with correct inteval:
S_R = interp1(Lambda_raw, S_R_raw, Lamda); 
S_G = interp1(Lambda_raw, S_G_raw, Lamda);
S_B = interp1(Lambda_raw, S_B_raw, Lamda);

%% Plotting the interpolated values
figure(1); plot(Lambda_raw, S_R_raw,'-ko'); hold on; plot(Lamda,S_R,'rx'); 
legend('S_{raw}','S_{interp}', 'Location','northwest');
xlabel('\lambda (nm)'); ylabel ('S');
title('Interpolated S_{R}(\lambda) for Canon D5 camera');
figure(2); plot(Lambda_raw, S_G_raw,'-ko'); hold on; plot(Lamda,S_G,'gx');
legend('S_{raw}','S_{interp}', 'Location','northwest');
xlabel('\lambda (nm)'); ylabel ('S');
title('Interpolated S_{G}(\lambda) for Canon D5 camera');
figure(3); plot(Lambda_raw, S_B_raw,'-ko'); hold on; plot(Lamda,S_B,'bx');
legend('S_{raw}','S_{interp}', 'Location','northeast');
xlabel('\lambda (nm)'); ylabel ('S');
title('Interpolated S_{B}(\lambda) for Canon D5 camera');

%% Reading the xls file
P_lambda = readmatrix('2_Spectral_Imaging/D65_and_A.xls'); %light source
%P_lambda = readmatrix('Fluorescents.xls'); %lightsource of fluorescent
R_lambda = readmatrix('2_Spectral_Imaging/MacbethColorChecker.xls'); %reflectivity of object

%The light source and macbeth color chart does not match the range of our
%desired wavelength - lightsource (300,830), macbeth(380,780).
%To match, we will only extract 400 to 700nm.
range = P_lambda(:, 1) > 399 & P_lambda(:, 1) < 701;
P_lambda = P_lambda(range, :);
range = R_lambda(:, 1) > 399 & R_lambda(:, 1) < 701;
R_lambda = R_lambda(range, :);

%% Applying equation 2
j = 2 %index of lightsource data, j=2 for A, j=8 for D65
t = tiledlayout(4,6,'TileSpacing','none','Padding','compact');
V_red_all = [];
for i = 2:25 %index of Macbeth patches, starts at i=2 
    %element-to-element multiplication and getting the integral/summation
    PRS_red = P_lambda(:,j).* R_lambda(:,i).*S_R(:,1); 
    PRS_red_sum = sum(PRS_red);
    PRS_green = P_lambda(:,j).* R_lambda(:,i).*S_G(:,1);
    PRS_green_sum = sum(PRS_green);
    PRS_blue = P_lambda(:,j).* R_lambda(:,i).*S_B(:,1);
    PRS_blue_sum = sum(PRS_blue);
    
    PS_red = P_lambda(:,j).* S_R(:,1);
    PS_red_sum = sum(PS_red);
    PS_green = P_lambda(:,j).* S_G(:,1);
    PS_green_sum = sum(PS_green);
    PS_blue = P_lambda(:,j).* S_B(:,1);
    PS_blue_sum = sum(PS_blue);
    
    V_red = (PRS_red_sum/PS_red_sum);
    V_green = (PRS_green_sum/PS_green_sum);
    V_blue = (PRS_blue_sum/PS_blue_sum);
    %V_red = PRS_red_sum;
    %V_red_all = [V_red_all, V_red];
    %V_green = PRS_green_sum;
    %V_blue = PRS_blue_sum;
    
    %creating the individual patches
    PATCH = ones(100,100,3);
    PATCH(:,:,1) = V_red;
    PATCH(:,:,2) = V_green;
    PATCH(:,:,3) = V_blue;
    nexttile; imshow(PATCH);
end
%% Rendering a scene/natural images
%[data] = textread('Spectral_Imaging/camera_8_spectra.txt','','delimiter', ' ');
[data] = textread('2_Spectral_Imaging/Nikon D1X.txt','','delimiter', ' ');
Lambda_raw = data(:,1); %1st column of data is the wavelength
S_R_raw = data(:,2); %2nd column contains the sensitivity on R channel
S_G_raw = data(:,3); %3rd column contains the sensitivity on G channel
S_B_raw = data(:,4); %4th column contains the sensitivity on B channel

Lamda = [400:10:700]; %400 to 700 nm in intervals of 5nm
Lamda = Lamda'; %transposing Lambda

%final spectral sensitivity with correct interval:
S_R = interp1(Lambda_raw, S_R_raw, Lamda); 
S_G = interp1(Lambda_raw, S_G_raw, Lamda);
S_B = interp1(Lambda_raw, S_B_raw, Lamda);

%load Spectral_Imaging/Scene7/ref_ribeira1bbb_reg1.mat
load 2_Spectral_Imaging/Scene4/ref_cyflower1bb_reg1.mat
R_lambda = reflectances;
R_lambda(:,:,33) = []; %removing 720nm
R_lambda(:,:,32) = []; %removing 710nm

%% 
for j = 2:6:6
    %P_lambda = readmatrix('Spectral_Imaging/Fluorescents.xls'); %lightsource of fluorescent
    P_lambda = readmatrix('2_Spectral_Imaging/D65_and_A.xls');
    %P_lambda(1,:)= []; %for Fluorescent light
    range = P_lambda(:, 1) > 399 & P_lambda(:, 1) < 701;
    P_lambda = interp1(P_lambda(:,1), P_lambda(:,j), Lamda); 
    
    V_red_all = [];
    V_green_all = [];
    V_blue_all = [];
    for x = 1:size(R_lambda,1)
        V_red_column = [];
        V_green_column = [];
        V_blue_column = [];
        for y = 1:size(R_lambda,2)
            R_lambda_1 = [];
            for z = 1:31
                R_lambda_1 = [R_lambda_1, R_lambda(x,y,z)];
            end
            R_lambda_1 = R_lambda_1.';
            PRS_red = P_lambda(:,1).* R_lambda_1(:,1).*S_R(:,1); 
            PRS_red_sum = sum(PRS_red);
            PRS_green = P_lambda(:,1).* R_lambda_1(:,1).*S_G(:,1);
            PRS_green_sum = sum(PRS_green);
            PRS_blue = P_lambda(:,1).* R_lambda_1(:,1).*S_B(:,1);
            PRS_blue_sum = sum(PRS_blue);
            
            PS_red = P_lambda(:,1).* S_R(:,1);
            PS_red_sum = sum(PS_red);
            PS_green = P_lambda(:,1).* S_G(:,1);
            PS_green_sum = sum(PS_green);
            PS_blue = P_lambda(:,1).* S_B(:,1);
            PS_blue_sum = sum(PS_blue);
    
            V_red = (PRS_red_sum/PS_red_sum);
            V_red_column = [V_red_column,V_red];
            V_green = (PRS_green_sum/PS_green_sum);
            V_green_column = [V_green_column,V_green];
            V_blue = (PRS_blue_sum/PS_blue_sum);
            V_blue_column = [V_blue_column,V_blue];
        end
        V_red_all = [V_red_all;V_red_column];
        V_green_all = [V_green_all;V_green_column];
        V_blue_all = [V_blue_all;V_blue_column];
    
    end
    
    figure(j)
    imshow(cat(3,V_red_all,V_green_all,V_blue_all))
end 