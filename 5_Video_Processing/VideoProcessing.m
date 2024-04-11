%% Tracker
%Maricor Soriano (c) 2022
cam =webcam;
preview(cam);
pause(5);
img = snapshot(cam);
roi = imcrop(img);
%closePreview(cam);

%% Get Color Model
meanR = mean2(roi(:,:,1)); stdR = std2(roi(:,:,1));
meanG = mean2(roi(:,:,2)); stdG = std2(roi(:,:,2));
meanB = mean2(roi(:,:,3)); stdB = std2(roi(:,:,3));

%% Track
for i = 1:100
    tic
    I = snapshot(cam);
    %threshold
    pR = (I(:,:,1)>(meanR- stdR))& (I(:,:,1)<(meanR + stdR));
    pG = (I(:,:,2)>(meanG- stdG))& (I(:,:,2)<(meanG + stdG));
    pB = (I(:,:,3)>(meanB- stdB))& (I(:,:,3)<(meanB + stdB));
    p = pR.*pG.*pB;
    imshow(p);
    toc
end

