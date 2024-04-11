imagefolder = dir('Portrait_painting/'); 
matrix_X = [];
for i=3:length(imagefolder) %we start at 3 since there are two extra files inside the directory
   filename = strcat('Portrait_painting/',imagefolder(i).name);
   test_face = im2gray(imread(filename));
   test_face = imresize(test_face, [50, 50]);
   test_face = im2gray(test_face); %make sure everything is grayscaled
   test_face = im2double(test_face);
   test_face = test_face(:)'; %turning the matrix into single horizontal array
   matrix_X = vertcat(matrix_X, test_face ); %vertically stacking each array
end
%% Applying PCA to matrix_X
[coeff, score, latent, tsquared, explained,mu] = pca(matrix_X); 
figure(1); plot(cumsum(explained));
xlabel("number of components"); ylabel("cumulative sum of explained")
%% Reconstructing faces using linear superpositions of weighted eigenfaces
j=35; %this is my face
test_face = matrix_X(j,:);
test_less_mu = test_face - mu;
PC = coeff(:,1:31); %remember that we only need a minimum of 31 principal components
a = PC'*test_less_mu';
image_compressed = PC*a;
image_compressed_recon = reshape(image_compressed(:,1),50,50);
figure, imshow(image_compressed_recon,[]); 
title('Compressed image with 31 components');
%% Quality check for varying number of principal components
j=50; %this is my face
test_face = matrix_X(j,:);
test_less_mu = test_face - mu;
t = tiledlayout(4,8,'TileSpacing','none','Padding','compact');
test_face_orig = test_face';
nexttile; imshow(reshape(test_face_orig(:,1),50,50),[]);
title("original");
for k=80:-3:1
    PC = coeff(:,1:k); 
    a = PC'*test_less_mu';
    image_compressed = PC*a;
    image_compressed_recon = reshape(image_compressed(:,1),50,50);
    nexttile; imshow(image_compressed_recon,[]);
    title(k);
end




