%loading the 4 images of different intensity functions
I1 = rgb2gray(imread('picture1.jpeg'));
I2 = rgb2gray(imread('picture2.jpeg'));
I3 = rgb2gray(imread('picture3.jpeg'));
I4 = rgb2gray(imread('picture4.jpeg'));
%figure(1);
%subplot(2,2,1); imshow(I1); title('I1, \theta = 0');
%subplot(2,2,2); imshow(I2); title('I2, \theta = \pi/2');
%subplot(2,2,3); imshow(I3); title('I3, \theta = \pi');
%subplot(2,2,4); imshow(I4); title('I4, \theta = 3\pi/2');
%% crop region of interest
[J1 RECT] = imcrop(I1);
J2 = imcrop(I2, RECT);
J3 = imcrop(I3, RECT);
J4 = imcrop(I4, RECT);
%figure(2);
%subplot(2,2,1); imshow(J1); title('\theta = 0');
%subplot(2,2,2); imshow(J2); title('\theta = \pi/2');
%subplot(2,2,3); imshow(J3); title('\theta = \pi');
%subplot(2,2,4); imshow(J4); title('\theta = 3\pi/2');

%% calculate phase
Num = double(J4) - double(J2);
Den = double(J1) - double(J3);
PHI = atan2(Num, Den);
%figure(3); imagesc(PHI);
%unwrap phase
Zobj = unwrap(PHI,[],2);
%figure(4); imagesc(Zobj);
%estimate reference
se = strel('disk',15);
background = imopen(Zobj,se);
%figure(5); imagesc(background);
%get phase difference
Obj = Zobj-background;
%figure(6); imagesc(Obj(:,1:200));
%figure(7); 
%s = surf(Obj(:,1:340));
%s.EdgeColor = 'none';

%% applying FFT to improve reconstructed image
Obj_fft = fftshift(fft2(Obj(:,1:200)));
size_fft = size(Obj_fft);
C=ones(size_fft);
%% 
for (y=round(size_fft(:,2)/2)+5 : size_fft(:,2)) %right mask
    for (x= round(size_fft(:,1)/2)-5: round(size_fft(:,1)/2)+5)
        C(x,y)=0;
    end
end
for (y= 1:round(size_fft(:,2)/2)-5) %left mask
    for (x= round(size_fft(:,1)/2)-5: round(size_fft(:,1)/2)+5)
        C(x,y)=0;
    end
end
%% 
for (y=round(size_fft(:,2)/2)-5 : round(size_fft(:,2)/2)+5) %right mask
    for (x= 1: round(size_fft(:,1)/2)-5)
        C(x,y)=0;
    end
end
for (y=round(size_fft(:,2)/2)-5 : round(size_fft(:,2)/2)+5) %left mask
    for (x= round(size_fft(:,1)/2)+5: round(size_fft(:,1)))
        C(x,y)=0;
    end
end
%% 
Obj_fft_filtered = Obj_fft.*C;
Obj_fft_filtered_inverse = ifft2(ifftshift(Obj_fft_filtered));

figure(8);
subplot(2,2,1); imagesc(log(abs(Obj_fft))); colormap(gray); title('FFT of Object'); 
subplot(2,2,2); imagesc(log(abs(Obj_fft_filtered))); colormap(gray);  title('Filtered FFT'); 
subplot(2,2,3); imagesc(Obj(:,1:200)); title('Reconstructed object');
subplot(2,2,4); imagesc(double(abs(Obj_fft_filtered_inverse))); title('Reconstructed Object with FFT');

