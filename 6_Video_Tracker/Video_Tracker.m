cam = webcam;
nnet = alexnet; %image size [227, 227]
%nnet = vgg19; %image size [224,224]
%nnet= resnet101; %image size [224,224]
%%
for i = 1:15
    vidFrame = cam.snapshot;
    pic = imresize(vidFrame,[224,224]);
    label =  classify(nnet,pic);imshow(pic)
    title(upper(char(label)))
end
%% object recognition from a video
nnet = resnet101;
vid = VideoReader('Video Tracker.mov');
% classify the frames
n=30; %number of frames we want to skip
t = tiledlayout(6,6);
while hasFrame(vid)
    vidFrame = readFrame(vid);
    vid.CurrentTime = vid.CurrentTime + n/vid.FrameRate; 
    vidFrame = imresize(vidFrame,[224,224]);
    label =  classify(nnet,vidFrame); 
    nexttile; imshow(vidFrame); title(upper(char(label)))
end
