% Initialize environment
close all
clear all
clc

% File and video processing setup
[Options] = Setup_Options();
[StackFilenames, DefaultPathname] = uigetfile('*.tif','Select .tif files to be processed');
CurrentFilename = StackFilenames;
CurrStackFilePath = strcat(DefaultPathname,CurrentFilename);
VideoFilePath=CurrStackFilePath;
VideoInfo = imfinfo(VideoFilePath);
NumFrames = length(VideoInfo);
ImageWidth = VideoInfo.Width; %in pixels
ImageHeight = VideoInfo.Height; %in pixels

% Parameters

NumberofPixels = 512;      % Total number of pixels
Focused_Frame = 17;        % Frame index known to be focused
Fs = 1;                    % Sampling frequency
threshold = 50;            % Global threshold for processing


% Thresholds for each frame
threshold1 = threshold;
threshold2 = threshold / 2;
threshold3 = threshold / 2;



% Preallocate video matrix
VideoMatrix = zeros(ImageHeight, ImageWidth, NumFrames, 'uint16');
framematrix = zeros(NumFrames, 1);

% Averaging over a 3x3 window around each pixel in video frames
windowSize = 3; 
kernel = ones(windowSize) / (windowSize^2); % Create averaging kernel

for b = 1:NumFrames
    % Read current frame from video
    CurrentFrameImage = imread(CurrStackFilePath, b);

    % Apply convolution for averaging
    AveragedFrame = conv2(double(CurrentFrameImage), kernel, 'same');

    % Store the processed frame in VideoMatrix
    VideoMatrix(:, :, b) = AveragedFrame;

    % Maintain a frame index
    framematrix(b) = b;     
end



index=zeros(NumFrames, NumberofPixels*2);

 for j=3:NumFrames-1
     % Calculate the difference between adjacent frames
    diff1=double(VideoMatrix(:,:,j)-VideoMatrix(:,:,j-1));
    diff2=double(VideoMatrix(:,:,j+1))-double(VideoMatrix(:,:,j));
    diff3=double(VideoMatrix(:,:,j-1))-double(VideoMatrix(:,:,j-2));
    % Create a logical mask of significant changes
    condition=(diff1>=threshold1).*(diff2>=-threshold2).*(diff3>=-threshold3).*diff1;
    % Find indices where condition is true
   indexcondition=find(condition>0);
 

for (i=1:numel(indexcondition))
    M = zeros(size(condition));
    M(indexcondition(i)) = 1; % location
b=condition(conv2(M,[1,1,1;1,0,1;1,1,1],'same')>0);
if condition(indexcondition(i))<=max(b)
    condition(indexcondition(i))=0;
end
end 
% Store indices of significant changes
    index(j-1,:)=[find(condition>0)', zeros(NumberofPixels*2-numel(find(condition>0)),1)'];
 end

%Definig new index to filter out indices that are adjacent  to indices already included in the finalindex from previous frames. 

 finalindex(1,:)=sort(index(1,:), 'descend');
for i=2:NumFrames
    c=finalindex(i-1,:);
    finalindex=[finalindex; ~ismember(index(i,:),c).*index(i,:).*~ismember(index(i,:),c+1).* ~ismember(index(i,:),c-1).* ~ismember(index(i,:),c-NumberofPixels).* ~ismember(index(i,:),c+NumberofPixels).*~ismember(index(i,:),c+2).*~ismember(index(i,:),c-2).*~ismember(index(i,:),c+2*NumberofPixels).*~ismember(index(i,:),c-2*NumberofPixels)];
    finalindex(i,:)=sort(finalindex(i,:), 'descend');
end

% Remove replicated elements from finalindex
C=finalindex;
D=zeros(size(C));
% Find unique elements across the entire array
[~,i] = unique(C);
% Place unique elements back into D at the same positions
D(i)=C(i);
% Flatten and sort the non-zero elements of D
infoarray=sort(nonzeros(reshape(D.',1,[])));

%Removing the pixels on sides, Removing 3 layers
%Sides=[SideWest SideEast SideBottom Sidetop];
% Define the pixels on the edges to be removed
Sides = [1:NumberofPixels*3, NumberofPixels^2-NumberofPixels*3+1:NumberofPixels^2, 1:NumberofPixels:NumberofPixels^2-NumberofPixels+1, ...
         2:NumberofPixels:NumberofPixels^2-NumberofPixels+2, 3:NumberofPixels:NumberofPixels^2-NumberofPixels+3, ...
         (NumberofPixels-2):NumberofPixels:(NumberofPixels^2-2), (NumberofPixels-1):NumberofPixels:(NumberofPixels^2-1), ...
         NumberofPixels:NumberofPixels:NumberofPixels^2];

Sides=unique(Sides);
% Remove the edge pixels from infoarray
infoarray=~ismember(infoarray,Sides).*infoarray;

% sort infoarray in descending order
infoarray=sort(infoarray,'descend');
% For removing 0 elements from infoarray (because it is problematic
%infoarray(i(j)) I will write:
for i=numel(infoarray):-1:1
    if infoarray(i)==0
        infoarray(i)=[];
    end
end

NumberofPixelsWithEvents=length(infoarray);

x=zeros(NumberofPixelsWithEvents,NumFrames);
xLowPassFilter=zeros(NumberofPixelsWithEvents,NumFrames);
count=zeros(NumberofPixelsWithEvents,1);

%The reason for choosing 50 is that I do not think I will see more than 50
%events in a pixel, it is arbitrary.

Ap=zeros(NumberofPixelsWithEvents,50);


for j=1:NumberofPixelsWithEvents
 for i=1:NumFrames
    sub=VideoMatrix(:,:,i);
    for k=1:9
    x(j,i)=x(j,i)+sub(infoarray(j)-5+k);
    end
    x(j,i)=x(j,i)/9;
 end 

   % Apply low pass filter
 xLowPassFilter(j,:)=lowpass(x(j,:),0.0001);

 A=xLowPassFilter(j,:);
% Find peaks in the filtered signal
 [~,App]=findpeaks(A(Focused_Frame:1000),Fs,'MinPeakProminence',200,'Annotate','extents');
 if numel(App)>0 & App>Focused_Frame
     App=App+Focused_Frame;
 Ap(j,:)=[App zeros(1,50-numel(App))];
 end
 count(j)=sum(Ap(j,:)>0);
end

%Although I think all counts==3 or some of them even higher are 
%correct, I will remove count>=3 to make sure i=everything is correct
%This is the code for removing pixels with equal or more than 3 events
%from count, Ap, and inforarray
i=find(count>=3);
i=unique(i);
i=sort(i,'descend');
for j=1:numel(i)
    count(i(j))=[];
    Ap(i(j),:)=[];
    infoarray(i(j))=[];
    x(i(j),:)=[];
    xLowPassFilter(i(j),:)=[];
end

%To make the frame matrix based on Pixels
%I do not think I will have more than 20 pixels for a frame. 
%50 is arbitrary
FramebyPixel=zeros(NumFrames,50);
for j=1:max(max(Ap))
    [i,l]=find(Ap==j);
    for k=1:numel(i)
        FramebyPixel(j,k)=infoarray(i(k));
    end
end

%I will start from Bigger pixels to smaller ones. Then, I will
FramebyPixel=sort(FramebyPixel,2,'descend');


%For removing the adjacent pixels in consecutive frames:
for i=1:NumFrames-1
    if numel(nonzeros(FramebyPixel(i,:)))>0
        for j=1:numel(nonzeros(FramebyPixel(i,:)))
            t=FramebyPixel(i,j);
            adj = [t, t+1, t-1, t+NumberofPixels, t-NumberofPixels, t-NumberofPixels+1, t-NumberofPixels-1, ...
       t+NumberofPixels+1, t+NumberofPixels-1, t+NumberofPixels-2, t+NumberofPixels+2, ...
       t-NumberofPixels+2, t-NumberofPixels-2, t+NumberofPixels*2-2, t-NumberofPixels*2+2, ...
       t+NumberofPixels*2-1, t-NumberofPixels*2+1, t+NumberofPixels*2, t-NumberofPixels*2, ...
       t+NumberofPixels*2+1, t-NumberofPixels*2-1, t+NumberofPixels*2+2, t-NumberofPixels*2-2, ...
       t+2, t-2, t+NumberofPixels*3, t-NumberofPixels*3, t-NumberofPixels*3+1, t+NumberofPixels*3-1, ...
       t+NumberofPixels*3-2, t-NumberofPixels*3+2, t+NumberofPixels*3+1, t-NumberofPixels*3-1, ...
       t+NumberofPixels+7, t-NumberofPixels+7, t+NumberofPixels*2-1, t-NumberofPixels*2+1, ...
       t+NumberofPixels*2-NumberofPixels*2+2];
            k=ismember(adj,FramebyPixel(i+1,:))+ismember(adj,FramebyPixel(i+2,:))+ismember(adj,FramebyPixel(i+3,:))+ismember(adj,FramebyPixel(i+4,:))+ismember(adj,FramebyPixel(i+5,:))+ismember(adj,FramebyPixel(i+6,:))+ismember(adj,FramebyPixel(i+7,:))++ismember(adj,FramebyPixel(i+8,:))++ismember(adj,FramebyPixel(i+9,:))+ismember(adj,FramebyPixel(i+10,:))+ismember(adj,FramebyPixel(i+11,:))+ismember(adj,FramebyPixel(i+12,:))+ismember(adj,FramebyPixel(i+13,:))+ismember(adj,FramebyPixel(i+14,:))+ismember(adj,FramebyPixel(i+15,:))+ismember(adj,FramebyPixel(i+16,:))+ismember(adj,FramebyPixel(i+17,:))++ismember(adj,FramebyPixel(i+18,:));
            if sum(k)>0
                FramebyPixel(i,j)=0;
            end  
            l=ismember(adj,FramebyPixel(i,j+1:end));
            if sum(l)>0
                FramebyPixel(i,j)=0;
            end 


        end
    end 
end


FramebyPixel=sort(FramebyPixel,2,'descend');
TotalNumberofEvents=sum(sum(FramebyPixel>0));
ContinuousCDF=zeros(1,NumFrames);
for i=2:NumFrames
    ContinuousCDF(i)=sum(FramebyPixel(i,:)>0)/TotalNumberofEvents+ContinuousCDF(i-1);
end
[DiscreteCDF,framediscrete,~]=unique(ContinuousCDF);


