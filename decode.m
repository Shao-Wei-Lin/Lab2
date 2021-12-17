warning('off','all');

v = VideoReader('0.mp4');

%% index argument
startingIndex = 1800 ;
frameIndex = 74;

totalNumberOfFrame = v.NumFrames;
centerHeight = v.Height/2;
thirdQuartile = v.Height*3/4;
centerWidth = v.Width/2;
selectedRow = thirdQuartile;

%% delimeters
Da = [-1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1];
Da = repmat(Da, 6, 1);
Da = Da(:)';

DoubleDa = [Da, Da];

Db = [1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1];
Db = repmat(Db, 6, 1);
Db = Db(:)';

Fa = [-1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1]; 
Fa = repmat(Fa, 6, 1);
Fa = Fa(:)';

Fb = [1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1];
Fb = repmat(Fb, 6, 1);
Fb = Fb(:)';

% filter = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7];

%% read frames
% for i = 1 : totalNumberOfFrame
frame = rgb2gray(read(v, frameIndex));
imshow(frame);

fprintf("Frame index : %d\n", frameIndex);

cluster = zeros(101, 1080);

% apply k-means on selected rows
for j = startingIndex : startingIndex+100
    targetRow = double(frame(j, :)');  
    receivedIndex = (kmeans(targetRow, 2, 'Replicates', 30) - 1)';
    small = -1;
    big = -1;
    for k = 1 : 1080
        if receivedIndex(k) == 0 && small == -1
            small = frame(j, k);
        end
        if receivedIndex(k) == 1 && big == -1
            big = frame(j, k);
        end
        if small ~= -1 && big ~= -1
            break
        end
    end

    if small > big
        cluster(j-startingIndex+1, :) = 1 - receivedIndex; 
    else
        cluster(j-startingIndex+1, :) = receivedIndex;
    end
end

majority = int8((sum(cluster, 1)/101) > 0.5);
majority(majority == 0) = -1;
% centerRow = majority;

% preprocess with low pass filter
% centerRow = xcorr(centerRow, filter);
% centerRow = centerRow(v.Width-3:2*v.Width-4);

% normalization
% normalizedInput = normalize(centerRow, 'range', [-1, 1]);

%% find max correlation delimeter
delimeters = [Da, Db, Fa, Fb];
delimeterNames = ['Da', 'Db', 'Fa', 'Fb'];
% for j = 1 : 4
%     corrValue = (xcorr(normalizedInput, delimeters(j)));
%     corrValue = corrValue(length(normalizedInput): (2*length(normalizedInput)-length(delimeters(j))));
%     [maxValue, pixelIndex] = max(corrValue);
%     fprintf("max correlation value for %s = %f, pixel index = %d\n", delimeterNames(j), maxValue, pixelIndex);
% end

autoCorrelationValue = [];

corrDoubleDa = (xcorr(majority, DoubleDa));
corrDoubleDa = corrDoubleDa(length(majority): (2*length(majority)-1));
[maxDoubleDa, pixelIndex] = max(corrDoubleDa);
fprintf("max correlation value for DoubleDa = %f, pixel index = %d to %d\n", maxDoubleDa, pixelIndex, pixelIndex+length(DoubleDa)-1);
autoCorrelationValue(end+1) = maxDoubleDa;

corrDa = (xcorr(majority, Da));
corrDa = corrDa(length(majority): (2*length(majority)-length(Da)));
[maxDa, pixelIndex] = max(corrDa);
fprintf("max correlation value for Da = %f, pixel index = %d to %d\n", maxDa, pixelIndex, pixelIndex+length(Da)-1);
autoCorrelationValue(end+1) = maxDa;

corrDb = (xcorr(majority, Db));
corrDb = corrDb(length(majority): (2*length(majority)-length(Db)));
[maxDb, pixelIndex] = max(corrDb);
fprintf("max correlation value for Db = %f, pixel index = %d to %d\n", maxDb, pixelIndex, pixelIndex+length(Db)-1);
autoCorrelationValue(end+1) = maxDb;

corrFa = (xcorr(majority, Fa));
corrFa = corrFa(length(majority): (2*length(majority)-length(Fa)));
[maxFa, pixelIndex] = max(corrFa);
fprintf("max correlation value for Fa = %f, pixel index = %d to %d\n", maxFa, pixelIndex, pixelIndex+length(Fa)-1);
autoCorrelationValue(end+1) = maxFa;

corrFb = (xcorr(majority, Fb));
corrFb = corrFb(length(majority): (2*length(majority)-length(Fb)));
[maxFb, pixelIndex] = max(corrFb);
fprintf("max correlation value for Fb = %f, pixel index = %d to %d\n", maxFb, pixelIndex, pixelIndex+length(Fb)-1);
autoCorrelationValue(end+1) = maxFb;

maxDelimiter = max([maxDa, maxDb, maxFa, maxFb]);

if maxDoubleDa >= 100
    fprintf("delimiter : DaDa\n");
elseif maxDelimiter == maxDa
    fprintf("delimiter : Da\n");
elseif maxDelimiter == maxDb
    fprintf("delimiter : Db\n");
elseif maxDelimiter == maxFa
    fprintf("delimiter : Fa\n");
elseif maxDelimiter == maxFb
    fprintf("delimiter :Fb\n");
end

% end
