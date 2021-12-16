v = VideoReader('0.mp4');

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

filter = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7];

%% read frames
% for i = 1 : totalNumberOfFrame
frame = rgb2gray(read(v, 10)); % current selected frame : 10, correct delimiter : Fa
% frame(selectedRow-3:selectedRow+3, 1:489) = 0;
imshow(frame);
centerRow = frame(selectedRow, :);

% preprocess with low pass filter
centerRow = xcorr(centerRow, filter);
centerRow = centerRow(v.Width-3:2*v.Width-4);

% normalization
normalizedInput = normalize(centerRow, 'range', [-1, 1]);

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

corrDoubleDa = (xcorr(normalizedInput, DoubleDa));
corrDoubleDa = corrDoubleDa(length(normalizedInput): (2*length(normalizedInput)-length(DoubleDa)));
[maxDoubleDa, pixelIndex] = max(corrDoubleDa);
fprintf("max correlation value for DoubleDa = %f, pixel index = %d to %d\n", maxDoubleDa, pixelIndex, pixelIndex+length(DoubleDa)-1);
autoCorrelationValue(end+1) = maxDoubleDa;

corrDa = (xcorr(normalizedInput, Da));
corrDa = corrDa(length(normalizedInput): (2*length(normalizedInput)-length(Da)));
[maxDa, pixelIndex] = max(corrDa);
fprintf("max correlation value for Da = %f, pixel index = %d to %d\n", maxDa, pixelIndex, pixelIndex+length(Da)-1);
autoCorrelationValue(end+1) = maxDa;

corrDb = (xcorr(normalizedInput, Db));
corrDb = corrDb(length(normalizedInput): (2*length(normalizedInput)-length(Db)));
[maxDb, pixelIndex] = max(corrDb);
fprintf("max correlation value for Db = %f, pixel index = %d to %d\n", maxDb, pixelIndex, pixelIndex+length(Db)-1);
autoCorrelationValue(end+1) = maxDb;

corrFa = (xcorr(normalizedInput, Fa));
corrFa = corrFa(length(normalizedInput): (2*length(normalizedInput)-length(Fa)));
[maxFa, pixelIndex] = max(corrFa);
fprintf("max correlation value for Fa = %f, pixel index = %d to %d\n", maxFa, pixelIndex, pixelIndex+length(Fa)-1);
autoCorrelationValue(end+1) = maxFa;

corrFb = (xcorr(normalizedInput, Fb));
corrFb = corrFb(length(normalizedInput): (2*length(normalizedInput)-length(Fb)));
[maxFb, pixelIndex] = max(corrFb);
fprintf("max correlation value for Fb = %f, pixel index = %d to %d\n", maxFb, pixelIndex, pixelIndex+length(Fb)-1);
autoCorrelationValue(end+1) = maxFb;
    
% end
