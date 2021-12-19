v = VideoReader('5.mp4');

%% index argument
startingIndex = 1800 ;
observationLength = 200;

totalNumberOfFrame = v.NumFrames;

%% delimeters
% DelimiterIndex : 1
Da = [-1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1];
Da = repmat(Da, 6, 1);
Da = Da(:)';

% DelimiterIndex : 0
DoubleDa = [Da, Da];

% DelimiterIndex : 2
Db = [1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1];
Db = repmat(Db, 6, 1);
Db = Db(:)';

% DelimiterIndex : 3
Fa = [-1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1]; 
Fa = repmat(Fa, 6, 1);
Fa = Fa(:)';

% DelimiterIndex : 4
Fb = [1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1];
Fb = repmat(Fb, 6, 1);
Fb = Fb(:)';

bitSequence = [];
bitFlag = 0;
prevDeli = -1;
missing = 0;
prevPixelIndex = 1080;

%% read frames
for i = 1 : totalNumberOfFrame
    frame = rgb2gray(read(v, i));  
    fprintf("Frame index : %d, ", i);
    
    %% apply k-means on selected rows
    cluster = zeros(101, 1080);
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
    
    %% auto correlation for max value delimiter and pixel index
    corrMaxValue = [];
    corrMaxIndex = [];
    
    corrDoubleDa = (xcorr(majority, DoubleDa));
    corrDoubleDa = corrDoubleDa(length(majority): (2*length(majority)-1));
    [maxDoubleDa, pixelIndexDoubleDa] = max(corrDoubleDa);
    
    corrDa = (xcorr(majority, Da));
    corrDa = corrDa(length(majority): (2*length(majority)-length(Da)));
    [maxDa, pixelIndexDa] = max(corrDa);
    
    corrDb = (xcorr(majority, Db));
    corrDb = corrDb(length(majority): (2*length(majority)-length(Db)));
    [maxDb, pixelIndexDb] = max(corrDb);
    
    corrFa = (xcorr(majority, Fa));
    corrFa = corrFa(length(majority): (2*length(majority)-length(Fa)));
    [maxFa, pixelIndexFa] = max(corrFa);
    
    corrFb = (xcorr(majority, Fb));
    corrFb = corrFb(length(majority): (2*length(majority)-length(Fb)));
    [maxFb, pixelIndexFb] = max(corrFb);
    
    maxDelimiter = max([maxDa, maxDb, maxFa, maxFb]);
    
    %% decode symbol
    if maxDoubleDa >= 100 % contains delimiter Double Da
        fprintf("delimiter : DaDa, index : %d\n", pixelIndexDoubleDa);
        if bitFlag == 4
            break;
        end
        bitFlag = bitFlag + 1; % start decoding
        prevPixelIndex = pixelIndexDoubleDa;
        prevDeli = 0;
    elseif maxDelimiter == maxDa
        fprintf("delimiter : Da, index : %d\n", pixelIndexDa);
        if bitFlag == 4
            break;
        end
        bitFlag = bitFlag + 1;
        prevPixelIndex = pixelIndexDa;
        prevDeli = 1;
    elseif maxDelimiter == maxDb
        fprintf("delimiter : Db, index : %d\n", pixelIndexDb);
        if bitFlag ~= 0
            if missing == 1 || pixelIndexDb > prevPixelIndex % need to amend a symbol
                target = majority(pixelIndexDb-observationLength : pixelIndexDb);
                if sum(target) >= 0
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                else
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                end  
                missing = 0; % set to no lag
                prevPixelIndex = 1080; % set to max width
            end
            if pixelIndexDb < v.Width-length(Db)-observationLength % Db is not too close to the end, use symbol after Db would be good
                target = majority(pixelIndexDb+length(Db) : pixelIndexDb+length(Db)+observationLength);
                if sum(target) >= 0
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                else
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                end
            else % Db is too close to the end, use symbols before
                missing = 1;
            end
        end
        prevDeli = 2;
    elseif maxDelimiter == maxFa
        fprintf("delimiter : Fa, index : %d\n", pixelIndexFa);
        if bitFlag ~= 0
            if pixelIndexFa < v.Width-length(Fa)-observationLength % Fa is not too close to the end, use symbol after Db would be good
                target = majority(pixelIndexFa+length(Fa) : pixelIndexFa+length(Fa)+observationLength);
                if sum(target) >= 0
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                else
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                end
            else % Fa is too close to the end
                target = majority(pixelIndexFa-observationLength : pixelIndexFa);
                if sum(target) >= 0
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                else
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                end  
            end
        end
        prevDeli = 3;        
    elseif maxDelimiter == maxFb
        fprintf("delimiter : Fb, index : %d\n", pixelIndexFb);
        if bitFlag ~= 0
            if pixelIndexFb < v.Width-length(Fb)-observationLength % Fb is not too close to the end, use symbol after Db would be good
                target = majority(pixelIndexFb+length(Fb) : pixelIndexFb+length(Fb)+observationLength);
                if sum(target) >= 0
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                else
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                end
            else % Fb is too close to the end
                target = majority(pixelIndexFb-observationLength : pixelIndexFb);
                if sum(target) >= 0
                    bitSequence(end+1) = 1;
%                     fprintf('decode value : 1\n');
                else
                    bitSequence(end+1) = 0;
%                     fprintf('decode value : 0\n');
                end
            end
            prevPixelIndex = pixelIndexFb; % record Fb's pixel index for Db reference
        end
        prevDeli = 4;
    end
end

%% display decode bit
% fprintf("100111011111101101001011000111\n"); % 0.mp4 answer

for j = 1 : 30
    fprintf("%d", bitSequence(j));
end
fprintf("\n");

