clear
clc
%
%  The script randomly moves half the number of files from one location to other  
%
%         source - Path for directory containing noisy speech files
%         Dest - Path for moving noisy speech files
%
%  Authors: Gautam Shreedhar Bhat
%
%  Copyright (c) 2019 by Gautam Shreedhar Bhat
%------------------------------------------------------------------------------------

source   = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\Part_1';
Dest     = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\Part_2';

for m = 1:1

    FileList = dir(fullfile(source, '*.wav'));
    %--- Select half the number of files
    N = length(FileList)/2;
   
    index    = randperm(numel(FileList), N);
    for k = 1:N
        file_source = fullfile(source, FileList(index(k)).name);
        disp(file_source);
        %--- Move files from source to destination
        movefile(file_source, Dest);
    end
end
disp('DONE!!!')
