clear
close all
clc

%
%  The script converts long speech files to short 30 secs files  
%
%         Speech_path - Path for directory containing speech files
%         Noise_path - Path for directory containing noise files
%         dirname - Path for saving noisy speech files
%         snr - signal to noise ratio
%
%  Authors: Gautam Shreedhar Bhat
%
%  Copyright (c) 2019 by Gautam Shreedhar Bhat
%------------------------------------------------------------------------------------

fs2 = 16000;
in_folder = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\1hrSpeech';
out_folder = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\Part_1\';
audio_files=dir(fullfile(in_folder,'*.wav'));

addpath(in_folder)

for n = 1:10
    filename=audio_files(n).name;
    [clean, fs] = audioread(filename);
    %--- Sampling rate of 16 KHz
    fs_16 = 16000;
    [p,q] = rat(fs_16/fs,0.001);
    resamp_Sig = resample(clean,p,q);
    %--- 30 secs each file
    len=floor(30*fs_16);
    Nframes = floor(length(resamp_Sig)/len)-1;
    k=1;
    for m = 1:Nframes
        final = resamp_Sig(k:k+len-1,:);
        name_out = strcat(out_folder, filename(1:end-4),'_',num2str(m),'.wav');
        %--- Write Audio
        audiowrite(name_out,final,fs_16);
        k = k+len;
    end
end
disp('DONE');