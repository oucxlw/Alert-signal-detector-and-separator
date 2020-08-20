clear
clc

%
%  The script adds alert signal from thge directory to the noisy speech files in a directory.
%
%         Speech_path - Path for directory containing noisy speech files
%         Siren_path - Path for directory containing alert signal files
%         dirname - Path for saving  alert signal combined with noisy speech files
%         snr - signal to noise ratio
%
%  Authors: Gautam Shreedhar Bhat
%
%  Copyright (c) 2019 by Gautam Shreedhar Bhat
%------------------------------------------------------------------------------------


Speech_path  = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\Part_1';
Siren_path   = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Siren Noise for database\pixel siren';
dirname      = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\siren_files';
snr  = 5;

speechFiles             = dir(Speech_path);
sirenFiles              = dir(Siren_path);

%--- Fetch noise and speech files
data.speechFiles        = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.sirenFiles         = sirenFiles(~ismember({sirenFiles.name},{'.' '..'}));

%--- Number of speech and noise files
data.nSpeechFiles       = numel(data.speechFiles);
data.nSirenFiles        = numel(data.sirenFiles);

Folder   = Speech_path;
FileList = dir(fullfile(Folder, '*.wav'));
N        = data.nSpeechFiles;
index    = randperm(numel(FileList), N);
M        = floor(data.nSpeechFiles/data.nSirenFiles);
cnt      = 1;
m = 1;

%--- For every noisy speech signal file in the directory, add alert signal
%signals based on the snr and number of alert signal files
for n    = 1:data.nSirenFiles
    sirenFile           = [data.sirenFiles(n).folder '\' data.sirenFiles(n).name];
    [sirenSig_2,fs_si]  = audioread(sirenFile);
    fs_res              = 16000;
    [p,q]               = rat(fs_res/fs_si,0.001);
    sirenSig_1          = resample(sirenSig_2,p,q);
    sirenSig            = [sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1];
    sirenSig_ch1  =  sirenSig(:,1);
    sirenSig_ch2  =  sirenSig(:,2);
    
        for m = cnt:cnt+M
            if (m <= data.nSpeechFiles)
            audiofile = [data.speechFiles(m).folder '\' data.speechFiles(m).name];
            %--- Read Audio
            [speech,fs] = audioread(audiofile);
            %--- Split two channels
            speech_ch1 = speech(:,1);
            speech_ch2 = speech(:,2);
            disp (m)
            
            [speech_noise1]         = addnoise(speech_ch1,sirenSig_ch1,snr,fs_res);
            [speech_noise2]         = addnoise(speech_ch2,sirenSig_ch2,snr,fs_res);
            
            speech_noise           = [speech_noise1,speech_noise2];
            speech_noise = speech_noise./max(abs(speech_noise));
            name_siren            = strcat(dirname,'\',num2str(n),'_',data.speechFiles(m).name);
            %--- Write Audio
            audiowrite(name_siren,speech_noise,16000);
            end
        end
        cnt=cnt+M+1;
    
end
