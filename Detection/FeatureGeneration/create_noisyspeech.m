clear
close all
%
%  The script adds noise files from noise directory to the speech files in speech directory.
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

rng(0);

Speech_path = 'C:\Users\gxs160730\Documents\Ph.D\noise and speech files to record\pixel recorded\two_clean';
Noise_path  = 'C:\Users\gxs160730\Documents\Ph.D\noise and speech files to record\pixel recorded\noise\Traffic';
dirname  = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\1hrSpeech\';
snr = 0;

speechFiles = dir(Speech_path);
noiseFiles  = dir(Noise_path);

%--- Fetch noise and speech files
data.speechFiles = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.noiseFiles = noiseFiles(~ismember({noiseFiles.name},{'.' '..'}));

%--- Number of speech and noise files
data.nSpeechFiles = numel(data.speechFiles);
data.nNoiseFiles  = numel(data.noiseFiles);
count=0;

%--- For every speech file in the directory, add noise signals based on the
%snr
for n = 1:data.nNoiseFiles
    noiseFile           = [data.noiseFiles(n).folder '\' data.noiseFiles(n).name];
    [noiseSig_2,fs_n]   = audioread(noiseFile);
    noiseSig            = [noiseSig_2;noiseSig_2];
    [p,q]              = rat(16000/fs_n,0.001);
    noiseSig            = resample(noiseSig,p,q);
    
    mkdir(dirname);
    cd(dirname);
    for s=1:data.nSpeechFiles
        speechFile             = [data.speechFiles(s).folder '\' data.speechFiles(s).name];
        [speechSig_2,fs_s]     = audioread(speechFile);
        [p1,q1]                = rat(16000/fs_n,0.001);
        speechSig              = resample(speechSig_2,p1,q1);
        
        [speech_noise1]         = addnoise(speechSig(:,1),noiseSig(:,1),snr,16000);
        speech_noise1           = speech_noise1/max(abs(speech_noise1));
        
        [speech_noise2]         = addnoise(speechSig(:,2),noiseSig(:,2),0,16000);
        speech_noise2           = speech_noise2/max(abs(speech_noise2));
        
        speech_noise           = [speech_noise1, speech_noise2];
        name_siren1             = strcat(dirname,'\',num2str(n),'_',data.speechFiles(s).name);
        
        %--- Write Audio
        audiowrite(name_siren1,speech_noise,16000)
    end
    cd('..')
    
end

