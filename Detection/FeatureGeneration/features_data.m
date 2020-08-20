clear
close all
clc

%
%  The script adds alert signal from thge directory to the noisy speech files in a directory.
%
%         Speech_path - Path for directory containing noisy speech files
%         Siren_path  - Path for directory containing alert signal combined with noisy speech files
%         dirname     - Path for saving  the .mat file
%         fileName    - Name of the .mat file
%
%  Authors: Gautam Shreedhar Bhat
%
%  Copyright (c) 2019 by Gautam Shreedhar Bhat
%------------------------------------------------------------------------------------

Speech_path               = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\Part_1_test';
Siren_path                = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Small_files\siren_files_test';
dirname                   = 'C:\Users\gxs160730\Documents\GitHub\Alert-signal-detector-and-separator\Detection\NoisySpeech\Traffic\Features';
fileName                  = 'Traffic_features';

speechFiles               = dir(Speech_path);
sirenFiles                = dir(Siren_path);

%--- Fetch noise and speech files
data.speechFiles          = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.sirenFiles           = sirenFiles(~ismember({sirenFiles.name},{'.' '..'}));

%--- Number of speech and noise files
data.nSpeechFiles         = numel(data.speechFiles);
data.nSirenFiles          = numel(data.sirenFiles);

%--- Initialize the parameters
params.overlap           = 256;
params.nfft              = 512;
params.window            = hanning(params.nfft);
params.n_mels            = 26;
params.fmin              = 300;
params.fmax              = 8000;
params.fs                = 16000;
params.filterbank        = buildFilterbank(params.fmin, params.fmax, params.n_mels, params.nfft, params.fs);

for n = 1:data.nSpeechFiles
    [speech,fs_si]         = audioread([data.speechFiles(n).folder '\' data.speechFiles(n).name]);
    speech_ch1             = [speech(:,1);ones(length(speech),1)];
    speech_ch2             = [speech(:,2);ones(length(speech),1)];
    
    for m = n:n
        
        [siren,fs_si]      = audioread([data.sirenFiles(m).folder '\' data.sirenFiles(m).name]);
        siren_ch1          = [siren(:,1);zeros(length(siren),1)];
        siren_ch2          = [siren(:,2);zeros(length(siren),1)];
        disp([m n])
        
        %--- Concatenate audio
        A = [siren_ch1 speech_ch1 speech_ch1 siren_ch1 speech_ch1 siren_ch1];
        B = A(:,randperm(size(A,2)));
        C = B(1:length(speech),:);
        D = B(length(speech)+1:end,:);
        C = reshape(C,[1,6*length(speech)]);
        
        %--- Generate features
        audioFrames     = reshape(C, params.overlap, []);
        procFrames      = params.window .* [audioFrames(:,1:end-1); audioFrames(:,2:end)];
        procFFT         = fft(procFrames, params.nfft);
        
        real_ch1 = real(procFFT(1:params.nfft/2+1,:));
        imag_ch1 = imag(procFFT(1:params.nfft/2+1,:));
        
        mag_procFFT     = abs(procFFT).^2;
        STFTMdb_tmp     = 10*log(mag_procFFT);
        STFTMdb         = STFTMdb_tmp(1:params.nfft/2+1,:);
        [len_rows, len_col] = size(procFrames);
        
        %--- Initialize matrix
        if (n==1)
            trainData  = zeros(len_col*data.nSirenFiles,(params.nfft/2+1) * 2 + 1);
        end
        
        %--- Generate labels
        labelFrames     = reshape(D, params.overlap, []);
        procLabels      = [labelFrames(:,1:end-1); labelFrames(:,2:end)];
        len_1           = length(procLabels(1,:));
        
        for i = 1:len_1
            if (sum (procLabels(:,i) ) > 0)
                Final_lables(:,i) = 0;
            else
                Final_lables(:,i) = 1;
            end
        end
        
        count = n;
        
        input = [real_ch1;imag_ch1];
        feature_vec = [input' Final_lables'];
        
        %--- Fill the matrix trainData
        if count == 1
            trainData(1:count*len_1,:) = feature_vec;
            old_feature_vec = feature_vec;
        else
            trainData((count-1)*len_1+1:count*len_1,:) = feature_vec;
        end
    end
    
end

%--- Write Features as .mat file
feature_name = strcat(dirname,'\',fileName);
save(feature_name,'trainData','-v7.3');
disp('DONE!!!')