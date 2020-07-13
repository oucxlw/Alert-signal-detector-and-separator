% Siren_Noisy_Speech

addpath('C:\Users\gxs160730\Documents\Ph.D\CNN-SE\CNN-VAD\Training Code\Functions\');
addpath('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Features');

Speech_path               = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\clean_part2';
Siren_path                = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\clean_siren';

speechFiles               = dir(Speech_path);
sirenFiles                = dir(Siren_path);

data.speechFiles          = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.sirenFiles           = sirenFiles(~ismember({sirenFiles.name},{'.' '..'}));

data.nSpeechFiles         = numel(data.speechFiles);
data.nSirenFiles          = numel(data.sirenFiles);

params.overlap           = 256;
params.nfft              = 512;
params.window            = hanning(params.nfft);
params.n_mels            = 26;
params.fmin              = 300;
params.fmax              = 8000;
params.fs                = 16000;
params.filterbank        = buildFilterbank(params.fmin, params.fmax, params.n_mels, params.nfft, params.fs);

cnt=1;cnt1=1;

%trainData               = zeros(9374*data.nSirenFiles,28); 
trainData                = zeros(11249*data.nSirenFiles,515); 

for n = 1:data.nSpeechFiles
    
    speech2                = [data.speechFiles(n).folder '\' data.speechFiles(n).name];
    [speech1,fs_si]        = audioread(speech2);
    
    speech_ch1                 = [speech1(:,1);ones(length(speech1),1)];
    speech_ch2                 = [speech1(:,2);ones(length(speech1),1)];    
    for m = n:n
        
        siren2             = [data.sirenFiles(m).folder '\' data.sirenFiles(m).name];
        [siren1,fs_si]     = audioread(siren2);
        siren_ch1          = [siren1(:,1);zeros(length(siren1),1)];
        siren_ch2          = [siren1(:,2);zeros(length(siren1),1)];
            
        disp([m n])
        
        A = [siren_ch1 speech_ch1 speech_ch1 siren_ch1 speech_ch1 siren_ch1 ];
        B = A(:,randperm(size(A,2)));
        C = B(1:length(speech1),:);
        D = B(length(speech1)+1:end,:);
        
        C = reshape(C,[1,6*length(speech1)]);
        plot(C);
        %% features
        audioFrames     = reshape(C, params.overlap, []);
        procFrames      = params.window .* [audioFrames(:,1:end-1); audioFrames(:,2:end)];
        procFFT         = fft(procFrames, params.nfft);
        
        real_ch1 = real(procFFT(1:params.nfft/2+1,:));
        imag_ch1 = imag(procFFT(1:params.nfft/2+1,:));
        
        mag_procFFT     = abs(procFFT).^2;
        STFTMdb1        = 10*log(mag_procFFT);        
        STFTMdb         = STFTMdb1(1:params.nfft/2+1,:);
                
        %% labels
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
        
%         feature_vec = [melFeatures' Z Final_lables'];
        input = [real_ch1;imag_ch1];
        feature_vec = [input' Final_lables'];
        
        if count == 1
            trainData(1:count*len_1,:) = feature_vec;
            old_feature_vec = feature_vec;
        else
            trainData((count-1)*len_1+1:count*len_1,:) = feature_vec;
        end
    end
   
end

filename1 = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\Features\','clean_','RI_features');
save(filename1,'trainData','-v7.3');
disp('DONE!!!')