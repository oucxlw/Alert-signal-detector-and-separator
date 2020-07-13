clear
clc

addpath('C:\Users\gxs160730\Documents\Ph.D\CNN-SE\CNN-VAD\Training Code\');
addpath('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Features');

Speech_path             = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Separation\NoisySpeech\Traffic\30sec_Noisy_Speech\605_part1';
Siren_path              = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Siren Noise for database\pixel siren';

speechFiles             = dir(Speech_path);
sirenFiles              = dir(Siren_path);

data.speechFiles        = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.sirenFiles         = sirenFiles(~ismember({sirenFiles.name},{'.' '..'}));

data.nSpeechFiles       = numel(data.speechFiles);
data.nSirenFiles        = numel(data.sirenFiles);



Folder   = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Separation\NoisySpeech\Traffic\30sec_Noisy_Speech\605_part1');
FileList = dir(fullfile(Folder, '*.wav'));
Cell = struct2cell(FileList);
Cell = Cell(1,:);
Cell = natsort(Cell);
N        = data.nSpeechFiles;

M        = floor(data.nSpeechFiles/data.nSirenFiles);

dirname  = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Small clean files\605_siren+clean';
cnt      = 1;
len_1 = 1874;
trainData = zeros(1874*N,1543); 
for n    = 1:data.nSirenFiles

    sirenFile           = [data.sirenFiles(n).folder '\' data.sirenFiles(n).name];
    [sirenSig_2,fs_si]  = audioread(sirenFile);
    [p,q]               = rat(16000/fs_si,0.001);
    sirenSig_1          = resample(sirenSig_2(:,1),p,q);
    sirenSig            = [sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1];
    
    for m = cnt:cnt+M
        if m<=N
            current_name = fullfile(Folder, Cell(m));
            current_name = string(current_name);
            [filepath,filename,fileext] = fileparts(current_name);

            [speech,fs] = audioread(current_name);
            speech1 = speech(:,1);
            speech = speech(:,2);
            
            sirenSig = sirenSig(1:length(speech));
            disp (m)

            speech_noise           = addnoise(speech,sirenSig,10,16000);
            speech_noise1          = addnoise(speech1,sirenSig,10,16000);
            speech_noise           = speech_noise/max(abs(speech_noise));
            name_siren1            = strcat(dirname,'\',num2str(n),'_',filename,fileext);
            
%             audiowrite(char(name_siren1),speech_noise,16000);
            %% features        
            params.overlap           = 256;
            params.nfft              = 512;
            params.window            = hanning(params.nfft);

            %noisy
            audioFrames_noi     = reshape(speech_noise, params.overlap, []);
            procFrames_noi      = params.window .* [audioFrames_noi(:,1:end-1); audioFrames_noi(:,2:end)];
            procFFT_noi         = fft(procFrames_noi, params.nfft);
            real_noi = real(procFFT_noi(1:params.nfft/2+1,:));
            imag_noi = imag(procFFT_noi(1:params.nfft/2+1,:));
            feature1 = [real_noi;imag_noi];
            
            audioFrames_noi1     = reshape(speech_noise1, params.overlap, []);
            procFrames_noi1      = params.window .* [audioFrames_noi1(:,1:end-1); audioFrames_noi1(:,2:end)];
            procFFT_noi1         = fft(procFrames_noi1, params.nfft);
            real_noi1 = real(procFFT_noi1(1:params.nfft/2+1,:));
            imag_noi1 = imag(procFFT_noi1(1:params.nfft/2+1,:));
            feature2 = [real_noi1;imag_noi1];
%             mag_procFFT_noi     = abs(procFFT_noi).^2;
%             STFTMdb1_noi        = 10*log(mag_procFFT_noi);
%             STFTMdb_noi         = STFTMdb1_noi(1:params.nfft/2+1,:);

            %siren
            audioFrames_si     = reshape(sirenSig, params.overlap, []);
            procFrames_si      = params.window .* [audioFrames_si(:,1:end-1); audioFrames_si(:,2:end)];
            procFFT_si         = fft(procFrames_si, params.nfft);
            real_si = real(procFFT_si(1:params.nfft/2+1,:));
            imag_si = imag(procFFT_si(1:params.nfft/2+1,:));
%             mag_procFFT_si     = abs(procFFT_si).^2;
%             STFTMdb1_si        = 10*log(mag_procFFT_si);
%             STFTMdb_si         = STFTMdb1_si(1:params.nfft/2+1,:);
            output2 = [real_si;imag_si];
            %clean
            audioFrames     = reshape(speech, params.overlap, []);
            procFrames      = params.window .* [audioFrames(:,1:end-1); audioFrames(:,2:end)];
            procFFT         = fft(procFrames, params.nfft);
            real_s = real(procFFT(1:params.nfft/2+1,:));
            imag_s = imag(procFFT(1:params.nfft/2+1,:));
            output1 = [real_s;imag_s];

            count = m;
            file_num = m*ones(1,1874);
            feature_vec = [feature1', feature2', output2',file_num'];
            if count == 1
                trainData(1:count*len_1,:) = feature_vec;
                old_feature_vec = feature_vec;
            else
                trainData((count-1)*len_1+1:count*len_1,:) = feature_vec;
            end
            
        end
    end
        cnt=cnt+M+1;
    
end
disp('DONE!!!')

filename1 = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Separation\Features\Traffic\Sep_T_0','_','RI_1_OP');
save(filename1,'trainData','-v7.3');
disp('DONE!!!')

%%%%%% TESTING FEATURES %%%%%%
% 
% trainData = zeros(len_1,2057); 
% 
% for n    = 1:data.nSirenFiles
% 
%     sirenFile           = [data.sirenFiles(n).folder '\' data.sirenFiles(n).name];
%     [sirenSig_2,fs_si]  = audioread(sirenFile);
%     [p,q]               = rat(16000/fs_si,0.001);
%     sirenSig_1          = resample(sirenSig_2(:,1),p,q);
%     sirenSig            = [sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1];
%     
%     for m = cnt:cnt+M
%         if m<=N
%             current_name = fullfile(Folder, Cell(m));
%             current_name = string(current_name);
%             [filepath,filename,fileext] = fileparts(current_name);
% 
%             [speech,fs] = audioread(current_name);
%             speech = speech(:,2);
%             speech1 = speech(:,1);
%             sirenSig = sirenSig(1:length(speech));
%             disp (m)
% 
%             speech_noise           = addnoise(speech,sirenSig,10,16000);
%             speech_noise1          = addnoise(speech1,sirenSig,10,16000);
%             speech_noise           = speech_noise/max(abs(speech_noise));
%             name_siren1            = strcat(dirname,'\',num2str(n),'_',filename,fileext);
% %             audiowrite(char(name_siren1),speech_noise,16000);
%             %% features        
%             params.overlap           = 256;
%             params.nfft              = 512;
%             params.window            = hanning(params.nfft);
% 
%             %noisy
%             audioFrames_noi     = reshape(speech_noise, params.overlap, []);
%             procFrames_noi      = params.window .* [audioFrames_noi(:,1:end-1); audioFrames_noi(:,2:end)];
%             procFFT_noi         = fft(procFrames_noi, params.nfft);
%             real_noi = real(procFFT_noi(1:params.nfft/2+1,:));
%             imag_noi = imag(procFFT_noi(1:params.nfft/2+1,:));
%             feature1 = [real_noi;imag_noi];
%             
%             audioFrames_noi1     = reshape(speech_noise1, params.overlap, []);
%             procFrames_noi1      = params.window .* [audioFrames_noi1(:,1:end-1); audioFrames_noi1(:,2:end)];
%             procFFT_noi1         = fft(procFrames_noi1, params.nfft);
%             real_noi1 = real(procFFT_noi1(1:params.nfft/2+1,:));
%             imag_noi1 = imag(procFFT_noi1(1:params.nfft/2+1,:));
%             feature2 = [real_noi1;imag_noi1];
% %             mag_procFFT_noi     = abs(procFFT_noi).^2;
% %             STFTMdb1_noi        = 10*log(mag_procFFT_noi);
% %             STFTMdb_noi         = STFTMdb1_noi(1:params.nfft/2+1,:);
% 
%             %siren
%             audioFrames_si     = reshape(sirenSig, params.overlap, []);
%             procFrames_si      = params.window .* [audioFrames_si(:,1:end-1); audioFrames_si(:,2:end)];
%             procFFT_si         = fft(procFrames_si, params.nfft);
%             real_si = real(procFFT_si(1:params.nfft/2+1,:));
%             imag_si = imag(procFFT_si(1:params.nfft/2+1,:));
% %             mag_procFFT_si     = abs(procFFT_si).^2;
% %             STFTMdb1_si        = 10*log(mag_procFFT_si);
% %             STFTMdb_si         = STFTMdb1_si(1:params.nfft/2+1,:);
%             output2 = [real_si;imag_si];
%             %clean
%             audioFrames     = reshape(speech, params.overlap, []);
%             procFrames      = params.window .* [audioFrames(:,1:end-1); audioFrames(:,2:end)];
%             procFFT         = fft(procFrames, params.nfft);
%             real_s = real(procFFT(1:params.nfft/2+1,:));
%             imag_s = imag(procFFT(1:params.nfft/2+1,:));
%             output1 = [real_s;imag_s];
% 
%             count = m;
%             file_num = m*ones(1,1874);
%             feature_vec = [feature1', feature2', output1', output2',file_num'];
%             if count == 1
%                 trainData(1:count*len_1,:) = feature_vec;
%                 old_feature_vec = feature_vec;
%             else
%                 trainData((count-1)*len_1+1:count*len_1,:) = feature_vec;
%             end
%             
%         end
%     end
%         cnt=cnt+M+1;
%     
% end
% disp('DONE!!!')
% 
% 
% filename1 = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\siren_features\Separation_M_10','_','Test_RI_features');
% save(filename1,'trainData','-v7.3');
% disp('DONE!!!')