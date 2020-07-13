clear
clc

addpath('C:\Users\gxs160730\Documents\Ph.D\CNN-SE\CNN-VAD\Training Code\Functions\');
addpath('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Features');

Speech_path             = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\clean_part1';
Siren_path              = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Siren Noise for database\pixel siren';

speechFiles             = dir(Speech_path);
sirenFiles              = dir(Siren_path);

data.speechFiles        = speechFiles(~ismember({speechFiles.name},{'.' '..'}));
data.sirenFiles         = sirenFiles(~ismember({sirenFiles.name},{'.' '..'}));

data.nSpeechFiles       = numel(data.speechFiles);
data.nSirenFiles        = numel(data.sirenFiles);

Folder   = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\clean_part1');
FileList = dir(fullfile(Folder, '*.wav'));
N        = data.nSpeechFiles;
index    = randperm(numel(FileList), N);

M        = floor(data.nSpeechFiles/data.nSirenFiles);

dirname  = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Data\Clean\clean_siren';
cnt      = 1;
for n    = 1:data.nSirenFiles
    sirenFile           = [data.sirenFiles(n).folder '\' data.sirenFiles(n).name];
    [sirenSig_2,fs_si]  = audioread(sirenFile);
    [p,q]               = rat(16000/fs_si,0.001);
    sirenSig_1          = resample(sirenSig_2,p,q);
    sirenSig            = [sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1;sirenSig_1];
    sirenSig_ch1  =  sirenSig(:,1);
    sirenSig_ch2  =  sirenSig(:,2);
    for m = cnt:cnt+M
        audiofile = [data.speechFiles(m).folder '\' data.speechFiles(m).name];
        
        [speech,fs] = audioread(audiofile);
        
        speech_ch1 = speech(:,1);
        speech_ch2 = speech(:,2);
        disp (m)
        [speech_noise1]         = addnoise(speech_ch1,sirenSig_ch1,5,16000);
        [speech_noise2]         = addnoise(speech_ch2,sirenSig_ch2,5,16000);
        speech_noise           = [speech_noise1,speech_noise2];
        speech_noise = speech_noise./max(abs(speech_noise));
        name_siren1            = strcat(dirname,'\',num2str(n),'_',data.speechFiles(m).name);
        audiowrite(name_siren1,speech_noise,16000);
    end
    cnt=cnt+M+1;
end
