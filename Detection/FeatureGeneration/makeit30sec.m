clear
close all
clc
fs2 = 16000;
folder = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Traffic\Noisy_Long_Files';
addpath('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Traffic\Noisy_Long_Files')
audio_files=dir(fullfile(folder,'*.wav'));

for n = 1:10
    filename=audio_files(n).name;
    [clean, fs] = audioread(filename);
    fs_16 = 16000;
    [p,q]              = rat(fs_16/fs,0.001);
    Sam_Sig            = resample(clean,p,q);
    len=floor(30*fs_16);
    Nframes = floor(length(Sam_Sig)/len)-1;
    k=1;
    for n1 = 1:Nframes
        final = Sam_Sig(k:k+len-1,:);
        filename_1 = strcat('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Traffic\Small_files\Part_1\', filename(1:end-4),'_',num2str(n1),'.wav');
        audiowrite(filename_1,final,fs_16);
        k = k+len;
    end
end
disp('DONE');