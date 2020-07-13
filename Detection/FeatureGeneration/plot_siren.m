%Plot signal
clear
clc
[x,fs] = audioread('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Machinery\M5\small_files\Test_features\M_5_A_5_MZ_Test_features_81.wav');
t = 1/16000*(1:length(x));
load ('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Babble\B0\B0A5\Outputs\B0A5_MZ1_28_pred.mat');
load ('C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Babble\B0\B0A5\Outputs\B0A5_MZ1_28_labels.mat');
% labels = feature_vec(:,15);
% plot(feature_vec)
% figure;
% plot(Y1')
% % plot(Y1)
% xlabel('time(s)')
% ylabel('Amplitude')

mul_ones = ones(1,256);
Y1 = double(Y1');
labels = double(labels);
k=1;
Y2=zeros(length(Y1)*256,1);
for i = 1: length(Y1)
    ans1 = Y1(i).*mul_ones';
    Y2(k:k+255,:) = ans1;
    k = k+255;
end

result = labels==Y1;
k2 = find(~result);

% plot(x)
% hold on;
% plot(Y2)
% xlabel('time(s)')
% ylabel('Amplitude')


adder = labels + Y1;
TP = length(find(adder == 2));
TN = length(find(adder == 0));
subtr = labels - Y1;
FP = length(find(subtr == -1));
FN = length(find(subtr == 1));


Acc = (TP+TN) /(TP+TN+FP+FN);
TPR = TP/(TP+FN);
FPR = FP /(FP + TN);
FNR = FN / (TP+FN);



% PLOT WAIL HI-LO ETC
% [x,fs] = audioread('C:\Users\gxs160730\Desktop\HiLo_Wail_Yelp\Pulsed_alarm_2.mp3');
% % plot(x);
% [noise,fs_n] = audioread('C:\Users\gxs160730\Documents\MATLAB\TOP 10 NOISES FROM DATABASE\Machinery\Machinery120.wav');
% % figure;
% % plot(noise);
% [p1,q1]                = rat(fs_n/fs,0.001);
% siren              = resample(x,p1,q1);
% addpath('C:\Users\gxs160730\Documents\Ph.D\CNN-SE\CNN-VAD\Training Code\');
% noisy =  addnoise(siren,noise,10,16000);
% plot(noisy);    
% % v = randn
% % for i=1:length(x)
% %     if(x(i)==0)
% %         x(i) = eps*10^6;
% %     end
% %     
% % end
% out = awgn(x,90);
% spectrogram(20*x(1:0.8*fs),512,500,512,fs,'yaxis'); xlabel('Time (in secs)'); ylabel('Frequency'); title('Pulsed alarms');colormap('parula');
