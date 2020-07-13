clear
clc

Dest      = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Machinery\Small_files\Part_2';
for m = 1:1
    
    Folder   = 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\Detection\Noisy_Speech\Machinery\Small_files\Part_1';
    FileList = dir(fullfile(Folder, '*.wav'));
    N = 605;
   
    index    = randperm(numel(FileList), N);
    for k = 1:N
        Source = fullfile(Folder, FileList(index(k)).name);
        disp(Source);
        movefile(Source, Dest);
    end
end
disp('DONE!!!')
%copyfile(Source, 'C:\Users\gxs160730\Documents\Ph.D\Siren Noise\speech+noise_M5\Test');