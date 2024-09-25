
%clear all

Compare_all_results = true;  % true if you want to compare all the results in ./output
                             
% if Compare_all_results = false, list the models you want to compare
Compared_list       = [string('model_name'), ...
                       string('model_name'), ...
                       string('model_name')];
                    %   Note: the output files are named by
                    %         ./output/{model_name}-STF.dat
                    %   and
                    %         ./output/{model_name}-spectrum.dat
                    


%-----------------------------------------------------
Work_dir = erase(which('Compare_spectrum.m'),'Compare_spectrum.m');
data_dir = [Work_dir '/output/'];

if(Compare_all_results)
    files     = dir(fullfile(data_dir,'*-spectrum.dat'));
else
    files = [];
    for i=1:length(Compared_list)
       model_name = strcat(Compared_list(i),'-spectrum.dat');
    files     = [files dir([char(data_dir) '/' char(model_name)])];
    end
end

model_num = length(files);
if(model_num==0)
    disp('There is no result in the ./output folder!')
end
for k = 1:model_num
files(k).name = erase(files(k).name,'-spectrum.dat');
end

clf
figure(1)
subplot(211)
for k = 1:model_num
    filename  = [data_dir files(k).name '-STF.dat'];
    data = importdata(filename);
    plot(data(:,1),data(:,2));
    hold on
end
legend(files(:).name)
title('Source time function')
xlabel('time (s)')
ylabel('STF')

subplot(212)
for k = 1:model_num
    filename  = [data_dir files(k).name '-spectrum.dat'];
    data = importdata(filename);
    loglog(data(:,1),data(:,2)/data(2,2));
    hold on
end
legend(files(:).name)
title('Spectrum of STF')
xlabel('Frequency (Hz)')
ylabel('Normalized spectrum')
xlim([0.005,10]);
ylim([10^-6,1])
% 
% loglog([0.04,40],[1,0.000001],'k');
% hold on
% loglog([0.04,40],[1,0.001],'k');
% hold off