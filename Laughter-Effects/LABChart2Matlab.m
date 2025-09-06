%% LABChart 2 Matlab

clear; close all; clc;

data_path = 'C:\Users\asnag\OneDrive - Università degli Studi di Milano\LT\LT MATLAB\Original\T0\';
d = dir(strcat(data_path,'*.mat'));

save_path = 'C:\Users\asnag\OneDrive - Università degli Studi di Milano\LT\LT MATLAB\Raw\T0\';

%%
for n = 1:length(d)
    
    %d(n) %Sanity Check
    data = load(strcat(data_path,d(n).name));
    event_marker = [data.com(:,3), data.com(:,5)];
    event_name = data.comtext;
    
    data_raw = data.data(data.datastart(1):data.datastart(2));
    for i=2:length(data.datastart)
        if i==6
            data_raw = [data_raw; data.data(data.datastart(i):end), data.data(end)];
        else
            data_raw = [data_raw; data.data(data.datastart(i):data.datastart(i+1))];
        end
    end
    data_raw = data_raw';
    writematrix(data_raw, strcat(save_path,d(n).name(1:end-4),'_raw.csv'))
    writematrix(event_marker, strcat(save_path,d(n).name(1:end-4),'_event.csv'))
    writematrix(event_name, strcat(save_path,d(n).name(1:end-4),'_event_name.csv'))

end