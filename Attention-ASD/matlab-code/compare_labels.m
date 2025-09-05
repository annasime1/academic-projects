clear all
%therapistdata
%patientdata
% File paths
%Zannino
%Dalfonso_Simeone
file1name = 'C:\Users\dalfo\Desktop\uni\magistrale\neuro\proj\Labeling\ESCS-R\Cavigliano_Lazarevski\Cavigliano\patientdata.txt';
file2name = 'C:\Users\dalfo\Desktop\uni\magistrale\neuro\proj\Labeling\ESCS-R\Lazarevski_Cavigliano\Lazarevski\patientdata.txt';
outputfile = 'C:\Users\dalfo\Desktop\uni\magistrale\neuro\proj\Labeling\ESCS-R\Cavigliano_Lazarevski\final\patientdata.txt';

% Import files
opts1 = detectImportOptions(file1name, 'Delimiter', ' '); 
file1 = readtable(file1name, opts1);

opts2 = detectImportOptions(file2name, 'Delimiter', ' '); 
file2 = readtable(file2name, opts2);

% Remove duplicates in both files
file1 = unique(file1, 'rows', 'stable');
file2 = unique(file2, 'rows', 'stable');

% Check for extra rows and print messages
if height(file1) > height(file2)
    [extraRows, iNON] = setdiff(file1(:, 1), file2(:, 1)); % Find rows in file1 not in file2
    if ~isempty(extraRows)
        disp('extra');
        if any(iNON == 1)
            disp('in first');
        end
    end
    file1(iNON, :) = []; % Remove extra rows from file1
elseif height(file2) > height(file1)
    [extraRows, iNON] = setdiff(file2(:, 1), file1(:, 1)); % Find rows in file2 not in file1
    if ~isempty(extraRows)
        disp('extra');
        if any(iNON == 1)
            disp('in first');
        end
    end
    file2(iNON, :) = []; % Remove extra rows from file2
end

% Combine data (ensure alignment)
[C,ia]=setdiff(file1,file2);
B=file2(ia,:);
final = file1; % Remove duplicates after combining
for i = 1:size(C, 1)
    % Find the corresponding index in file2
    idxFile2 = ia(i);
    
    % Check if the second column of the row from file2 contains the word 'nowhere'
    if string(file2{idxFile2, 2})=='nowhere'
        % If the row in file2 contains 'nowhere', overwrite the row in file1
        final(ia(i), :) = file2(idxFile2, :);
    end
end
% Write final table with full precision
fid = fopen(outputfile, 'w');
for i = 1:height(final)
    fprintf(fid, '%.7f %s\n', final{i, 1}, final{i, 2}{1}); % Adjust format as needed
end
fclose(fid);

disp('File saved with full precision.');
