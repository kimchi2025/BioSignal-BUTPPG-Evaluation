% Read Csv file
% all 48 original files

T = readtable('subj_hr_info.csv');

% Notes:
% - T.Ear_finger % 0 (ear), 1 (finger)
% 
% - T.Motion 
%    0 (rest),  1 (higher pressure of finger/ear on camera),
%    2 (moving the finger/ear on the lens), 3 (walking), 4 (coughing), 
%    5 (laughing), 6 (changing the light), 7 (talking)
%    {8 (2 and 3 together). personally added 8 because it would return NaN
%    otherwise and they were marked as 2;3 in the file}
% 
% - T.Quality %1=good, 0=poor

%% Make a New Table with Our Data
% may seem repetitive but this is because reading samples may take a while
% therefore this is to preserve the data at multiple steps as a "backup"

% number of singals to be read in currently
% values used for paper: 51, 3888, and 48
% this was done multiple times, but the most recent one is what the code is set up for
n = 48; 

% preallocate arrays
ID = zeros(n,1);
RefHR = zeros(n,1);
Location = zeros(n, 1);
Motion = zeros(n,1);
Quality = zeros(n,1);
signal = zeros(300, n);
time = zeros(300,n);

% for file reading
str2 = "butppg";
str3 = "/";
str4 = "_PPG";

for i=1:n
    str1 = T.ID(i);
    fileName = str2 + str3 + num2str(str1) + str3 + num2str(str1) + str4;

    % check if fileName file doesn't exist and skip if true
    if ~isfile((fileName + ".dat"))
        continue
    end
    % assign variables
    ID(i) = str1;
    RefHR(i) = T.HR(i);
    Location(i) = T.Ear_finger(i);
    Motion(i) = T.Motion(i);
    Quality(i) = T.Quality(i);

    [sig,fs,tm] = rdsamp(char(fileName));

    % old files have different formatting than new files. I make them the
    % same here. 300 x 3 (new). old were 1 x 300
    if mean(tm) == 0
        sig = sig';
        tm = linspace(0,length(sig)/fs,length(sig));
        tm = tm';
    end
    if size(sig,2) > 1
        sig = sig(:,1); % take channel 1
    end
    signal(:,i) = sig;
    time(:,i) = tm;
    disp(['Loaded ', num2str(i),' to array. Out of ', num2str(n)]) %Progress update.
    disp([num2str((i/n)*100),'% done.'])

end


T48 = table(ID,RefHR,Location,Motion,Quality);
T48.Signal = signal';
T48.Time = time';
to_delete = (T48.ID == 0);
T48(to_delete,:) = [];
% to access signal do new_T.Signal(1,:) to access the row with the signal

%% Save as CSV file for backup (exclude Signal and Time depending on need)
backup = T48;
% separator = ' ';
% backup.Signal = cellfun(@(x) strjoin(arrayfun(@num2str, x, 'UniformOutput', false), ' '), num2cell(backup.Signal, 2), 'UniformOutput', false);
% backup.Time = cellfun(@(x) strjoin(arrayfun(@num2str, x, 'UniformOutput', false), ' '), num2cell(backup.Time, 2), 'UniformOutput', false);
writetable(backup,'Backup48Data.csv','Delimiter',' ');

%% How you would read it from CSV
% and turn the strings back to numbers

PPG_T = readtable('myData.csv');

% call like this to get the signal vector. replace 1 with the desired row
% more efficient to save the table as a .mat file than do this however. which is what we did instead
example_sig2 = str2num(cell2mat(hi.Signal(1)));