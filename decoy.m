function decoy()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               decoy.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. DESCRIPTION
% decoy.m generates decoy protein sequences for target-decoy based FDR
% estimation. The input file contains target protein sequences in fasta
% format. The sequences are read by the fastaread function and reversed.
% All decoy sequences will have a "REV_" suffix.
%
% 2. INPUT SYNTAX
% decoy()
% 
% 3. OUTPUT
% The output is a .fasta file containing target and decoy sequences
%
% 4. DEPENDENCIES
% - Bioinformatics Toolbox for the fastaread function
% - Statistics & Machine Learning toolbox for randsample function
%
% Created: 18/7/2019
% Last update: 22/7/2019
% (c) Joris Meurs, MSc

clc
% Choose decoy type
decoyInput = questdlg('Choose a decoy type','Decoy type','Reverse','Random','Reverse');
switch decoyInput
    case 'Reverse'
       decoyType = 1; 
       fprintf('Decoy type: reverse\n'); 
    case 'Random'
       decoyType = 2; 
       fprintf('Decoy type: random\n'); 
    otherwise
        return
end

% Locate .fasta file
[fastaFile,fastaPath] = uigetfile('.fasta',...
    'Select .FASTA File',...
    'C:\Users\paxjm9\Documents\FASTA Files');
if isequal(fastaFile,0)
    return
end
fastaLocation = fullfile(fastaPath,fastaFile);

% Read .fasta file
[headers,sequences] = fastaread(fastaLocation);

% Generate reverse sequences
if decoyType == 1
    reverseHeaders = [];
    reverseSequences = [];
    if ~iscell(headers) % Only one target sequence
        temp_header = []; temp_sequence = [];
        temp_header = headers;
        temp_sequence = sequences;
        entryLoc = [];
        entryLoc = find(temp_header=='|');
        reverseHeaders = [temp_header(1:entryLoc(2)-1) '_REVERSED' temp_header(entryLoc(2):end)];
        reverseSequences = fliplr(temp_sequence);
    else % Multiple sequences
        wb = waitbar(0,'Generating decoy sequences...');
        for j = 1:length(headers)
            waitbar(j/length(headers),wb,'Generating decoy sequences...');
            temp_header = []; temp_sequence = [];
            temp_header = headers{j};
            temp_sequence = sequences{j};
            entryLoc = [];
            entryLoc = find(temp_header=='|');
            reverseHeaders{j} = [temp_header(1:entryLoc(2)-1) '_REVERSED' temp_header(entryLoc(2):end)];
            reverseSequences{j} = fliplr(temp_sequence);
        end
    end
end

% Randomise protein sequences
if decoyType == 2
    randomHeaders = [];
    randomSequences = [];
    if ~iscell(headers) % Only one target sequence
        temp_header = []; temp_sequence = [];
        temp_header = headers;
        temp_sequence = sequences;
        entryLoc = [];
        entryLoc = find(temp_header=='|');
        randomHeaders = [temp_header(1:entryLoc(2)-1) '_RANDOM' temp_header(entryLoc(2):end)];
        randomSequences = randsample(temp_sequence,length(temp_sequence));
    else % Multiple sequences
        wb = waitbar(0,'Generating decoy sequences...');
        for j = 1:length(headers)
            waitbar(j/length(headers),wb,'Generating decoy sequences...');
            temp_header = []; temp_sequence = [];
            temp_header = headers{j};
            temp_sequence = sequences{j};
            entryLoc = [];
            entryLoc = find(temp_header=='|');
            randomHeaders{j} = [temp_header(1:entryLoc(2)-1) '_RANDOM' temp_header(entryLoc(2):end)];
            randomSequences{j} = randsample(temp_sequence,length(temp_sequence));
        end
    end 
end

% Write all sequences to new .fasta file
extensionLoc = find(fastaFile=='.');
strippedFileName = fastaFile(1:extensionLoc-1);
fileID = fopen([strippedFileName '_decoy.fasta'],'w');

if ~iscell(headers) % Handle single target sequence
    % Header line
    fprintf(fileID,'>%s\n',headers);
    % Sequence
    fprintf(fileID,'%s\n',sequences);
    
    % Deal with decoy type
    if decoyType == 1
        % Reverse header line
        fprintf(fileID,'>%s\n',reverseHeaders);
        % Reverse sequence
        fprintf(fileID,'%s\n',reverseSequences);
    end
    if decoyType == 2
        % Reverse header line
        fprintf(fileID,'>%s\n',randomHeaders);
        % Reverse sequence
        fprintf(fileID,'%s\n',randomSequences);        
    end
else  
    % Write first actual sequences
    waitbar(0,wb,'Writing target sequences to file...');
    for j = 1:length(sequences)
       waitbar(j/length(sequences),wb,'Writing target sequences to file...'); 
       % Header line
       fprintf(fileID,'>%s\n',headers{j});
       % Sequence
       fprintf(fileID,'%s\n',sequences{j});
    end

    % Write decoy sequences
    waitbar(0,wb,'Writing decoy sequences to file...');
    % Deal with decoy type
    if decoyType == 1
        for j = 1:length(reverseSequences)
           waitbar(j/length(sequences),wb,'Writing decoy sequences to file...');
           % Reverse header line
           fprintf(fileID,'>%s\n',reverseHeaders{j});
           % Reverse sequence
           fprintf(fileID,'%s\n',reverseSequences{j}); 
        end
    end
    if decoyType == 2
        for j = 1:length(randomSequences)
           waitbar(j/length(sequences),wb,'Writing decoy sequences to file...');
           % Reverse header line
           fprintf(fileID,'>%s\n',randomHeaders{j});
           % Reverse sequence
           fprintf(fileID,'%s\n',randomSequences{j}); 
        end        
    end
    delete(wb);
end
fclose(fileID);
fclose('all');
fprintf('Decoy(s) generated \n');
fileLoc = fullfile(cd,[strippedFileName '_decoy.fasta']);
fprintf('File location: %s \n',fileLoc);

end
