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
% 2. INPUT
% To call the function type "decoy()" in the command window. No input 
% arguments required 
% 
% 3. OUTPUT
% The output is a .fasta file containing target and decoy sequences
%
% 4. DEPENDENCIES
% - Bioinformatics Toolbox for the fastaread function
%
% Created: 18/7/2019
% Joris Meurs, MSc

% Deal with input arguments
if nargin > 0
    error('No input allowed!');
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
reverseHeaders = [];
reverseSequences = [];
if ~iscell(headers) % Only one target sequence
    temp_header = []; temp_sequence = [];
    temp_header = headers;
    temp_sequence = sequences;
    reverseHeaders = ['REV_' temp_header];
    reverseSequences = fliplr(temp_sequence);
else % Multiple sequences
    wb = waitbar(0,'Generating decoy sequences...');
    for j = 1:length(headers)
        waitbar(j/length(headers),wb,'Generating decoy sequences...');
        temp_header = []; temp_sequence = [];
        temp_header = headers{j};
        temp_sequence = sequences{j};
        reverseHeaders{j} = ['REV_' temp_header];
        reverseSequences{j} = fliplr(temp_sequence);
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
    
    % Reverse header line
    fprintf(fileID,'>%s\n',reverseHeaders);
    % Reverse sequence
    fprintf(fileID,'%s\n',reverseSequences);
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
    for j = 1:length(reverseSequences)
       waitbar(j/length(sequences),wb,'Writing decoy sequences to file...');
       % Reverse header line
       fprintf(fileID,'>%s\n',reverseHeaders{j});
       % Reverse sequence
       fprintf(fileID,'%s\n',reverseSequences{j}); 
    end
    delete(wb);
end
fclose(fileID);
fclose('all');
fprintf('Decoy(s) generated \n');
fileLoc = fullfile(cd,[strippedFileName '_decoy.fasta']);
fprintf('File location: %s \n',fileLoc);

end
