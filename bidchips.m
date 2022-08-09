% 
% **************************** BIDCHIPS ver 1.0 ******************************
%
% This MATLAB function file implements the BIDCHIPS package, which can be used
% for the quantification and removal of biases in ChIP-seq datasets resulting in
% a purified ChIP-seq binding signal.
%
% The package accompanies the paper titled "BIDCHIPS: Bias decomposition and
% removal from ChIP-seq data clarifies true binding signal and its functional
% correlates" by P Ramachandran, GA Palidwor, and TJ Perkins, Epigenetics &
% Chromatin, 2015, 8:33
%
% Copyright (c) 2015 by the authors of the above paper.  All rights reserved.
% 
% ---------------------------------------------------------------------------------
% Citation: If you use this package, please cite: BIDCHIPS: bias decomposition
% and removal from ChIP-seq data clarifies true binding signal and its
% functional correlates, Parameswaran Ramachandran, Gareth A. Palidwor, and
% Theodore J. Perkins, Epigenetics & Chromatin 2015, 8:33
% ---------------------------------------------------------------------------------
% 
% Licence: This code is free to use for non-commerical purposes.
%
% FEEDBACK: For questions, comments, e-mail rpara26@gmail.com, tperkins@ohri.ca
% 
% Posted on www.perkinslab.ca in Sept. 2015
%
% USAGE: 
% --------- 
% 
% The software takes in a list of inputs including the names of the ChIP-seq and
% other short-read sequencing data files (in BED format), preprocesses them to
% remove duplicate reads, and then builds a genome-wide background model to
% quantify the contributions of different bias signals. Then, if the user inputs
% a set of genomic intervals, the software predicts the background levels within
% these intervals in terms of read counts. Finally, the true binding signal
% levels are estimated in these intervals in terms of a set of purified read
% counts. 
% 
% NOTES: 
% --------
%   (1) Due to significantly large memory requirements (RAM), we recommend
%       running this software in a cloud-based or cluster computing environment
% 
%   (2) The software stores the regression coefficients of the background model
%       in a MAT file. If the user had input a filename containing a set of
%       prediction intervals, the software carries out predictions. Else, it
%       terminates.
% 
%   (3) All subfunctions necessary to use this code are contained in this file. 
% 
% The basic call to the main function is as follows:
% 
% bidchips(DataFldr, MappabFldr, GenmFldr, WinSz, ReadLength, ChrListFileName, PredictIntvlsFileName, TrtFileName, DNaseIFileName, iDNAFileName, IgGFileName)
% 
% INPUTs:
% ----------
% The input parameters are defined as:
% 
% DataFldr              -> the full path to the folder/directory containing
%                          all input data files, i.e., the files corresponding
%                          to the last six input arguments
%                          
% MappabFldr            -> the full path to the folder/directory containing the
%                          mappability intervals for the chosen read length and
%                          organism. Only base pairs that are perfectly mappable
%                          should be listed (i.e., with a mappability data value
%                          of 1). Each chromosome should correspond to a
%                          separate file (with filename extension *.map), and
%                          the file content should have the format:
%                          ---
%                          track type=wiggle_0 name="CRG Align 36" description="Alignability of 36mers by GEM from ENCODE/CRG(Guigo)"
%                          chr22	16054663	16054692	1
%                          chr22	16054704	16054714	1
%                          chr22	16055096	16055100	1
%                          ---
%                          These *.map files (e.g., chr1.map) can be downloaded
%                          as follows:
%                          Go to the UCSC Table Browser at https://genome.ucsc.edu/cgi-bin/hgTables
%                          Choose organism, genome, and assembly
%                          Choose group: "mapping and sequencing"; and track: "Mappability"
%                          Choose table: x-mer mappability (e.g., 36-mer)
%                          Choose region: genome; and position: a corresponding chr name
%                          filter by data value = 1
%                          output format: "data points"
%                          Choose an output filename: e.g., "chr1.map"
%                        
% GenmFldr              -> the full path to the folder/directory containing the
%                          genome sequence (one file per chr) as FASTA files
%                          (e.g., chr1.fa)
%                          Genomes can be downloaded from: http://hgdownload.soe.ucsc.edu/downloads.html
%                          
% WinSz                 -> the user-defined window size (in base pairs) to be used for building
%                          the background model (e.g., 129) 
%                          NOTE: Although the peaks/intervals would be of
%                          varying lengths, the background model will be built
%                          using a single window size input here! So, we suggest
%                          that the user choose this window size carefully,
%                          representing an average size that will capture the
%                          background reasonably well in all prediction
%                          intervals (if given)
% 
% ReadLength            -> the most common read length found in the BED files;
%                          this should also match with the downloaded
%                          mappability data; For example, if the ChIP-seq and
%                          the control datasets contain primarily (or wholly)
%                          36-bp reads, then this argument should be 36.
%                          Accordingly, the 36-bp mappability track should be
%                          downloaded 
%                          NOTE: In our experience, we have noticed that
%                          frequently, due to a number of reasons, all reads in
%                          a given ChIP-seq dataset would not be of a single
%                          read length. This is common in ENCODE data. So,
%                          although an approximation in some sense, it should be
%                          fair enough to take the most common read length found
%                          in all datasets involved as the value for this
%                          argument
% 
% ChrListFileName       -> the name of the text file containing the list of
%                          all the chromosomes and their sizes in base pairs
%                          corresponding to a given organism (e.g.,
%                          chrList.txt). This list should be exhaustive, i.e.,
%                          it should list all possible chromosomes for the
%                          organism. The software will then take a subset of
%                          this list for analysis, depending on whether all
%                          datasets contain all chromosomes or not                         
%                          For instance, the file should contain the following
%                          lines for the human genome (tab delimited):
%                          chr1  249250621  
%                          chr2  243199373
%                          chr3  198022430
%                          chr4  191154276
%                          chr5  180915260
%                          chr6  171115067
%                          chr7  159138663
%                          chr8  146364022
%                          chr9  141213431
%                          chr10 135534747
%                          chr11 135006516
%                          chr12 133851895
%                          chr13 115169878
%                          chr14 107349540
%                          chr15 102531392
%                          chr16 90354753
%                          chr17 81195210
%                          chr18 78077248
%                          chr19 59128983
%                          chr20 63025520
%                          chr21 48129895
%                          chr22 51304566
%                          chrX  155270560
%                          chrY  59373566 
% 
% PredictIntvlsFileName -> [OPTIONAL] a BED formatted file containing the intervals
%                                     where the purified ChIP-seq signal (in
%                                     terms of the read counts) needs to be
%                                     estimated. This file can be the output
%                                     from a peak-calling software (such as
%                                     MACS) or user-defined.
%                                     See
%                                     https://genome.ucsc.edu/FAQ/FAQformat.html#format1
%                                     for more about the BED format 
%                                     NOTE: If no prediction is required, then
%                                     simply input an empty vector, [].  PLEASE
%                                     DON'T SKIP THE INPUT!
% 
% TrtFileName           -> the main ChIP-seq treatment file name (e.g., ATF3.bed)
% 
% DNaseIFileName        -> the DNaseI hypersensitivity data file name (*.bed)
% 
% iDNAFileName          -> the input DNA control data file name (*.bed)
% 
% IgGFileName           -> the IgG control data file name (*.bed)
% 
% 
% OUTPUTs:
% ----------
% 
% .bed file   -> The output BED file containing the observed ChIP-seq signal
%                levels (in terms of read counts), the predicted background signal
%                levels (in terms of read counts), and the difference of these two
%                values (which is our "Purified" signal level).  The exact columns
%                are as follows:
%                | "Chr Name" | "0-based Start" | "1-based End" | "Observed NRds" | "Background NRds" | "Purified Signal" |
% 
% .mat & .txt -> A MAT file, as well as a TXT file with the computed linear
%                model coefficients, that can be used in future if desired
% 
% 

%---------- PART 1 - LOADING ALL DATA AND BUILDING BACKGROUND MODEL ----------------

function bidchips(DataFldr, MappabFldr, GenmFldr, WinSz, ReadLength, ChrListFileName, PredictIntvlsFileName, TrtFileName, DNaseIFileName, iDNAFileName, IgGFileName)
tic;
% Checking if prediction is necessary
if isempty(PredictIntvlsFileName)
    disp('No prediction intervals given; so only the background model will be built and saved!')
    YesforPredict = false;
else
    YesforPredict = true;
end

[~,otptnam,~] = fileparts(PredictIntvlsFileName);

disp(' ')
disp('***************************************')
disp('       --- BIDCHIPS  ver 1.0 ---       ')
disp('***************************************')
disp(' ')
disp('Inputs Given: ')
disp('---------------------')
disp(['Main Treatment       : ' TrtFileName])
disp(['DNaseI               : ' DNaseIFileName])
disp(['Input DNA Ctrl       : ' iDNAFileName])
disp(['IgG Ctrl             : ' IgGFileName])
if YesforPredict
    disp(['Prediction Intervals : ' PredictIntvlsFileName])
else
    disp('Prediction Intervals : NONE')
end
disp(['Chromosome List File : ' ChrListFileName])
disp(['Read Length          : ' num2str(ReadLength)])
disp(['Window Size          : ' num2str(WinSz)])
disp(['Data Folder          : ' DataFldr])
disp(['Mappability Folder   : ' MappabFldr])
disp(['Genome Folder        : ' GenmFldr])
disp(' ')

% Check whether DataFldr, MappabFldr, and GenmFldr end with a file separator; Add file separator if needed
if DataFldr(end) == '/' || DataFldr(end) == '\'
    DataFldr(end)= filesep;
else
    DataFldr(end+1)= filesep;
end
if MappabFldr(end) == '/' || MappabFldr(end) == '\'
    MappabFldr(end)= filesep;
else
    MappabFldr(end+1)= filesep;
end
if GenmFldr(end) == '/' || GenmFldr(end) == '\'
    GenmFldr(end)= filesep;
else
    GenmFldr(end+1)= filesep;
end
[~, justTrtname, ~] = fileparts(TrtFileName);
% Process Chr List & Sizes File
fid = fopen([DataFldr ChrListFileName]);
ListofchrsAndSzs = textscan(fid, '%s %f %*[^\n]','collectoutput',1);
% Cols that have been Read: {'Chr' 'size'}
status = fclose(fid);
if status ~= 0
    disp(['The ' DataFldr ChrListFileName ' text file did not close properly!'])
end
ListofchrsAndSzs = [ ListofchrsAndSzs{1} num2cell(ListofchrsAndSzs{2}) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading ALL BED files 
% 
% NOTE: The read start locations are incremented by 1 in the "ReadBED" function
% file to transform them to the one-based co-ordinate system. By default, as per
% the BED format, the start locations are in zero-based co-ordinate system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Treatment Data (Transcription factor or Histone mark ChIP-seq)
disp(['Reading Main TREATMENT data : ' TrtFileName])
[TrtReadsStruct, TrtchrList] = ReadBED(DataFldr, TrtFileName);

% DNaseI Data (Chromatin accessibility)
disp(['Reading DNaseI data : ' DNaseIFileName])
[DNaseIReadsStruct, DNaseIchrList] = ReadBED(DataFldr, DNaseIFileName);

% iDNA Data (Input DNA control)
disp(['Reading iDNA data : ' iDNAFileName])
[iDNAReadsStruct, iDNAchrList] = ReadBED(DataFldr, iDNAFileName);

% IgG Data (IgG control)
disp(['Reading IgG data : ' IgGFileName])
[IgGReadsStruct, IgGchrList] = ReadBED(DataFldr, IgGFileName);

%---------------------------------------------------------------------------------
% Getting a consensus chr list based on intersection of chrs from each dataset
% for which reads are available
% 
% We assume that if Plus reads are present, then Minus reads would also be
% present. We don't check for the existence of both Plus and Minus strand reads
% for each chr separately.
%---------------------------------------------------------------------------------
disp('Checking to pick only those chrs for which all datasets have reads ...')
chrIntersect = intersect(TrtchrList,DNaseIchrList);
chrIntersect = intersect(chrIntersect,iDNAchrList);
chrIntersect = intersect(chrIntersect,IgGchrList);
[lia,~] = ismember(ListofchrsAndSzs(:,1), chrIntersect);
ListofchrsAndSzs = ListofchrsAndSzs(lia,:);
disp('Final chromosome set to proceed with : ')
disp([{'Chromosome' 'Size'; '--------------------' '---------'} ; ListofchrsAndSzs]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading GENOME Sequence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading the genome sequence ...')
for jis = 1:size(ListofchrsAndSzs,1)
    curchr = ListofchrsAndSzs{jis,1};
    fprintf('%s ',curchr)
    fid = fopen([GenmFldr curchr '.fa']);
    ftext = textscan(fid,'%s','delimiter','\n','commentstyle','t');
    status = fclose(fid);
    if status ~= 0
        disp(['The ' GenmFldr curchr ' file did not close properly!'])
    end
    ftext = ftext{:}; % Getting rid of the extra cell layer
    Genm.(curchr) = upper(sprintf('%s', ftext{1:end}));
    % The above is the line that joins all the lines into a single line
    % Code idea to use sprintf is from strjoin function in File Exchange submission
end
fprintf('%s\n',' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Mappability files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading the mappability files ...')
for jis = 1:size(ListofchrsAndSzs,1)
    curchr = ListofchrsAndSzs{jis,1};
    fprintf('%s ',curchr)
    fid = fopen([MappabFldr curchr '.map']);
    MainCell = textscan(fid,'%*s %f %f %f %*[^\n]', 'commentstyle', 't', 'collectoutput', 1); 
    % Reading only the range starts and ends, as they all should correspond to
    % the current chr (going by the downloaded '.map' file name)
    status = fclose(fid);
    if status ~= 0
        disp(['The ' MappabFldr curchr ' file did not close properly!'])
    end
    Mapab.(curchr).RngStrts = MainCell{1}(:,1)+1; % The "+1" is to make the start locations 1-based
    Mapab.(curchr).RngEnds  = MainCell{1}(:,2);
end
fprintf('%s\n',' ')

%--------------- PART 1 - Building the GENOME-WIDE BACKGROUND MODEL ---------------
% Step 1: Make consecutive windows genome-wide; Then use the read locations,
% genome sequence, and mappability intervals that have been loaded into memory
% to CREATE feature TRACKS needed to build the background regression model
%
% Step 2: Build the background model and save coefficients into memory
%
% NOTE: These are specific for a window size, and they will change if the window
% size changes
%----------------------------------------------------------------------------------
disp('*** BUILDING THE BACKGROUND MODEL ***')
disp('Making windows and processing all tracks for: ')
for jis = 1:size(ListofchrsAndSzs,1)    
    curchr = ListofchrsAndSzs{jis,1};
    curchrSiz = ListofchrsAndSzs{jis,2};
    fprintf('%s ',curchr)
    
    % ***** Making consecutive chromosome-wide windows *****
    FrstRd = 1;
    LstRd = WinSz*floor(curchrSiz/WinSz);
    histEdges = FrstRd:WinSz:LstRd+1;
        
    %***** MAIN CHIP-SEQ TREATMENT *****  
    CmBndRds = [TrtReadsStruct.([curchr '_Plus']).RdStrts ; TrtReadsStruct.([curchr '_Minus']).RdEnds]; 
    RdsInWndws = histc(CmBndRds, histEdges); % Getting read counts in each window
    RdsInWndws(end) = []; % Removing the last element
    TrtStruct.(curchr) = RdsInWndws;
    % Making binary vectors for use with peaks (if needed)
    if YesforPredict
        BIN_Trt.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_Trt.(curchr)(CmBndRds) = 1;
    end
    
    %***** DNaseI *****
    CmBndRds = [DNaseIReadsStruct.([curchr '_Plus']).RdStrts ;  DNaseIReadsStruct.([curchr '_Minus']).RdEnds]; 
    RdsInWndws = histc(CmBndRds, histEdges); % Getting read counts in each window
    RdsInWndws(end) = []; % Removing the last element
    DNaseIStruct.(curchr) = RdsInWndws;
    % Making binary vectors for use with peaks (if needed)
    if YesforPredict
        BIN_DNaseI.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_DNaseI.(curchr)(CmBndRds) = 1;
    end
    
    %***** iDNA *****
    CmBndRds = [iDNAReadsStruct.([curchr '_Plus']).RdStrts  ; iDNAReadsStruct.([curchr '_Minus']).RdEnds]; 
    RdsInWndws = histc(CmBndRds, histEdges); % Getting read counts in each window
    RdsInWndws(end) = []; % Removing the last element
    iDNAStruct.(curchr) = RdsInWndws;
    % Making binary vectors for use with peaks (if needed)
    if YesforPredict
        BIN_iDNA.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_iDNA.(curchr)(CmBndRds) = 1;
    end
    
    %***** IgG *****
    CmBndRds = [IgGReadsStruct.([curchr '_Plus']).RdStrts ;  IgGReadsStruct.([curchr '_Minus']).RdEnds]; 
    RdsInWndws = histc(CmBndRds, histEdges); % Getting read counts in each window
    RdsInWndws(end) = []; % Removing the last element
    IgGStruct.(curchr) = RdsInWndws;
    % Making binary vectors for use with peaks (if needed)
    if YesforPredict
        BIN_IgG.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_IgG.(curchr)(CmBndRds) = 1;
    end
    
    %***** MAPPABILITY *****
    % What we have is positive strand mappability data (we can get -ve strand mappability by appropriate shifting)
    MapabStrts = Mapab.(curchr).RngStrts(:);
    MapabEnds  = Mapab.(curchr).RngEnds(:);
    idx = arrayfun(@(x,y) x:y, MapabStrts, MapabEnds, 'uniformoutput',0);
    % idx contains the indices of all mappable base pairs locations, i.e., every element within the mappable intervals
    idx = horzcat(idx{:});
    %
    % Counting Mappable base pairs in windows
    % +ve strand mappability
    MappabinWndwsP = histc(idx, histEdges);
    MappabinWndwsP(end) = []; % Ignoring the last (superfluous) bin 
    % This bin corresponds to only the (LstRd+1)th position, which is extra, and
    % whose purpose was just to trick the histc function to include all required
    % counts in the penultimate window
    %
    % -ve strand mappability
    % We need to get -ve strand mappability counts too! 
    MappabinWndwsM = histc(idx,max(1,( histEdges-(ReadLength-1) )));
    MappabinWndwsM(end) = []; % Ignoring the last (superfluous) bin
    MappabinWndws = (MappabinWndwsP + MappabinWndwsM)/(2*WinSz); % Computing the fraction of basepairs that are mappable! NOTE: the "2" in denominator
    MappabStruct.(curchr) = MappabinWndws;
    % Making binary vectors for use with peaks (if needed)
    % + strand mappability
    if YesforPredict
        BIN_map.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_map.(curchr)(idx) = 1;
    end
    
    %***** GC CONTENT *****
    % Finding GC locations (G & C indices)
    idx = sort([strfind(upper(Genm.(curchr)),'G') strfind(upper(Genm.(curchr)),'C')]);
    %
    % Counting GC in windows
    GCinWndws = histc(idx, histEdges);
    GCinWndws(end) = []; % Ignoring the last (superfluous) bin
    GCinWndws = GCinWndws/WinSz; % Computing the fractions
    GCStruct.(curchr) = GCinWndws;
    % Making binary vectors for use with peaks (if needed)
    if YesforPredict
        BIN_gc.(curchr) = zeros(1,curchrSiz,'uint8');
        BIN_gc.(curchr)(idx) = 1;
    end      
end
fprintf('%s\n',' ')

%-----------------------------------------------------
% Linear Regression & Saving all coefficients
%-----------------------------------------------------
disp('Doing linear regression ...')
TrtNRds         = [];
DNaseINRds      = [];
iDNANRds        = [];
IgGNRds         = [];
MappabinWndws   = [];
GCinWndws       = [];
for jis = 1:size(ListofchrsAndSzs,1)
    curchr = ListofchrsAndSzs{jis,1};
    % Joining all windows for all chrs together
    TrtNRds = [TrtNRds ; TrtStruct.(curchr)(:)];
    DNaseINRds = [DNaseINRds ; DNaseIStruct.(curchr)(:)];
    iDNANRds = [iDNANRds ; iDNAStruct.(curchr)(:)];
    IgGNRds = [IgGNRds ; IgGStruct.(curchr)(:)];
    MappabinWndws = [MappabinWndws ; MappabStruct.(curchr)(:)];
    GCinWndws = [GCinWndws ; GCStruct.(curchr)(:)];
end
% DOING REGRESSION
X = [ones(size(GCinWndws)) MappabinWndws GCinWndws DNaseINRds iDNANRds IgGNRds];
y = TrtNRds;
tm1 = X'*X;
tm2 = X'*y;
b = tm1\tm2; % Closed form: inv(X'*X)*X'*y - Using the \ operator is more efficient and accurate according to MATLAB documentation
% SAVING COEFFICIENTS
Coefs.(['Win_' num2str(WinSz)]) = b;
disp('Saving the regression coefficients for future use ...')
save(['Saved_Coefs_' justTrtname ], 'Coefs');
% in txt file
fid = fopen(['Saved_Coefs_' justTrtname '.txt'],'wt');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n', 'Constant', 'Mappab', 'GC', 'DNaseI', 'iDNA', 'IgG');
cellfun(@(x,y,z,a,b,c) fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',x,y,z,a,b,c), ...
    num2cell(b(1)),     num2cell(b(2)),     num2cell(b(3)),     num2cell(b(4)),     num2cell(b(5)),     num2cell(b(6)) );
status = fclose(fid);
if status ~= 0
    disp('One of the text files did not close properly!')
end

%---------- PART 2 - PREDICT IN PEAKS USING BACKGROUND MODEL --------------------------
%
% If the user has provided a "prediction intervals" file, then in this part, we
% use the regression coefficients computed earlier representing the background
% factors to predict the background effects inside the given intervals
%
% Then, we subtract these effects from the observed number of reads in Trt track
% to estimate the true binding portion of the signal
%
%--------------------------------------------------------------------------------------

if YesforPredict
    
    %*** Read the Prediction Intervals file (these could also be peaks) ***
    fid = fopen([DataFldr PredictIntvlsFileName]);
    pksMAINCell = textscan(fid, '%s %f %f %*[^\n]','collectoutput',1);
    % Cols that have been Read: { 'Chr' 'start' 'end' } | Ignore all columns beyond the first three
    status = fclose(fid);
    if status ~= 0
        disp('The text file did not close properly!')
    end
    % Consolidating pksMAINCell into a single cell array
    pksMAINCell = [pksMAINCell{1} num2cell(pksMAINCell{2})];
    numofpks = size(pksMAINCell,1);
    disp(['No. of prediction intervals found: ' num2str(numofpks)]);
    % Finding the unique list of chrs found in the peak (or) intervals data
    chrsFoundinpks = unique(pksMAINCell(:,1));
    [lia,~] = ismember(chrsFoundinpks, ListofchrsAndSzs(:,1));
    % Now, we will ignore all those chrs from the prediction intervals (or peaks)
    % list that are NOT in the common consensus list! This is because, in order to
    % predict background accurately within an interval in a particular chr, we need
    % to have the values from all the background tracks. This is not guaranteed
    % unless the chr is on the consensus list
    chrsToIgnore = chrsFoundinpks(~lia);
    % Ignoring all chromosomes in chrsToIgnore
    disp('Ignoring the following chromosome entries from intervals list that are NOT in the consensus chromosome list: ')
    disp('( chr   | no. of entries left out)')
    for chrcounter = 1:length(chrsToIgnore)
        cchrTI = chrsToIgnore{chrcounter};
        ignrIndx = strcmpi(pksMAINCell(:,1), cchrTI);
        pksMAINCell(ignrIndx,:) = [];
        fprintf('%7s | %d\n', cchrTI, sum(ignrIndx));
    end
    fprintf('%s\n',' ');
    numofpks = size(pksMAINCell,1);
    disp(['No. of REMAINING prediction intervals: ' num2str(numofpks)]);
    
    % *** Processing all intervals (or peaks) WITHOUT Looping ***
    %
    % NOTE: Although the peaks/intervals would be of varying sizes, the
    % background model has been built using a single window size that was given
    % by the user! So, we suggest that the user choose this window size
    % carefully, representing an average size that will capture the background
    % reasonably well in all the given prediction intervals
    %
    pkschrs   = pksMAINCell(:,1);
    peakStrts = pksMAINCell(:,2);
    peakEnds  = pksMAINCell(:,3);
    peakLens  = cell2mat(peakEnds) - cell2mat(peakStrts) + 1;
    ShftdpeakStrts = num2cell( cell2mat(peakStrts)-(ReadLength-1) );
    ShftdpeakEnds  = num2cell( cell2mat(peakEnds)-(ReadLength-1) );
    
    % Preparing tracks inside peak intervals
    TrtNRdsinPks       = cellfun(@(x,L,U) sum(BIN_Trt.(x)(L:U)), pkschrs, peakStrts, peakEnds);
    IgGNRdsinPks       = cellfun(@(x,L,U) sum(BIN_IgG.(x)(L:U)), pkschrs, peakStrts, peakEnds);
    iDNANRdsinPks      = cellfun(@(x,L,U) sum(BIN_iDNA.(x)(L:U)), pkschrs, peakStrts, peakEnds);
    DNaseINRdsinPks    = cellfun(@(x,L,U) sum(BIN_DNaseI.(x)(L:U)), pkschrs, peakStrts, peakEnds);
    GCinWndwsinPks     = cellfun(@(x,L,U) sum(BIN_gc.(x)(L:U)), pkschrs, peakStrts, peakEnds)./peakLens;
    mapinOrigPks       = cellfun(@(x,L,U) sum(BIN_map.(x)(L:U)), pkschrs, peakStrts, peakEnds);            % + strand
    mapinShftdPks      = cellfun(@(x,L,U) sum(BIN_map.(x)(L:U)), pkschrs, ShftdpeakStrts, ShftdpeakEnds);  % - strand
    MappabinWndwsinPks = (mapinOrigPks + mapinShftdPks)./(2*peakLens);
    % Using background regression coefficients to predict IN PEAKS
    X = [ones(size(GCinWndwsinPks)) MappabinWndwsinPks GCinWndwsinPks DNaseINRdsinPks iDNANRdsinPks IgGNRdsinPks];
    fitline = X*b(:);  % the variable 'b' did not change, and was not erased from workspace! So, we can directly use it here!
   
    % Computing PURIFIED ChIP-seq binding signal based on our prediction
    TrtNRdsinPks = TrtNRdsinPks(:);
    fitline = fitline(:);
    PureSignal = TrtNRdsinPks - fitline;
    disp('Creating the OUTPUT file with columns: | "Chr Name" | "0-based Start" | "1-based End" | "Observed NRds" | "Background NRds" | "Purified Signal" | ')
    fid = fopen([otptnam '_ModelWinSz_' num2str(WinSz) '_OtPt.bed'],'wt');
    cellfun(@(x,y,z,a,b,c) fprintf(fid,'%s\t%d\t%d\t%f\t%f\t%f\n',x,y,z,a,b,c), pkschrs, peakStrts, peakEnds, num2cell(TrtNRdsinPks), num2cell(fitline), num2cell(PureSignal) );
    % Cols Written: { | "Chr Name" | "0-based Start" | "1-based End" | "Observed NRds" | "Background NRds" | "Purified Signal" | }
    status = fclose(fid);
    if status ~= 0
        disp('OUTPUT text files did not close properly!')
    end    
    disp('------- DONE Predicting!  Check OUTPUT file! -------')
    disp(' ')
else
    disp('No "Prediction Intervals" file name given!  So, terminating after building the background model and saving it!')
end
toc
exit

% ****************************************************************
%                BED File Processing Function
% **************************************************************** 
% 
% Read a BED file - separate into positive- and negative-strand reads,
% de-duplicate the reads, sort, and then return a structure containing the reads
% for each chr, and the list of chrs for which reads were found in the dataset
function [ReadsStruct, chrList] = ReadBED(DataFldr, fullbedname)
fid = fopen([DataFldr fullbedname]);
% tTemp = fgetl(fid);  % Uncomment to reject the first line if this line has header info
MainCell = textscan(fid, '%s %f %f %*s %*s %s %*[^\n]', 'commentstyle', 't');  % Read only 4 columns and ignore the rest
status = fclose(fid);
if status ~= 0
    disp(['The BED file ' DataFldr fullbedname ' did not close properly!'])
end

% Just getting the total number of reads read
szMain = size(MainCell{2});
TotNoRds = szMain(1);
disp(['Total number of reads found in the BED file : ' num2str(TotNoRds)]);

% Obtaining the list of chromosomes present using "unique"
[chrList,~,chrIdx] = unique(MainCell{1});
disp('Separating reads, removing duplicates, and sorting ...')
AccumRds = struct;
numUnqRds = [];
for jis = 1:length(chrList)
    % Plus
    Indx = (chrIdx == jis) & strcmpi(MainCell{4}, '+');
    [uPlusRdStrts,uIndx,~] = unique(MainCell{2}(Indx));
    noNuPlusRdEnds = MainCell{3}(Indx);
    uPlusRdEnds = noNuPlusRdEnds(uIndx);
    AccumRds.([chrList{jis} '_Plus']).RdStrts = uPlusRdStrts(:)+1; % The "+1" is to make the starts 1-based, as they are 0-based in the BED file
    AccumRds.([chrList{jis} '_Plus']).RdEnds  = uPlusRdEnds(:);
    %     Use the two lines below if the first colmn in the BED file has only nos. like "1"
    %     instead of "chr1"
    %     AccumRds.(['chr' chrLst{jis} '_Plus']).RdStrts = uPlusRdStrts;
    %     AccumRds.(['chr' chrLst{jis} '_Plus']).RdEnds  = uPlusRdEnds;
    
    % Minus
    Indx = (chrIdx == jis) & strcmpi(MainCell{4}, '-');
    [uMinusRdStrts,uIndx,~] = unique(MainCell{2}(Indx));
    noNuMinusRdEnds = MainCell{3}(Indx);
    uMinusRdEnds = noNuMinusRdEnds(uIndx);
    AccumRds.([chrList{jis} '_Minus']).RdStrts = uMinusRdStrts(:)+1; % The "+1" is to make the starts 1-based, as they are 0-based in the BED file
    AccumRds.([chrList{jis} '_Minus']).RdEnds  = uMinusRdEnds(:);
    %     Use the two lines below if the first colmn in BED file has only nos. like "1"
    %     instead of "chr1"
    %     AccumRds.(['chr' chrLst{jis} '_Minus']).RdStrts = uMinusRdStrts;
    %     AccumRds.(['chr' chrLst{jis} '_Minus']).RdEnds  = uMinusRdEnds;
    numUnqRds = [numUnqRds length(uPlusRdStrts)+length(uMinusRdStrts)];
end
ReadsStruct = orderfields(AccumRds);
disp(['Total number of unique reads found : ' num2str(sum(numUnqRds))]);
disp('... DONE ')
disp(' ')



