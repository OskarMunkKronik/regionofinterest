function[mzroi,MSroi,PeaksMat,runTime,sizeMZRoI,Rt,Dt]=ROIpeaks_ACN(FileName,Options)
%ROIpeaks_ACN finds regions of interest in the chromatographic and mass
%spectral dimension. It performs the binning of the mass spectral
%measurements by organizing the cell array with all mass spectral
%measurements into a matrix. Then it sorts the rows according to m/z value.
%Afterwich it performs binning, by binning the mass spectral measurements by taking the first value + mzerror
%all m/z which falls in this range are binned together. At the end the
%intensity weighted mean (wmean = 1) or median (wmean = 0)
%m/z is calculated for the mass bin.

%Input:
%FileName: The name of the CDF file to be processed.
%Options:
%RtInt: RtInt(1):   Is the beginning of the retention time interval in scan number wanting to be
%                   processed. RtInt(2): End of retention time interval in scan numbers
%mzerror:           is the user specified allowed mass error for the ROI
%                   procedure to use. Typically set to a multiple of the mass accuracy of
%                   the mass spectrometer
%ppm:               Defines whether the mzerror is defined in ppm (ppm==1) or Da (ppm=0)
%thresh:            Is the intensity threshold below which data points are
%                   exluded.
%minroi:            filters away ROIs with less consecutive scans than specified in
%minroi
%GapAllowed:        Allows for a number of missing scan points inside the ROI
%IMS:               true if you have IMS data. The minroi filter will be done on the
%                   mobilograms
%wmean:             i                   ntensity weighted mean (wmean = 1) or median (wmean = 0) m/z is calculated for the mass bin.
%CollapseRoIs:      If true: than m/z traces are created with ROIs
%                   of the same m/z, if false ROIs confined in m/z and Rt and Dt are
%                   stored in MSroi.

%Output:
%mzroi:         Median m/z of the m/z's in each region of interest.
%MSroi:         Sparse Data matrix with summed intensities. Dimensions (m/z x Scan number)
%Rawdata:       Sorted raw data with m/z, Intensity, scan number and m/z groups.
%runTime:       processing time after loading of the data.
%sizeMZRoI:     Dimensions of the MSroi matrix
%Rt:            Chromatographic retention time
%Dt:            Drift time

%Authors: Oskar Munk Kronik/Giorgio Tomasi
%email:   oskarmunkkronik@gmail.com


%% Initialize

if (exist(FileName,'file') ~= 2), error('Invalid file path'); end
if (nargin < 2 || isempty(Options)), Options = struct(); end
sInfo   = ncinfo(FileName,'point_count');
nScans  = sInfo.Size;
Options = checkOptions(Options,nScans);

%%
% tic
%Organize data

%Load cdf and create PeaksMat
pts                = ncread(FileName,'point_count');
varNam             = ncinfo(FileName);
if ismember({varNam.Variables.Name},'scan_index')
    scanIndex          = ncread(FileName,'scan_index');
else
    scanIndex = [0;cumsum(pts(1:end-1))];
end
Rt                 = single(ncread(FileName,'scan_acquisition_time'));
scanIndex(end + 1) = scanIndex(end) + pts(end);
%Read drift times
if Options.IMS
    Dt = single(ncread(FileName,'drift_acquistisiton_time'));
else
    Dt =[];
end
if (isequal(Options.RtInt,[1, nScans]))
    readNpoints = scanIndex(end);
    readPars    = {};
else
    readFirstPosition    = scanIndex(Options.RtInt(1)) + 1;
    readNPoints          = scanIndex(Options.RtInt(2)) - scanIndex(Options.RtInt(1));
    pts                  = pts(Options.RtInt(1):Options.RtInt(2));
    Rt                   = Rt(Options.RtInt(1):Options.RtInt(2));
    readPars             = {readFirstPosition,readNPoints};
    if Options.IMS
        Dt = Dt(Options.RtInt(1):Options.RtInt(2));
    end
end
PeaksMat      = zeros(readNpoints,4,'double');
PeaksMat(:,1) = ncread(FileName,'mass_values',readPars{:});
PeaksMat(:,2) = ncread(FileName,'intensity_values',readPars{:});
if Options.IMS
    scan_number = ncread(FileName,'scan_acquisition_number');
    PeaksMat(:,3) = repelem(scan_number,pts)';
%     [~,~,uRt]=unique(Rt);
%      PeaksMat(:,5) = uRt(PeaksMat(:,3));
      
else
    PeaksMat(:,3) = repelem(Options.RtInt(1):Options.RtInt(2),pts)';
end

if Options.RefMass
    % Load RefMass
    Rt_RefMass                 =    single(ncread([Options.RefMass_dir,'\',FileName],'scan_acquisition_time'));
    mass_values_RefMass        =    ncread([Options.RefMass_dir,'\',FileName],'mass_values');
    intensity_values_RefMass   =    ncread([Options.RefMass_dir,'\',FileName],'intensity_values');

    pts_RefMass                        =    ncread([Options.RefMass_dir,'\',FileName],'point_count');
    %     scan_index_RefMass        =    = [0;cumsum(pts(1:end-1))];
    scan_index_RefMass        =     repelem( [1:length(pts_RefMass)], ncread([Options.RefMass_dir,'\',FileName],'point_count'));
    Ind_RefMass                =    intensity_values_RefMass > Options.RefMass_thresh & abs(mass_values_RefMass-Options.RefMass_mz) < Options.RefMass_mzDev;
    %     Rt_RefMass = Rt_RefMass(Ind_RefMass);


    Rt_Sample = Rt;%(PeaksMat(:,3));
    RefMass_closest = zeros(length(Rt_Sample),1);
    nSam = 1;
    lenRefRt = length(Rt_RefMass);
    RefMass_mz_diff = nan(lenRefRt,1);
    RefMass_closest(1) = 1;
    nRef = 1;
    start = find(Rt_Sample>Rt_RefMass(nRef),1,'first');
    for nSam = start : length(Rt_Sample)

        Rt_diff =  Rt_Sample(nSam) - Rt_RefMass(nRef);
        if nRef < lenRefRt & Rt_diff > Rt_RefMass(nRef+1) - Rt_Sample(nSam)
            RefMass_closest(nSam) = 1;
            sub_mz  = mass_values_RefMass(scan_index_RefMass == nRef);
            sub_int = intensity_values_RefMass(scan_index_RefMass == nRef);
            %            [max_intensity,max_ind] = max(sub_int(Options.RefMass_mz-sub_mz(sub_int > Options.RefMass_thresh & abs(sub_mz-Options.RefMass_mz) < Options.RefMass_mzDev)));
            if ~isempty(Options.RefMass_mz-sub_mz(sub_int > Options.RefMass_thresh & abs(sub_mz-Options.RefMass_mz) < Options.RefMass_mzDev))
                RefMass_mz_diff(nRef) =  Options.RefMass_mz-sub_mz(sub_int > Options.RefMass_thresh & abs(sub_mz-Options.RefMass_mz) < Options.RefMass_mzDev);
                %            RefMass_mz_diff(nRef) = sub_mz(Ind);
            end
            nRef = nRef + 1;

        end
    end
    RefMass_mz_diff(isnan(RefMass_mz_diff)) = median(RefMass_mz_diff,'omitnan');
    RefMass_closest =  repelem( [cumsum(RefMass_closest)], pts);
    RefMass_mz_diff =  RefMass_mz_diff(RefMass_closest);
    LockMass_Corr = RefMass_mz_diff./Options.RefMass_mz;
    PeaksMat(:,1) = (1+LockMass_Corr).*PeaksMat(:,1);
end
%         nSam
% Read driftTimes

aa = tic;
% AboveThresh=SampleStruct.intensity_values.data>Options.thresh;
% PeaksMat=PeaksMat(SampleStruct.intensity_values.data>Options.thresh,:);
flag        = ':';
if (Options.prefilter)

    PeaksMat = PeaksMat(PeaksMat(:,2) >= Options.thresh,:);

elseif (Options.thresh > 0) % Necessary? && Options.GapAllowed > 0
    flag = false(readNpoints,1);
end


%Sorts the data so that the m/z's are sorted in ascending order:
PeaksMat = sortrows(PeaksMat,1,'ascend');
mzVec    = PeaksMat(:,1);
%Performs the grouping of the individual m/z's:
mzDev = Options.mzerror;
if (Options.ppm), mzDev = mzDev *10^-6; end

    function wU = calculateULppm(m,mzDev), wU = m * (1 + mzDev); end % NOTE: the correct value is m * (1 + mzDev)/(1 - 0.5 * mzDev)
    function wU = calculateULDa(m,mzDev), wU = m + mzDev; end

if (Options.ppm), updateUL = @calculateULppm;
else, updateUL = @calculateULDa;
end

%Performs the binning
Repeat   = true;
RepeatNo = 1;
nPts     = size(PeaksMat,1);

while Repeat

    if (~Options.prefilter && Options.thresh > 0)

        if (RepeatNo == 1)

            kStart = 1;
            while (PeaksMat(kStart,2) <= Options.thresh && kStart <= nPts)
                kStart = kStart + 1;
            end
            flag(kStart)       = true;
            PeaksMat(kStart,4) = 1;
            wU                 = updateUL(PeaksMat(kStart,1),0.5 * mzDev);

        end

        for (k = kStart + 1:nPts)

            flag(k) = PeaksMat(k,2) >= Options.thresh;
            if (flag(k) && PeaksMat(k,1) > wU)
                wU            = updateUL(PeaksMat(k,1),mzDev);
                PeaksMat(k,4) = 1;
            end

        end

    else

        PeaksMat(1,4) = 1;
        wU            = updateUL(PeaksMat(1,1),0.5 * mzDev);
        for (k = 2:nPts)
            if (PeaksMat(k,1) > wU)
                wU            = updateUL(PeaksMat(k,1),mzDev);
                PeaksMat(k,4) = 1;
            end
        end

    end
    PeaksMat(:,4) = cumsum(PeaksMat(:,4));
    if Options.wmean
        w     = accumarray(PeaksMat(flag,4),PeaksMat(flag,2));
        w     = PeaksMat(flag,2)./w(PeaksMat(flag,4));
        mzroi = accumarray(PeaksMat(flag,4),w.*PeaksMat(flag,1));
    else
        mzroi = accumarray(PeaksMat(flag,4),PeaksMat(flag,1),[],@median);
    end
    if Options.ppm
        Repeat = or(sum(diff(mzroi(2:end))./mzroi(2:end-1)<mzDev)>0 , (diff(mzroi(1:2))./mzroi(1)<mzDev*0.5)>0);
    else
        % Repeat=sum(diff(mzroi)<mzDev)>0;
        Repeat = or(sum(diff(mzroi(2:end))<mzDev)>0 , (diff(mzroi(1:2))<mzDev*0.5)>0);

    end
    PeaksMat(:,1)=mzroi(PeaksMat(:,4));
    if Repeat

        RepeatNo = RepeatNo + 1;
        PeaksMat(:,4) = 0;

    end

end


% end
%Sorts the data so that the rows of the final matrix are in ascending order
%as the subsequent proces should be faster then.
if (~Options.prefilter)
    [PeaksMat,order] = sortrows(PeaksMat,[4,3]);
    flag             = flag(order);
else
    [PeaksMat,order] = sortrows(PeaksMat(flag,:),[4,3]);
    mzVec    = mzVec(order);
    mzVec    = mzVec(flag);

end
% Raw data prior to the chromatographic filter
% Chromatographic filter: excludes measurements which has <
% chromatographically consecutive points than specified by minroi
% PeaksMat1 = PeaksMat;
PeaksMat(:,1)=round(mzroi(PeaksMat(:,4)),Options.nDecimals);
if Options.minroi > 1
    [PeaksMat,gaps,~,flagMZ] = minroiFilter2(PeaksMat,flag,Options);
    mzVec = repelem(mzVec(flagMZ(1:end-1)),gaps+1); %Removes exluded m/z due to minroiFilter2 and replicates m/z of gaps
    if (Options.GapAllowed && any(gaps)), PeaksMat = fillGaps(PeaksMat,gaps); end
    %Calculates new number of data points in PeaksMat after minroiFilter2
    nPts = height(PeaksMat);
end

%Collapse regions with the same m/z
if Options.CollapseRoIs
    PeaksMat(:,4) = 0;
    PeaksMat(1,4) = 1;

    currentMZ = PeaksMat(1,1);

    for (k = 2:nPts)
        if (PeaksMat(k,1)-currentMZ) > 10^-Options.nDecimals
            currentMZ      = PeaksMat(k,1);
            PeaksMat(k,4)  = 1;
        end
    end
    PeaksMat(:,4) = cumsum(PeaksMat(:,4));
end

% if Options.IMS
%     if Options.minroi > 0
% %         [~,~,uRt]=unique(Rt);
%             
% %         PeaksMat(:,5) = uRt(PeaksMat(:,3));
%         PeaksMat(:,[1:2,5,4,3]) = PeaksMat;
%         [~,~,uRt]=unique(PeaksMat(:,3));
%         if (~Options.prefilter)
%             [PeaksMat,order] = sortrows(PeaksMat,[4,3]);
%             flag             = flag(order);
%         else
%             [PeaksMat,order] = sortrows(PeaksMat(flag,:),[4,3]);
%             mzVec    = mzVec(order);
%             mzVec    = mzVec(flag);
% 
%         end
%         
%         if Options.NumTrace >1 
%             ind = 1:2:max(PeaksMat(:,3));
%             PeaksMat1 = PeaksMat(ismember(uRt,ind),:);
%             mzVec1    = mzVec(ismember(uRt,ind));
% %             flag1 = flag(find(ismember(uRt,ind)));
%             PeaksMat2 = PeaksMat(~ismember(uRt,ind),:);
%             mzVec2    = mzVec(~ismember(uRt,ind));
% %             flag2 = flag(~ismember(uRt,ind));
%         [PeaksMat1,gaps,~,flagMZ] = minroiFilter2(PeaksMat1,flag,struct('GapAllowed',Options.GapAllowed,'minroi',Options.minroiChrom,'prefilter',Options.prefilter,'IMS_chrom_interpolation',true));
%         mzVec1 = repelem(mzVec1(flagMZ(1:end-1)),gaps+1); %Removes exluded m/z due to minroiFilter2 and replicates m/z of gaps
%         if (Options.GapAllowed && any(gaps)), PeaksMat1 = fillGaps(PeaksMat1,gaps); end
%           
%         
%         [PeaksMat2,gaps,~,flagMZ] = minroiFilter2(PeaksMat2,flag,struct('GapAllowed',Options.GapAllowed,'minroi',Options.minroiChrom,'prefilter',Options.prefilter,'IMS_chrom_interpolation',true));
%         mzVec2 = repelem(mzVec2(flagMZ(1:end-1)),gaps+1); %Removes exluded m/z due to minroiFilter2 and replicates m/z of gaps
%         if (Options.GapAllowed && any(gaps)), PeaksMat2 = fillGaps(PeaksMat2,gaps); end
%         PeaksMat = cat(1,PeaksMat1,PeaksMat2);
%         mzVec = cat(1,mzVec1,mzVec2);
%         
%         end 
%         PeaksMat(:,[1:2,5,4,3]) = PeaksMat;
%         [~,~,PeaksMat(:,4)] = unique( PeaksMat(:,4));
%         
%         %Calculates new number of data points in PeaksMat after minroiFilter2
%         nPts = height(PeaksMat);
%     end
% end
% Gives the remaining groups after the chromatographic filter new group
% number so the group numbers go from 1-number of groups remaining.
sizeMZRoI = [PeaksMat(end,4),Options.RtInt(2)]; % sorted by construction

% Sums the intensities belonging to each individual measurement which occurs
% at the same scan point and m/z group.
MSroi = accumarray(PeaksMat(:,[4,3]),double(PeaksMat(:,2)),sizeMZRoI,[],Options.fillIn,~isnan(Options.fillIn));
if (any(MSroi < 0,'all')), warning('Some RoI contain negative values'); end
% mzroi = mzroi(1:sizeMZRoI(1));
% Calculate median m/z's of the individual m/z values found by the ROI procedure.
PeaksMat(:,1) = mzVec;
if Options.wmean == 1
    w=accumarray(PeaksMat(:,4),PeaksMat(:,2));
    w=PeaksMat(:,2)./w(PeaksMat(:,4));
    mzroi = accumarray(PeaksMat(:,4),w.*PeaksMat(:,1));
else
    mzroi = accumarray(PeaksMat(:,4),PeaksMat(:,1),[],@median);
end
runTime = toc(aa);
end

function PeaksMat = fillGaps(PeaksMat,gaps)
[gapsPosition,~,nPoints] = find(gaps);
if (isempty(gapsPosition)), return; end
nPoints      = double(nPoints);
gapsFirst    = PeaksMat(gapsPosition - 1,:);
gapsLast     = PeaksMat(gapsPosition,:);
maxGap       = max(nPoints);
interpCoeff  = zeros(maxGap);
PeaksMat     = repelem(PeaksMat,gaps + 1,1);
gapsPosition = gapsPosition + [0;cumsum(nPoints(1:end - 1))];
for (i = 1:maxGap), interpCoeff(1:i,i) = (1:i)/(i + 1); end
for (i = 1:length(gapsPosition))
    coeff                  = interpCoeff(1:nPoints(i),nPoints(i)); % For legibility
    rowIndex               = gapsPosition(i) + (0:nPoints(i) - 1); % For legibility
    PeaksMat(rowIndex,1:2) = (1 - coeff) * gapsFirst(i,1:2) + coeff * gapsLast(i,1:2);
    PeaksMat(rowIndex,3)   = (gapsFirst(i,3) + 1):(gapsLast(i,3) - 1);
end

end


function  [PeaksMat,gaps,flagRoI,flag] = minroiFilter2(PeaksMat,flag,Options)
%Filters out MS measurements,which is not doesn't have
switch find([8 16 32 64] > log2(PeaksMat(end,4)),1,'first') + 1 % Reduces memory footprint
    case 2, precision = 'uint16';
    case 3, precision = 'uint32';
    otherwise, precision = 'uint64';
end
% bb       = tic();
flagRoI  = zeros(size(PeaksMat,1),1,precision);
gaps     = zeros(size(PeaksMat,1),1,'uint8'); % Gaps < 255 scans, which seems reasonable
maxGap   = Options.GapAllowed + 1;
lastScan = size(PeaksMat,1) + 1;
s                   = 1;
currentRegion       = PeaksMat(1,4);
regionCounter       = 1;
regionLengthScans   = zeros(size(PeaksMat,1),1,precision);
regionLengthPoints  = zeros(size(PeaksMat,1),1,precision);
PeaksMat(end + 1,4) = PeaksMat(end,4) + 1;
lastRT              = PeaksMat(1,3) - 1;
lastRTflag          = lastRT;
for k = 1:lastScan

    dRT = PeaksMat(k,3) - lastRT;
    if (PeaksMat(k,4) ~= currentRegion || dRT > maxGap) % if (dRT <= maxGap), if (flag(k)), else, end

        regionLengthScans(regionCounter)  = (PeaksMat(k-1,3) - PeaksMat(s,3) + 1); % This is be used for additional filtering based on the region's length. k == 1 will never get in here because of initialisation
        regionLengthPoints(regionCounter) = k-s;
        flagRoI(s)                        = 1;
        s                                 = k;
        regionCounter                     = regionCounter + 1;
        currentRegion                     = PeaksMat(s,4);
        if (~Options.prefilter), flag(lastRTflag + 1:k) = false; end

    else, gaps(k) = dRT - 1;
    end
    lastRT = PeaksMat(k,3);
    if (~Options.prefilter && flag(k)), lastRTflag = k; end

end
regionLengthPoints(regionCounter) = 1;
regionLengthScans(regionCounter)  = 1;
tmp                               = regionLengthPoints(1:regionCounter) >= feval(precision,Options.minroi) & regionLengthScans(1:regionCounter) >= feval(precision,Options.minroi); %#ok<FVAL> % filters based on number of points in a region and the length of the region
flag                              = repelem(tmp,regionLengthPoints(1:regionCounter));
PeaksMat                          = PeaksMat(flag,:);
gaps                              = gaps(flag);
flagRoI                           = flagRoI(flag);
% if ~isfield(Options,"IMS_chrom_interpolation")
    PeaksMat(:,4)                     = cumsum(flagRoI);
% else

% end
end

function outOptions = checkOptions(inOptions,len)
outOptions = struct('minroi',10,...
    'minroiChrom',12,...
    'mzerror',10,...
    'ppm',true,...
    'thresh',0,...
    'wmean',true,...
    'RtInt',[1,len],...
    'GapAllowed',0,...
    'verbose',true,...
    'prefilter',true,...
    'fillIn',0, ...
    'CollapseRoIs',true, ...
    'nDecimals',4, ...
    'IMS',false,...
    'RefMass',false,...
    'RefMass_mz',    [],...
    'RefMass_mzDev', 0.1,...
    'RefMass_thresh', 10^3,...
    'RefMass_dir','' ...
    );

if (~isstruct(inOptions)), error('Invalid options'); end
optionNames = fieldnames(outOptions);
if (isfield(inOptions,'verbose')), outOptions.verbose = inOptions.verbose; end
optionNames(strcmp(optionNames,'verbose')) = [];
for (i = 1:length(optionNames))

    if (isfield(inOptions,optionNames{i}))

        outOptions.(optionNames{i}) = inOptions.(optionNames{i});
        if (outOptions.verbose)
            patt = '\n [\b%s]\b: %g';
            if (length(outOptions.(optionNames{i})) > 1), patt = strcat(patt,' - %g'); end
            fprintf(1,patt,optionNames{i},outOptions.(optionNames{i}));
        end

    else

        if (outOptions.verbose)

            fprintf(1,'\n [\b%s]\b was not specified.',optionNames{i})
            switch optionNames{i}
                case 'ppm',       fprintf(1,' Use TRUE (default) for Dalton and true for ppm')
                case 'thresh',    fprintf(1,' It is defined as an absolute intensity')
                case 'wmean',     fprintf(1,' If 1/TRUE -> use a weighted mean to calculate the m/z value, 0/FALSE = Non-weighted median')
                case 'RtInt',     fprintf(1,' The whole file is processed')
                case 'GapAllowed',fprintf(1,' n allows for gaps of at most n scans in RoI''s (0 - default)')
                case 'prefilter', fprintf(1,' If TRUE (default), the values below in absolute intensity are filtered before the ROI definition.')
                case 'fillIn',    fprintf(1,' 0 by default')
                case 'CollapseRoIs',    fprintf(1,' TRUE by default')
                case 'nDecimals',    fprintf(1,' 4 by default')
                case 'IMS',        fprintf(1,'  false by default')
                case 'RefMass',        fprintf(1,'  false by default')
                case 'RefMass_mz',        fprintf(1,'  empty by default')
                case 'RefMass_mzDev',        fprintf(1,'  +/- 0.1 Da by default')
                case 'RefMass_thresh',        fprintf(1,'  10^3 by default')
                case 'RefMass_dir',        fprintf(1,'Empty by default')
%                 case 'NumTrace',        fprintf(1,'Empty by default')
                otherwise, fprintf(1,' Default is used')

            end

        end

    end

end
if (outOptions.verbose), fprintf(1,'\n\n'); end

end

