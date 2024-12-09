function [mzroi_aug,MSroi_aug,mzroiAll]=MSroiaug_ACN2(mzroi,MSroi,Options)%,mzerror,ppm,wmean,RtInt,OneDim)
%Input:
%mzroi: a cell array (number of samples x 1) containing the mzroi's of each sample (number of
%m/z's x 1)

%MSroi: a cell array (number of samples x 1) containing the intensities of each sample for all m/z's and Scan points(number of
%m/z's x number of scans)

%mzerror: Sets the allowed difference between two subsequent m/z values in
%either ppm or Da. The mz values are sorted in ascending order

%ppm: defines whether the mzerror is defined in ppm (ppm==1) or Da (ppm=0)


%Note:
%If MSe data is used the cell array should be (number of samples x 2)
%first coloumn is the low energy trace, second coloumn is the high
%energy trace

%Output:
%mzroi_aug: a vector where with th mzroi's for all the samples have
%been binned and the mean has been calculated
%MSroi_aug: a three-way array (number of mzroi_aug's x number of scans x number of samples) containing the summed intensities of the
%augmented data matrices for the individual samples.

%Author: Oskar Munk Kronik
%email: oskarmunkkronik@gmail.com
%% Initialize
if ~isfield(Options,'minroi')
    fprintf(1,'minroi was not specified, please write this into the Options struct\n')
end    %Options.minroi=4;
if ~isfield(Options,'mzerror')
    fprintf(1,'mzerror was not specified, please write this into the Options struct\n')
end
if ~isfield(Options,'ppm')
    fprintf(1,'whether to use ppm or Da was not specified, please write this into the Options struct (1 = ppm, 0 = Dalton)\n')
end
if ~isfield(Options,'thresh')
    fprintf(1,'The noise threshold (thresh) was not specified, please write this into the Options struct. It is defined as an absolute intensity\n')
end
if ~isfield(Options,'wmean')
    fprintf(1,'wmean was not specified, please write this into the Options (1 = using a weighted mean to calculate the m/z value (recommended), 0 = Non weighted median)  struct\n')
end
if ~isfield(Options,'RtInt')
    fprintf(1,'RtInt was not specified, please write this into the Options struct. Otherwise the entire CDF will be processed\n')
    
end
% if ~isfield(Options,'OneDim')
%     fprintf(1,'It is not specified if the data is one dimensional (OneDim = 1) or two dimensional (OneDim ~= 1), please write this into the Options struct\n')
% end
% if ~isfield(Options,'Unfolded2DTensor')
%     fprintf(1,'It is not specified if the output should be an unfolded 2D tensor (NO = 0, Yes = 1 )\n')
%     Options.Unfolded2DTensor=1;
% end
%%
close all force
% Waitbar's Message Field
Msg = 'Augmenting  Progress...!';
NbrePts = size(MSroi,1)*size(MSroi,2) ;
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);

mzroiAll=[];
OneNorm=[];
%dimensions([m/z, Sample number, order of m/z in sample, new pos in
%augmented matrix]

trace=[];
%dimensions([number of all m/z's for all samples in both high and low energy trace, number of scans]
%Organize the data into matrix

%Find the sample with the minimum number of scan points. Then the scan
%points from 1- eScanMin is used for the other samples.
Options.sz=nan(size(mzroi,1),size(mzroi,2));
for T=1:size(mzroi,2)
    for k=1:size(mzroi,1)
        if ~isempty(MSroi{k,T})
            Options.sz(k,T)=size(MSroi{k,T},2);
        end
    end
end
% if or(Options.OneDim==1 , Options.Unfolded2DTensor==1)
    Options.eScanMin=min(Options.sz,[],'all');
% end
% Options.RtInt = length(Options.RtInt(1):Options.RtInt(2));
%Concatonates the mzroi's of all samples and high and low energy trace.
%THIS could be optimized and parallelized.
% for T=1:size(mzroi,2)

if ~isfield(Options,'RtInt')
    for T=1:size(mzroi,2)
        for k=1:size(mzroi,1)
            Options.RtInt(k,T,:)=[1,Options.sz(k,T)];
        end
    end
end

for T=1:size(mzroi,2)
    for k=1:size(mzroi,1)
        if isempty(mzroi{k,T})
        else
            
            %
            mzInclude=sum(MSroi{k,T}(:,Options.RtInt(k,T,1):Options.RtInt(k,T,2))>0,2)>0;
            MSroi{k,T}=MSroi{k,T}(mzInclude,:);
            %                OneNorm=sum(MSroi{k,T},2);
            OneNorm=cat(1,OneNorm,sum(MSroi{k,T},2));
            mzroiAll=cat(1,mzroiAll,cat(2,mzroi{k,T}(mzInclude),repelem(k,length(mzroi{k,T}(mzInclude)))',(1:length(mzroi{k,T}(mzInclude)))'));
            trace=cat(1,trace,repelem(T,length(mzroi{k,T}(mzInclude)))');
        end
    end
end
%Sort rows
mzroiAll=cat(2,mzroiAll,zeros(size(mzroiAll,1),1));
mzroiAll=cat(2,mzroiAll,trace);
if issparse(mzroiAll)
    mzroiAll=full(mzroiAll);
end
[mzroiAll,Ind]=sortrows(mzroiAll,1);
OneNorm=OneNorm(Ind);

%Performs the binning
Repeat=true;
while Repeat
    mz=1;k=1;
    while k<size(mzroiAll,1)

        if k == 1
        if Options.ppm==1
            wU=mzroiAll(k,1)+mzroiAll(k,1)*Options.mzerror*0.5*10^-6;
        else
            wU=mzroiAll(k,1)+Options.mzerror*0.5;
        end
        else 
             if Options.ppm==1
            wU=mzroiAll(k,1)+mzroiAll(k,1)*Options.mzerror*10^-6;
        else
            wU=mzroiAll(k,1)+Options.mzerror;
             end
        end 



        while mzroiAll(k,1)<wU && k<size(mzroiAll,1)
            mzroiAll(k,4)=mz;
            k=k+1;
        end
        mz=mz+1;
    end
    
    if mzroiAll(k,1)<wU
        mzroiAll(k,4)=mz-1;
        
    else
        mzroiAll(k,4)=mz;
    end
    
    %Calculate mzroi_aug
    if Options.wmean==1
        w=accumarray(mzroiAll(:,4),OneNorm);
        w=OneNorm./w(mzroiAll(:,4));
        mzroi_aug = accumarray(mzroiAll(:,4),w.*mzroiAll(:,1));
    else
        mzroi_aug=accumarray(mzroiAll(:,4),mzroiAll(:,1),[],@median);
    end
    %Calculates whether to repeat
    if Options.ppm
        Repeat = or(sum(diff(mzroi_aug(2:end))./mzroi_aug(2:end-1)<Options.mzerror*10^-6)>0 , (diff(mzroi_aug(1:2))./mzroi_aug(1)<Options.mzerror*0.5*10^-6)>0);
    else
         Repeat = or(sum(diff(mzroi_aug(2:end))<Options.mzerror)>0 , (diff(mzroi_aug(1:2))<Options.mzerror*0.5)>0);

    end
    if Repeat
        mzroiAll(:,1)=mzroi_aug(mzroiAll(:,4));
        
    end
end




%Organizes the data into a three-way array. If the data was acquired in
%MSe mode, the high and low energy trace are summed. THIS could be
%optimized with accumarray


% if Options.OneDim==1
    Options.szMSroi=[max(mzroiAll(:,4)),Options.eScanMin];
    [MSroi_aug]=SumMSlevel(MSroi,mzroi,mzroiAll,Options);
% else
%     Options.sz=  [max(mzroiAll(:,4)),sum(Options.sz)];
%     [MSroi_aug]=SumMZgroups4Pulsed2D(MSroi,mzroi,mzroiAll,Options);
%     if Options.Unfolded2DTensor==1
%         [MSroi_aug]=Unfolded2DTensor(MSroi,mzroi,mzroiAll,Options);
%     end
% 
% end

%Calculates the mean m/z values of the individual m/z values from the different samples:
if Options.wmean==1
    w=accumarray(mzroiAll(:,4),OneNorm);
    w=OneNorm./w(mzroiAll(:,4));
    mzroi_aug = accumarray(mzroiAll(:,4),w.*mzroiAll(:,1));
else
    mzroi_aug=accumarray(mzroiAll(:,4),mzroiAll(:,1),[],@median);
end
mzroiAll=cat(2,mzroiAll,trace);
end

%%  Sum m/z across groups and Mass level 1 and 2:
function [MSroi_aug]=SumMSlevel(MSroi,mzroi,mzroiAll,Options)
%Pre-allocate space for Augmented and summed intensities
MSroi_aug=cell(max(mzroiAll(:,2)),size(mzroi,2));
Options.szMSroi=repmat(Options.szMSroi,[max(mzroiAll(:,2)),1]);
% MSroi_aug=zeros(max(mzroiAll(:,4)),Options.eScanMin,max(mzroiAll(:,2)),size(MSroi,2));
for T=1:size(mzroi,2)
    MSroiT=MSroi(:,T);
    MSroi_augT=MSroi_aug(:,T);
    % RtInt=squeeze(Options.RtInt(k,T,:))';

    parfor k=1:max(mzroiAll(:,2))
        %         hWaitbarMsgQueue.send(0);
        % MSroi_augT{k}=SumXAug(mzroiAll,MSroiT{k}(:,[RtInt(1):RtInt(2)]),k,T,Options);
       
                MSroi_augT{k}=SumXAug(mzroiAll,MSroiT{k}(:,Options.RtInt(k,T,1):min(Options.RtInt(k,:,2))),k,T,Options);
    end
    MSroi_aug(:,T)=MSroi_augT;
end
if Options.SumMSLevels 
if size(mzroi,2)>1
    for k=1:max(mzroiAll(:,2))
        minScan = min(size(MSroi_aug{k,1},2),size(MSroi_aug{k,2},2));
        MSroi_aug{k,1}=MSroi_aug{k,1}(:,1:minScan)+MSroi_aug{k,2}(:,1:minScan);
        %         MSroi_aug(:,:,k,T)=SumXAug(mzroiAll,MSroiT{k}(:,Options.RtInt(k,T,1):Options.RtInt(k,T,2)),k,T,Options);
    end
    MSroi_aug(:,2)=[];
end
end
end 
% %%  Unfolded 2D tensor:
% function [MSroi_aug]=Unfolded2DTensor(MSroi,mzroi,mzroiAll,Options)
% %Pre-allocate space for Augmented and summed intensities
% MSroi_aug=cell(max(mzroiAll(:,2)),size(MSroi,2));
% for T=1:size(mzroi,2)
%     MSroiT=MSroi(:,T);
%     %     [s,e]=deal(max(mzroiAll(:,2)),1);
%     %     s(1)=1;e(1)=0;
%     for k=1:max(mzroiAll(:,2))
%         %         e(k)=e+size(MSroiT{k},2);
%         Options.szMSroi(k,:)=[max(mzroiAll(:,4)),size(MSroiT{k},2)];
%         %         s(k)=s+size(MSroiT{k},2);
%     end
% % end
% Options.szMSroi(:,2)=max( Options.szMSroi(:,2));
% % MSroi_aug=zeros(max(mzroiAll(:,4)),Options.szMSroi(1,2),max(mzroiAll(:,2)),size(MSroi,2));
% 
% % for T=1:size(mzroi,2)
% %     MSroiT=MSroi(:,T);
%     parfor k=1:max(mzroiAll(:,2))
%         %         hWaitbarMsgQueue.send(0);
%         MSroiT{k}=SumXAug(mzroiAll,MSroiT{k}(:,Options.RtInt(k,T,1):Options.RtInt(k,T,2)),k,T,Options);
%         %         MSroi_aug(:,:,k,T)=SumXAug(mzroiAll,MSroiT{k}(:,Options.RtInt(k,T,1):Options.RtInt(k,T,2)),k,T,Options);
%     end
%     MSroi_aug(:,T)=MSroiT;
% end
% % MSroi_aug=sum(MSroi_aug,4);
% end
% 
% %%
% %Sum m/z across groups and Mass level 1 and 2 for pulsed Elution:
% function [MSroi_aug]=SumMZgroups4Pulsed2D(MSroi,mzroi,mzroiAll,Options)
% % szMSroi=[max(mzroiAll(:,4)),eScanMin];
% % MSroi_aug=sparse(max(mzroiAll(:,4)),sum(Options.sz));
% 
% for T=1:size(mzroi,2)
%     MSroiT=MSroi(:,T);
%     %     [s,e]=deal(max(mzroiAll(:,2)),1);
%     %     s(1)=1;e(1)=0;
%     for k=1:max(mzroiAll(:,2))
%         %         e(k)=e+size(MSroiT{k},2);
%         Options.szMSroi(k,:)=[max(mzroiAll(:,4)),size(MSroiT{k},2)];
%         %         s(k)=s+size(MSroiT{k},2);
%     end
%     parfor  k=1:max(mzroiAll(:,2))%par
%         %         hWaitbarMsgQueue.send(0);
%         MSroiT{k}=SumXAug(mzroiAll,MSroiT{k},k,T,Options);%szMSroi(k,:));
%     end
%     MSroi_aug=MSroiT{1};
%     if size(MSroiT,1)>1
%     for k=2:max(mzroiAll(:,2))
%         MSroi_aug=cat(2,MSroi_aug,MSroiT{k});
%     end
%     end 
% end
% end
