

%% Load Options file 
path2OptionsFile = 'D:\data\Oskar\ShortCommunication\CodeGit\roi-paper\ACN\regionofinterest\';
cd(path2OptionsFile)
S_OptionsStruct  

%% Get list of CDF files that needs to be subjected to ROI 
cd(Options.Paths.CDF)
fileList = dir('*.CDF');
Options.ROI.mzOptimization.SamplesForOptimization = randi([1,length(fileList)],Options.ROI.mzOptimization.nSamples,1);
save([Options.Paths.save2mat,'OptionsFile.mat'],"Options","path2OptionsFile")
%% Optimize 
%Pre-allocate space 
[mzerror_Sample,nmz] = deal(zeros(length(Options.ROI.mzOptimization.MZmultiplyFactors),length(Options.ROI.mzOptimization.SamplesForOptimization)));
parpool(3)
for n = 1:Options.ROI.mzOptimization.nSamples
    k = Options.ROI.mzOptimization.SamplesForOptimization(n);
    [mzerror_Sample(:,n),~,nmz(:,n)]=OptParamRoi(fileList(k).name,Options.ROI.mzOptimization.mzerror,Options.ROI.mzOptimization.MZmultiplyFactors,Options.ROI.minroi,Options.ROI);
end

%Calculate the median mzerrror of all the samples used in the optimizationn
%       
Options.ROI.mzerror = median(mzerror_Sample,'all');

%% Plot for visualisation
figure
p = plot( Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,nmz(:,n));
hold on
ps = scatter(Options.ROI.mzerror,nmz(Options.ROI.mzerror == Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,n),'red','filled');
legend([p,ps],{fileList(Options.ROI.mzOptimization.SamplesForOptimization).name,'selected m/z error'})
xlabel('\delta_m_/_z','FontAngle','italic')
ylabel('Number of m/z traces / ROIs')



%% ROI processing 
[mzroi,MSroi,Rt] = deal(cell(length(fileList),Options.ROI.NumTrace));
NbrePts = length(fileList)*Options.ROI.NumTrace;

close all
% Waitbar's Message Field
Msg = ['Sample ',num2str(1),' - ROI Progress...!'];
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);



for T = Options.ROI.NumTrace
    fprintf(1,'Trace: %i/%i\n',T,Options.ROI.NumTrace)

    for k = 1:length(fileList)
        if T == 1
            FileName = fileList(k).name;
        else 

            FileName = dir([Options.Paths.CDF_MS2,fileList(k).name(1:end-Options.ROI.RemoveSuffix),'.*']);
        end 
        hWaitbarMsgQueue.send(0);
              
        % ROI
        [mzroi{k,T},MSroi{k,T},~,     ~,      ~,     Rt,~,~] = ROIpeaks_ACN(FileName,Options.ROI);
    end

end
%Save to mat
save([Options.Paths.save2mat,'mzroi.mat'],'mzroi',"Options")
save([Options.Paths.save2mat,'MSroi.mat'],'mzroi','MSroi',"Rt","Options")



%% Augment
[mzroi_aug,MSroi_aug]=MSroiaug_ACN2(mzroi,MSroi,Options.ROI);

%Save to mat
save([Options.Paths.save2mat,'mzroi_aug.mat'],'mzroi_aug',"Options")
save([Options.Paths.save2mat,'MSroi_aug.mat'],'MSroi_aug',"Rt","Options")
    