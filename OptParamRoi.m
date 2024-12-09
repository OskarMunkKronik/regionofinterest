%%ROI Optimization
function [mzDev,minroi,nmz]=OptParamRoi(FileName,MassAccuracy,MZmultiplyFactors,Vecminroi,Options)%,thresh,ppm,wmean,sScan,eScan)

%waitbar
  close all force
    NbrePts =length(MZmultiplyFactors*Vecminroi);
        % Waitbar's Message Field
    Msg = 'Determination of the mass accuracy...!';
    % Create ParFor Waitbar
    [~,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);
Vecmzerror=repelem(MassAccuracy,length(MZmultiplyFactors)).*MZmultiplyFactors;

% %% Optimization
% OptionsROI = struct(    
ppm = true;
    thresh = Options.thresh;
 wmean = Options.wmean;
 GapAllowed = Options.GapAllowed;
 prefilter = Options.prefilter;
 verbose   = Options.verbose;
 fillIn    = Options.fillIn;

nmz=zeros(length(Vecminroi),length(Vecmzerror));

for mr=1:length(Vecminroi)
    minroi=Vecminroi(mr);
    lenMZ=length(Vecmzerror);
    fprintf(1,'minroi: %i/%i\n',mr,length(Vecminroi))

    parfor mass=1:length(Vecmzerror)
         hWaitbarMsgQueue.send(0);
           
        % Options.mzerror=Vecmzerror(mass);
      
        mzroi=ROIpeaks_ACN(FileName,struct('mzerror',Vecmzerror(mass),'minroi',minroi,'ppm',ppm,'thresh',thresh,'wmean',wmean,'GapAllowed',GapAllowed,'prefilter',prefilter,'verbose',verbose,'fillIn',fillIn));%,thresh,mzerror,minroi,sScan,eScan,ppm,wmean);
     
        nmz(mr,mass)=length(mzroi);
          end
end
[minroi,mzDev]=find(nmz(:,:,1)==max(nmz(:,:,1),[],'all'));
if length(Vecminroi)>1
    figure
    [XCo,YCo]=meshgrid(Vecminroi,Vecmzerror);
    surf(XCo,YCo,nmz')
    colorbar
    axis tight
    
    view(0,90)
else
    figure
    plot(Vecmzerror,nmz,'-o','Color','b','LineWidth',3)
    if Options.ppm==1
        xlabel('Mass Deviation (ppm)')
    else
        xlabel('Mass Deviation (Da)')
    end
    ylabel('number of ROIs')
    hold on
    scatter(Vecmzerror(mzDev),nmz(minroi,mzDev),[],'filled','r')
    axis tight
end
minroi=Vecminroi(minroi);
mzDev=Vecmzerror(mzDev);

end