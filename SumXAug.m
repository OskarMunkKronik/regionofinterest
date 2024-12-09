function [X,posXrNew]=SumXAug(mzroiAll,MSroi,k,T,Options)%szMSroi)

pos=mzroiAll(:,2)==k & mzroiAll(:,5)==T;
[posXr(:,1),posXc(:,1),ValX(:,1)]=find(MSroi);

posX= [posXr,posXc,ValX];
posX=sortrows(posX,1);



mzroiSample=mzroiAll(pos,:);
posXrNew=zeros(length(posXr),1);
h=1;
szPosX=size(posX);
posX=cat(1,posX,repelem(nan,szPosX(2)));
szPosX=size(posX);
for mz=1:size(mzroiSample,1)
    while mzroiSample(mz,3)==posX(h,1) && h<=szPosX(1)
        posXrNew(h)=mzroiSample(mz,4);
        
        h=h+1;
    end
end

X=accumarray([posXrNew,posX(1:end-1,2)],posX(1:end-1,3),[Options.szMSroi(k,1),Options.RtInt(k,T,2)],[],[],1);
end