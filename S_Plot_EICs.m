%Plotting the results from the augmented matrix. k is the sample index, T
%is the MS level 

mzTarget = 268.191; % Target m/z
sz = size(MSroi_aug);
col = 'rb';
sign = [1,-1];
h = 1;
clf
for k = 1:sz(1)
    
    for T = 1:sz(2)
        subplot(sz(1),2,h)
        [~,a] = min(abs(mzroi_aug-mzTarget));
        p(T) = plot(Rt{k,T}./60,MSroi_aug{k,T}(a,:),Color=col(T));

        hold on
    
        subplot(sz(1),2,h+1)
        [~,max_ind] = max(MSroi_aug{k,T}(a,:));
        stem(mzroi_aug,MSroi_aug{k,T}(:,max_ind).*sign(T),'Color',col(T),'Marker','none')
        hold on 

    end
    subplot(sz(1),2,h)
    legend(p,{'MS1','MS2'});
    xlabel('Rt (min)')
    ylabel('Counts')
    title(fileList(k).name)
    
    subplot(sz(1),2,h+1)
    legend(p,{'MS1','MS2'});
    xlabel('m/z')
    ylabel('Counts')
    h = h+2;
    title(fileList(k).name)
    
end
