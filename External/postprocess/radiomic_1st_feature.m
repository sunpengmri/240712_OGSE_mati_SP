% compute the 1st order feature 
% refer:  
% Aerts HJ, Velazquez ER, Leijenaar RT, et al. Decoding tumour 
% phenotype by noninvasive imaging using a quantitative radiomics approach. Nature communications 2014;5:4006.
function fea = radiomic_1st_feature(t2img,maskimg)


if size(t2img) ~= size(maskimg)
    error('please check your mask img and t2 img, for they do not have the same size..\n');
else
    maskimg(maskimg~=0)=1;
    roi_img = maskimg.* t2img;
end

%%%% check the mask size also be feature volume size
temp = sum(maskimg(:));
if temp<=50
    fprintf('mask size is too small\n');
end

%%  %% 1st order feature base on 
% fea.vol = temp;
roi_vec = roi_img(find(maskimg~=0));
fea.size = length(roi_vec);
fea.energy = sum(roi_vec.^2,'omitnan')./length(roi_vec(roi_vec>0));
fea.kurtosis = kurtosis(roi_vec);
fea.maximum = max(roi_vec,[],'all');
fea.mean = mean(roi_vec,'omitnan');
fea.mad = mad(roi_vec);
fea.median = median(roi_vec,'omitnan');
fea.minimum = min(roi_vec,[],'all');
fea.range = range(roi_vec);
fea.rms = rms(roi_vec,'omitnan');
fea.skewness = skewness(roi_vec);
fea.std = std(roi_vec,'omitnan');
fea.var = var(roi_vec,'omitnan');
fea.pct10 =  prctile(roi_vec,10);
fea.pct25 =  prctile(roi_vec,25);
fea.pct50 =  prctile(roi_vec,50);
fea.pct75 =  prctile(roi_vec,75);
fea.pct90 =  prctile(roi_vec,90);
[fea.entropy fea.uniformity] = jt_entropy(roi_vec,100);