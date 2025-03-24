%% This code is used for post OGSE analysis using the MatiGUI app
% 1. calculate ADC maps and rADC
% 2. calculate cellularity
% 3. calculate histogram analysis 

%% set path and parameters
function [] = post_mati(folder_path,result_path)
% nii data obtained from mrtrix after diffusion preprocess
project_dir = folder_path;
% subject=subject_name;
result_path=result_path;
% mkdir(result_path)

pgse_name='pgse.nii.gz';
ogsen1_name= 'ogsen1.nii.gz';
ogsen2_name= 'ogsen2.nii.gz';
ogse_name='ogse.nii.gz';
jinyinhu = 1;

ogsen2_nii=[project_dir  '/' 'ogsen2.nii.gz'];
ogse_mc_nii=[project_dir  '/' 'ogse_mc.nii.gz'];

if isfile(ogsen2_nii)
    nii_combined = 0;
else
    nii_combined = 1;
end

if isfile(ogse_mc_nii)
    isMC = 1;
else
    isMC = 0;
end

if (nii_combined==1)
    if (isMC==0)
        ogse_name='ogse.nii.gz';
    else
        ogse_name='ogse_mc.nii.gz';
    end
else
    pgse_name='pgse.nii.gz';
    ogsen1_name= 'ogsen1.nii.gz';
    ogsen2_name= 'ogsen2.nii.gz';
end
if jinyinhu
    pgse_bvalue_raw = [0,250,500,750,1000,1400,1800];
    ogseN1_bvalue_raw = [250,500,650,800];
    ogseN2_bvalue_raw = [100,200,250];
    
    pgse_bvalue_fitting = [0,250,500,750,1000,1400,1800];
    ogseN1_bvalue_fitting = [0,250,500,650,800];
    ogseN2_bvalue_fitting = [0,100,200,250];
else
    pgse_bvalue_raw = [0,250,500,1000,1500];
    ogseN1_bvalue_raw = [0,250,500,800,1200];
    ogseN2_bvalue_raw = [0,100,200,300];
    
    pgse_bvalue_fitting = [0,250,500,1000,1500];
    ogseN1_bvalue_fitting = [0,250,500,800,1200];
    ogseN2_bvalue_fitting = [0,100,200,300];

end

[~, pgse_idx] = ismember(pgse_bvalue_fitting, pgse_bvalue_raw);
[~, ogsen1_idx] = ismember(ogseN1_bvalue_fitting, ogseN1_bvalue_raw);
[~, ogsen2_idx] = ismember(ogseN2_bvalue_fitting, ogseN2_bvalue_raw);


% impulsed results
vin_nii=[result_path '/' 'vin.nii.gz']; 
d_nii=[result_path '/' 'd.nii.gz']; 
Dex_nii=[result_path '/' 'Dex.nii.gz']; 

% mask
mask_whole_nii='mask_tumor.nii.gz';
mask_tumor_nii='mask_tumor.nii.gz';
mask_other_nii='mask_tumor.nii.gz';


mask_whole_nii=[project_dir '/'  mask_whole_nii]; 
mask_tumor_nii=[project_dir '/'  mask_tumor_nii]; 
mask_other_nii=[project_dir '/'  mask_other_nii]; 

info_mask_whole=niftiinfo(mask_whole_nii);
mask = double(niftiread(info_mask_whole));
info_mask_tumor=niftiinfo(mask_tumor_nii);
mask_tumor = double(niftiread(info_mask_tumor));
info_mask_other=niftiinfo(mask_other_nii);
mask_other = double(niftiread(info_mask_other));

% 计算mask_tumor中非零元素的层面，并选出中间层面
non_zero_slices = find(squeeze(sum(sum(mask_tumor, 1), 2)) > 0);
slicei = non_zero_slices(round(length(non_zero_slices) / 2));

%% 01. loda and reorganize data
if (nii_combined == 1)
    % readin combined ogse data
    ogse_nii=[project_dir '/' ogse_name]; 
    info_ogse=niftiinfo(ogse_nii);
    ogse_raw = double(niftiread(info_ogse));
    ogse_raw=rot90(flipud(ogse_raw),3); % this for mrtrix convert nii
    % ogse_raw=rot90(ogse_raw,1); % this is for mrcrongl
    % imshow(ogse_raw(:,:,5,1),[])
    mask=rot90(flipud(mask),3); 
    mask_tumor=rot90(flipud(mask_tumor),3); 
    mask_other=rot90(flipud(mask_other),3);
else % 
    pgse_nii=[project_dir '/'  'pgse.nii.gz']; 
    info_pgse=niftiinfo(pgse_nii);
    pgse_raw = double(niftiread(info_pgse));
    % pgse_raw=rot90(flipud(pgse_raw),3); % breast
    pgse_raw=rot90((pgse_raw),1); % head nii from mricron
    % imshow(pgse_raw(:,:,5,1),[])
    ogsen1_nii=[project_dir '/' 'ogsen1.nii.gz']; 
    info_ogsen1=niftiinfo(ogsen1_nii);
    ogsen1_raw = double(niftiread(info_ogsen1));
    % ogsen1_raw=rot90(flipud(ogsen1_raw),3); % breast
    ogsen1_raw=rot90((ogsen1_raw),1); % head nii from mricron
    % imshow(ogsen1_raw(:,:,5,1),[])

    ogsen2_nii=[project_dir '/' 'ogsen2.nii.gz']; 
    info_ogsen2=niftiinfo(ogsen2_nii);
    ogsen2_raw = double(niftiread(info_ogsen2));
    % ogsen2_raw=rot90(flipud(ogsen2_raw),3); % breast
    ogsen2_raw=rot90((ogsen2_raw),1); % head nii from mricron
    % imshow(ogsen2_raw(:,:,5,1),[])
    % mask=rot90(flipud(mask),3); % breast
    % mask_tumor=rot90(flipud(mask_tumor),3); % breast
    % mask_other=rot90(flipud(mask_other),3); % breast

    mask=rot90((mask),1); % head nii from mricron
    mask_tumor=rot90((mask_tumor),1); % head nii from mricron
    mask_other=rot90((mask_other),1); % head nii from mricron
end

info_vin=niftiinfo(vin_nii);
vin_map = double(niftiread(info_vin));
% ogsen2_raw=rot90(flipud(ogsen2_raw),3); % breast
vin_map=rot90((vin_map),1); % head nii from mricron

info_d=niftiinfo(d_nii);
d_map = double(niftiread(info_d));
% ogsen2_raw=rot90(flipud(ogsen2_raw),3); % breast
d_map=rot90((d_map),1); % head nii from mricron

info_Dex=niftiinfo(Dex_nii);
Dex_map = double(niftiread(info_Dex));
% ogsen2_raw=rot90(flipud(ogsen2_raw),3); % breast
Dex_map=rot90((Dex_map),1); % head nii from mricron

%%%% check Data and Mask
figure(1) % raw data and mask check
if (nii_combined == 1)
    subplot(3,1,1)
    imshow(squeeze(ogse_raw(:,:,slicei,1)),[]);
    title("ogse")
    
    subplot(3,1,2)
    imshow(squeeze(mask_tumor(:,:,slicei)),[]); % multi slice
    % imshow(mask_tumor,[]); % multi slice
    title("maskTumor")

    subplot(3,1,3)
    imshow(squeeze(mask_other(:,:,slicei)),[]); % multi slice
    % imshow(mask_tumor,[]); % multi slice
    title("maskOther")

else
    subplot(2,2,1)
    imshow(squeeze(pgse_raw(:,:,slicei,1)),[]);
    title("pgse")
    
    subplot(2,2,2)
    imshow(squeeze(ogsen1_raw(:,:,slicei,1)),[]);
    title("ogseN1")
    
    subplot(2,2,3)
    imshow(squeeze(ogsen2_raw(:,:,slicei,1)),[]);
    title("ogseN2")
    
    subplot(2,2,4)
    imshow(mask_tumor(:,:,slicei),[]);
    title("mask3")
end


%%
Nacq_ogsen1 = length(ogseN1_bvalue_raw);% total number of acquisition points
Nacq_ogsen2 = length(ogseN2_bvalue_raw);% 
Nacq_pgse = length(pgse_bvalue_raw);
if nii_combined==1 & ~jinyinhu
    ogsen1_raw=ogse_raw(:,:,:,1:Nacq_ogsen1);
    ogsen2_raw=ogse_raw(:,:,:,1+Nacq_ogsen1:Nacq_ogsen1+Nacq_ogsen2);
    pgse_raw=ogse_raw(:,:,:,1+Nacq_ogsen1+Nacq_ogsen2:end);
end

if jinyinhu
   pgse_raw=ogse_raw(:,:,:,1:Nacq_pgse);
   ogsen1_raw(:,:,:,2:Nacq_ogsen1+1)=ogse_raw(:,:,:,1+Nacq_pgse:Nacq_pgse+Nacq_ogsen1);
   ogsen1_raw(:,:,:,1)=ogse_raw(:,:,:,1);
   ogsen2_raw(:,:,:,2:Nacq_ogsen2+1)=ogse_raw(:,:,:,1+Nacq_ogsen1+Nacq_pgse:end);
   ogsen2_raw(:,:,:,1)=ogse_raw(:,:,:,1);
   ogsen1_idx=ogsen1_idx+1;
   ogsen2_idx=ogsen2_idx+1;
end

% extract diff images for ADC fitting
pgse_adc_fitting=pgse_raw(:,:,:,pgse_idx);
ogsen1_adc_fitting=ogsen1_raw(:,:,:,ogsen1_idx);
ogsen2_adc_fitting=ogsen2_raw(:,:,:,ogsen2_idx);

%% fitting
warning('off', 'all');
[PGSE_ADC_Map] = ADCMap(pgse_adc_fitting, pgse_bvalue_fitting, 0);
PGSE_ADC_Map(PGSE_ADC_Map<0 | PGSE_ADC_Map>0.01) = nan;
[OGSEN1_ADC_Map] = ADCMap(ogsen1_adc_fitting, ogseN1_bvalue_fitting, 0);
OGSEN1_ADC_Map(OGSEN1_ADC_Map<0 | OGSEN1_ADC_Map>0.01) = nan;
[OGSEN2_ADC_Map] = ADCMap(ogsen2_adc_fitting, ogseN2_bvalue_fitting, 0);
OGSEN2_ADC_Map(OGSEN2_ADC_Map<0 | OGSEN2_ADC_Map>0.01) = nan;
warning('on', 'all');
% rADC
delta_ADC_Map_N1=(OGSEN1_ADC_Map-PGSE_ADC_Map)./(OGSEN1_ADC_Map+eps).*100;
delta_ADC_Map_N2=(OGSEN2_ADC_Map-PGSE_ADC_Map)./(OGSEN2_ADC_Map+eps).*100;

PGSE_ADC_Map=PGSE_ADC_Map.*mask;
OGSEN1_ADC_Map=OGSEN1_ADC_Map.*mask;
OGSEN2_ADC_Map=OGSEN2_ADC_Map.*mask;
delta_ADC_Map_N1=delta_ADC_Map_N1.*mask;
delta_ADC_Map_N2=delta_ADC_Map_N2.*mask;

% cellularity calculation
cellularity_map=vin_map./d_map .* 100;


%% ADC map
figure(10)
subtitle('Different ADC maps from different OGSE pulse');
subplot(2,4,1)
imshow(PGSE_ADC_Map(:,:,slicei),[0,0.003]);colormap jet;
title("PGSE")
colorbar
subplot(2,4,2)
imshow(OGSEN1_ADC_Map(:,:,slicei),[0,0.003]);colormap jet;
title("OGSEN1")
colorbar
subplot(2,4,3)
imshow(OGSEN2_ADC_Map(:,:,slicei),[0,0.003]);colormap jet;  
title("OGSEN2")
colorbar

colorbar
subplot(2,4,4)
imshow(delta_ADC_Map_N1(:,:,slicei),[0,50]);colormap jet;  
title("rADC1")
colorbar

subplot(2,4,5)
imshow(delta_ADC_Map_N2(:,:,slicei),[0,50]);colormap jet;  
title("rADC2")
colorbar
saveas(gcf, [result_path 'ADC.png'])

subplot(2,4,6)
imshow(vin_map(:,:,slicei),[0,0.4]);colormap jet;  
title("vin")
colorbar
saveas(gcf, [result_path 'ADC.png'])

subplot(2,4,7)
imshow(d_map(:,:,slicei),[0,30]);colormap jet;  
title("d")
colorbar
saveas(gcf, [result_path 'ADC.png'])

subplot(2,4,8)
imshow(cellularity_map(:,:,slicei),[0,2]);colormap jet;  
title("cellularity")
colorbar
saveas(gcf, [result_path 'ADC.png'])

%% histogram analysis
%% statistics saving
map_name = ["Vin" "Diameter" "Cellularity" "Dex" ...
    "PGSE_ADC" "OGSEN1_ADC" "OGSEN2_ADC" "delta_ADC_N1" "delta_ADC_N2"]; 

maps = zeros(size(mask,1),size(mask,2),size(mask,3),length(map_name));
maps(:,:,:,1) = vin_map;
maps(:,:,:,2) = d_map;
maps(:,:,:,3)= cellularity_map;
maps(:,:,:,4)=Dex_map;
maps(:,:,:,5)=PGSE_ADC_Map;
maps(:,:,:,6)=OGSEN1_ADC_Map;
maps(:,:,:,7)=OGSEN2_ADC_Map;
maps(:,:,:,8)=delta_ADC_Map_N1;
maps(:,:,:,9)=delta_ADC_Map_N2;


% map_name = ["Vin" "Diameter" "Cellularity" "Dex" ...
    % "PGSE_ADC" "OGSEN1_ADC" "OGSEN2_ADC" "delta_ADC_N1" "delta_ADC_N2"]; 
stats = zeros(length(map_name),20); 
for i = 1:length(map_name)
        % relocate data for statistics
        map_data = squeeze(maps(:,:,:,i));
        map_tag = char(map_name(i));
        fea_tumor=radiomic_1st_feature(map_data,mask_tumor);
        fea_control=radiomic_1st_feature(map_data,mask_other);
        % statistic 
        stats(i, 1) = fea_tumor.mean;
        stats(i, 2) = fea_tumor.median;
        stats(i, 3)  = fea_tumor.pct10;
        stats(i, 4) = fea_tumor.pct25;
        stats(i,5) =  fea_tumor.pct50;
        stats(i,6) =  fea_tumor.pct75;
        stats(i,7) =  fea_tumor.pct90;
        stats(i, 8) = fea_tumor.skewness;
        stats(i, 9) = fea_tumor.kurtosis;
        stats(i, 10) = fea_tumor.uniformity;
        stats(i, 11) = fea_control.mean;
        stats(i, 12) = fea_control.median;
        stats(i, 13)  = fea_control.pct10;
        stats(i, 14) = fea_control.pct25;
        stats(i,15) =  fea_control.pct50;
        stats(i,16) =  fea_control.pct75;
        stats(i,17) =  fea_control.pct90;
        stats(i, 18) = fea_control.skewness;
        stats(i, 19) = fea_control.kurtosis;
        stats(i, 20) = fea_control.uniformity;
        % histgram
        is_histgram_figure=0;
        if is_histgram_figure
            tumor_index=find(mask_tumor==1);
            normal_index=find(mask_other==1);
            tumor_roi = map_data(tumor_index);
            normal_roi = map_data(normal_index);
            figure(i+20)
            set(figure(i+20),'Position',[100,100,800,600]);
            %     subplot(1,2,1);
            h1 = histogram(tumor_roi);
            h1.Normalization = 'probability';h1.BinWidth = 0.0025;legend('TumorRaw')
            hold on;
            h2 = histogram(normal_roi);    
            h2.Normalization = 'probability';h2.BinWidth = 0.0025;legend('Tumor','Normal');
            title(map_tag);
            figure_path = result_path;figure_name = [map_tag '.png'];
            saveas(gcf,[figure_path figure_name])
            hold off;
        end
end
    T = array2table([map_name',stats],'VariableNames',{'parameters' 'fea_tumor.mean' 'fea_tumor.median' 'fea_tumor.pct10' 'fea_tumor.pct25' 'fea_tumor.pct50' 'fea_tumor.pct75' 'fea_tumor.pct90' ...
        'fea_tumor.skewness' 'fea_tumor.kurtosis' ' fea_tumor.uniformity' 'fea_control.mean' 'fea_control.median' 'fea_control.pct10' 'fea_control.pct25' 'fea_control.pct50' 'fea_control.pct75' 'fea_control.pct90' ...
         'fea_control.skewness' 'fea_control.kurtosis' ' fea_control.uniformity'});
    writetable(T, [result_path 'result.xls']);


%% data saving
%ADCMap
PGSE_ADC_slice_nii=PGSE_ADC_Map;
temp(:,:,:) = PGSE_ADC_slice_nii;
PGSE_ADC=temp;
PGSE_ADC_nii=[result_path '/' 'PGSE_ADC']; 

OGSEN1_ADC_slice_nii=OGSEN1_ADC_Map;
temp(:,:,:) = OGSEN1_ADC_slice_nii;
OGSEN1_ADC=temp;
OGSEN1_ADC_nii=[result_path '/' 'OGSEN1_ADC']; 

OGSEN2_ADC_slice_nii=OGSEN2_ADC_Map;
temp(:,:,:) = OGSEN2_ADC_slice_nii;
OGSEN2_ADC=temp;
OGSEN2_ADC_nii=[result_path '/' 'OGSEN2_ADC'];

delta_ADC_N1_slice_nii=delta_ADC_Map_N1;
temp(:,:,:) = delta_ADC_N1_slice_nii;
delta_ADC_N1=temp;
delta_ADC_N1_nii=[result_path '/' 'delta_ADC_N1']; 

delta_ADC_N2_slice_nii=delta_ADC_Map_N2;
temp(:,:,:) = delta_ADC_N2_slice_nii;
delta_ADC_N2=temp;
delta_ADC_N2_nii=[result_path '/' 'delta_ADC_N2']; 

cellularity_map_nii=cellularity_map;
% temp(:,:,:) = cellularity_map_nii;
% cellularity_nii=temp;
cellularity_nii=[result_path '/' 'cellularity']; 


s_export_nii(mask_whole_nii,PGSE_ADC, PGSE_ADC_nii);
s_export_nii(mask_whole_nii,OGSEN1_ADC, OGSEN1_ADC_nii);
s_export_nii(mask_whole_nii,OGSEN2_ADC, OGSEN2_ADC_nii);
s_export_nii(mask_whole_nii,delta_ADC_N1, delta_ADC_N1_nii);
s_export_nii(mask_whole_nii,delta_ADC_N2, delta_ADC_N2_nii);
s_export_nii(mask_whole_nii,cellularity_map_nii, cellularity_nii);
% s_export_nii_mask(mask_whole_nii,tumor, tumor_nii);
