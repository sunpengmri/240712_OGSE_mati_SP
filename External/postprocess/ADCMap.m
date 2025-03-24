%for M0 fitting
%angus 2019.5.19
function [Map] = ADCMap(indata, bvalues, minimum)
 % tes is a vector of the te values in msec.
 % indata is assumed to be dimensioned as X * Y * te - i.e., same slice, multiple te
 % ADCmap
indata = indata(:,:,:,2:end)./(indata(:,:,:,1)+eps);
bvalues = bvalues(2:end); 
rows = size(indata,1);
cols = size(indata,2);
N    = size(indata,3);
ly   = length(bvalues);
bvalues = reshape(bvalues, ly, 1);  % must be a columns vector

Map = zeros(rows,cols);       
indata = abs(indata);
 % Threshold images before fitting
thresh = sum( (abs(indata(:,:,:,:)) > minimum), 4) > (ly-1);
inData = abs(indata);   % in case user has input complex data
for s=1:N
    for r=1:rows
        for c=1:cols
            if thresh(r,c,s) 
                ydata = inData(r,c,s,:);
                ydata = log(reshape(ydata,ly,1));
                lfit = ydata\-bvalues;        
                Map(r,c,s) = lfit(1);   
            end
        end
    end
end
Map = 1./(Map+eps);
 
