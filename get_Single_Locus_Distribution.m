function [clineWidthSamples, afdiffs] = get_Single_Locus_Distribution(indexes)

if nargin < 1
     indexes = 1:2000; % for m = 0.05, s = 0.02
%    indexes = 2001:3000; % for m = 0.01, s = 0.02
%    indexes = 3001:4000; % for m = 0.02, s = 0.02
%    indexes = 4001:7000; % for m = 0.1, s = 0.02
end

n = numel(indexes);

orig = pwd;

clineWidthSamples = zeros(n * 1000,1);

startRow = 1;

for i = 1:n
    
    wi = indexes(i);
    
    data = importdata(['../sle_data/sle_test' num2str(wi) '/AlleleFreqTS.txt']);
    data = data.data;
    
    qGlobal = data(:,7);
    nr = numel(qGlobal);
    x = qGlobal(1:(nr-1));
    y = qGlobal(2:nr);
    x = x - 0.5;
    y = y - 0.5;
    prod = x .* y; % this will be negative when x and y are on opposite sides of 0.5
    lastCrossing = find(prod < 0, 1, 'last') + 1;
    
    if numel(lastCrossing) > 0
        
        afdiffs = (data(1:lastCrossing,9) - data(1:lastCrossing,8));
        
        clineWidthSamples(startRow:(startRow + lastCrossing - 1)) = afdiffs;
        
        startRow = startRow + lastCrossing;
        
    end
end

clineWidthSamples(startRow:(numel(clineWidthSamples))) = [];

afdiffs = clineWidthSamples;

clineWidthSamples = 1 ./ clineWidthSamples;
%clineWidthSamples = clineWidthSamples(isfinite(clineWidthSamples));


cd(orig);


nsamps = numel(clineWidthSamples);

fpt = fopen('ClineWidthSamples.txt', 'w');
for i = 1:nsamps
    fprintf(fpt, '%E %E\n', clineWidthSamples(i), afdiffs(i));
end
fclose(fpt);
