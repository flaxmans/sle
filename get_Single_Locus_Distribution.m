function clineWidthSamples = get_Single_Locus_Distribution(indexes)

if nargin < 1
    indexes = 1:2000;
end

n = numel(indexes);

orig = pwd;

clineWidthSamples = zeros(n * 1000,1);

startRow = 1;

for i = 1:n
    
    wi = indexes(i);
    
    data = importdata(['../sle_tests/sle_test' num2str(wi) '/AlleleFreqTS.txt']);
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
        
        clineData = data(1:lastCrossing,10);
        
        clineWidthSamples(startRow:(startRow + lastCrossing - 1)) = clineData;
        
        startRow = startRow + lastCrossing;
        
    end
end

clineWidthSamples(startRow:(numel(clineWidthSamples))) = [];

cd(orig);

clineWidthSamples = clineWidthSamples(isfinite(clineWidthSamples));

nsamps = numel(clineWidthSamples);

fpt = fopen('ClineWidthSamples.txt', 'w');
for i = 1:nsamps
    fprintf(fpt, '%E\n', clineWidthSamples(i));
end
fclose(fpt);
