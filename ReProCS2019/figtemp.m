thresh_pt = 5e-1;

transN = zeros(10, 10);
transA = zeros(10, 10);

for ii = 1 : 10
    for jj = 1  : 10
        for kk = 1 : 100
            if(PhaseTransNORST(ii, jj, kk) <= thresh_pt)
                transN(ii, jj) = transN(ii, jj) + 1;
            end
            
            if(PhaseTransAltProj(ii, jj, kk) <= thresh_pt)
                transA(ii, jj) = transA(ii, jj) + 1;
            end
        end
    end
end

transN
transA


figure
subplot(211)
imagesc(transN)
colormap('gray');
subplot(212)
imagesc(transA)
colormap('gray');
       
        
X = transN;
[numi, numj] = size(X);

%co-ordinate axis labels
% rrange = [1 : 10];
% nrange = unique(ceil(linspace(50, 1000, 10)));


FileName = 'data_TIT/PhaseTransNORST.dat';
fileID = fopen(FileName, 'w'); %can vary read/write/overwrite if needed

%note here that all entries are integers
for ii = 1 : numi
    for jj = 1 : numj
        fprintf(fileID, '%.2f %d %d\n', outfracrowrange(ii), rrange(jj), X(ii , jj));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);
