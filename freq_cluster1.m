
FST_all_records_cell = FST;


Temp = [];
for i = 1:length(FST)
    temp = cell2mat(FST(i));
    Temp = [Temp;temp(:)]; % make cell to vector
end

FST = [Temp, zeros(length(Temp),1)];% make 2 domensional
% clustering within 0.01 distance
[clustersCentroids,clustersGeoMedians,clustersXY] = ...
    clusterXYpoints(FST,0.01,1,'point','merge');
% succuss rate > 20%
clustersXY_slt = [];
for i = 1:length(clustersXY)
    if length(cell2mat(clustersXY(i))) > (max(n_order)-min(n_order))*0.2
        clustersXY_slt = [ clustersXY_slt; clustersXY(i)];
    end
end
% mean of each cluster
for i = 1:length(clustersXY_slt)
    temp_mean = mean(cell2mat(clustersXY_slt(i)));
    clustersXY_mean = [clustersXY_mean; temp_mean];
end

clustersXY_mean = clustersXY_mean(:,1);

    
