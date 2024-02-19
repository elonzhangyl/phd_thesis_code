function FREQ_cluster_mean = freq_extract_oma(data1,n_order,...
    freq_bot,freq_top,maxdist,sucRate)
% sucRate=0.2;maxdist=0.01;freq_bot = 1.45;freq_top = 1.7;



fs = 20;
FREQ_cluster_mean = zeros(size(data1,[2,3]))';
for nwd = 1:size(data1,2)
    nwd
    for nsd = 1:size(data1,3)
        [w, frf] = myspectrum(data1(:,nwd,nsd), fs);
        % selection of the FRF range betwen 0.02 and 0.38 Hz
        [w1, frf1] = select_frf(w, frf, 0.02, 2.5) ;
        % direct identification at order 12
%         for n_order = 20:30
        [~, ~, ~, FST] = lscf(w1, frf1, n_order);


        %% clustering
        Temp = [];
        for i = 1:length(FST)
            temp = cell2mat(FST(i));
            Temp = [Temp;temp(:)]; % make cell to vector
        end

        FST = [Temp, zeros(length(Temp),1)];% make 2 domensional
        % clustering within 0.01 distance
        clustersXY = clusterXYpoints(FST,maxdist,1,'centroid','merge');
        % succuss rate > 20%
        clustersXY_slt = [];
        for i = 1:length(clustersXY)
            if length(cell2mat(clustersXY(i))) > (max(n_order)-min(n_order))*sucRate
                clustersXY_slt = [ clustersXY_slt; clustersXY(i)];
            end
        end
        % mean of each cluster
        clustersXY_mean = [];
        for i = 1:length(clustersXY_slt)
            temp_mean = mean(cell2mat(clustersXY_slt(i)));
            clustersXY_mean = [clustersXY_mean; temp_mean];
        end

        clustersXY_mean = clustersXY_mean(:,1);
        %% end of clustering
        % select SS2
        for i = 1:length(clustersXY_mean)
            if clustersXY_mean(i)> freq_bot && clustersXY_mean(i)<freq_top
                freq_cluster_mean = clustersXY_mean(i);
            end
        end

        FREQ_cluster_mean(nsd,nwd) = freq_cluster_mean;
    end
end

end

