%% freq extract
node = 11;
wdow = 2048;
fs = 20;
overlap = [];
nfft = 2^nextpow2(wdow)*8;
    
load('parked_dense_node_seed1To210and1001To1300_y.mat')
data1 = squeeze(dataset1(:,node,:,:));
% freq = freq_extract(data1,lower,upper)

pk_f = zeros(size(data1,[2,3]));
for nwd = 1:size(data1,2)
    for nsd = 1:size(data1,3)
        [pxx,f] = pwelch(data1(:,nwd,nsd),wdow,overlap,nfft,fs);
        [pk,loc] = findpeaks(...
            pxx(floor(1.45/(fs/2)*(nfft/2 + 1)):...
            floor(1.8/(fs/2)*(nfft/2 + 1))),...
            'SortStr','descend','NPeaks',1);
        pk_f(nwd,nsd) = f(floor(1.45/(fs/2)*(nfft/2 + 1))+loc-1);
    end
end
freq_sim_owt1 = pk_f';% ref
freq_sim_owt2 = freq_sim_owt1;% ref

load('parked_water25m_seed211To420and1301To1600_y.mat')
% [pxx1,f1] = pwelch(data1(:,1,1),wdow,overlap,nfft,fs);
% [pxx2,f2] = pwelch(data1(:,51,1),wdow,overlap,nfft,fs);
% plot(f1,log10(pxx1),f2,log10(pxx2))
data1 = squeeze(dataset1(:,node,:,:));
pk_f = zeros(size(data1,[2,3]));
for nwd = 1:size(data1,2)
    for nsd = 1:size(data1,3)
        [pxx,f] = pwelch(data1(:,nwd,nsd),wdow,overlap,nfft,fs);
        [pk,loc] = findpeaks(...
            pxx(floor(1.4/(fs/2)*(nfft/2 + 1)):...
            floor(1.8/(fs/2)*(nfft/2 + 1))),...
            'SortStr','descend','NPeaks',1);
        pk_f(nwd,nsd) = f(floor(1.4/(fs/2)*(nfft/2 + 1))+loc-1);
    end
end
freq_sim_owt3 = pk_f';%ref 

% plot(freq_train_owt1(:),'.');hold on
% plot(freq_train_owt3(:),'.');hold off
figure
plot(freq_sim_owt1(:),freq_sim_owt3(:),'.')

load('parked_water20m_dam0.1_seed1891To1900_y.mat')
data1 = squeeze(dataset1(:,node,:,:));
pk_f = zeros(size(data1,[2,3]));
for nwd = 1:size(data1,2)
    for nsd = 1:size(data1,3)
        [pxx,f] = pwelch(data1(:,nwd,nsd),wdow,overlap,nfft,fs);
        [pk,loc] = findpeaks(...
            pxx(floor(1.45/(fs/2)*(nfft/2 + 1)):...
            floor(1.7/(fs/2)*(nfft/2 + 1))),...
            'SortStr','descend','NPeaks',1);
        pk_f(nwd,nsd) = f(floor(1.45/(fs/2)*(nfft/2 + 1))+loc-1);
    end
end
freq_sim_owt1_dam = pk_f';%test_dam

plot(1:510,freq_sim_owt1(:),'.',1:510,freq_sim_owt3(:),'.',...
    1:510,freq_sim_owt1_dam(:),'.')
legend('undamaged 1','undamaged 2','damaged 1')
xlabel('water level')
ylabel('feature')