clear
%% acceleration
name = 'water23m';
wd = [11];
NWD = length(wd);
dataset1 = zeros(1200*10,26,NWD,10);

for nwd = 1:NWD
    sd = (1:10)+floor((wd(nwd)-3)/2)*10;
    NSD = length(sd);
    for nsd = 1:NSD
        fileID = fopen(['C:\loong\BLADED WORKPLACE\',...
            '1.OC3_Monopile\11.windfarm\',name,...
            '\',name,'_wind',num2str(wd(nwd)),'_seed',num2str(sd(nsd)),...
            '\',name,'_wind',num2str(wd(nwd)),'_seed',num2str(sd(nsd)),...
            '.$58']);
        a = fread(fileID,[1,6*26*(1200*10)],'float32');%6 DOF; 8 sensor;
        % 6 rows; 5 colomns; 1200*10 pages (can be cut off at any time)
        fclose(fileID);
        b = reshape(a,6,26,(1200*10));
        c = squeeze(b(1,:,:))';% select 1st DOF; sensor 1:8.
        dataset1(:,:,nwd,nsd) = c;
    end
end
save([name,'_seed',...
    num2str(sd(1)),'To',num2str(sd(end)),'.mat'],'dataset1');





