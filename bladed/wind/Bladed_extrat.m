clear

sd = 101:102;
NSD = length(sd);
%% acceleration
dataset1 = zeros(65536*NSD,1);
for nsd = 1:NSD
    fileID = fopen(['C:\loong\BLADED WORKPLACE\wind',...
        '\wind_seed',num2str(sd(nsd)),...
        '\wind.wnd']);
    a = fread(fileID,[1,3*25*25*65536],'float32');
    fclose(fileID);
    b = reshape(a,3,25,25,65536);
    c = squeeze(b(1,12,12,:))';

    dataset1(65536*(nsd-1)+1:65536*nsd,:) = c;
end
save(['wind',...
    num2str(sd(1)),'To',num2str(sd(end)),'.mat'],'dataset1');



