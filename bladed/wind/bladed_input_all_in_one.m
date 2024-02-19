clear

filepath = 'C:\loong\BLADED WORKPLACE\1.OC3_Monopile\0.pool2004\';
% template = 'und.$PJ';
template = 'wind.$PJ';
expression1 = '101.wnd';
expression2 = 'SEED	101';
expression3 = 'UBAR	 11.4';

wd = 3:4;
NWD = length(wd);
for nwd = 1:NWD
    sd = (1:10)+wd(nwd)*100;
    NSD = length(sd);
    for nsd = 1:NSD
        % creat folder
        mkdir([filepath,template(1:end-4),'_',num2str(wd(nwd)),...
            '_seed',num2str(sd(nsd))]);
    end
end

for nwd = 1:NWD
    sd = (1)+wd(nwd)*100;
    NSD = length(sd);
    for nsd = 1:NSD
        %read
        Str = strsplit(fileread(template), '\n');
    %     fclose(template);
        %edit
        replace1 = [template(1:end-4),'_',num2str(wd(nwd)),'_seed',num2str(sd(nsd)),'.wnd'];
        newStr1 = regexprep(Str,expression1,replace1);
%         replace2 = ['SEED	',num2str(sd(nsd))];
%         newStr2 = regexprep(newStr1,expression2,replace2);
        replace3 = ['UBAR	 ',num2str(wd(nwd))];
        newStr3 = regexprep(newStr1,expression3,replace3);
        %write
        filePh = fopen([filepath,template(1:end-4),'_',...
            num2str(wd(nwd)),'_seed',num2str(sd(nsd)),...    
            '\',template(1:end-4),'_',num2str(wd(nwd)),...
            '_seed',num2str(sd(nsd)),'.txt'],'w');
%         fprintf(filePh,'%s\n',newStr2{:});
        fprintf(filePh,'%s\n',newStr3{:});
        fclose(filePh);
        %change file extention
        directory = [filepath,template(1:end-4),'_',...
            num2str(wd(nwd)),'_seed',...
            num2str(sd(nsd))];    
        file = dir(directory); 
        oldName = cell(length(file)-2, 1);
        for ii = 3:length(file)
           oldName{ii-2} = file(ii).name;
        end
        newname = [template(1:end-4),'_',...
            num2str(wd(nwd)),'_seed',num2str(sd(nsd)),'.$PJ']; 
        movefile([directory '\' oldName{1}],...
            [directory '\' newname]);
    end
end