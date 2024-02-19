clear

filepath = 'C:\loong\BLADED WORKPLACE\1.OC3_Monopile\11.windfarm\2.pool.und\';
expression1 = 'wind3_seed1\\wind3_seed1.wnd';
% expression1 = '.wnd';
% expression2 = 'IDUM	101';

wd = [11];
NWD = length(wd);
name = 'water23m';

for nwd = 1:NWD
%     template = ['und_wind',num2str(wd(nwd)),'.$PJ'];
    template = [name,'_wind',num2str(wd(nwd)),'.$PJ'];
    sd = (1:10)+floor((wd(nwd)-3)/2)*10;
    NSD = length(sd);
    for nsd = 1:NSD
        % creat folder
        mkdir([filepath,template(1:end-4),...
            '_seed',num2str(sd(nsd))]);

        %read
        Str = strsplit(fileread(template), '\n');
        replace1 = ['wind',num2str(wd(nwd)),'_seed',num2str(sd(nsd)),...
            '\\wind',num2str(wd(nwd)),'_seed',num2str(sd(nsd)),'.wnd'];
        newStr1 = regexprep(Str,expression1,replace1);
%         replace2 = ['IDUM	',num2str(sd(nsd))];
%         newStr2 = regexprep(newStr1,expression2,replace2);
        %write
        filePh = fopen([filepath,template(1:end-4),'_seed',num2str(sd(nsd)),...    
            '\',template(1:end-4),'_seed',num2str(sd(nsd)),'.txt'],'w');
        fprintf(filePh,'%s\n',newStr1{:});
        fclose(filePh);
        %change file extention
        directory = [filepath,template(1:end-4),'_seed',...
            num2str(sd(nsd))];    
        file = dir(directory); 
        oldName = cell(length(file)-2, 1);
        for ii = 3:length(file)
           oldName{ii-2} = file(ii).name;
        end
        newname = [template(1:end-4),'_seed',num2str(sd(nsd)),'.$PJ']; 
        movefile([directory '\' oldName{1}],...
            [directory '\' newname]);
    end
end

