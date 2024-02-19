clear

filepath = 'C:\loong\BLADED WORKPLACE\1.OC3_Monopile\0.pool2008\';
% template = 'und.$PJ';
template = 'wind.$PJ';
expression1 = '101.wnd';
expression2 = 'SEED	101';
expression3 = 'UBAR	 11.4';

wd = [3 11];
NWD = length(wd);

for nwd = 1:NWD
    template = ['wind',num2str(wd(nwd)),'.$PJ'];
    sd = (1:10)+floor((wd(nwd)-3)/2)*10+400;
    NSD = length(sd);
    for nsd = 1:NSD
        % creat folder
        mkdir([filepath,template(1:end-4),...
            '_seed',num2str(sd(nsd))]);
        
        %read
        Str = strsplit(fileread(template), '\n');
    %     fclose(template);
        %edit
        replace1 = [template(1:end-4),'_seed',num2str(sd(nsd)),'.wnd'];
        newStr1 = regexprep(Str,expression1,replace1);
        replace2 = ['SEED	',num2str(sd(nsd))];
        newStr2 = regexprep(newStr1,expression2,replace2);
        %write
        filePh = fopen([filepath,template(1:end-4),...
            '_seed',num2str(sd(nsd)),...    
            '\',template(1:end-4),...
            '_seed',num2str(sd(nsd)),'.txt'],'w');
        fprintf(filePh,'%s\n',newStr2{:});
        fclose(filePh);
        %change file extention
        directory = [filepath,template(1:end-4),...
            '_seed',...
            num2str(sd(nsd))];    
        file = dir(directory); 
        oldName = cell(length(file)-2, 1);
        for ii = 3:length(file)
           oldName{ii-2} = file(ii).name;
        end
        newname = [template(1:end-4),...
            '_seed',num2str(sd(nsd)),'.$PJ']; 
        movefile([directory '\' oldName{1}],...
            [directory '\' newname]);
    end
end