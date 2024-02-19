t = linspace(0,1210,1210*30)';
x = chirp(TimeStep,0.2,1200,0.3,[],90)*1e5; 


% information
% get(0,'ScreenSize') and = [1 1 1536 864]
% screen physical size 20.75 inch * 11.67 inch
% screen: 74 pixels per inch; 29.1 pixels per centimeter
% 3 columns 6cm; 2 columns 9cm; 1 column 18cm; 1.5 column 14cm

figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.78]*29.1);
p = plot(t(1:100*30),x(1:100*30));
xlabel('Time (s)');
ylabel('Amplitude (m)')
% p.LineWidth = 1.5;
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'input_chirp_1e5.eps','Resolution',1000)

%% gwn
t = linspace(0,1210,1210*30)';
seed = 101;
x_gwn = wgn(1210*30,1,0,1,seed)*1e6; 
figurewidth = 9; %cm
f = figure('Position',[10 10 figurewidth figurewidth*0.78]*29.1);
p = plot(t(1:100*30),x_gwn(1:100*30));
xlabel('Time (s)');
ylabel('Amplitude (m)')
% p.LineWidth = 1.5;
set(findall(f,'-property','FontSize'),'FontSize',7)
exportgraphics(f,'input_gwn_1e6.eps','Resolution',1000)

