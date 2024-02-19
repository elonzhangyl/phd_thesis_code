% ref:Damage detection in uncertain nonlinear systems based on
% stochastic Volterra series
clear
clc
m = 0.26;c = 1.36;k1 = 5.49e3;k2 = 3.24e4;k3 = 4.68e7;
amp = 1; 
alpha = 0.98;
fs = 512;
fs_noise = fs;
tstart = 0;
tfinal = 20;
u0 = [0; 0];
refine = 4;

% Options for ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@events);

% pt = @(t) chirp(t,15,tfinal,30);
% pt = @(t) cos(2*pi*5*t);
dataset2 = zeros(fs*(tfinal-10),10);
seed_list = 99001:99120;
for sd = 1:length(seed_list)
    seed = seed_list(sd)

    noise0 = wgn(tfinal*fs_noise+1,1,0,1,seed);
    noise = lowpass(noise0,2,fs_noise,'Steepness',0.99);
%     noise = noise0;
%     bandpassSpecs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
%         10,15,30,35,60,0.5,60,fs);
%     bandpassFilt  = design(bandpassSpecs);
%     noise = filter(bandpassFilt,noise0);
    noise(1) = 0;
%     pwelch(noise,[],[],[],fs_noise)
    t_noise = tstart:1/fs:tfinal;

odefun1 = @(t,u,m,c,k1,k2,k3,amp,t_noise,noise) ...
    [u(2); (-c*u(2)-k1*u(1)-k2*u(1)^2-k3*u(1)^3+...
    amp*interp1(t_noise,noise,t))/m];
odefun2 = @(t,u,m,c,k1,k2,k3,amp,alpha,t_noise,noise) ...
    [u(2); (-c*u(2)-alpha*k1*u(1)-k2*u(1)^2-k3*u(1)^3+...
    amp*interp1(t_noise,noise,t))/m];

% Accumulators for ODE solve output
tout = tstart;  % global solution time
uout = u0';     % global solution position
teout = [];     % time at which events occur
ueout = [];     % position at which events occur
ieout = [];     % flags which event trigger the switch

% Choose starting equation
if u0(1) > 0
    flag = 0;
else
    flag = 1;
end
% Main loop
tout = 0;
while tout(end) < tfinal
%     tout(end)
    tspan = [tout(end) tfinal];
    if length(tspan)<2
        break
    end
    if flag == 0
        [t, u] = ode45(@(t,u) ...
            odefun1(t,u,m,c,k1,k2,k3,amp,t_noise,noise),...
            tspan,u0,options);
        flag = 1;
    else
        [t, u] = ode45(@(t,u) ...
            odefun2(t,u,m,c,k1,k2,k3,amp,alpha,t_noise,noise),...
            tspan,u0,options);
        flag = 0;
    end
    
        nt = length(t);
        tout = [tout; t(2:nt)];
        uout = [uout; u(2:nt,:)];
%         t1eout = [t1eout; te];% Events at tstart are never reported.
%         u1eout = [u1eout; ue];
%         i1eout = [i1eout; ie];
        u0 = [u(nt,1); u(nt,2)];
%         options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
%           'MaxStep',t(nt)-t(1));

end

    tout_inter = tstart:1/fs:tfinal-1/fs;
    uout_inter = interp1(tout,uout,tout_inter);
%     pwelch(uout_inter(:,1),[],[],[],fs)
    
    ddu = zeros(size(uout));%acceleration
    for i = 1:length(tout)
        if uout(i,1) > 0
            ddu(i,:) = odefun1(tout(i),uout(i,:),...
                m,c,k1,k2,k3,amp,t_noise,noise);
        else
            ddu(i,:) = odefun2(tout(i),uout(i,:),...
                m,c,k1,k2,k3,amp,alpha,t_noise,noise);
        end
    end
    
    ddu_inter = interp1(tout,ddu,tout_inter);

    dataset2(:,sd) = ddu_inter(fs*10:end-1,2);
%     plot(tout,uout(:,2));grid on
%     plot(tout,ddu(:,2),tout_inter,ddu_inter(:,2));grid on
%     pwelch(dataset2(1:end,2),[],[],[],fs)
end

save(['duffingCrack_dam',num2str(alpha),'_amp',num2str(amp),...
    '_seed',num2str(min(seed_list)),'_',num2str(max(seed_list)),'.mat'],...
    'dataset2');

function [value,isterminal,direction] = events(t,u)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = u(1);     % detect velocitu = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
end


