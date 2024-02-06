%% A ca2+ puff model based on integrodifferential equations (Hawker et al., 2024)
% Code is based on the ca2+ puff model by Cao et al., (2013,2014) and used the
% method by Brady (1972) to write the ODEs for the h42 gating variable in
% the Cao et al. model as integrodifferential equations

% Reduced 6-state model 
clear all
format long

c0=100; % resting calcium concentration

Jrelease=200000;  % calcium released via a single IPR channel

% constant IPR parameters
q45=11;q54=3330;q12=1240;q21=88;q23=3;q32=69;q26=10500;q62=4010;

Vs=Jrelease*20;
Ks=12000;
% sl=500;  % step length

Jleak=Vs*100/(100+Ks);

% dye buffer
BT=20000; % total calcium buffer concentration (nM)
k_on=0.15;   % on-rate or binding rate (/nM/s)
Kd=2000;   % dissociation constant (nM)
k_off=Kd*k_on;  % off-rate or dissociation rate
B=k_off*BT/(k_on*c0+k_off);  % resting free buffer concentration
B_new=B;

Y=[c0,B];
pars=[Jrelease, 0, Jleak, BT, Vs, Ks, k_on, k_off];

Ca_model(0, Y,pars);

% construct vectors and set initial conditions

c=c0;   % initial calcium concentration
Num_IPR=10;  % number of IPRs
state=4*ones(Num_IPR,1); % number indicates the state in order (C1 C2 C3 C4 O5 O6)

dt0=1e-4; % initial time stepsize in second
dt=dt0;

% Delay
tau=3; % second(s)
nPast=tau/dt0; % no. elements 

cm_nPast=c(end)*ones(Num_IPR,nPast);
cm=c(end)*ones(Num_IPR,1);

time=0; % time series
time_c=0:dt0:tau; % continuous time series
time_c(end)=[];

g=zeros(Num_IPR,1);
r1=rand(Num_IPR,1);
zc=zeros(Num_IPR,1);
%%

IP3=0.1; % muM 

V24=60+253*1.2^3./(IP3.^3+1.2^3);
k24=479+70*1.2^2./(IP3.^2+1.2^2);
n24=6.3+1.72*IP3.^2./(IP3.^2+1.2^2);
kn24=79750+17368*1.2.^2./(IP3.^2+1.2^2);
nn24=8.2*IP3.^2./(IP3.^2+1.5^2);
a24=1+30*0.5.^2./(IP3.^2+0.5^2);

V42=100;
k42=398+257*IP3.^4./(IP3.^4+3.6^4);
n42=5.9+5.3*1.2.^2./(IP3.^2+1.2^2);
kn42=170+70000*IP3.^3./(IP3.^3+6.5^3);
nn42=3.2+4.88*IP3.^2./(IP3.^2+1.3^2);
a42=1.8*IP3.^2./(IP3.^2+0.58^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h24_inf=kn24^nn24./(kn24^nn24+c0.^nn24);
h24=h24_inf*ones(Num_IPR,1);
h24_new=h24;
h24_track=h24;

ah24=40;
Vh24=0; 
Kh24=10000;

m24_inf=c0.^n24./(k24^n24+c0.^n24);
m24=m24_inf*ones(Num_IPR,1);
m24_new=m24;
m24_track=m24;

am24=100;
Vm24=0; 
Km24=10000;

m42_inf=c0.^n42./(k42.^n42+c0.^n42);
m42=m42_inf*ones(Num_IPR,1);
m42_new=m42;
m42_track=m42;

am42=100;
Vm42=0;
Km42=10000;

h42_inf=kn42^nn42./(kn42^nn42+c0.^nn42);
h42=h42_inf*ones(Num_IPR,1);
h42_new=h42;
h42_track=h42;

ah42=0.5;
Vh42=100;
Kh42=20000;

%%
h42r=ah42+Vh42*(cm_nPast.^7)./(Km42^7+cm_nPast.^7);
h42gInf=kn42^nn42./(kn42^nn42 + c0.^nn42);
alpha_h42=h42r.*h42gInf;
%%
IT=1; % incremental time per iteration (s)
Numtimes=10; % number of iterations
dis=100;


for full=1:Numtimes

tic

Tmax=IT*full;

while time(end)<Tmax
        if abs(imag(dt))>0 || real(dt)<0
            break
        end

        No=0;
        for openind=1:Num_IPR
            if state(openind,end)==5 || state(openind,end)==6
                No=No+1;
            end
        end

        Y0=[c(end), B(end)];

        % here use RK4 to solve the ODE
        pars(2)=No;
        Y_new=RK4(@(t,Y) Ca_model(t,Y, pars), time, Y0, dt);

        c_new=Y_new(1);
        B_new=Y_new(2);

        dt1=dt*ones(Num_IPR,1);

    for i=1:Num_IPR

        cm(i)=c_new+120000*heaviside(state(i,end)-4.5);

        %% Gating variable calculation

        % Steady state values
        m42_inf=cm(i).^n42./(k42.^n42+cm(i).^n42);
        m42_new(i)=m42_inf;

        m24_inf=cm(i).^n24./(k24^n24+cm(i).^n24);
        m24_new(i)=m24_inf;

        h24_inf=kn24^nn24./(kn24^nn24+cm(i).^nn24);
        h24_new(i)=h24_inf;

        % Brady (1972) integrodifferential equation
        h42r(i,end)=ah42+Vh42*(cm(i).^7)./(Kh42^7+cm(i).^7); 
        h42gInf(i)=kn42^nn42./(kn42^nn42 + cm(i).^nn42);
        alpha_h42(i,end)=h42r(i,end)*h42gInf(i);
        h42_new(i)=gatingSolutionMH(h42r(i,end-nPast+1:dis:end),alpha_h42(i,end-nPast+1:dis:end),diff(time_c(end-nPast+1:dis:end)),h42_inf);

        %% calculating q24 and q42 rates

        q24=a24+V24*(1-m24(i).*h24(i));
        q42=a42+V42*m42(i).*h42(i);

        Ta_old=[0 q12 0 0 0 0 ;q21 0 q23 q24 0 q26;0 q32 0 0 0 0;0 q42 0 0 q45 0;...
            0 0 0 q54 0 0;0 q62 0 0 0 0];

        q24=a24+V24*(1-m24_new(i).*h24_new(i));
        q42=a42+V42*m42_new(i).*h42_new(i);

        Ta_new=[0 q12 0 0 0 0 ;q21 0 q23 q24 0 q26;0 q32 0 0 0 0;0 q42 0 0 q45 0;...
            0 0 0 q54 0 0;0 q62 0 0 0 0];

        g_old=g(i);
        g_new=g_old+(sum(Ta_old(state(i,end),:))+sum(Ta_new(state(i,end),:)))/2*dt;
        epsilon=log(1/r1(i));

        if g_new>=epsilon  % transition occurs
            dt1(i)=(epsilon-g_old)/(g_new-g_old)*dt;
        end
    end
    [dt_min,index]=min(abs(dt1));
    if index==1 && dt_min/dt==1
        index=0;
    end
    dt=dt_min;

    % here use RK4 to solve the ODE
    pars(2)=No;
    Y_new=RK4(@(t,Y) Ca_model(t,Y, pars), time, [c(end), B(end)], dt);

    c_new=Y_new(1);
    B_new=Y_new(2);

    %% calculating calcium concentration and gating variables
    for i=1:Num_IPR

        cm(i)=c_new+120000*heaviside(state(i,end)-4.5);
        %% Gating variable calculation

        % Steady state values
        m42_inf=cm(i).^n42./(k42.^n42+cm(i).^n42);
        m42_new(i)=m42_inf;

        m24_inf=cm(i).^n24./(k24^n24+cm(i).^n24);
        m24_new(i)=m24_inf;

        h24_inf=kn24^nn24./(kn24^nn24+cm(i).^nn24);
        h24_new(i)=h24_inf;

        % Brady (1972) integrodifferential equation
        h42r(i,end)=ah42+Vh42*(cm(i).^7)./(Kh42^7+cm(i).^7); 
        h42gInf(i)=kn42^nn42./(kn42^nn42 + cm(i).^nn42);
        alpha_h42(i,end)=h42r(i,end)*h42gInf(i);
        h42_new(i)=gatingSolutionMH(h42r(i,end-nPast+1:dis:end),alpha_h42(i,end-nPast+1:dis:end),diff(time_c(end-nPast+1:dis:end)),h42_inf);

        %% calculating q24 and q42 rates
        q24=a24+V24*(1-m24(i).*h24(i));
        q42=a42+V42*m42(i).*h42(i);

        Ta_old=[0 q12 0 0 0 0;q21 0 q23 q24 0 q26;0 q32 0 0 0 0;0 q42 0 0 q45 0;...
            0 0 0 q54 0 0;0 q62 0 0 0 0];

        q24=a24+V24*(1-m24_new(i).*h24_new(i));
        q42=a42+V42*m42_new(i).*h42_new(i);

        Ta_new=[0 q12 0 0 0 0;q21 0 q23 q24 0 q26;0 q32 0 0 0 0;0 q42 0 0 q45 0;...
            0 0 0 q54 0 0;0 q62 0 0 0 0];

        g(i)=g(i)+(sum(Ta_old(state(i,end),:))+sum(Ta_new(state(i,end),:)))/2*dt;

    end

    c=[c,c_new];
    time=[time,time(end)+dt];
    time_c=[time_c, time_c(end)+dt];
 
    %% tracking the gating variable solutions
    h24=h24_new;
    h24_track=[h24_track,h24_new];

    m24=m24_new;
    m24_track=[m24_track,m24_new];

    h42=h42_new;
    h42_track=[h42_track,h42_new];
   
    m42=m42_new;
    m42_track=[m42_track,m42_new];

    B=[B,B_new];

    if index==0
        state=[state,state(:,end)];
        dt=dt0;
    else

        q24=a24+V24*(1-m24(index).*h24(index));
        q42=a42+V42*m42(index).*h42(index);

        Ta_new=[0 q12 0 0 0 0 ;q21 0 q23 q24 0 q26;0 q32 0 0 0 0;0 q42 0 0 q45 0;...
            0 0 0 q54 0 0;0 q62 0 0 0 0];
        r2=rand;
        previous=Ta_new(state(index,end),:);
        kk=1;
        while sum(previous(1:kk))/sum(previous)<r2
            kk=kk+1;
        end
        current=state(:,end);
        current(index)=kk;
        state=[state,current];
        r1(index)=rand;
        g(index)=0;

    end


end

toc

save(['puff_10IPRs_IP3_01uM_Brady_reduced6s_',num2str(full)])

time=time(end); % time series
c=c(end);   % resting calcium concentration
state=state(:,end); % number indicates the state in order (C1 C2 C3 C4 O5 O6)
h42_track=h42_track(:,end);
m42_track=m42_track(:,end);
h24_track=h24_track(:,end);
m24_track=m24_track(:,end);
B=B(end);
% 
h42r=h42r(:,end-nPast+1:end);
alpha_h42=alpha_h42(:,end-nPast+1:end);

time_c=time_c(end-nPast+1:end);
end


%%
totaltime=[];
Ca=[];
h42_total=[];
m42_total=[];
h24_total=[];
m24_total=[];
state_total=[];
cm_total=[];
fluo=[];

for ldf=1:Numtimes

load(['puff_10IPRs_IP3_01uM_Brady_reduced6s_',num2str(ldf)])

totaltime=[totaltime,time];
Ca=[Ca,c];

h42_total=[h42_total,h42_track];
m42_total=[m42_total,m42_track];
h24_total=[h24_total,h24_track];
m24_total=[m24_total,m24_track];
state_total=[state_total,state];
fluo=[fluo,B];
end

save('puff_10IPRs_IP3_01uM_Brady_reduced6s')









