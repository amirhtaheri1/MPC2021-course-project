clc
clear all
close all
format bank
tic

plot_flag = 1;
curtailment = 0;
global PBpb
PBpb=0;
up_sample_ratio = 4
Data_Ts_120 = xlsread('Data_from_paper.xls');
Psunny = interp(Data_Ts_120(:,2),up_sample_ratio);
Pcloudy = interp(Data_Ts_120(:,3),up_sample_ratio);
Pwind = interp(Data_Ts_120(:,4),up_sample_ratio);
Load_up = interp(Data_Ts_120(:,5),up_sample_ratio);
Load_curtailed = interp(Data_Ts_120(:,6),up_sample_ratio);

if curtailment
   Load_up = Load_curtailed; 
end
% % testing (down sampling for faster simulation)
% down_sample_ratio = 5
% Psunny = downsample(Data_Ts_120(:,2),down_sample_ratio);
% Pcloudy = downsample(Data_Ts_120(:,3),down_sample_ratio);
% Pwind = downsample(Data_Ts_120(:,4),down_sample_ratio);
% Load_up = downsample(Data_Ts_120(:,5),down_sample_ratio);
% Load_curtailed = downsample(Data_Ts_120(:,6),down_sample_ratio);

% tariff
mid=214.7/120000;
peak=215.7/120000;
base=205.5/120000;

tariff=zeros(1,2880);

tariff(1:719)=mid;
tariff(720:959)=peak;
tariff(960:2039)=base;
tariff(2040:2519)=peak;
tariff(2520:2880)=mid;

%% Simulation details
if plot_flag
    % Load
    figure; xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,Load_up,'k','LineWidth',2); xlim([0 24]); ylim([0 1.3]); xticks([0:2:24]);
    xlabel('Time (h)'); ylabel('Power (kW)');
    
    % Curtailed load
    figure; xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,Load_up,'k','LineWidth',2); xlim([0 24]); ylim([0 1.3]); xticks([0:2:24]);
    hold on; plot(xaxis,Load_curtailed,'LineWidth',2);
    xlabel('Time (h)'); ylabel('Power (kW)');
    legend('Load','Curtailed load','Location','northwest');
    
    figure;
    subplot(2,1,1);
    xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,Psunny,'LineWidth',2); xlim([0 24]); ylim([0 2]); xticks([0:2:24]);
    xlabel('Time (h)'); ylabel('Power (kW)');
    legend('Psunny','Location','northwest');
    
    subplot(2,1,2);
    xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,Pcloudy,'r','LineWidth',2); xlim([0 24]); ylim([0 1.5]); xticks([0:2:24]);
    xlabel('Time (h)'); ylabel('Power (kW)');
    legend('Pcloudy','Location','northwest');
    
    figure; xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,Pwind,'g','LineWidth',2); xlim([0 24]); ylim([0 1.3]); xticks([0:2:24]);
    xlabel('Time (h)'); ylabel('Power (kW)');
    
    figure; xaxis=linspace(0,24,length(Load_up));
    plot(xaxis,12e4*tariff,'m','LineWidth',2); xlim([0 24]); ylim([0 250]); xticks([0:2:24]);
    xlabel('Time (h)'); ylabel('Tariff (R$/MW)');
    ylim([])
return;
end

%% Create MPC object
nlobj = nlmpc(2,2,'MV',[1,2,3],'MD',[4,5]);
Ts=30;
nlobj.Ts = 30;
nlobj.PredictionHorizon = 5;
nlobj.ControlHorizon = 5;

nlobj.Model.StateFcn = 'SOC_state';
nlobj.Model.IsContinuousTime = false;
nlobj.Model.OutputFcn = @(x,u,Ts) [x(1);x(2)];
nlobj.Model.NumberOfParameters =1;

x0 = [0; 0];
u0 = [0; 0; 0; 0; 0];
EKF = extendedKalmanFilter(@SOC_state,@SOC_output);
EKF.State = x0;
validateFcns(nlobj,x0,u0(1:3)',u0(4:5)',{Ts});

xref = x0';

nlobj.ManipulatedVariables(1).Min = -6;
nlobj.ManipulatedVariables(2).Min = -3.072;
nlobj.ManipulatedVariables(3).Min = 0.06;

nlobj.ManipulatedVariables(1).Max = 75;
nlobj.ManipulatedVariables(2).Max = 3.072;
nlobj.ManipulatedVariables(3).Max = 1.5;

nlobj.ManipulatedVariables(1).RateMin = -5;
nlobj.ManipulatedVariables(2).RateMin = -3.072;
nlobj.ManipulatedVariables(3).RateMin = -0.06;

nlobj.ManipulatedVariables(1).RateMax = 5;
nlobj.ManipulatedVariables(2).RateMax = 3.072;
nlobj.ManipulatedVariables(3).RateMax = 1.5;

nlobj.States(1).Min = -20;
nlobj.States(2).Min = -40;

nlobj.States(1).Max = 15;
nlobj.States(2).Max = 35;


% nlobj.MeasuredDisturbances.ScaleFactor=1; 

nlobj.Optimization.CustomCostFcn = @CostFunc;
% nlobj.Optimization.ReplaceStandardCost = true;

% nlobj.Model.NumberOfParameters = 1;
Nt=length(Load_up);
xHistory = x0;
x=x0;
uk = u0;
xk=x;
% p_renw=zeros(1,Nt);
p_renw=Pwind;
rng default
p_noise=0.2*(rand(1,Nt)-0.5);
Balance = zeros(1,Nt);
nloptions = nlmpcmoveopt;
PBpb_vector = zeros(1,Nt);
saveControl=zeros(5,Nt);
saveState=zeros(2,Nt);
for i = 1:Nt
    i
    saveControl(:,i) = uk;
%     uk
    saveState(:,i) = xk;
    nloptions.Parameters = {{Load_up(i),p_renw(i),tariff(i),saveControl(:,i),saveState(:,i)}};
    % Correct previous prediction
    xk = correct(EKF,x);
    [uk(1:3),info,info2] = nlmpcmove(nlobj,xk,uk(1:3),xref,uk(4:5)',nloptions);
    % uK=[uk;uk(4)]
    % Predict state for next iteration
    uk(4)=p_renw(i);%+p_noise(i);
    uk(5)=p_noise(i);
    saveControl(:,i) = uk;
    PBpb_vector(i) = uk(1)-uk(2)-uk(3)+uk(4);
%     PBpb_vector(i)=PBpb;
    
    % Energy balance testing
    Balance(i)=saveControl(1,i)-saveControl(2,i)-saveControl(3,i)+saveControl(4,i)-PBpb_vector(i);
    
    predict(EKF,uk,Ts);
    % Implement optimal control move
    x = SOC_state(xk,uk);
    y = x;    
    xHistory = [xHistory x];
%     i/Nt*100
end
%SOC
figure
% subplot(2,1,1)
plot(xHistory(1,:)+60,'LineWidth',2)
xlabel('Sample')
% ylabel('SOC Pb')
% title('SOC Pb')
% subplot(2,1,2)
hold on
plot(xHistory(2,:)+60,'LineWidth',2)
ylim([20 90])
% xlabel('time')
% ylabel('SOC Li')
% title('SOC Li')
legend('SOC Pb','SOC Li');
xlim([1 Nt]);

%Control
figure()
plot(saveControl(1,:),'LineWidth',2)
hold on
plot(saveControl(2,:),'LineWidth',2)
hold on
plot(saveControl(3,:),'LineWidth',2)
hold on
plot(saveControl(4,:),'LineWidth',2)
hold on
plot(PBpb_vector,'LineWidth',2)
hold on
plot(Balance,'LineWidth',2)
legend('Pgird','PBli','Pload','Pgen','Ppb','Balance')
xlim([1 Nt]);
figure()
plot(saveControl(3,:))
hold on
plot(Load_up)
xlabel('Sample');
xlim([1 Nt]);
toc
format