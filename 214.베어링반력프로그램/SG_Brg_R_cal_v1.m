% - 화면을 보고 사용할 데이터 범위 선정 필요

close all; clear all; clc
global t y CO 

%       1  2  3  4  5  6  7  8  9
%  Time C1 C2 Ta A1 A2 A3 B1 B2 B3.    % PM_CCW_Launch_Final
f = 128; % Sampling frequency

%% Import the raw data ( by Manual )
%  -- time 제외, Tacho. 포함(마지막 ch)
load Raw_data_CCW_PM_Launch_final.mat;
data_all=Raw_data_CCW_PM_Launch_final(:,2:end); % time removal
t=Raw_data_CCW_PM_Launch_final(:,1);

% ----- start time set ----->  Need to revise 
data_all=data_all(0.5*10^5:end,:); 
t=Raw_data_CCW_PM_Launch_final(0.5*10^5:end,1);

ch=size(data_all,2);

%% 채널별 계측 시간이 다른경우 데이터 수 일치
%%% ------ NAN data removal
NAN_sum=sum(isnan(data_all),2);

for i=1:length(data_all) % find NAN start no.
    k=i;
    if NAN_sum(i)==0
        break
    end
end

for i=length(data_all):-1:1  % find NAN end no.
    kk=i;
    if NAN_sum(i)==0
        break
    end
end

data_all_cut=data_all(k:kk,:);
data_all_cut=data_all;

data = data_all_cut;

%% 데이터 트림

[m,n]=size(data);

data_nf=data(1,1:8);
data_normal=data;
data_normal(:,1:8)=data(:,1:8)-data_nf(:,1:8); %타코 데이터 추출
t=t-t(1,1); %노말라이징
t_max=round(max(t(:,1)));
idx_t=floor(m/t_max); %T최댓값에 따른 샘플링
%변화량(diff)
for idx=t_max:t_max:m-t_max
    slp(idx/t_max,1:8)=data_normal(idx+t_max,1:8)-data_normal(idx,1:8);
end


for slp_i=1:8
    slp_lim(:,slp_i)={find(abs(slp(:,slp_i))>=max(abs(slp(:,slp_i)))/10)};
end

slp_min=cellfun(@min, slp_lim);
slp_max=cellfun(@max, slp_lim);

%데이터 변화 구간 최대/최소
data_min=max(slp_min)*t_max;
data_max=min(slp_max)*t_max;

%데이터 트림
for i=1:9
data_trim(:,i)=data_normal([data_min:data_max],i);
end

%data_sm=smoothdata(data_trim,'sgolay');
%% --- fill up empty data with previous data
% -- filter, mean, median 계산시 NAN 제거 필요

size_data = size(data_trim,1);
count=0;
for i=1:size(data_trim,2)
    for j=1:size_data
        if isnan(data_trim(j,i)) == 1
            data_trim(j,i) = data_trim(j-1,i);
            count=count+1;
        end
    end
end

% ---- Tacho amp normalization  &  make positive 
Tacho=data_trim(:,end);
% Tacho=Tacho * -1 - 1000000;
Tacho_med=median(Tacho(:,end));
Tacho_amp=max(Tacho) - min(Tacho);
Tacho=Tacho-Tacho_med;
Tacho_factor=10;
Tacho=Tacho/Tacho_amp * Tacho_factor;
if mean(Tacho) < 0
    Tacho = -1 * Tacho; % Normalized Tacho 0 ~ 1
end

% Tacho_med=median(data(:,end))
% Tacho_amp=max(data(:,end)) - min(data(:,end))
% Tacho=data(:,end);
% Tacho=Tacho-Tacho_med;
% Tacho=Tacho/Tacho_amp
% data(:,end) = data(:,end) / Tacho_amp; 
% Tacho_amp=10 - -1000
% data(:,end) = data(:,end) / 4000;  % Tacho signal mag. normalize
data_trim(:,end)=Tacho;
figure(100);
plot(data_trim);
axis tight;

%% Low pass filtering & cutting unnecessary data
% -- filter setting
n = 5;            % Order
Wn = 0.005;       % Cutoff freq (2000 Hz)
ftype = 'low';    % Low pass
[b,a] = butter(n,Wn,ftype);

% -- filter process & filter front data removal
data_fil=filter(b,a,data_trim(:,1:end-1)); % Tacho 제외
t=t(1:size_data,1);
figure(103);plot(t,data_trim,t,data_fil);

x_shift=1000; % Filtered signal shift to x-dir. 
data_fil=data_fil(x_shift:end,:); % 앞으로 이동한 만큼 뒤쪽 데이터 감소
data_fil(:,end+1)=data_trim(x_shift:end,end);  % Tacho 포함
t_fil=t(1:end-x_shift+1);
tacho_fil=Tacho(1:end-x_shift+1);
figure(104);
plot(t,data_trim,t_fil,data_fil);
axis tight;
 
% Tach=data(:,end);
% Tach_fil=Tach(x_shift:end);
% figure(104);plot(t,Tach,t_fil,Tach_fil)
% data_fil(:,size(data_fil,2)+1)=Tach_fil;

%% Find the start & end points
for i = 10:floor(size(data_fil,1)/f) - 1 % convert into sec unit
    p = (i-1) * f + 1;
    q = i * f;
    for j=1:ch
        mean_table(i,j) = mean(data_fil(p:q,j));
    end
    mean_table(i,ch+1) = mean(abs(mean_table(i,1:ch)));
end

% data parsing
mean_diff = diff(mean_table(:,6));
temp_s = []; temp_e = []; k_s = 0; k_e = 0;
for i = 1:size(mean_diff)
    if abs(mean_diff(i,1)) > 0.1 %                                                              %%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.1 to ?
        k_s = k_s + 1;
        temp_s(k_s,1) = i;
    end
    if abs(mean_diff(i,1)) > 0.02 %                                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.02 to ?
        k_e = k_e + 1;
        temp_e(k_e,1) = i;
    end
end
point_s = (temp_s(2) + 1) * f;   % Start point                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%% + 1 to ?
point_e = (temp_e(end) - 5) * f; % End point                                               %%%%%%%%%%%%%%%%%%%%%%%%%%%% - 1 to ?

data_fil_cut = data_fil(point_s:point_e,:);
t_fil_cut=t_fil(point_s:point_e,:);
tacho_fil_cut=tacho_fil(point_s:point_e,:);
figure(200);
plot(t_fil_cut,data_fil_cut);
axis tight;

%% Downsampling (128 Hz to 128/dwn Hz)
dwn = 60; %                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%% 60 to ?
data_fil_cut_dwn = [];

for i = 1:floor(size(data_fil_cut)/dwn)
    data_fil_cut_dwn(i,:) = data_fil_cut((i-1)*dwn+1,:);
    t_fil_cut_dwn(i,1)=t_fil((i-1)*dwn+1);
    tacho_fil_cut_dwn(i,1)=tacho_fil((i-1)*dwn+1);
end

figure(300);
plot(t_fil_cut_dwn,data_fil_cut_dwn);
axis tight;

%% Sine Fitting
x=t_fil_cut_dwn(:,1);
y=data_fil_cut_dwn; %Sine + noise
[d_m,d_n]=size(t_fil_cut_dwn);
xstart=x(1);
xend=x(end);
xLength = d_m;
tic;

for ch=1:8
    [SineP]=sineFit(x,y(:,ch),ch);
    SinP(:,ch)=SineP;
    toc;
    xSstep=xend/min(xLength/100,1/SinP(3,ch)*0.1);
    xS=xstart:xSstep:xend;
    xFFT=xstart-xLength*.02:xLength/100:xend+xLength*.02;%FFTcurve longer
    y_f_r(:,ch)=SinP(1,ch)+SinP(2,ch)*sin(2*pi*SinP(3,ch)*x+SinP(4,ch));%result
end

data=horzcat(y_f_r,tacho_fil_cut_dwn);
figure(400);
plot(x,data);
axis tight;

%% Strain calculation at Tacho.
% Tacho(Peak) detection
peak = []; p = 0;
for i=1:size(y_f_r)-1
    if y_f_r(i,end) > Tacho_factor*.2 %                                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%% > 20 to ?
        p = p + 1;
        peak(p,1) = i;
    end
end
peak(size(peak,1)+1,1) = peak(size(peak,1),1) + 10000;

P_seg = []; seg = 0; col = 1;
for i=1:size(peak,1)-1
    if peak(i,1)+1 == peak(i+1,1)
        seg = seg + 1;
        P_seg(seg,col) = peak(i,1);
    end
    if peak(i,1)+1 < peak(i+1,1)
        seg = seg + 1;
        P_seg(seg,col) = peak(i,1);
        seg = 0;
        col = col + 1;
    end
    if i == size(peak)
        P_seg(seg,col) = peak(i,1);
    end
end
P_no = [];
for i=1:size(P_seg,2)
    P_no(1,i) = nnz(P_seg(:,i)); % count the number of nonzero component
end

% peak_calibration
P_st=P_seg(1,:);
P_del=diff(P_st);
P_del_avg=mean(P_del);

