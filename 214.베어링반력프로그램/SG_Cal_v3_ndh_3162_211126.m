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


%% --- fill up empty data with previous data
% -- filter, mean, median 계산시 NAN 제거 필요

size_data = size(data,1);
count=0;
for i=1:size(data,2)
    for j=1:size_data
        if isnan(data(j,i)) == 1
            data(j,i) = data(j-1,i);
            count=count+1;
        end
    end
end

% ---- Tacho amp normalization  &  make positive 
Tacho=data(:,end);
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
data(:,end)=Tacho;
figure(100); plot(data);

%% Low pass filtering & cutting unnecessary data
% -- filter setting
n = 5;            % Order
Wn = 0.005;       % Cutoff freq (2000 Hz)
ftype = 'low';    % Low pass
[b,a] = butter(n,Wn,ftype);

% -- filter process & filter front data removal
data_fil=filter(b,a,data(:,1:end-1)); % Tacho 제외
figure(103);plot(t,data,t,data_fil);

x_shift=1000; % Filtered signal shift to x-dir. 
data_fil=data_fil(x_shift:end,:); % 앞으로 이동한 만큼 뒤쪽 데이터 감소
data_fil(:,end+1)=data(x_shift:end,end);  % Tacho 포함
t_fil=t(1:end-x_shift+1);
figure(104);plot(t,data,t_fil,data_fil)
 
% Tach=data(:,end);
% Tach_fil=Tach(x_shift:end);
% figure(104);plot(t,Tach,t_fil,Tach_fil)
% data_fil(:,size(data_fil,2)+1)=Tach_fil;

%%
% Find the start & end points
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
% t_fil_cut=t_fil(point_s:point_e,:);
figure(200); plot(data_fil_cut);

%% Downsampling (128 Hz to 128/dwn Hz)
dwn = 60; %                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%% 60 to ?
data_fil_cut_dwn = [];

for i = 1:floor(size(data_fil_cut)/dwn)
    data_fil_cut_dwn(i,:) = data_fil_cut((i-1)*dwn+1,:);
    t_fil_cut_dwn(i,1)=t_fil((i-1)*dwn+1);
end

figure(300); plot(data_fil_cut_dwn);

%% Fitting the data into sine wave
N_iter = 20; %                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%% 60 to ?
SG_fil = data_fil_cut_dwn(:,1:end-1);
% t = data_fil_cut_dwn(:,1);
t = t_fil_cut_dwn;
N_ch = size(SG_fil,2);

% ------ normalization of SG signal
med=median(SG_fil);
for i=1:size(med,2)
    SG_fil(:,i)=SG_fil(:,i)-med(i);
end

temp_amp = 0.5*(max(SG_fil)-min(SG_fil));

for ch=1:N_ch;
    r=0;
    for i=1:N_iter
        % --------- Need to improve 3 ~ 5 initial value
        % x_iv=[20 .03 6 0.25 0.21]; % initial value for guessing
%         x_iv=[20 .03 1 0.25 0.21]; % initial value for guessing
%         x_iv=x_iv.*rand(1,length(x_iv));
        x_iv(1)=temp_amp(ch); % initial value for guessing
        x_iv(2)=0.003; % initial value for guessing
        x_iv(3:5)=rand(1,3); 
        
        y=SG_fil(:,ch);  % SG signal 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --  optimization
        
        % option=optimset('Display','iter','TolFun',1e-4, 'TolX', 1e-4);
        option=optimset('TolFun',1e-3, 'TolX', 1e-3);
        
        % xx=fminsearch(@Test_opt3, x_iv, option);
        xx=fminsearch(@Test_opt4, x_iv, option);
        
        Error(i)=CO;
        X_est(i,:)=xx;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ----- results plot
        if rem(i,10)==1
            r=r+1;
            t_iter=floor((N_iter-1)/10)+1;
            figure('Name',sprintf(' Fitting Ch %d (%d / %d)', ch,r,t_iter),'NumberTitle','off');
        end
        k=rem(i,10);
        if k==0;
            k=10;
        end
        f_est=xx(1)*sin(xx(2)*pi*t-xx(5))+xx(3)+xx(4)*t;
        subplot(5,2,k);plot(t,y,t,f_est,'.');grid on;
        if Error(i) < 0.02
            break
        end
    end
    
    %%%% -- minimum error result
    ind=find(Error==min(Error));
    % F_est=X_est(ind,1)*sin(X_est(ind,2)*pi*t-X_est(ind,5))+X_est(ind,3)+X_est(ind,4)*t;
    % figure('Name',sprintf('Fiting Result : Ch %d (Case : %d)', ch, ind),'NumberTitle','off');
    % plot(t,y,t,F_est,'.');grid on;
    
    %%% -- save Results
    % eval( ['R_X_est_ch' num2str(ch) ' = X_est;'] );
    % eval( ['R_Ind_ch' num2str(ch) ' = ind;'] );
    % eval( ['R_X_est_ind_ch' num2str(ch) ' = X_est(ind,:);'] );
    % eval( ['R_F_est_ch' num2str(ch) ' = F_est;'] );
    
    R_X_est(ch,:)=X_est(ind,:);
    X_ind(ch,:)=ind;
    R_X_error(ch,:)=Error(ind);
end

% F_est : sine fitted signal
for ch=1:N_ch;
    F_est(:,ch)=R_X_est(ch,1)*sin(R_X_est(ch,2)*pi*t-R_X_est(ch,5))+R_X_est(ch,3)+R_X_est(ch,4)*t;
    y=SG_fil(:,ch);
    figure('Name',sprintf('Fiting Result : Ch %d (Case : %d)', ch, X_ind(ch)),'NumberTitle','off'); 
    plot(t,y,t,F_est(:,ch),'.');grid on;
end

%% DeTrend
for ch=1:N_ch;
    SG_fil_det(:,ch)=SG_fil(:,ch)-R_X_est(ch,4)*t;
    y=SG_fil(:,ch);
    figure('Name',sprintf('Detrend Result : Ch %d', ch),'NumberTitle','off');
    plot(t,y,t,SG_fil_det(:,ch));legend('Original','Detrended');grid on;
end

%% Strain calculation at Tacho.
% Tacho(Peak) detection
peak = []; p = 0;
for i=1:size(SG_fil)
    if data_fil_cut_dwn(i,end) > Tacho_factor*.2 %                                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%% > 20 to ?
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
P_st=P_seg(1,:)
P_del=diff(P_st)
P_del_avg=mean(P_del)

%% Peak calibration

% P_calib=[142 416 416+271 958 1231 1504 1504+271 2046 2319]; % AM_CW
P_calib=P_st; % AM_CCW, PM_CCW
% P_calib=[[51 325 599 599+271 1140 1414 1688 1688+271 2229 2503 2776]];% AM_CCW_Launch

P_phase=[1 3 5 7 9];% AM_CW 
% P_phase=[2 4 6 8];% AM_CCW PMCCW

% read strain at tach.
ST_Tacho_raw=data_fil_cut_dwn(P_calib(P_phase),:); % data_fil_cut_dwn == SG_fil ;
ST_Tacho_det=SG_fil_det(P_calib(P_phase),:); %

ST_Tacho_raw_diff=diff(ST_Tacho_raw);
ST_Tacho_det_diff=diff(ST_Tacho_det);

ST_Tacho_raw_AVG=mean(abs(ST_Tacho_raw_diff));
ST_Tacho_det_AVG=mean(abs(ST_Tacho_det_diff));

%% Max-Min Strain

p_start1=1;p_end1=6;
temp_data1=data_fil_cut_dwn(P_calib(p_start1):P_calib(p_end1),:);
% temp_data11=SG_fil_det(peak_calib(p_start1):peak_calib(p_end1),:);
figure;
plot(temp_data1);
% hold on;plot(temp_data11,'r');
ST_Peak1(1,:)=max(data_fil_cut_dwn(P_calib(p_start1):P_calib(p_end1),:));
ST_Peak1(2,:)=min(data_fil_cut_dwn(P_calib(p_start1):P_calib(p_end1),:));
ST_Peak1(3,:)=ST_Peak1(1,:)-ST_Peak1(2,:);

ST_Peak1_Det(1,:)=max(SG_fil_det(P_calib(p_start1):P_calib(p_end1),:));
ST_Peak1_Det(2,:)=min(SG_fil_det(P_calib(p_start1):P_calib(p_end1),:));
ST_Peak1_Det(3,:)=ST_Peak1_Det(1,:)-ST_Peak1_Det(2,:);

p_start2=4;p_end2=9;
temp_data2=data_fil_cut_dwn(P_calib(p_start2):P_calib(p_end2),:);
figure
plot(temp_data2)
ST_Peak2(1,:)=max(data_fil_cut_dwn(P_calib(p_start2):P_calib(p_end2),:));
ST_Peak2(2,:)=min(data_fil_cut_dwn(P_calib(p_start2):P_calib(p_end2),:));
ST_Peak2(3,:)=ST_Peak2(1,:)-ST_Peak2(2,:);

ST_Peak2_Det(1,:)=max(SG_fil_det(P_calib(p_start2):P_calib(p_end2),:));
ST_Peak2_Det(2,:)=min(SG_fil_det(P_calib(p_start2):P_calib(p_end2),:));
ST_Peak2_Det(3,:)=ST_Peak2_Det(1,:)-ST_Peak2_Det(2,:);

 
% for i=1:size(ST_mag_raw,1)-1;
%     ST_mag_det_diff_AVG2(i,:)=ST_mag_det_diff_AVG;
%     ST_mag_raw_diff_AVG2(i,:)=ST_mag_raw_diff_AVG;
% end
% ST_mag_det_diff_AVG2-abs(ST_mag_det_diff)

% strain_temp = zeros(size(peak_seg,2),ch); strain_sine_temp = zeros(size(peak_seg,2),ch); strain_sine_det_temp = zeros(size(peak_seg,2),ch); % "size(peak_seg,2)" means the number of tacho peak
% for i=1:ch
%     for j=1:size(peak_seg,2)
%         for k=1:peak_no(j)
%             strain_temp(j,i) = strain_temp(j,i) + data_fil_cut_dwn(peak_seg(k,j),i+1);
%             strain_sine_temp(j,i) = strain_sine_temp(j,i) + F_est(peak_seg(k,j),i);
%             strain_sine_det_temp(j,i) = strain_sine_det_temp(j,i) + SG_fil_det(peak_seg(k,j),i);
%         end
%     end
%     strain_temp(:,i) = (strain_temp(:,i) ./ peak_no');
%     strain_sine_temp(:,i) = (strain_sine_temp(:,i) ./ peak_no');
%     strain_sine_det_temp(:,i) = (strain_sine_det_temp(:,i) ./ peak_no');
% end
% 
% 
% if size(peak_seg,2) == 2
%     strain = strain_temp(2,:) - strain_temp(1,:);
%     strain_sine = strain_sine_temp(2,:) - strain_sine_temp(1,:);
%     strain_sine_det = strain_sine_det_temp(2,:) - strain_sine_det_temp(1,:);
% end
% if size(peak_seg,2) == 3
%     strain = ((strain_temp(2,:) - strain_temp(1,:)) + (strain_temp(2,:) - strain_temp(3,:)))/2;
%     strain_sine = ((strain_sine_temp(2,:) - strain_sine_temp(1,:)) + (strain_sine_temp(2,:) - strain_sine_temp(3,:)))/2;
%     strain_sine_det = ((strain_sine_det_temp(2,:) - strain_sine_det_temp(1,:)) + (strain_sine_det_temp(2,:) - strain_sine_det_temp(3,:)))/2;
% end
% if size(peak_seg,2) == 4
%     strain = ((strain_temp(2,:) - strain_temp(1,:)) + (strain_temp(2,:) - strain_temp(3,:)) + (strain_temp(3,:) - strain_temp(4,:)))/3;
%     strain_sine = ((strain_sine_temp(2,:) - strain_sine_temp(1,:)) + (strain_sine_temp(2,:) - strain_sine_temp(3,:)) + (strain_sine_temp(3,:) - strain_sine_temp(4,:)))/3;
%     strain_sine_det = ((strain_sine_det_temp(2,:) - strain_sine_det_temp(1,:)) + (strain_sine_det_temp(2,:) - strain_sine_det_temp(3,:)) + (strain_sine_det_temp(3,:) - strain_sine_det_temp(4,:)))/3;
% end
% 
% F_est(:,7) = data_fil_cut_dwn(:,ch+2);
% figure(400); plot(data_fil_cut_dwn(:,2:ch+2));
% figure(500); plot(F_est);
% figure(600); plot(SG_fil_det);
% RESULT = [strain; strain_sine; strain_sine_det]' * 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Manual calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Tacho=data_fil_cut_dwn(:,7);
% N_tach=find(Tacho>40);
% N_tach_read=[111 419 729 1041];
% ST_tacho=data_fil_cut_dwn(N_tach_read',:)';
% ST_tacho_amp1=(ST_tacho(:,1)-ST_tacho(:,3))/2
% ST_tacho_amp2=(ST_tacho(:,3)-ST_tacho(:,1))/2


