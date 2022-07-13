function SineParams=sineFit(x,y,ch,varargin)
if nargin>3 && varargin{1}==0
  boolGraphic=false;
else
  boolGraphic=true;
end
%% FFT
pi2=2*pi;
NumSamples=length(x);
T=x(2)-x(1);
fNy=1/(2*T);%나이키스트 주파수
offs=mean(y);%오프셋 벨류
y_m=y-offs;%오프셋 소거
n = 128*2^nextpow2(NumSamples);
Y = fft(y_m,n);
n2=floor(n/2);
P2 = abs(Y/NumSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = (0:n2)/n/T;% 주파수 스케일
%FFT 파라미터
[maxFFT,maxFFTindx]=max(P1);%최대값 크기,위치
fPeak=fs(maxFFTindx);
Phip=angle(Y(maxFFTindx))+pi/2;
Phip=Phip-x(1)*fPeak*pi2;
if Phip<0;Phip=2*pi+Phip;end
%% Fitting
if numel(x)<12 || x(end)-x(1)<5/fPeak
  offs=(max(y)+min(y))/2;
elseif abs(offs)<0.1 && maxFFT>0.9
  offs=0;
end
paramsFFTp=[offs,maxFFT,fPeak,Phip];
P1P=P1(1:maxFFTindx);
if maxFFTindx>0.99*n2 && ~(maxFFTindx+2<n2 && P1(n2-2)<0.7*maxFFT)
  fIndxExtra=[round(maxFFTindx*0.995);find(P1P<0.7*maxFFT,1,'last');find(P1P<0.5*maxFFT,1,'last')];
  AExtra=(max(y)-min(y))/2;
  Numf=3; %주파수 상수
  PeakInds = find(diff(sign(diff(P1))))+1;%FFT peack 위치
  PeakVals=sortrows([PeakInds(:),fs(PeakInds)',P1(PeakInds)'],3);%솔팅
  if numel(PeakInds)>1 && PeakVals(end,3)*.95<PeakVals(end-1,3)
    fIndxExtra=[fIndxExtra;PeakVals(end-1,1)];
    Numf=4;
  end
elseif x(end)-x(1)<1/fPeak
  fIndxExtra=[round(maxFFTindx*0.9);find(P1P<0.6*maxFFT,1,'last')];
  AExtra=(max(y)-min(y));
  Numf=2;
else
  Numf=1;
  paramsFFT=paramsFFTp;
end
if Numf>1
  fExtra=fs(fIndxExtra);
  PhiExtra=(angle(Y(fIndxExtra))+pi/2-x(1)*fExtra*pi2);
  AExtra=repelem(AExtra,Numf);
  offExtra=repelem(offs,Numf);
  paramsFFT=[offExtra',AExtra',fExtra',PhiExtra'];
end
paramsOut=zeros(Numf,6);
%fmin탐색
for i=1:Numf
  fun = @(SineParams)sseval(SineParams,x,y);
  [SineParams,SE] =fminsearch(fun,paramsFFT(i,:),...%SE:에러 제곱
    optimset('MaxFunEvals',200000,'MaxIter',200000));
  %주파수 변환
  if SineParams(3)<0
    SineParams(3)=-SineParams(3);
    SineParams(4)=pi-SineParams(4);%sin(2*pi*-f-phi)=sin(2*pi*f+phi+pi)
  end
  %크기 변환
  if SineParams(2)<0
    SineParams(2)=-SineParams(2);
    SineParams(4)=SineParams(4)+pi;
  end
  MSE=SE/numel(x);
  paramsOut(i,:)=[SineParams,MSE,MSE];
  if SineParams(3)>fNy
    paramsOut(i,5)=Inf;
  end
end
%% MSE에 따른 최적값 선정
[MSEmin,MSEminIndx]=min(paramsOut(:,5));
SineParams=paramsOut(MSEminIndx,1:4);

if MSEmin<=0.00001 || ...
    NumSamples<5 || ... 
    (NumSamples==5 && SineParams(3)<0.8*paramsFFT(1,3)) ||... 
    (MSEmin<1 && x(end)-x(1)<0.5/SineParams(3))
  maxAmp=66*maxFFT;
elseif MSEmin>0.3
  maxAmp=4*maxFFT;
elseif MSEmin>0.01
  maxAmp=6*maxFFT;
elseif MSEmin>0.001
  maxAmp=18*maxFFT;
else
  maxAmp=33*maxFFT;
end
if SineParams(2)>maxAmp || SineParams(3)>fNy
  SineParams=paramsFFTp;
  MSE=(sum((y - (SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*x+SineParams(4)))).^2))/numel(x);
  SineParams(5)=-MSE;
else
  SineParams(5)=paramsOut(MSEminIndx,6);
end

SineParams(4)=rem(SineParams(4),pi2);
if SineParams(4)<0
  SineParams(4)=SineParams(4)+pi2;
end
if boolGraphic
  PlotResults(x,y,SineParams,paramsFFTp,fs,P1,maxFFTindx,maxFFT,ch);
end

function sse = sseval(SineParams,x,y)
offs = SineParams(1);
A =SineParams(2);
f=SineParams(3);
Phi=SineParams(4);
sse = sum((y - (offs+A*sin(2*pi*f*x+Phi))).^2);

%% Plot results (optional)
function PlotResults(x,y,SineParams,paramsFFTp,fs,P1,maxFFTindx,maxFFT,ch)
xstart=x(1);
xend=x(end);
xLength=xend-xstart;
xSstep=min(xLength/100,1/SineParams(3)*0.1);
xS=xstart:xSstep:xend;
xFFT=xstart-xLength*.02:xLength/100:xend+xLength*.02;%FFTcurve longer
y5=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*xS+SineParams(4));%result
yFFT=paramsFFTp(1)+paramsFFTp(2)*sin(2*pi*paramsFFTp(3)*xFFT+paramsFFTp(4));
ch=ch;
%time plot:
%hFigPlotSin = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotSin' );
%if isempty(hFigPlotSin)
%  screensize=get(0, 'MonitorPositions');
%  hFigPlotSin=figure('Tag','Fig$PlotSin','Name','Sinus ',...
%    'OuterPosition',[960,screensize(1,4)/2,screensize(1,3)-960,screensize(1,4)/2]);
drawnow
%end
%time series plot
figure('Name',sprintf('Sinefitting Result : Ch %d', ch),'NumberTitle','off');
cla reset;
plot(x,y,'k.');%time series as dots
xlabel('Time [s]');
hold on;
pIn=plot(x,y,'r-');%time series as line
pFFT=plot(xFFT,yFFT,'color',[0.7 0.7 0.7]);
pResult=plot(xS,y5,'b-');%result
legend([pIn,pResult,pFFT],'Input','Result', 'FFT peak');
hold off;
grid on;
axis tight;


%FFT plot:
%hFigPlotFFT = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotFFT' );
%if isempty(hFigPlotFFT)
%  hFigPlotFFT=figure('Tag','Fig$PlotFFT','Name','FFT',...
%    'OuterPosition',[960,40,screensize(1,3)-960,screensize(1,4)/2-45]);
drawnow
%end

figure('Name',sprintf('FFT Result : Ch %d', ch),'NumberTitle','off');
cla reset;
pFFTin=plot(fs,P1,'r-');
xlabel('Frequency [Hz]');
ylabel('Amplitude')
hold on;
pFFTmax=plot(fs(maxFFTindx),maxFFT,'r+','MarkerSize',12);%max FFT
pFFTresult=plot(SineParams(3),SineParams(2),'b+','LineWidth',2);
plot([SineParams(3),SineParams(3)],[0,max(max(P1)*1.01,SineParams(2))],'b-');
hLeg=legend([pFFTin,pFFTresult,pFFTmax],'Input',...
  ['Result:     ' num2str(SineParams(2),3) ', ' num2str(SineParams(3),3) ' Hz'],...
  ['max FFT:  ' num2str(maxFFT,3) ', ' num2str(fs(maxFFTindx),3) ' Hz'],...
  'Location','best');
title(hLeg,'        amplitude | frequency','FontSize',8);
hold off;
grid on;
disp(['Result:        y= ' num2str(SineParams(1)) ' + ' num2str(SineParams(2)) ...
  ' * sin(2*pi*' num2str(SineParams(3)) '+' num2str(SineParams(4)) ')   MSE: ' num2str(SineParams(5))]);
disp(['FFT:           y= ' num2str(paramsFFTp(1)) ' + ' num2str(paramsFFTp(2)) ...
  ' * sin(2*pi*' num2str(paramsFFTp(3)) '+' num2str(paramsFFTp(4)) ')' ]);