    x=t_fil_cut_dwn(:,1);
    y=data_fil_cut_dwn(:,1);%Sine + noise
    tic;
    [SineP]=sineFit(x,y) %Offset, amplitude, frquency, phase, MSE
    toc
%uncomment following line if you want to save y=f(x) and run it sineFitDemo
% save('xy.mat','x','y');