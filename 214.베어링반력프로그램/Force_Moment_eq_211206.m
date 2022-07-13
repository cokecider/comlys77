% Location [mm]
L_R  = [4045   13495  21276.2 29151.2];
L_SG = [16583  20310  22152   30027  ];
L_BNNT = 40;
L_N    = 227.5;
L_Prop = 1470;

% weight [kgf]
W_BNNT = 25656; % 17983;
W_N    = 22315;
W_Prop = 772181.2; % 720113.8
W_AFT  = W_BNNT + W_N + W_Prop;
g=9.81;

temp_no_max = size(Shaft_element,2);
Shaft_element_ext(:,1:7) = Shaft_element(:,temp_no_max-6:end);
S_ext=Shaft_element_ext;

% vloume cal [unit : mm3]
for i=1:size(Shaft_element_ext,1)
    temp_vol= pi() * (1/3/4 * (S_ext(i,4)^2 + S_ext(i,5)^2 + S_ext(i,4)*S_ext(i,5)) - S_ext(i,6)^2/4) * S_ext(i,3);
%     temp_vol= (S_ext(i,4)^2-S_ext(i,6)^2)/4 * pi() * S_ext(i,3);
    S_ext(i,8)=temp_vol;
    S_ext(i,9)=S_ext(i,7)/S_ext(i,8);
end

% CG of SG location [mm]
for i=1:length(L_SG);
    index_k(i) = min( find ( S_ext(:,2) > L_SG(i) ) ) - 1;
end
for j=1:length(L_SG)
    M_sum(j)=0;
    W_sum(j)=0;
    L_SG_CG(j)=0;
    temp_W=0;
    temp_M=0;
    
    for i=1:index_k(j)
        if i==index_k(j)
            temp_W = ( L_SG(j) - S_ext(i,2) ) / S_ext(i,3) * S_ext(i,7) * g;
            temp_M = ( S_ext(i,2) + 0.5 * (L_SG(j) - S_ext(i,2) )) * temp_W;
        else
            temp_W = S_ext(i,7) * g; 
            temp_M = ( S_ext(i,2) + 0.5* S_ext(i,3) ) * temp_W;
        end
        M_sum(j) = M_sum(j) + temp_M;
        W_sum(j) = W_sum(j) + temp_W;
    end
    
    L_SG_CG(j)=M_sum(j)/W_sum(j);
end

% Force eq.
Eq_R1 = W_AFT + W_sum(1);
Eq_R2 = W_AFT + W_sum(2);
Eq_R3 = W_AFT + W_sum(3);
Eq_R4 = W_AFT + W_sum(4);

% Moment eq.
% M_SG=[10000 20000 30000 40000];
M_SG=[146308.633 -218598.352 -345914.767 -495346.505];
M_AFT = L_BNNT*W_BNNT + L_N*W_N + L_Prop*W_Prop;

M_AFT = M_AFT / 10^3 ;
L_SG_CG = L_SG_CG / 10^3 ;

Eq_M1    = M_AFT - M_SG(1) + L_SG_CG(1) * W_sum(1);
Eq_M2    = M_AFT + M_SG(2) + L_SG_CG(2) * W_sum(2);
Eq_M3    = M_AFT + M_SG(3) + L_SG_CG(3) * W_sum(3);
Eq_M4    = M_AFT + M_SG(4) + L_SG_CG(4) * W_sum(4);


% MAtrix calculation

A=[      1      1      0      0       1       0       0       0;   
    L_R(1) L_R(2)      0      0 L_SG(1)       0       0       0;
         1      1      0      0       0       1       0       0;   
    L_R(1) L_R(2)      0      0       0 L_SG(2)       0       0;
         1      1      1      0       0       0       1       0;   
    L_R(1) L_R(2) L_R(3)      0       0       0 L_SG(3)       0;    
         1      1      1      1       0       0       0       1;   
    L_R(1) L_R(2) L_R(3) L_R(4)       0       0       0 L_SG(4) ];

Eq=[Eq_R1; Eq_M1; Eq_R2; Eq_M2;  Eq_R3;   Eq_M3;  Eq_R4;  Eq_M4];


R = inv(A) * Eq;




