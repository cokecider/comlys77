function A_line_fit = line_fit(A)
%% Line min/max elimination

[m,n]=size(A);

A_nf=A(1,2);
A_normal=A;
A_normal(:,2)=A(:,2)-A_nf;

t_max=max(A(:,1));

idx_t=round(m/t_max);

for idx=t_max:t_max:m-t_max
    slp(idx/t_max,1)=A_normal(idx+t_max,2)-A_normal(idx,2);
end

slp_lim=find(abs(slp(:,1))>=max(abs(slp(:,1)))/10);
slp_min=min(slp_lim);
slp_max=max(slp_lim);

A_min=slp_min*t_max;
A_max=slp_max*t_max;

A_trim(:,1)=A_normal([A_min:A_max],1);
A_trim(:,2)=A_normal([A_min:A_max],2);

 A_sm=smoothdata(A_trim,'sgolay');

%   자세한 설명 위치

end