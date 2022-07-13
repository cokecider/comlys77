function A_wsvd_inv = wsvd_inv(A, scale_factor)

[U, S, V] = svd(A);

S_T = zeros(size(S));
[m,n] = size(S);
D = min(m,n);

for dd = 1:D
  
  if dd == D
    S_T(dd,dd) = 1/S(dd,dd)*scale_factor;
  else
    S_T(dd,dd) = 1/S(dd,dd);
  end
  
end

A_wsvd_inv = V * S_T' * U';

end
% function