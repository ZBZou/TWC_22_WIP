function X = HOGMT(s, h, ep)
h_r = real(h);
h_i = imag(h);
H = [h_r, -h_i;h_i, h_r];
[U,temp_sig,V] = svd(H);
Sig = diag(temp_sig);
for i = 1:length(Sig)                
   if Sig(i) >= ep
      N = i; % N most contributed eigenfunctions
   end
end 

% construct X
sr = real(s);
si = imag(s);

S = [sr;si];

X = zeros(length(S),1);

for i = 1:N
   xn = dot(S, U(:,i))/Sig(i);
   X = X + xn*V(:,i);
end

Xr = X(1:length(s),:);
Xi = X(length(s)+1:end,:)*1i;
X = Xr+Xi;