function DmFy = chebdifft(Fy,M_dev)

% The function Dmf = chebdifft(f,M_dev) computes the M'th
% approximate Chebyshev derivatives of the data vector y.
% A Fast Fourier Transform is used compute the Chebyshev cofficients
% of the data vector. A recursion formula is used to compute the
% Chebyshev coefficients for each derivative. A FFT is then used again
% to compute the derivatives in physical space.
% 
%  Input:
%  f:        Vector containing function values at the Chebyshev points
%            x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  M:        Derivative required (positive integer)
%
%  Output:
%  ym:       Vector containing approximate M'th derivative
%            
%  J.A.C. Weideman, S.C. Reddy 2000.  
[N,n_modes]=size(Fy);
DmFy=zeros(N,n_modes);
for j=1:n_modes
f=Fy(:,j);                                      % Make sure f is a vector     
%N=length(f);      
a0=real(fft([f; flipud(f(2:N-1))]));               % Extend and compute fft
a0=a0(1:N).*[0.5; ones(N-2,1); 0.5]/(N-1);   % a0 contains Chebyshev 

if M_dev==0
    back=[2*a0(1);a0(2:N-1);2*a0(N);flipud(a0(2:N-1))];
else
    a=[a0 zeros(N,M_dev)];                           % Recursion formula
    for ell=1:M_dev                                  % for computing coefficients
        a(N-ell,ell+1)=2*(N-ell)*a(N-ell+1,ell);   % of ell'th derivative 
        for k=N-ell-2:-1:1
            a(k+1,ell+1)=a(k+3,ell+1)+2*(k+1)*a(k+2,ell);
        end;
        a(1,ell+1)=a(2,ell)+a(3,ell+1)/2;
    end;
    back=[2*a(1,M_dev+1); a(2:N-1,M_dev+1); 2*a(N,M_dev+1); flipud(a(2:N-1,M_dev+1))];
end
Dmf=0.5*fft(back);                            % Transform back to
Dmf=Dmf(1:N);                                 % physical space

if max(abs(imag(f))) == 0; Dmf = real(Dmf); end    % Real data in,
                                                    %real derivative out
DmFy(:,j)=Dmf;
end
                                                    
end
