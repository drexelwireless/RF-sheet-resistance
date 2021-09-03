%This function uses polyfit to approximate Z0L
% function [Ro, Lo, Co, Go, gamma_t, gamma, beta, s11r, s12r, s11d, s12d, ...
%     vp,z_opt,Z0Lr,Z0L,x] = RLCG_Z(data1, dl, ub, lb, Lra, Cra)
function [Ro, Lo, Co, Go,gamma, gamma_t, Z0L,s11m,s21m,s11r,s21r,z_opt] = RLCG_func(data1, ub, lb, dl)
% function [Ro, Lo, Co, Go,gamma, gamma_t, Z0L, s11r, s21r, z_opt] = RLCG_Z2(data1, ub, lb)
freq =1e-9.* data1(:,1);
Z0 = 50; %Source impedance
n=length(data1);
L = 0.041; %length in meters
% L = 0.01; %length in meters

s11 = 10.^(data1(:,2)/20).*exp(1i*(data1(:,3)*pi/180)); %formula confirmed by VNA
s21 = 10.^(data1(:,4)/20).*exp(1i*(data1(:,5)*pi/180));
s12 = 10.^(data1(:,6)/20).*exp(1i*(data1(:,7)*pi/180)); %formula confirmed by VNA
s22 = 10.^(data1(:,8)/20).*exp(1i*(data1(:,9)*pi/180));
% s11 = 10.^(data1(:,2)/10).*exp(1i*(data1(:,3)*pi/180)); 
% s21 = 10.^(data1(:,4)/10).*exp(1i*(data1(:,5)*pi/180));
% s12=s21;
% s22=s11;
s11m = s11;
s21m = s21;
s12m = s12;
s22m = s22;

%extraction of beta
A = ((1+s11).*(1-s22) + s12.*s21)./(2.*s21);
gamma = (1/L).*acosh(A); %total gamma in Np/m
beta = imag(gamma);
beta2 = ((beta(2)-beta(1))./(freq(2)-freq(1)).*(freq-freq(2)) + beta(2));
%beta done
s_DUT = zeros(n,4);
%reference plane shift correction
for m=1:n
   c1 = [exp(-1i*beta2(m)*dl), 0; 0, exp(-1i*beta2(m)*dl)];
%    s_DUT = inv(c1)*[s11(m), s12(m); s21(m), s22(m)]*inv(c1);

   sdut = inv(c1)*[s11(m), s12(m); s21(m), s22(m)]*inv(c1);
   s_DUT(m,1) = sdut(1,1);
   s_DUT(m,2) = sdut(1,2);
   s_DUT(m,3) = s_DUT(m,2);
   s_DUT(m,4) = s_DUT(m,1);
end

s11 = s_DUT(:,1);
s12 = s_DUT(:,2);
s21 = s12;
s22 = s11;
%ABCD needed to get gamma and Z0L
A = ((1+s11).*(1-s22) + s12.*s21)./(2.*s21);
B = Z0* ((1+s11).*(1+s22) - s12.*s21)./(2.*s21);
C = (1/Z0).*((1-s11).*(1-s22) - s12.*s21)./(2.*s21);
D = ((1-s11).*(1+s22) + s12.*s21)./(2.*s21);
gamma = (1/L).*acosh(A); %total gamma in Np/m
alpha = real(gamma);
%smooth gamma
% gamma = movmean(gamma,20);


Z0L = sqrt(B./C);
beta = imag(gamma);
beta2 = ((beta(2)-beta(1))./(freq(2)-freq(1)).*(freq-freq(2)) + beta(2));
p = polyfit(freq,alpha,2);
y = polyval(p,freq);
% y = smoothdata(alpha);
gamma = real(gamma) + 1i.*beta2;
% gamma = y + 1i.*beta2; %polyfit
%optimization of Z = Zr + 1i*Zi
x0 = [0, 0];
A = [0,0;0,0];
b = [0,0]';
Aeq = [];
beq = [];
% e3 = @(x)objf(x,gamma,s11,s21); %should use the reconstructed s-params
z1 = real(Z0L);
z2 = imag(Z0L);
p1 = polyfit(freq,z1,2);
y1 = polyval(p1,freq);
p2 = polyfit(freq,z2,2);
y2 = polyval(p2,freq);
e3 = @(x)objf(x,gamma,s11,s21,Z0L);
x = fmincon(e3,x0,A,b,Aeq,beq,lb,ub);
%polyfit real and imaginary part of Z0L

% z_opt = y1 + 1i.*y2 + ones(16,1).* x(1) + ones(16,1).* 1i*x(2);
z_opt = y1 + ones(n,1).* x(1) + ones(n,1).* 1i*x(2);

% z_opt = x(1) + 1i*x(2);

Ro = real(gamma.*z_opt);
Lo = imag(gamma.*z_opt)./(2*pi.*freq);
Co = imag(gamma./z_opt)./(2*pi.*freq);
Go = real(gamma./z_opt);
gamma_t = sqrt((Ro+1i*2*pi*freq.*Lo).*(Go + 1i*2*pi.*freq.*Co));
zr = z_opt;
dsr = 2*zr.*50.*cosh(gamma.*L) + (zr.^2 + 50^2) .* sinh(gamma.*L);
s11r = (1./dsr).*(zr.^2 - 50^2).*sinh(gamma.*L); %s11 reconstructed
s21r = (1./dsr).*(2.*zr.*50);

function f = objf(x, g, s11, s21, Z0L)
z0 = 50;
L = 0.041;
% L = 0.01;

z1f = real(Z0L);
z2f = imag(Z0L);
p1f = polyfit(freq,z1f,2);
y1f = polyval(p1f,freq);
p2f = polyfit(freq,z2f,2);
y2f = polyval(p2f,freq);

% z = y1f + 1i*y2f + ones(1,n).*x(1) + 1i.*ones(1,n).*x(2);
z = y1f + ones(1,n).*x(1) + 1i.*ones(1,n).*x(2);

% z = x(1) + 1i*x(2);
ds = 2*z.*50.*cosh(g.*L) + (z.^2 + 50^2) .* sinh(g.*L);
s11x = (1./ds).*(z.^2-z0^2).*sinh(g.*L);
s21x = (1./ds).*2.*z*z0;
f = sum(sum((real(s11x+s21x-s11-s21)).^2)+sum((imag(s11x+s21x-s11-s21)).^2));
end
end