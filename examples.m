%Some examples illustrating the use of the functions
%   cfrac_interpolate
%   evalcfrac
%   prz_cfrac
%Examples are based on 
%Y. Nakatsukasa, O. Sete, and L. N. Trefethen. 
%The AAA algorithm for rational approximation, 
%SIAM J. Sci. Comput., 40(3):A1494--A1522, 2018.
%

close all;
clear all;

%example 0: complex square example
npts = 2000;
Z = rand(npts,1) + 1i*rand(npts,1);
ff = @(z) sqrt(z.*(1-z));
[aa,xx] = cfrac_interpolate(Z,ff(Z));
[pol,zer, res] = prz_cfrac(aa,xx);
%plotting
figure()
plot(Z,'o'); hold on;
plot(pol,'*', 'Color', 'r');
xlim([-1 2])
ylim([-1.5 1.5])
grid on
pbaspect([1 1 1])
xlabel('Real')
ylabel('Imag')
legend('data location','poles')
title('Clustering of poles near branch cuts')


%example 1: gamma function
Z = linspace(-1.5,1.5);
ff = @(z) gamma(z);
[aa,xx] = cfrac_interpolate(Z,ff(Z));
zz = linspace(-3.5,4.5,1000).';
figure()
plot(zz,evalcfrac(aa,xx,zz))
xlim([-3.5 4.5])
ylim([-8 8])
grid on
xlabel('x')
ylabel('r(x)')
title('Rational approximation of the Gamma function')
[pol,zer, res] = prz_cfrac(aa,xx);
disp(['poles on the negative axis:' num2str(pol(find(real(pol)<0))')])

%example 4: Analytic functions in the unit disk
ff = @(z) tan(z);
Z = exp(linspace(-1,1, 128)*1i*2*pi); %equidistant points on the unit circle
alpha = rand(1000,1)*2*pi; %random angles
ZZ = rand(1000,1).*(cos(alpha)+1i*sin(alpha)); %random points inside the unit circle
nn = 30;
errors = inf(nn,1);
for k=1:nn
    [aa,xx] = cfrac_interpolate(Z,ff(Z), 10^-13,k);
    errors(k) = max(abs(ff(ZZ)-evalcfrac(aa,xx,ZZ)));
end
figure()
plot(log10(errors),'o')
grid on
xlabel('number of points')
ylabel('max error inside the unit circle')
title('rational approximation to tan(z) sampled on the unit circle ')

%example 5: Analytic functions with nearby branch points
ff = @(z) log(1.1-z);
Z = exp(linspace(-1,1, 256)*1i*2*pi); %equidistant points on the unit circle
alpha = rand(1000,1)*2*pi; %random angles
ZZ = rand(1000,1).*(cos(alpha)+1i*sin(alpha));%random points inside the unit circle
nn = 40;
errors = inf(nn,1);
for k=1:nn
    [aa,xx] = cfrac_interpolate(Z,ff(Z), 10^-13,k);
    errors(k) = max(abs(ff(ZZ)-evalcfrac(aa,xx,ZZ)));
end
figure()
plot(log10(errors),'o')
grid on
xlabel('number of points')
ylabel('max error inside the unit circle')
title('rational approximation to log(1.1-z) sampled on the unit circle ')


%example 7:Approximation in other connected domains.
npts = 2000;
Z = rand(npts,1).*1i.*(-1).^(randi(2,npts,1)-1) + rand(npts,1)*10 ;
ff = @(z) 1./besselj(0,z);

[aa,xx] = cfrac_interpolate(Z,ff(Z), 1e-13);
[pol,zer, res] = prz_cfrac(aa, xx);
disp(['poles in rectangle:'])
pol(real(pol)>=0 & real(pol)<=10 & abs(imag(pol))>=0 & imag(pol)<1)

figure()
plot(Z,'o'); hold on;
plot(pol,'*', 'Color', 'r');
xlim([-5 15])
ylim([-3 3])
grid on
xlabel('Real')
ylabel('Imag')
legend('data location','poles')
title('Approximation of 1/J0(z) in 2000 random points in a rectangle in the complex plane')
