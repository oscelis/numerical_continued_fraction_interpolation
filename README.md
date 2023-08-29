# Numerical Continued Fraction Interpolation.

This repository contains MATLAB codes for performing [Numerical Continued Fraction Interpolation](https://doi.org/10.48550/arXiv.2109.10529).

## Examples.
Some of the below examples are inspired by [AAA](https://doi.org/10.1137/16M1106122).

### Approximation in a rectangular domain.
We approximate the reciprocal Bessel function $1/J0(z)$ in $2000$ uniformly distributed random points in the rectangle defined by the corners $±i$ and $10±i$.
 The interpolation data is setup as below. 
```matlab
npts = 2000;
Z = rand(npts,1).*1i.*(-1).^(randi(2,npts,1)-1) + rand(npts,1)*10 ;
ff = @(z) 1./besselj(0,z);
```
Then we construct the continued fraction rational interpolant with (default) tollerance $10^{-13}$. The adaptive procedure stops after chosing $n=25$ points.

```matlab
[aa,xx] = cfrac_interpolate(Z,ff(Z), 1e-13);
```
The poles of this rational function can be obtained as follows.
```matlab
[pol,zer, res] = prz_cfrac(aa, xx);
```
It can be checked that the obtained poles inside the rectangle agree well with the [roots of the Bessel function](https://mathworld.wolfram.com/BesselFunctionZeros.html) in that domain.
```matlab
pol(real(pol)>=0 & real(pol)<=10 & abs(imag(pol))>=0 & imag(pol)<1)
```
```matlab
ans =
   2.4048 + 0.0000i
   5.5201 + 0.0000i
   8.6537 - 0.0000i
```
A figure of the data and the location of the obtained poles is shown below.
```matlab
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
```
![bessel](https://github.com/oscelis/numerical_continued_fraction_interpolation/assets/7952417/d0cc58af-923e-49e6-b8cf-e664cb64effb)

### Analytic functions in the unit disk.
We approximate the function $tan(z)$ sampled in 128 equidistant points on the unit circle 
```matlab
ff = @(z) tan(z);
Z = exp(linspace(-1,1, 128)*1i*2*pi); %equidistant points on the unit circle
```
and check how well it approximates inside the unit disk.
```matlab
alpha = rand(1000,1)*2*pi; %random angles
ZZ = rand(1000,1).*(cos(alpha)+1i*sin(alpha)); %random points inside the unit circle
```
It is observed that, as more interpolation points are allowed to be taken, then 10-digit accuracy is quickly achieved.  
```matlab
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
```
![tanz](https://github.com/oscelis/numerical_continued_fraction_interpolation/assets/7952417/763c0548-63f9-4b74-9fb4-3fc7f42873fa)
