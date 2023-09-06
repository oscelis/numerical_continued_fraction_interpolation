# Numerical Continued Fraction Interpolation.

This repository contains MATLAB codes for performing [Numerical Continued Fraction Interpolation](https://doi.org/10.48550/arXiv.2109.10529).

If you use these codes in any published research, please cite the following publication(s):
* O. Salazar Celis, **Numerical continued fraction interpolation**, *Ukrainian Mathematical Journal*, 2023, (accepted - in press) [arXiv](https://doi.org/10.48550/arXiv.2109.10529).
* O. Salazar Celis, **Adaptive Thiele interpolation**, *ACM Communications in Computer Algebra*, 56(3):125--132, 2022 ([doi](https://doi.org/10.1145/3594252.3594254))


# Examples.
Some of the below examples are inspired by [AAA](https://doi.org/10.1137/16M1106122).

## Approximating the Gamma function.
We approximate the Gamma function $\Gamma(z)$ on the real line using 100 linearly spaced points between $-1.5$ and $1.5$
```matlab
Z = linspace(-1.5,1.5);
ff = @(z) gamma(z);
```
One obtains the continued fraction rational interpolant as follows. The adaptive procedure stops after chosing $n=19$ points when the default tollerance of $10^{-13}$ is reached on all given points.
```matlab
[aa,xx] = cfrac_interpolate(Z,ff(Z));
```
One can plot the approximation on a wider domain such as ($-3.5, 4.5)$
```matlab
zz = linspace(-3.5,4.5,1000).';
figure()
plot(zz,evalcfrac(aa,xx,zz))
xlim([-3.5 4.5])
ylim([-8 8])
grid on
xlabel('x')
ylabel('r(x)')
title('Rational approximation of the Gamma function')
```
![gamma](https://github.com/oscelis/numerical_continued_fraction_interpolation/assets/7952417/b46fd555-8d59-4c49-868e-506236de9868)

The poles, zeros and residues of this rational function are obtained as follows.
```matlab
[pol,zer, res] = prz_cfrac(aa,xx);
```
It can be checked that the obtained poles on the negative axis and their associated residues agree well with the first ones of $\Gamma(z)$ although their accuracy decreases as one moves further from the interpolation interval. Recall that $\Gamma(z)$ has poles $0,-1,-2,-3,\ldots$ and associated [residues](https://en.wikipedia.org/wiki/Gamma_function#Residues) $Res(\Gamma,-k) = (-1)^k /k!$.
```matlab
[pol(find(real(pol)<0)) res(find(real(pol)<0))]
```
```matlab
ans =
  -0.0000 + 0.0000i   1.0000 - 0.0000i
  -1.0000 + 0.0000i  -1.0000 + 0.0000i
  -2.0000 + 0.0000i   0.5000 - 0.0000i
  -3.0033 + 0.0000i  -0.1697 + 0.0000i
  -3.7948 + 0.0000i   0.0368 - 0.0000i```
```

## Approximation in a rectangular domain.
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

## Analytic functions in the unit disk with nearby poles.
We approximate the function $tan(z)$ sampled in 128 equidistant points on the unit circle 
```matlab
ff = @(z) tan(z);
Z = exp(linspace(-1,1, 128)*1i*2*pi); %equidistant points on the unit circle
```
and check how well it approximates inside the unit disk.
```matlab
alpha = rand(1000,1)*2*pi; %random angles
ZZ = rand(1000,1).*(cos(alpha)+1i*sin(alpha)); %random points inside the unit disk
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
ylabel('max error inside the unit disk')
title('rational approximation to tan(z) sampled on the unit circle ')
```
![tanz](https://github.com/oscelis/numerical_continued_fraction_interpolation/assets/7952417/4b0cc71c-c320-4b9e-930c-1631ba3238c4)

The function $tan(z)$ has poles at $(k+1/2)\pi$ for $k \in \mathbb{Z}$. As expected from a rational approximation method, the poles of the rational function obtained with $n=17$ capture those of tan(z) closest to the unit circle $\pm \pi/2 \approx 1.5708$ quite well.
```matlab
[pol,zer, res] = prz_cfrac(aa,xx);
```
```matlab
pol(abs(pol)< 2)
```
```matlab
ans =
  -1.5708 + 0.0000i
   1.5708 + 0.0000i
```
## Analytic functions in the unit disk with nearby branch point singularities.
Rational approximations are also usefull when singularities are not poles but branch point singularities. 
In a similar fashion as before we approximate $log(1.1-z)$ sampled in 256 equidistant points on the unit circle
```matlab
ff = @(z) log(1.1-z);
Z = exp(linspace(-1,1, 256)*1i*2*pi); %equidistant points on the unit circle
```
Also here 10-digit accuracy is quickly achieved inside the unit disk. 
![log](https://github.com/oscelis/numerical_continued_fraction_interpolation/assets/7952417/588b39ec-fe57-4aae-b557-1f1217a45882)
