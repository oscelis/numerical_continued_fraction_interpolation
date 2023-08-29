% FUNCTION NAME:
%   cfrac_interpolate
%
% DESCRIPTION:
%   Returns coefficients of a continued fraction 
%       aa1 + (x-zz1)/aa2+... 
%   interpolating given real or complex data (xx,ff)
%
%
% INPUT:
%   xx  - (vector double) distinct locations of interpolation points
%   ff  - (vector double) function values at xx
%   tol - optional (double) stop adaptive sampling when tollerance is 
%         reached on remaining points, by default 1e-13
%   NN  - optional (integer) max. number of interpolation points to use, 
%         by default all points are used unless tol is reached earlier
%
%
% OUTPUT:
%   aa - (vector double) inverse differences
%   zz = (vector double) associated locations used
%
% ASSUMPTIONS AND LIMITATIONS:
%   xx and ff should be finite
%
%
% REFERENCE:
%   Numerical Continued Fraction Interpolation
%   Oliver Salazar Celis
%   https://doi.org/10.48550/arXiv.2109.10529
%
%
% See also evalcfrac, prz_cfrac
function [aa,zz] = cfrac_interpolate(xx,ff,tol,NN)
%Returns coefficients of continued fraction aa1 + (x-zz1)/aa2+...
% interpolating data (xx,ff)
if nargin<3, tol = 1e-13; end       %default tol
if nargin<4, NN = length(xx); end   %how many points to use maximun
xx = xx(:); ff= ff(:);
aa=nan(NN,1); zz=nan(NN,1); rr=nan(NN,1);
    for k=1:NN %main loop
        if k==1 %init
            rr = ff; %inverse differences
            [~,i] = min(abs(ff)); %smallest value
        else
            [~,i]= max(abs(evalcfrac(aa,zz,xx)-ff)); %adaptive choice
            rr=(xx-zz(k-1))./(rr-aa(k-1)); %inverse differences
        end
        aa(k)=rr(i);zz(k)=xx(i); %store cfrac coef
        ff(i) = [] ;rr(i) =[]; xx(i)=[]; %reduce data
        if k < NN 
            if max(abs(evalcfrac(aa,zz,xx)-ff)) < tol*max(abs(ff))
                disp(['Target precision reached earlier at n=' num2str(k)]);
                break
            end
        end
    end
    aa = aa(~isnan(aa)); zz = zz(~isnan(zz));
end
