% BISECTION Finds function zeros.
%
% function by Quarteroni et al. (Quarteroni et al.,
%   Scientific computing with MatLab and Octave
%
% ZERO=BISECTION(FUN,A,B,TOL,NMAX) tries to find a zero
% ZERO of the continuous function FUN in the interval
% [A,B] using the bisection method. If
% the search fails an error message is displayed.
% FUN is a function handle associated with an anonymous
% function or a Matlab function.
% ZERO=BISECTION(FUN,A,B,TOL,NMAX,P1,P2,...) passes
% parameters P1,P2,... to the function FUN(X,P1,P2,...)
% [ZERO,RES,NITER]=BISECTION(FUN,...) returns the value
% of the residual in ZERO and the iteration number at
% which ZERO was computed.
function [zero,res,niter]=bisection(fun,a,b,tol, nmax,varargin)
x = [a, (a+b)*0.5, b];
fx = fun(x,varargin{:});
if fx(1)*fx(3) > 0
    error([' The sign of the function at the ',...
        'endpoints of the interval must be different\n']);
elseif fx(1) == 0
    zero = a; res = 0; niter = 0; return
elseif fx(3) == 0
    zero = b; res = 0; niter = 0; return
end
niter = 0;
I = (b - a)*0.5;
while I >= tol & niter < nmax
    niter = niter + 1;
    if fx(1)*fx(2) <  0
        x(3) = x(2);
        x(2) = x(1)+(x(3)-x(1))*0.5;
        fx = fun(x,varargin{:});
        I = (x(3)-x(1))*0.5;
    elseif fx(2)*fx(3) < 0
        x(1) = x(2);
        x(2) = x(1)+(x(3)-x(1))*0.5;
        fx = fun(x,varargin{:});
        I = (x(3)-x(1))*0.5;
    else
        x(2) = x(find(fx==0)); I = 0;
    end
end
if  (niter==nmax && I > tol)
    fprintf(['Bisection stopped without converging ',...
        'to the desired tolerance because the \n',...
        'maximum number of iterations was reached\n']);
end
zero = x(2);
x = x(2);
res = fun(x,varargin{:});
end