function [x,nw,np,nt] = osc_wold(x,y,nocomp,iter,tol)
% [x,nw,np,nt] = osc_wold(x,y,nocomp,iter,tol)
% Orthogonal singal correction method developed by Fearn et. al
%
%Input
%x:            The matrix of predictor variables (x)
%y:            Predicted variable(s) (y), scaled as desired, 
%nocomp: The number of OSC components to calculate.
%iter:         Optional input variables are the maximum number of
%               iterations used in attempting to maximize the variance
%               captured by orthogonal component (default = 0),
%tol:          the tolerance on percent of x variance to consider
%               in formation of the final w vector (tol, default = 99.9).
%
%Output
%x:          OSC filtered signals;
%nw:       OSC corrected weights, 
%np:       OSC corrected loads vector
%nt:        OSC corrected scores that were used in making the correction. 
% 
%  Once the calibration is done, new (scaled) x data can be corrected by 
%  newx = x - x*nw*inv(np'*nw)*np';
%
%I/O: [nx,nw,np,nt] = osccalc(x,y,nocomp,iter,tol);
%
%See also: CROSSVAL
%
%Edit by Da Chen, Janurary 2, 2009

[m,n] = size(x);
nw = zeros(n,nocomp);
np = zeros(n,nocomp);
nt = zeros(m,nocomp);
tx=x;
if nargin < 4 | isempty(iter)
  iter = 1000;
end
if nargin < 5 | isempty(tol)
  tol = 99.9;
end

for i=1:nocomp
  [u,s,v]=svds(x,1);
  t = u(:,1)*s(1);
  %p=v(:,1);
  %t=x*p/(p'*p);
  
   dif=1;
   k=0;
   told=t;
   
  while dif>1e-3
    k=k+1;
    t=t-y*inv(y'*y)*y'*t;
    ty=t;
    w = somesimpls(tx,ty,tol);
    w = w/norm(w);
    ty=[];
    %t=tt*inv(tt'*tt)*tt'*t;
    t=x*w;
    dif=norm(t-told)/norm(t);
    told=t;
    if k > iter
       dif = 0;
    end
  end
 
 %w=t'*inv(x*x')*x;
 p=t'*x/(t'*t);
 x=x-t*p;
 p=p';
 np(:,i) = p;
 nw(:,i) = w;
 nt(:,i) = t;
end

%---------------------------------------------------
function [reg,i] = somesimpls(x,y,tol)
%return SIMPLS weights which correspond to a total variance of "tol"
%I/O: [weights,i] = sim(x,y,tol)

s      = x'*y;
totvar = sum(sum(x.^2));
total  = 0;
i      = 0;

while total < tol
  
  i  = i+1;

  rr = s;           %weights from covar.
  tt = x*rr;
  normtt = norm(tt);
  tt = tt/normtt;
  rr = rr/normtt;
  pp = (tt'*x)';
  
  qq = y'*tt;
  uu = y*qq;
  vv = pp;
  if i > 1
    vv = vv - basis*(basis'*pp);
    uu = uu - loads{1,1}*(loads{1,1}'*uu);
  end
  vv = vv/norm(vv);
  s  = s - vv*(vv'*s);
   
  total = total + (pp'*pp)/totvar*100;
  
  wts(:,i)        = rr;           % x-block weights
  loads{1,1}(:,i) = tt;           % x-block scores
  loads{2,2}(:,i) = qq;           % y-block loadings
  basis(:,i)      = vv;           % basis of x-loadings
  
end

if i > 1
  reg = sum((wts*diag(loads{2,2}))');
else
  reg = (wts*loads{2,2})';
end
reg = reg'; 
