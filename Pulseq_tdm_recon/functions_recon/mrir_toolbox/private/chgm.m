function [a,b,x,hg]=chgm(a,b,x,hg,varargin);
%     ===================================================
%     Purpose: Compute confluent hypergeometric function
%     M(a,b,x)
%     Input  : a  --- Parameter
%     b  --- Parameter(b <> 0,-1,-2,...)
%     x  --- Argument
%     Output:  HG --- M(a,b,x)
%     Routine called: GAMMA for computing â(x)
%     ===================================================

%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

  ta=[];tb=[];xg=[];tba=[];
  a0=a;
  a1=a;
  x0=x;
  hg=0.0d0;

  if(b == 0.0d0|b == -abs(fix(b)));
    hg=1.0d+300;
  elseif(a == 0.0d0|x == 0.0d0);
    hg=1.0d0;
  elseif(a == -1.0d0);
    hg=1.0d0-x./b;
  elseif(a == b);
    hg=exp(x);
  elseif(a-b == 1.0d0);
    hg=(1.0d0+x./b).*exp(x);
  elseif(a == 1.0d0&b == 2.0d0);
    hg=(exp(x)-1.0d0)./x;
  elseif(a == fix(a)&a < 0.0d0);
    m=fix(-a);
    r=1.0d0;
    hg=1.0d0;
    for  k=1:m;
      r=r.*(a+k-1.0d0)./k./(b+k-1.0d0).*x;
      hg=hg+r;
    end;  k=m+1;
  end;

  if(hg ~= 0.0d0)return; end;

  if(x < 0.0d0);
    a=b-a;
    a0=a;
    x=abs(x);
  end;
  if(a < 2.0d0)nl=0; end;
  if(a >= 2.0d0);
    nl=1;
    la=fix(a);
    a=a-la-1.0d0;
  end;
  for  n=0:nl;
    if(a0 >= 2.0d0)a=a+1.0d0; end;
    if(x <= 30.0d0+abs(b)|a < 0.0d0);
      hg=1.0d0;
      rg=1.0d0;
      for  j=1:500;
	rg=rg.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*x;
	hg=hg+rg;
	if(abs(rg./hg)< 1.0d-15)break; end;
      end;
    else;
      [ta]=gamma(a);
      [tb]=gamma(b);
      xg=b-a;
      [tba]=gamma(xg);
      sum1=1.0d0;
      sum2=1.0d0;
      r1=1.0d0;
      r2=1.0d0;
      for  i=1:8;
	r1=-r1.*(a+i-1.0d0).*(a-b+i)./(x.*i);
	r2=-r2.*(b-a+i-1.0d0).*(a-i)./(x.*i);
	sum1=sum1+r1;
	sum2=sum2+r2;
      end;  i=8+1;
      hg1=tb./tba.*x.^(-a).*cos(pi.*a).*sum1;
      hg2=tb./ta.*exp(x).*x.^(a-b).*sum2;
      hg=hg1+hg2;
    end;
    if(n == 0)y0=hg; end;
    if(n == 1)y1=hg; end;
  end;
  if(a0 >= 2.0d0);
    for  i=1:la-1;
      hg=((2.0d0.*a-b+x).*y1+(b-a).*y0)./a;
      y0=y1;
      y1=hg;
      a=a+1.0d0;
    end;  i=la-1+1;
  end;
  if(x0 < 0.0d0)hg=hg.*exp(x0); end;
  a=a1;
  x=x0;
  return;
end
