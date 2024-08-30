function psi0 = ic(x, k, lx, epsilon)
% PSI0 defines initial condition  psi0 = psi(x,0)
nmode=k;        % one mode=1; two modes=2; three modes=3;
                % gaussian=4; two gaussians=5; gaussian wavepacket=6;
                % standing=7; sawtooth=8; compact traveling wave=11
                % travelling wave = 12;
                % multiple period travelling wave = 13
                %nmode=10 for saved output with compact support

if nmode==0
    psi0 = exp(1i.*x);
elseif nmode==1.5
    psi0 = cos(x);
elseif nmode==1
   amp1 = .2;
   psi0 = amp1.*exp(1i.*x);
elseif nmode==2
   amp1 = 1.;  amp2 = .5*amp1*1i;
   %psi0 = amp1.*exp(1i.*x) + amp2.*exp(2i.*x);
   psi0 = -exp(1i.*x) + 0.5*exp(2i.*(x +2*pi^2));
elseif nmode==3
    amp1 = .001;  amp2 = .5*amp1*1i;  amp3 = 0.5i*amp1;
    psi0 = amp1.*exp(1i.*x) + amp2.*exp(2i.*(x-2.8)) + amp3*exp(-1i.*x);
elseif nmode==4
    amp0 = 1.; 
    psi0 = amp0*exp(-8.*(x-pi/2).^2);
elseif nmode==5
    amp0 = .6; amp1 = -.4*1i;
    psi0 = amp0.*exp(-5.*(x-pi/2).^2) + amp1.*exp(-5.*(x-3*pi/2).^2);    
elseif nmode==6
    amp0 = .5;
    psi0 = amp0.*exp(-3.*(x-pi).^2).*exp(64i.*x);
elseif nmode==7
    amp0 = .2; amp1 = -amp0^2;
    psi0 = amp0*exp(1i*x) + amp0*exp(-1i*x);
    psi0 = psi0 + amp1*exp(3i*x) + amp1*exp(-3i*x);
elseif nmode==8
    amp0=0.2;
    psi0(1:lx/4) = x(1:lx/4)/pi;                           % 0 <= x <= pi/2
    psi0((lx/4)+1:3*lx/4) = (pi - x((lx/4)+1:3*lx/4))/pi;  % pi/2 < x <= 3*pi/2
    psi0((3*lx/4)+1:lx) = (x((3*lx/4)+1:lx) - 2*pi)/pi;    % 3*pi/2 < x < 2*pi
    psi0(1:lx/2) = amp0*psi0(1:lx/2);
    psi0(lx/2+1:lx) = 1i*amp0*psi0(lx/2+1:lx);
    %psi0 = amp0*psi0;
%{
elseif nmode==10
    %alpha = 1.; amp0=1.;
    alpha = 2.; amp0 = 0.4;
    amp0=amp0*(pi^2/4)^alpha; %1
    psi0(1:lx/4) = 0.;                           % 0 <= x <= pi/2
    psi0((lx/4)+1:3*lx/4) = (4/pi^2)^alpha*amp0.*(pi^2/4 - (x((lx/4)+1:3*lx/4)-pi).^2).^alpha;  % pi/2 < x <= 3*pi/2
    psi0((3*lx/4)+1:lx) = 0.;    % 3*pi/2 < x < 2*pi
    %psi0 = psi0.*exp(1i.*10.*x);
%}
elseif nmode==10
    %alpha = 1.; amp0=1.;
    alpha = 2.; amp0 = 0.4;
    amp0=amp0*(pi^2/4)^alpha; %1
    psi0(1:lx/4) = 0.;                           % 0 <= x <= pi/2
    psi0((lx/4)+1:3*lx/4) = (4/pi^2)^alpha*amp0.*(pi^2/4 - (x((lx/4)+1:3*lx/4)-pi).^2).^alpha;  % pi/2 < x <= 3*pi/2
    psi0((3*lx/4)+1:lx) = 0.;    % 3*pi/2 < x < 2*pi
    %psi0 = psi0.*exp(1i.*10.*x);  
elseif nmode==101
    amp0=.08; %1
    alpha = 4.; %1.
    psi0(1:lx/4) = 0.;                           % 0 <= x <= pi/2
    psi0((lx/4)+1:3*lx/4) = amp0.*(pi^2/4 - (x((lx/4)+1:3*lx/4)-pi).^2).^alpha;  % pi/2 < x <= 3*pi/2
    psi0((3*lx/4)+1:lx) = 0.;    % 3*pi/2 < x < 2*pi

elseif nmode==11
    s = 1/sqrt(2); c = sqrt(1-s^2);
    s = s/sqrt(2); c=c/sqrt(2);
    psi0(1:lx/4+1) = 0.;                    % 0 <= x <= pi/2
    psi0((3*lx/4)+1:lx) = 0.;               % 3*pi/2 < x < 2*pi   
    nic= lx/2;
    theta = 1:nic;
    theta = pi*theta/length(theta);
    xi = 0.*theta;
    for ii = 1:nic
        xi(ii) = ftw(theta(ii),c,s);
    end
    u = c*cos(xi)-s;
    v = c*sin(xi);
    psi0((lx/4)+1:3*lx/4) = u + 1i.*v;  % pi/2 < x <= 3*pi/2
elseif nmode==12 %traveling wave solution
    s = 1./sqrt(2);
    %s = 0.8; %0\le s\le 1
    c=sqrt(1-s.^2);
    xi=0.*x;
    for ii = 1:lx
        xi(ii) = ftw(x(ii),c,s);
    end
    u = c*cos(xi)-s;
    v = c*sin(xi);
    psi0 = u + 1i.*v;
elseif nmode==13
    nper = 8;
    %s = .4; %1./sqrt(2);
    s=.8;
    c=sqrt(1-s.^2);
    s = s./sqrt(nper);
    c = c./sqrt(nper);
    xi=0.*x;
    for ii = 1:lx
        xi(ii) = ftw(x(ii),c,s);
    end
    u = c*cos(xi)-s;
    v = c*sin(xi);
    psi0 = u + 1i.*v;    
    psi0 = (1.+.25*cos(x)).*psi0;
elseif nmode==9
    amp0 = .4;
    theta = .4*exp(-10.*(x-pi).^2);
    psi0 = amp0*exp(1i*theta);
elseif nmode==14 %For wkb comparison
    rr(1:lx/4+1) = 0.;
    pp(1:lx/4+1) = -pi/2;%-pi^3/8 + (pi./2).^3./3;                     % 0 <= x <= pi/2
    rr((lx/4)+2:3*lx/4) = pi^2/4 - (x((lx/4)+2:3*lx/4)-pi).^2;  % pi/2 < x < 3*pi/2
    pp((lx/4)+2:3*lx/4) =  x((lx/4)+2:3*lx/4)-pi; %pi^2/4.*(x((lx/4)+2:3*lx/4)-pi)...
                            %- (x((lx/4)+2:3*lx/4)-pi).^3./3;
    rr((3*lx/4)+1:lx) = 0.;                   
    pp((3*lx/4)+1:lx) = pi/2; %pi^3/8 - (pi./2).^3./3;                 % 3*pi/2 <= x < 2*pi
    pp = .1*pp;
    psi0 = rr.*exp(1i*pp);
    figure(101); clf;
    k = make_k(lx);
    rr = abs(psi0);
    r2k = imag(conj(psi0).*deriv(psi0,k));
    plot(x, rr, 'b', x, r2k, 'g');
    xlabel('x'); ylabel('Blue = r, Green = r^2 k');
    axis tight;  
elseif nmode==15
    amp0= .5;
    alpha = 1./3.;
    psi0(1:lx/4+1) = 0.;                           % 0 <= x <= pi/2
    yy = x((lx/4)+2:3*lx/4)-pi/2;
    zz = 3*pi/2-x((lx/4)+2:3*lx/4);
    a = alpha; % + 0.2*yy.*zz;
    rr0 = amp0.*(yy.*zz/4).^a;  % pi/2 < x < 3*pi/2
    phi0=0.*yy;
    %phi0 = -log(yy.*zz);
    %phi0 = (yy.*zz).^(1.-3*alpha);
    psi0((lx/4)+2:3*lx/4) = rr0.*exp(i*phi0);
    psi0((3*lx/4)+1:lx) = 0.;    % 3*pi/2 <=x < 2*pi
elseif nmode==16
    amp0= .5;
    alpha = 1/3;
    psi0(1:lx/4+1) = 0.;                           % 0 <= x <= pi/2
    yy = x((lx/4)+2:3*lx/4)-pi/2;
    zz = 3*pi/2-x((lx/4)+2:3*lx/4);
    a = alpha; % + 0.2*yy.*zz;
    rr0 = amp0.*(yy.*zz/4).^a;  % pi/2 < x < 3*pi/2
    phi0 = 0; %log(yy.*zz);
    psi0((lx/4)+2:3*lx/4) = rr0.*exp(2i*phi0);
    psi0((3*lx/4)+1:lx) = 0.;    % 3*pi/2 <=x < 2*pi    
elseif nmode==17
    psi0(1:lx) = 0.;               
    yy = x((lx/4)+2:3*lx/4)-pi/2;
    zz = yy.^2- pi^2;
    psi0((lx/4)+2:3*lx/4) = yy.*sqrt(zz);
    psi0 = psi0./max(psi0);
elseif nmode==18
    psi0(1:lx) = 0.;
    yy = 3*pi/2-x((lx/4)+2:3*lx/4);
    alpha = 1./2;
    %beta = 1 - 2*alpha; %r^2 k = 1
    %beta = 1 - 3*alpha; %r^3 k = 1
    %phi0 = yy.^beta.*(pi/2-yy); %log(yy);
    phi0 = 0.;
    rr0 = yy.^alpha; %sqrt(yy)
    rr0 = rr0.*exp(-1./(pi-yy))/2;
    %phi0 = (pi/2-yy)./rr0.^2; %log(yy);
    %ww = 0.*yy + 0.2;
    %rr0 = bsxfun(@min,rr0,ww);
    psi0((lx/4)+2:3*lx/4) = rr0.*exp(1i*phi0);
elseif nmode==19
    amp1 = 1.;
    %psi0 = amp1.*exp(1i.*x) + amp2.*exp(2i.*x);
    psi0 = amp1*cos(x);
elseif nmode==20
    c = 10^6;
    psi0 = exp(-c * x.^2);
elseif nmode==21
    c = 5;
    nu = 1;
    % define waveNum, r, and phi
    waveNum = zeros(1, lx);
    r = zeros(1, lx);
    phi = zeros(1, lx);
    waveNum(lx/2+1:3*lx/4) = sqrt(c/2/nu) .* ((x(3*lx/4) - x((lx/2)+1:3*lx/4)) .^(-1/2));
    waveNum(lx/4+1:lx/2) = flip(waveNum(lx/2+1:3*lx/4));
    r(lx/4+1:3*lx/4) = sqrt(c ./ waveNum(lx/4+1:3*lx/4));
    phi(lx/2+1:3*lx/4) = -2 * sqrt(c/2/nu) .* ((x(3*lx/4) - x((lx/2)+1:3*lx/4)) .^(1/2));
    %phi(lx/2+1:3*lx/4) = waveNum(lx/2+1:3*lx/4) .* x(lx/2+1:3*lx/4);
    phi(lx/4+1:lx/2) = flip(phi(lx/2+1:3*lx/4));
    psi0(1:lx/4) = 0.;                         
    psi0((lx/4)+1:3*lx/4) = epsilon .* r(lx/4+1:3*lx/4) .* exp(1i/(epsilon^2) .* phi(lx/4+1:3*lx/4));
    psi0((3*lx/4)+1:lx) = 0;
elseif nmode==22 % another traveling wave solution
    psi0 = zeros(size(x));
    for j=1:size(psi0, 2)
        psi0(j) = NLStravelingSolu(x(j), 0);
    end
elseif nmode == 23 % the asymptotic solution of MRS initial value
    ic_mode = 1;
    k = 0; % placeholder for frequency cutoff
    uv = mrs_uzero(x, k, length(x), ic_mode, 1, 1);
    u0 = uv(1 : length(uv)/2);
    v0 = uv(length(uv)/2+1 : length(uv));
    psi0 = 1/sqrt(6) .* v0 - 1i/sqrt(6) .* u0;
    %psi0 = 1/sqrt(6) .* u0 + 1i/sqrt(6) .* v0;
else
    amp0=.2;
    psi0 = amp0*(...
           rand()*exp(1i*(x - 2*pi*rand())) + ...
           rand()*exp(-1i*(x - 2*pi*rand() )) + ...
           rand()*exp(2*1i*(x - 2*pi*rand())) + ...
           rand()*exp(-2*1i*(x - 2*pi*rand() ))     );
end
end

function xi = ftw(theta,c,s)
% XI finds implicit variable xi from theta in the traveling wave solution
myfun = @(xi,theta,c,s) (c^2+s^2).*xi - 2*c*s*sin(xi)-theta; % function
fun = @(xi) myfun(xi,theta,c,s);
xi = fzero(fun,theta);
end
