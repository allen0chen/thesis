function uv = mrs_uzero(x, k, lx, ic_mode, ampu0, ampv0)
% MRS_UZERO defines initial condition  u0 = u(x,0), v0 = v(x,0)

if ic_mode == 1
    % u=v=0 at x=pi/2
    delta = 0.1;
    u = ampu0*cos(x).*exp(sin(x));
    v = ampv0*(exp(sin(x-pi/4))-exp(1/sqrt(2)));
    u = u + delta*ampu0*sin(10*x);
elseif ic_mode == 2
    % doesn't have common zero
    delta = 0.1;
    u = -ampu0*cos(x).*exp(sin(x));
    v = ampv0*(exp(sin(x-pi/4))-exp(1/sqrt(2)));
    u = u + delta*ampu0*cos(10*x);
elseif ic_mode == 3
    % u is sum of two cosine functions, v is 0
    u = 2*cos(x) + cos(2*(x +2*pi^2));
    u = ampu0*u;
    v = zeros(1,lx);
elseif ic_mode == 4
    % straight line: u=v=x
    u = x;
    v = x;
elseif ic_mode == 5
    % correspond to a DQS initial data
    delta = 1;
    nmode = 10;
    A = ic(x, nmode, lx, delta);
    u = ampu0 * -sqrt(6) * imag(A);
    v = ampv0 * sqrt(6) * real(A);
elseif ic_mode == 6
    % two shocks at x=pi/2 and x=3*pi/2
    % u = -1 + (x-pi/2)/pi when pi/2<x<3*pi/2
    u = zeros(1, lx);
    for j=(lx/4+1):(3*lx/4)
        u(1, j) = ampu0 * (-1 + 4*(j-lx/4-1)/lx);
    end
    v = 0 .* x;
elseif ic_mode == 7
    % two shocks at x=pi/2 and x=3*pi/2
    % u = 1 - (x-pi/2)/pi when pi/2<x<3*pi/2
    u = zeros(1, lx);
    for j=(lx/4+1):(3*lx/4)
        u(1, j) = ampu0 * (1 - 4*(j-lx/4-1)/lx);
    end
    v = 0 .* x;
elseif ic_mode == 8
    % similar to case 6&7, but piecewise continuous
    u = zeros(1, lx);
    for j=(lx/4+1):(3*lx/4)
        u(1, j) = ampu0 * (1 - 4*(j-lx/4-1)/lx);
    end
    for j=(lx/8+1):(lx/4)
        u(1, j) = ampu0 * (j-lx/8-1) * 8/lx;
    end
    for j=(3*lx/4+1):(7*lx/8)
        u(1, j) = ampu0 * (-1 + (j-3*lx/4-1) * 8/lx);
    end
    v = 0 .* x;
elseif ic_mode == 9
    % stationary shock, has a shock at x = pi, u=eps when x<pi, u=-eps when x>pi
    u = zeros(1, lx);
    for j=1:(lx/2)
        u(1, j) = ampu0;
    end
    for j=(lx/2+1):lx
        u(1, j) = -1*ampu0;
    end
    v = 0 .* x;
elseif ic_mode == 10
    % has a shock at x = pi, but not symmetric. u=eps when x<pi, u=0 when x>pi
    u = zeros(1, lx);
    for j=1:(lx/2)
        u(1, j) = ampu0;
    end
    v = 0 .* x;
elseif ic_mode == 11
    % has a shock at x = pi, but piecewise continuous, i.e. decreases
    % to zero
    u = zeros(1, lx);
    for j=(lx/4+1):(lx/2)
        u(1, j) = ampu0 * (j-lx/4) * 4/lx;
    end
    for j=(lx/2+1):(3*lx/4)
        u(1, j) = ampu0 * (-1 + (j-lx/2-1) * 4/lx);
    end
    v = 0 .* x;
elseif ic_mode == 12
    % same as case 11, but with u, v swapped
    v = zeros(1, lx);
    for j=(lx/4+1):(lx/2)
        v(1, j) = ampv0 * (j-lx/4) * 4/lx;
    end
    for j=(lx/2+1):(3*lx/4)
        v(1, j) = ampv0 * (-1 + (j-lx/2-1) * 4/lx);
    end
    u = 0 .* x;
elseif ic_mode == 13
    % mollified version of case 11
    alpha = 4;
    temp_piece1 = @(y) 2*(y-pi/2)/pi;
    temp_piece2 = @(y) 2*(y-3*pi/2)/pi;
    u = ampu0.* (mollification(temp_piece1, pi/2, pi, x, lx, pi/alpha) + mollification(temp_piece2, pi, 3*pi/2, x, lx, pi/alpha));
    v = 0 .* x;
elseif ic_mode == 14
    % mollified version of case 6
    alpha = 16;
    epsilon = pi/alpha;
    u_smooth = zeros(1, lx);
    for j=1:lx
        if (x(j) <= pi/2-epsilon) || (x(j) >= pi+epsilon)
            u_smooth(j) = 0;
        elseif (pi/2-epsilon < x(j)) && (x(j) <= pi/2+epsilon)
            mullification = @(y) 1/2/epsilon .* exp(-1 ./ (1-abs((x(j)-y)./epsilon).^2)) .* sin(4*y);
            u_smooth(j) = integral(mullification, pi/2, x(j)+epsilon);
        elseif (pi-epsilon <= x(j)) && (x(j) < pi+epsilon)
            mullification = @(y) 1/2/epsilon .* exp(-1 ./ (1-abs((x(j)-y)./epsilon).^2)) .* sin(4*y);
            u_smooth(j) = integral(mullification, x(j)-epsilon, pi);
        else
            mullification = @(y) 1/2/epsilon .* exp(-1 ./ (1-abs((x(j)-y)./epsilon).^2)) .* sin(4*y);
            u_smooth(j) = integral(mullification, x(j)-epsilon, x(j)+epsilon);
        end
    end
    u = ampu0 * u_smooth;
    v = 0 .* x;
elseif ic_mode == 15
    % testing with local perturabation on u, u is continuous, but not C
    % infinity
    delta = 0.1;
    u = ampu0*cos(x).*exp(sin(x));
    v = ampv0*(exp(sin(x-pi/4))-exp(1/sqrt(2)));
    % Perturbation is centeralized around zero of u0 and v0, which is
    % located at x=pi/2, or x(lx/4+1)
    alpha = 4;
    localizer = zeros(1, lx);
    for j=(lx/4+1-lx/2/alpha):(lx/4+1+lx/2/alpha)
        localizer(j) = 1;
    end
    pert = delta/alpha * localizer .* cos(10*alpha*x);
    u = u + ampu0*pert;
elseif ic_mode == 16
    % mollified version case 15
    delta = 0.1;
    u = ampu0*cos(x).*exp(sin(x));
    v = ampv0*(exp(sin(x-pi/4))-exp(1/sqrt(2)));
    % Perturbation is centeralized around zero of u0 and v0, which is
    % located at x=pi/2, or x(lx/4+1)
    alpha = 4;
    temp = @(y) cos(10*alpha*y);
    pert = delta/alpha .* mollification(temp, pi/2-pi/alpha, pi/2+pi/alpha, x, lx, pi/alpha/4);
    figure(333)
    plot(x, u)
    u = u + ampu0*pert;
elseif ic_mode == 17
    u = ampu0*cos(x);
    v = zeros(1,lx);
elseif ic_mode == 18
    u = zeros(1,lx);
    v = ampu0*cos(x);
elseif ic_mode == 19
    % correspond to a DQS initial data, but does not use dqs ic file
    u = zeros(1, lx);
    u((lx/4)+1:3*lx/4) = ampu0 .* (pi^2/4 - (x((lx/4)+1:3*lx/4)-pi).^2).^2;
    v = zeros(1, lx);
elseif ic_mode == 20
    % correspond to a DQS initial data, but does not use dqs ic file
    u = zeros(1, lx);
    v = zeros(1, lx);
    v((lx/4)+1:3*lx/4) = ampv0 * sqrt(6) .* (pi^2/4 - (x((lx/4)+1:3*lx/4)-pi).^2).^2;
else
    disp('No IC specified')
end

figure(222)
plot(x, u, 'blue')
hold on
plot(x, v, 'red')
axis('tight')

uv = [u v];

end