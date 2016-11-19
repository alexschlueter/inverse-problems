function antideriv()
    y = @(x) (1-cos(2*pi*x))/2/pi;
    u = @(x) sin(2*pi*x);
    alphas = [0.2, 0.1, 0.01];
    delta = 0.2;
    k = 30;
    ydelta = @(x) y(x) + sqrt(2)*delta*sin(2*pi*k*x);
    x = linspace(0, 1, 100);
    plot(x, y(x));
    figure, plot(x, ydelta(x));
    figure, subplot(2, 2, 1), plot(x, u(x));
    for i = 1:length(alphas)
        subplot(2, 2, i+1), plot(x, inv(alphas(i), ydelta, x));
    end
end

function inv = inv(alpha, y, x)
    inv = zeros(1, length(x));
    n = 1;
    while sigma(n) > alpha
        inv = inv + 1/sigma(n)*integral(@(t) y(t).*vn(n, t), 0, 1)*un(n, x);
        n = n+1;
    end
end

function sigma = sigma(n)
    sigma = 2/(2*n-1)/pi;
end

function un = un(n, x)
    un = sqrt(2)*cos(x/sigma(n));
end

function vn = vn(n, x)
    vn = sqrt(2)*sin(x/sigma(n));
end