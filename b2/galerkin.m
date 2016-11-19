function galerkin()
    global n
    n = 20;
    r = @(t) sin(3*pi*t);
    %r = @(t) (t-0.5).^2+0.2;
    y = yFactory(r);
    
    Mphi = zeros(n, n);
    %Mpsi = zeros(n, n);
    disp('Assembling matrices...')
    for i = 1:n
        for j = 1:i
            fprintf('i = %d j = %d\n', i, j);
            Mphi(i,j) = integral(@(x) phi(i, x).*phi(j, x), 0, 1);
            %Mpsi(i,j) = integral(@(x) psi(i, x).*psi(j, x), 0, 1, 'ArrayValued', true);
        end
    end
    Mphi = Mphi + Mphi' - diag(diag(Mphi));
    %Mpsi = Mpsi + Mpsi' - diag(diag(Mpsi));
    Mpsi = eye(n, n)/n;
    
    f = zeros(n, 1);
    
    disp('Assembling right side...')
    for i = 1:n
        fprintf('i = %d\n', i);
        f(i) = integral(@(x) y(x).*phi(i, x), 0, 1);
    end
    
    disp('Solving linear equations...');
    c = Mphi\f;
    lambda = Mpsi\c;
    
    x = linspace(0, 0.999, 100);
    disp('Plotting results...')
    plot(x, r(x))   
    figure, plot(x, y(x))
    figure, plot(x, rsol(lambda, x))
end

function k = k(x, t)
    k = zeros(size(x));
    idx = x < t;
    k(idx) = x(idx).*(1-t);
    k(~idx) = t.*(1-x(~idx));
end

function y = yFactory(r)
    y = @(x) integral(@(t) k(x, t).*r(t), 0, 1, 'ArrayValued', true);
end

function phi = phi(i, x)
global n
    phi = n*integral(@(t) k(x, t), (i-1)/n, i/n, 'ArrayValued', true);
end

function psi = psi(i, t)
    global n
    psi = (i-1)/n <= t && t < i/n;
end

function rsol = rsol(lambda, t)
    global n
    Apsi = zeros(size(t, 2), n);
    for i = 1:size(t, 2)
        for j = 1:n
            Apsi(i, j) = psi(j, t(i));
        end
    end
    rsol = Apsi*lambda;
end