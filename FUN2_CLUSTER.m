label = zeros(1000, 1);
mu_new = mu;
eps = 1e-6;
delta = 1;
while (delta > eps)
    mu = mu_new;
    for i =1:1000
        y = repmat (X(i, :), k, 1);
        dist = y - mu;
        d = sum(dist.*dist,2);
        j = find(d==min(d));
        label(i) = j;
    end
    for j = 1 : k
        order = find(label == j);
        mu_new(j, :) = mean(X(order, :), 1);
    end
    delta = sqrt(sum(sum((mu-mu_new).*(mu-mu_new))));
end
label = zeros(1000, 1);
for i = 1 : 1000
    R = repmat(X(i,:),k,1) - mu;
    Residual = sum(R.*R,2);
    j = find(Residual == min(Residual));
    label(i) = j;
end
% Construct map function
s = zeros(k, 1);
for j =1 : k
    order = find(label==j);
    Y = X_label(order);
    s(j) = mode(Y);
end
map_label =zeros(1000, 1);
for j = 1 : k
    map_label(label==j) = s(j);
end
figure;
hold on;
for i =1:1000
    if map_label(i)==1
        plot(X(i,1),X(i,2),'r.');
    elseif map_label(i)==2
        plot(X(i,1),X(i,2),'b.');
    elseif map_label(i)==3
        plot(X(i,1),X(i,2),'k.');
    elseif map_label(i)==4
        plot(X(i,1),X(i,2),'g.');
    else
        plot(X(i,1),X(i,2),'m.');
    end
end
% show the cluster center
for i = 1 : 5
    plot(mu(i,1),mu(i,2),'yo','LineWidth',3);
end
