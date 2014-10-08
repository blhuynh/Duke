clear; clc
t = round(rand(1,30));

amp = 1;
tau = 0.5;
f = zeros(size(t));

for k = 1:length(t)
    if t(k) == 1
        expvec = exp(-[1:length(t)]/tau);
        veclength = length(t)-k+1;
        expvec = expvec(1:veclength);
        f(k:end) = f(k:end)+expvec;
    end
end

figure(1); clf; hold on
plot(1:length(t),f);
vline(find(t==1))