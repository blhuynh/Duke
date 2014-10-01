clear; clc

a = fileload('part6.txt');
a = reshape(a,2,numel(a)/2)';

figure(1); clf; hold on
plot(log10(a(:,1)),a(:,2),'k-x','MarkerSize',5);
xlabel('log(Distance) (\mum)')
ylabel('Threshold (mA)')
title('Threshold-Distance Relationship')
setfont

print -dpng part6