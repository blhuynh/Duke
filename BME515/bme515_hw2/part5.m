clear; clc

a = fileload('part5.txt');
a = reshape(a,2,numel(a)/2)';

figure(1); clf; hold on
plot(a(:,1),a(:,2),'k-x','MarkerSize',5);
xlabel('Fiber Diameter (\mum)')
ylabel('Threshold (mA)')
title('Threshold-Fiber Diameter Relationship')
setfont

print -dpng part5