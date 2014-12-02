D1 = normrnd(9, 2, 100, 1);
D2 = normrnd(3, 1, 100, 1);
Dvalues = sort([D1' D2']');
bincentres = 0:0.5:15;
Dvalues_hist = hist(Dvalues,bincentres);
figure
plot(bincentres,Dvalues_hist,'LineWidth',3)
set(gca,'FontSize',14)
xlabel('D (um)','FontSize',14)
ylabel('Count','FontSize',14)
title('Histogram of D values for compound nerve','FontSize',14)

save('Dvalues.mat','Dvalues','D1','D2')