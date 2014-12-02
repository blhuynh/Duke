Fvalues = 10 + (20-10).*rand(100,1);

bincentres = 8:0.5:22;
Fvalues_hist = hist(Fvalues(1:100),bincentres);
figure
plot(bincentres,Fvalues_hist,'LineWidth',3)
set(gca,'FontSize',14)
xlabel('FR (Hz)','FontSize',14)
ylabel('Count','FontSize',14)
title('Histogram of firing rates for compound nerve','FontSize',14)

save('Fvalues.mat','Fvalues')