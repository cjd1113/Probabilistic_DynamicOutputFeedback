
clear all

% AverageElapsed = Elapsed/10

% AverageElapsed = 9.4;
AverageElapsed = 0.57;
% AverageElapsed = Elapsed_musyn/4

Nmcepotential = [1:1:1e6];

CPUHours = Nmcepotential*AverageElapsed*(1/3600);

d = 0.05;
Accuracy = sqrt(1./(2*Nmcepotential./log(2/d)));

figure(1)
yyaxis left
loglog(Nmcepotential,CPUHours,'linewidth',2)
ylabel('CPU Hours','interpreter','latex','fontsize',14)
yyaxis right
loglog(Nmcepotential,Accuracy,'linewidth',2)
ylim([0 1])
ylabel('Probability Estimate Accuracy, 95\% Confidence','interpreter','latex','fontsize',14)
xlabel('Number of Monte Carlo Samples','interpreter','latex','fontsize',14)
title('Raw CPU Time vs. Monte Carlo Samples with Chernoff Bound Accuracy','interpreter','latex','fontsize',16)

numCores = [12 24 48 96];
for i = 1:length(numCores)
    ActualTime(i,:) = CPUHours./numCores(i);
end



figure(2)
hold on
for i = 1:length(numCores)
    yyaxis left
    plot(Nmcepotential,ActualTime(i,:),'-','linewidth',2)
end
hold off
yax = gca;
set(yax,'YScale','log')
% set(gca,'loglog')
ylabel('Approximated Actual Run Time, Hours','interpreter','latex','fontsize',16)
yyaxis right
loglog(Nmcepotential,Accuracy,'linewidth',2)
ylim([0 1])
ylabel('Probability Estimate Accuracy, 95\% Confidence','interpreter','latex','fontsize',16)
xlabel('Number of Monte Carlo Samples','interpreter','latex','fontsize',16)
str1 = sprintf('Enhanced Cost Function Execution Time vs. Monte Carlo Samples');
str2 = sprintf('with Chernoff Bound Accuracy, Effect of Parallelization');
title('\begin{tabular}{c} Enhanced Cost Function Execution Time vs. Monte Carlo Samples \\ with Chernoff Bound Accuracy, Effect of Parallelization \end{tabular}','interpreter','latex','fontsize',18)
% title({str1,str2},'interpreter','latex','fontsize',18)
% title('Enhanced Cost Function Execution Time vs. Monte Carlo Samples with Chernoff Bound Accuracy, Effect of Parallelization','interpreter','latex','fontsize',18)

x1 = 300;
y1 = 0.01;
txt1 = '\leftarrow 96 Cores';
text(x1,y1,txt1)

x2 = 3;
y2 = 0.01;
txt2 = '12 cores \rightarrow';
text(x2,y2,txt2)


