% After running Fig 3 create schematic - contains panels of Fig 1

figure; hold on;
plot(sim2.tday, sim2.Iday, 'k', 'LineWidth', 2); 
plot(sim2.tday, sim2.Lday, 'r', 'LineWidth', 2); 
plot(sim2.tday, sim2.Ipmean, 'b', 'LineWidth', 2); 
box off; grid off; xlim(t2); hold off;
%ylabel('$I_t$', 'FontSize', fnt);
%xlabel('$t$ (d)', 'FontSize', fnt);


figure; hold on;
plot(sim2.tday, sim2.wch(1,:), 'Color', grey2, 'LineWidth', 2); 
plot(sim2.tday, sim2.wch(4,:), 'Color', grey1, 'LineWidth', 2); 
box off; grid off; xlim([1 30]); hold off;
%ylabel('$w_u$', 'FontSize', fnt);
%xlabel('$u$ (d)', 'FontSize', fnt);


figure; hold on;
plot(sim2.tday, est2.Rmean-1, 'r', 'LineWidth', 2); 
plot(sim2.tday, est2.Ommean-1, 'b', 'LineWidth', 2); 
box off; grid off; xlim(t2); hold off;
