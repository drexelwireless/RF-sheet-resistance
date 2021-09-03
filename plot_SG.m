function [] = plot_SG(freq,gamma_t, gammar, s11r, s12r, s11d, s12d)
figure
% plot(freq, real(gamma_t));
alpha = real(gamma_t);
% alpha= movmean(alpha, 20);
% alpha = smooth(freq,alpha,1,'rloess');
plot(freq, alpha);
hold on;

ylabel('Attenuation Constant (Np/m)');
hold on;
plot(freq(1:10:end), real(gammar(1:10:end)),'ms','MarkerSize',5,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);
% plot(freq, real(gammar),'m--');
% legend('\alpha from RLCG', '\alpha measured','location','best');

yyaxis right;
plot(freq, imag(gamma_t),'r-.');
hold on;
plot(freq(1:10:end), imag(gammar(1:10:end)),'kd','MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[1 .6 .6]);

xlabel('Frequency (GHz)'), ylabel('Phase Constant (rad/m)')
legend('\alpha from RLCG', '\alpha measured','\beta from RLCG', ...
    '\beta measured','location','best');

figure
subplot(2,2,1)
plot(freq, 20*log10(abs(s12d)),'k--');
hold on, grid on;
plot(freq, 20*log10(abs(s12r)));
legend('S21 measured','S21 reconstructed','Location','Best');
xlabel('Frequency (GHz)'), ylabel('S12 (dB)');

subplot(2,2,2)
plot(freq, angle(s12d),'k--');
hold on, grid on;
xlabel('Freq (GHz)'), ylabel('S21 (angle)')
plot(freq, angle(s12r));
ylim([-2*pi, 2*pi]);
legend('S12 ang measured', 'S12 ang reconst','location', 'best');


subplot(2,2,3)
plot(freq, 20*log10(abs(s11d)),'k--');
hold on;
xlabel('Frequency (GHz)'), ylabel('S11 (dB)');
plot(freq, 20*log10(abs(s11r)));
legend('S11 measured', 'S11 reconst','location', 'best');

subplot(2,2,4)
plot(freq, angle(s11d),'k--');
hold on, grid on;
ylim([-2*pi, 2*pi]);
xlabel('Frequency (GHz)'), ylabel('S11 (angle)')
plot(freq, angle(s11r));
legend('S11 ang measured', 'S11 ang reconst','location', 'best');
end
