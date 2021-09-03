function [] = plot_RLCG(data1,Ro, Lo, Co, Go)
freq =1e-9.* data1(:,1);
figure
subplot(2,2,1)
plot(freq, Ro);
ylabel('Ro (Ohm/m)'), xlabel('Freq (GHz)');
grid on,axis tight;
subplot(2,2,2)
plot(freq, Lo);
ylabel('Lo (nH/m)'), xlabel('Freq (GHz)');
grid on,axis tight;
subplot(2,2,3)
plot(freq, (10^3)*Co);
ylabel('Co (pF/m)'), xlabel('Freq (GHz)');
grid on,axis tight;
subplot(2,2,4)
plot(freq, Go);
ylabel('Go (Mho/m)'), xlabel('Freq (GHz)');
grid on,axis tight;
end