for j=1:4096,
plot(z(:,1),t_focus(j,1:151));
axis([min(z(:,1)) max(z(:,1)) -300 300]);
hold on;
plot(-(j-150)/fs*1500*[1 1],[-300 300],'r');
plot(-(j-150-38*2/1500*fs)/fs*1500*[1 1],[-300 300],'r');
plot((j-150-(117.5*2-27.5*2)/1500*fs)/fs*1500*[1 1],[-300 300],'r');
hold off;
drawnow;
end
