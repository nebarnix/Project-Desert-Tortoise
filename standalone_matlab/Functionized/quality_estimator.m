%Calc a quality point every 100 samples
clear phase_estimatori;
alpha = 0.00005;
phase_estimator = 1.5708;
idx2 = 1;
for idx=1:numel(d_phasei)
    
    phase_estimator = phase_estimator * (1.0 - alpha) + alpha * abs(d_phasei(idx));
    %phase_estimatori(idx) = 10 * log10( (1.5708 - phase_estimator) ^ 2);
    if(mod(idx,100) == 0)
        
        phase_estimatori(idx2) = 10 * log10( (1.5708 - phase_estimator) ^ 2);
        idx2 = idx2 + 1;
    end
end
%% Plot

figure(20);
axisHandles = plotyy(BasebandRawTime(1:100:end-100),phase_estimatori,frameTime(:,1),minorFrameID);

axis(axisHandles(1),[0 max(BasebandRawTime) -15 0]);
axis(axisHandles(2),[0 max(BasebandRawTime) 0 600]);
%get(gca,'Children')
%get(axisHandles(1),'Children')
%set(axisHandles(2), 'tick', 1);
set(axisHandles(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

set(get(axisHandles(2),'Children'), 'MarkerSize', 1);
set(get(axisHandles(2),'Children'), 'Marker', '.');
set(get(axisHandles(2),'Children'), 'lineStyle',  'none');
set(axisHandles(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
%axis([0 900 0 2]);