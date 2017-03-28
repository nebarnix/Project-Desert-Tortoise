 
% 	double y1 = abs(dataStreamIn(idx))*abs(dataStreamIn(idx));
% 	d_y1 = alpha*y1 + beta*d_y1;
% 
% 	double y2 = abs(dataStreamIn(idx))*abs(dataStreamIn(idx))*abs(dataStreamIn(idx))*abs(dataStreamIn(idx));
% 	d_y2 = alpha*y2 + beta*d_y2;
%       }
%       return noutput_items;
%     }
% 
%     double
%     mpsk_snr_est_m2m4::snr()
%     {
%       double y1_2 = d_y1*d_y1;
%       d_signal = sqrt(2*y1_2 - d_y2);
%       d_noise = d_y1 - sqrt(2*y1_2 - d_y2);
%       return 10.0*log10(d_signal / d_noise);
%     }
%     
%     
clear SNR d_noise d_signal;
alpha = 0.0001;
beta = 1.0-alpha;
idx2 = 1;
d_y1 = 0;
d_y2 = 0;
for idx=1:numel(dataStreamLPF)
    y1 = abs(dataStreamLPF(idx)) .^ 2;
	d_y1 = alpha*y1 + beta*d_y1;
    
    y2 = abs(dataStreamLPF(idx)) .^ 4;
	d_y2 = alpha*y2 + beta*d_y2;
    
    %phase_estimator = phase_estimator * (1.0 - alpha) + alpha * abs(d_phasei(idx));
    %phase_estimatori(idx) = 10 * log10( (1.5708 - phase_estimator) ^ 2);
    if(mod(idx,100) == 0)
        
       y1_2 = d_y1^2;
       d_signal = sqrt(2*y1_2 - d_y2);
       d_noise = d_y1 - sqrt(2*y1_2 - d_y2);
       SNR(idx2) = 10.0*log10(d_signal / d_noise);
       idx2 = idx2 +1;
    end
end

%%
figure(13);
axisHandles = plotyy(BasebandRawTime(1:100:end-10),SNR,frameTime(:,1),minorFrameID);

axis(axisHandles(1),[0 max(BasebandRawTime) -5 18]);
axis(axisHandles(2),[0 max(BasebandRawTime) 0 600]);
%get(gca,'Children')
%get(axisHandles(1),'Children')
%set(axisHandles(2), 'tick', 1);
set(axisHandles(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')

set(get(axisHandles(2),'Children'), 'MarkerSize', 1);
set(get(axisHandles(2),'Children'), 'Marker', '.');
set(get(axisHandles(2),'Children'), 'lineStyle',  'none');
set(axisHandles(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')