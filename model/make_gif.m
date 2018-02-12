function make_gif()

input_signal = csvread('input_signal.csv');
zI = input_signal(1,:);
tI = input_signal(2,:);
f_tI = input_signal(3,:);
disp(f_tI)

FT_input_signal = csvread('FT_input_signal.csv');
kI = FT_input_signal(1,:);
f_kI = FT_input_signal(2,:);

FT_output_signal = csvread('FT_output_signal.csv');
kO = FT_output_signal(1,:);
f_kO = FT_output_signal(2,:);
kappa = FT_output_signal(3,:);

output_signal = csvread('radio_output.csv');
zO = output_signal(1,:);
tO = output_signal(2,:);
f_tO = output_signal(3,:);

% Useful extracted variables
maxZ = max(zI);
minAmp = min([min(f_tI),min(f_tO)]);
maxAmp = max([max(f_tI),max(f_tO)]);

% Plotting
figure(1)
plot(zI, f_tI, 'o', zO, f_tO, 'o')
title('Compare Input and Output in Time Domain')
xlabel('Time (s)')
ylabel('Amplitude (V/m)')
hold on

figure(2)
stem(kI, abs(f_kI))
hold on
stem(kO, abs(f_kO))
title('Compare Input and Output in Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude (V/m)')


% Plot the first frame:
t = exp(-kappa.*zI);

figure(3)
input = plot(zI,f_tI*t(1));
axis([0,maxZ,minAmp*3, maxAmp*3])

gif('results.gif','DelayTime',2,'frame',gcf)

for k = 2:length(f_tI)
   set(input,'Ydata',f_tI*t(k))
   gif
end

end