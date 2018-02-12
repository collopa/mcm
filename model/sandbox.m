function sandbox()

% Variables
keys = {'A1', 'omega1', 'delta1', 'B1', ...
        'A2', 'omega2', 'delta2', 'B2'};
values = [1, 1, 0, 0, ...
          0.5, 10, pi/3, 0];
dict = containers.Map(keys,values);
max_time = 5*pi; % how long is sim?
ns = 100; % how many samples are we taking?
Fs = max_time/ns; % sample frequency
L = ns + 1;
distance = 5*pi; % how far our wave travels in the frame

% Time domain components
t = (0:Fs:max_time);
z = (0:2*pi/100:2*pi);
phi1 = dict('A1') * sin(dict('omega1') * t + dict('delta1') -z) + dict('B1');
phi2 = dict('A2') * sin(dict('omega2') * t + dict('delta2') -z) + dict('B2');

% Signals
signal = phi1 + phi2;
omega = -ns/2:ns/L:ns/2-ns/L;
FTsig = fft(signal);
IFTsig = ifft(FTsig);

% Plotting
figure(1)
plot(t, signal, 'o', t, real(IFTsig), 'o');
hold on
figure(2)
stem(omega, abs(FTsig));
figure(3)
plot(t, signal, ...
     t, exp(-1)*real(IFTsig), ...
     t, exp(-2)*real(IFTsig), ...
     t, exp(-3)*real(IFTsig));
figure(4)
plot(z, exp(-0.5*z).*real(IFTsig), 'o');


figure(5)
% Plot the first frame:
carla = sin(t-z);
disp(carla)
joy = exp(-z);
h = plot(t,carla*joy(1));

axis([0,distance,-dict('A1')*5,dict('A1')*5])

gif('test.gif','DelayTime',0.2,'frame',gcf)

for k = 2:ns
   set(h,'Ydata',carla*joy(k))
   gif
end


end