clear, close all;
filename = 'three_sinusoids.mp3'; % Путь к аудиофайлу

[y, Fs] = audioread(filename); % Загрузка аудиоданных и частоты дискретизации

% Вывод информации о звуковом файле
fprintf('File: %s\n', filename);
fprintf('Duration: %.2f seconds\n', length(y) / Fs);
fprintf('Sampling rate: %d Hz\n', Fs);
fprintf('Number of channels: %d\n', size(y, 2));

% Построение графика аудиоданных
t = (0:length(y)-1) / Fs; % Временная шкала
figure;
subplot(2,1,1)
plot(t, y); grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Input');
xlim([0 2.5]); % Ограничение оси X


filename = 'output.wav'; % Путь к аудиофайлу

[x, fs] = audioread(filename); % Загрузка аудиоданных и частоты дискретизации

% Вывод информации о звуковом файле
fprintf('File: %s\n', filename);
fprintf('Duration: %.2f seconds\n', length(y) / fs);
fprintf('Sampling rate: %d Hz\n', fs);
fprintf('Number of channels: %d\n', size(y, 2));

t = (0:length(x)-1) / fs; % Временная шкала

subplot(2,1,2)
plot(t, x); grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Output');
xlim([0 2.5]); % Ограничение оси X