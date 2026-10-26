function make_video(M, ylim_vals, nomeFile, fps)

%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRI DA VARIARE SE NECESSARIO %%%%%%%%%%%%%%%

t = linspace(0, 10, 201).';
x = linspace(0, 1, 102).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Gestione dei parametri opzionali
if nargin < 3 || isempty(nomeFile)
    nomeFile = 'video_output.mp4';
end
if nargin < 4 || isempty(fps)
    fps = 30;
end

% Estrazione dimensioni della matrice
[~, numColonne] = size(M);

% Configurazione VideoWriter
v = VideoWriter(nomeFile, 'MPEG-4');
v.FrameRate = fps;
open(v);

% Creazione della figura (Visible off per renderla più veloce in background)
hFig = figure('Visible', 'off');
hPlot = plot(x, M(:, 1), 'LineWidth', 2, 'Color', 'b');
grid on;
xlabel('');
% ylabel('YLABEL');
xlim([x(1) x(end)]);
ylim(ylim_vals);

% Ciclo di generazione frame
for j = 1:numColonne
    set(hPlot, 'YData', M(:, j));
    title(sprintf('Tempo: t = %.2f s', t(j)), Interpreter="latex");  % DA CAMBIARE SE SERVE

    frame = getframe(hFig);
    writeVideo(v, frame);
end

% Chiusura e salvataggio
close(v);
close(hFig);

fprintf('File salvato con successo: %s\n', nomeFile);
end