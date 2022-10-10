clear; clc;

%% Plot Solution
% Read 'Xminus' and 'Xplus' returned by program from .txt files
Xminus = readmatrix('xminus.txt');%idx  = isnan(Xminus); Xminus(idx,:) = [];
Xplus = readmatrix('xplus.txt');%idx  = isnan(Xplus); Xplus(idx,:) = [];


% Draw the 'Xplus' region red
for i = 1:size(Xplus, 1)
    rectangle('Position', Xplus(i,:), 'FaceColor', [1 0.5 0.5],'EdgeColor', [0 0 0]);
end

% Draw the 'Xminus' region green (overlap Xplus region)
for i = 1:size(Xminus, 1)
    rectangle('Position', Xminus(i,:), 'FaceColor', [0.5 1 0.5],'EdgeColor', [0 0 0]);
end


% Limit the x- and y-axes
xlim([min([Xminus(:,1); Xplus(:,1)]) max([Xminus(:,1); Xplus(:,1)])]); 
ylim([min([Xminus(:,2); Xplus(:,2)]) max([Xminus(:,2); Xplus(:,2)])]);
grid on;

%% Plot Example
t = linspace(0.1, 5.1);
y = 5 * exp(-1 * t);
figure(2); plot(t, y, 'k'); hold on; plot(t, y * 0.5, 'g'); plot(t, y * 1.5, 'r'); hold off; grid on;
