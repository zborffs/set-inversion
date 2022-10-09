clear; clc;

% Read 'Xminus' and 'Xplus' returned by program from .txt files
Xminus = readmatrix('xminus.txt');
Xplus = readmatrix('xplus.txt');

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