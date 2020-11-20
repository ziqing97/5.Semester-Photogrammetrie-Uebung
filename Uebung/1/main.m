% Photogrametrie Uebung 1
% Ziqing Yu 3218051  Roman Germann
% 22/12/2019

clc;
clear all;
close all;


%% Daten lesen
% Punktkoordinaten in globale Koordinatensystem
fname = 'signals_loc.txt';
fileid = fopen(fname,'r');
punktmenge10 = [2563, 2583, 2593, 2893, 2913, 3004, 3006, 3007, 3008, 3009, 3026, 309019, 310019, 311019, 410019, 411029, 509019, 510019, 510029, 511019, 9001, 9002];
punktmenge11 = [2583, 2593, 2923, 3006, 3007, 3009, 310019, 311019, 312019, 410019, 411029, 411039, 412019, 413029, 510019, 510029, 511019, 512019, 9001, 9002];

Cell = {};
i=1;
while ~feof(fileid)
    str=fgetl(fileid);
    if ~isempty(str)
            Cell{i} = str;
            data = split(str);
            name(i, :) = str2double(data(1));
            x(i, :) = str2double(data(2));
            y(i, :) = str2double(data(3));
            z(i, :) = str2double(data(4));
            i=i+1;
    end
end

% Pixelkoordinaten auf dem Bild
size10 = length(punktmenge10);
size11 = length(punktmenge11);
koordinaten10 = zeros(size10,3);
koordinaten11 = zeros(size11,3);

for i = 1 : size10
    k = find(name == punktmenge10(i));
    koordinaten10(i,1) = x(k);
    koordinaten10(i,2) = y(k);
    koordinaten10(i,3) = z(k); 
end

for i = 1 : size11
    k = find(name == punktmenge11(i));
    koordinaten11(i,1) = x(k);
    koordinaten11(i,2) = y(k);
    koordinaten11(i,3) = z(k);
end

fclose(fileid);

i = 1;
pixelkoordinaten10 = zeros(size10, 2);
fname10 = 'bildkoordinaten10.txt';
fileid10 = fopen(fname10,'r');
while ~feof(fileid10)
    str10=fgetl(fileid10);
    if ~isempty(str10)
            Cell{i} = str10;
            data10 = split(str10);
            pixelkoordinaten10(i,1) = str2double(data10(2));
            pixelkoordinaten10(i,2) = str2double(data10(3));
            i=i+1;
    end
end

i = 1;
pixelkoordinaten11 = zeros(size11, 2);
fname11 = 'bildkoordinaten11.txt';
fileid11 = fopen(fname11,'r');
while ~feof(fileid11)
    str11 = fgetl(fileid11);
    if ~isempty(str11)
            Cell{i} = str11;
            data11 = split(str11);
            pixelkoordinaten11(i,1) = str2double(data11(2));
            pixelkoordinaten11(i,2) = str2double(data11(3));
            i=i+1;
    end
end
fclose('all');        % alle geoeffnete txt schliessen
clear x y z           % die brauche ich nicht mehr

%% Rechnung
c = 0.12;
 % Pixelkoordinaten zu Bildkoordinaten
R = [83 + 1 / 3, 0; 0, -(83 + 1 / 3 )];
T = [3839.5;6911.5];
bildkoordinaten10 = zeros(size10,2);
bildkoordinaten11 = zeros(size11,2);
for i = 1 : size10
    parameter = R \ (pixelkoordinaten10(i,:)' - T);
    bildkoordinaten10(i,1) = parameter(1);
    bildkoordinaten10(i,2) = parameter(2);
end
bildkoordinaten10 = bildkoordinaten10 / 1000;  % mm zu m

for i = 1 : size11
    parameter = R \ (pixelkoordinaten11(i,:)' - T);
    bildkoordinaten11(i,1) = parameter(1);
    bildkoordinaten11(i,2) = parameter(2);
end
bildkoordinaten11 = bildkoordinaten11 / 1000;  % mm zu m

% Bestimmung der Naeherungswerte fuer den RRS 
y10 = [bildkoordinaten10(:,1); bildkoordinaten10(:,2)];
A10 = [ones(size10,1),zeros(size10,1),koordinaten10(:,1),koordinaten10(:,2);
       zeros(size10,1),ones(size10,1),koordinaten10(:,2),-koordinaten10(:,1)];
x10 = (A10' * A10) \ A10' * y10;
kappa10 = atan2(x10(4),x10(3));
M10 = sqrt(x10(3)^2 + x10(4)^2);
X0_10 = (x10(2) * sin(kappa10) - x10(1) * cos(kappa10)) / M10;
Y0_10 = - (x10(1) * sin(kappa10) + x10(2) * cos(kappa10)) / M10;
z_10 = koordinaten10(:,3);
Z0_10 = mean(z_10)  + c / M10;


y11 = [bildkoordinaten11(:,1);bildkoordinaten11(:,2)];
A11 = [ones(size11,1),zeros(size11,1),koordinaten11(:,1),koordinaten11(:,2);
       zeros(size11,1),ones(size11,1),koordinaten11(:,2),-koordinaten11(:,1)];
x11 = (A11' * A11) \ A11' * y11;
kappa11 = atan2(x11(4),x11(3));
M11 = sqrt(x11(4)^2 + x11(3)^2);
X0_11 = (x11(2) * sin(kappa11) - x11(1) * cos(kappa11)) / M11;
Y0_11 = - (x11(1) * sin(kappa11) + x11(2) * cos(kappa11)) / M11;
z_11 = koordinaten11(:,3);
Z0_11 = sum(z_11) / size11 + c / M11;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

xi = [X0_10; Y0_10; Z0_10; 0; 0; kappa10];
% Berechnung des RRS durch Parameterschaetzung (alle Punkte) (Q ist Inversion der Normalgleichugnsmatrix)
[X0_10, Y0_10, Z0_10, omega10, phi10, kappa10, Q10] = RRS(c, X0_10, Y0_10, Z0_10, kappa10, bildkoordinaten10, koordinaten10);
[X0_11, Y0_11, Z0_11, omega11, phi11, kappa11, Q11] = RRS(c, X0_11, Y0_11, Z0_11, kappa11, bildkoordinaten11, koordinaten11);

% Berechnung des RRS durch Parameterschaetzung (4 Punkte in den Ecken)
% hier nehme ich Punkte 509019 511019 3026 311019 fuer Bild 20010010 und
% 510019 512019 310019 312019 fuer Bild 20010011
bildkoordinaten10_eck = [bildkoordinaten10(17,1),bildkoordinaten10(17,2);
                         bildkoordinaten10(20,1),bildkoordinaten10(20,2);
                         bildkoordinaten10(11,1),bildkoordinaten10(11,2);
                         bildkoordinaten10(14,1),bildkoordinaten10(14,2)];
koordinaten10_eck = [koordinaten10(17,1),koordinaten10(17,2),koordinaten10(17,3);
                     koordinaten10(20,1),koordinaten10(20,2),koordinaten10(20,3);
                     koordinaten10(11,1),koordinaten10(11,2),koordinaten10(11,3);
                     koordinaten10(14,1),koordinaten10(14,2),koordinaten10(14,3)];
bildkoordinaten11_eck = [bildkoordinaten11(15,1),bildkoordinaten11(15,2);
                         bildkoordinaten11(18,1),bildkoordinaten11(18,2);
                         bildkoordinaten11(7,1),bildkoordinaten11(7,2);
                         bildkoordinaten11(9,1),bildkoordinaten11(9,2)];
koordinaten11_eck = [koordinaten11(15,1),koordinaten11(15,2),koordinaten11(15,3);
                     koordinaten11(18,1),koordinaten11(18,2),koordinaten11(18,3);
                     koordinaten11(7,1),koordinaten11(7,2),koordinaten11(7,3);
                     koordinaten11(9,1),koordinaten11(9,2),koordinaten11(9,3)];         
                 
[X0_10_eck, Y0_10_eck, Z0_10_eck, omega10_eck, phi10_eck, kappa10_eck, Q10_eck] = RRS(c, X0_10, Y0_10, Z0_10, kappa10, bildkoordinaten10_eck, koordinaten10_eck);
[X0_11_eck, Y0_11_eck, Z0_11_eck, omega11_eck, phi11_eck, kappa11_eck, Q11_eck] = RRS(c, X0_11, Y0_11, Z0_11, kappa11, bildkoordinaten11_eck, koordinaten11_eck);

% Berechnung des RRS durch Parameterschaetzung (4 Punkte dicht liegen)
% hier nehme ich Punkte 3008 3007 9002 9001 fuer Bild 20010010 und 
% 2923 413029 411039 2583 fuer Bild 20010011
bildkoordinaten10_dicht = [bildkoordinaten10(8,1),bildkoordinaten10(8,2);
                         bildkoordinaten10(9,1),bildkoordinaten10(9,2);
                         bildkoordinaten10(21,1),bildkoordinaten10(21,2);
                         bildkoordinaten10(22,1),bildkoordinaten10(22,2)];
koordinaten10_dicht = [koordinaten10(8,1),koordinaten10(8,2),koordinaten10(8,3);
                     koordinaten10(9,1),koordinaten10(9,2),koordinaten10(9,3);
                     koordinaten10(21,1),koordinaten10(21,2),koordinaten10(21,3);
                     koordinaten10(22,1),koordinaten10(22,2),koordinaten10(22,3)];
[X0_10_dicht, Y0_10_dicht, Z0_10_dicht, omega10_dicht, phi10_dicht, kappa10_dicht, Q10_dicht] = RRS(c, X0_10, Y0_10, Z0_10, kappa10, bildkoordinaten10_dicht, koordinaten10_dicht);
bildkoordinaten11_dicht = [bildkoordinaten11(3,1),bildkoordinaten11(3,2);
                         bildkoordinaten11(14,1),bildkoordinaten11(14,2);
                         bildkoordinaten11(13,1),bildkoordinaten11(13,2);
                         bildkoordinaten11(1,1),bildkoordinaten11(1,2)];
koordinaten11_dicht = [koordinaten11(3,1),koordinaten11(3,2),koordinaten11(3,3);
                     koordinaten11(14,1),koordinaten11(14,2),koordinaten11(14,3);
                     koordinaten11(13,1),koordinaten11(13,2),koordinaten11(13,3);
                     koordinaten11(1,1),koordinaten11(1,2),koordinaten11(1,3)];
[X0_11_dicht, Y0_11_dicht, Z0_11_dicht, omega11_dicht, phi11_dicht, kappa11_dicht, Q11_dicht] = RRS(c, X0_11, Y0_11, Z0_11, kappa11, bildkoordinaten11_dicht, koordinaten11_dicht);
% Genauigkeit
genau_10 = zeros(6,1);
for i = 1 : 6
    genau_10(i) = sqrt(Q10(i,i));
end
genau_11 = zeros(6,1);
for i = 1 : 6
    genau_11(i) = sqrt(Q11(i,i));
end

genau_10_eck = zeros(6,1);
for i = 1 : 6
    genau_10_eck(i) = sqrt(Q10_eck(i,i));
end
genau_11_eck = zeros(6,1);
for i = 1 : 6
    genau_11_eck(i) = sqrt(Q11_eck(i,i));
end

genau_10_dicht = zeros(6,1);
for i = 1 : 6
    genau_10_dicht(i) = sqrt(Q10_dicht(i,i));
end
genau_11_dicht = zeros(6,1);
for i = 1 : 6
    genau_11_dicht(i) = sqrt(Q11_dicht(i,i));
end

% Differenz
diff_10_AB = [X0_10, Y0_10, Z0_10, omega10, phi10, kappa10] - [X0_10_eck, Y0_10_eck, Z0_10_eck, omega10_eck, phi10_eck, kappa10_eck];
diff_10_AC = [X0_10, Y0_10, Z0_10, omega10, phi10, kappa10] - [X0_10_dicht, Y0_10_dicht, Z0_10_dicht, omega10_dicht, phi10_dicht, kappa10_dicht];

diff_11_AB = [X0_11, Y0_11, Z0_11, omega11, phi11, kappa11] - [X0_11_eck, Y0_11_eck, Z0_11_eck, omega11_eck, phi11_eck, kappa11_eck];
diff_11_AC = [X0_11, Y0_11, Z0_11, omega11, phi11, kappa11] - [X0_11_dicht, Y0_11_dicht, Z0_11_dicht, omega11_dicht, phi11_dicht, kappa11_dicht];