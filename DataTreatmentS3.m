clc
clear variables
close all

% Definir la carpeta donde están los archivos
folderPath = 'D:\Doctorado\MaquinaMBM\Datos MBM'; % Cambia esto a tu ruta real

% Número total de archivos
numFiles = 160;

% % Parámetros de DBSCAN
% epsilon = 0.5;  % Radio de vecindad
% minPts = 5;     % Mínima cantidad de puntos para formar un cluster

% Inicializar una celda para almacenar los datos
DataXarm = cell(numFiles, 1);
DataSTM = cell(numFiles, 1);

% Leer cada archivo y almacenarlo en la celda
for i = 1:numFiles
    varName = sprintf('daniel_xarm%d', i);
    filename = fullfile(folderPath, sprintf('daniel_xarm%d.csv', i));
    if isfile(filename)
        dataStruct.(varName) = readmatrix(filename); % Almacenar en estructura
        DataXarm{i} = dataStruct.(varName);
    else
        warning('Archivo no encontrado: %s', filename);
    end
end

% Elimina celdas con datos en cero
DataXarm(any(cellfun(@(x) isequal(x, 0), DataXarm), 2), :) = [];

% Elimina celdas vacias o en cero
DataXarm(cellfun(@isempty, DataXarm) | cellfun(@(x) isequal(x, 0), DataXarm)) = [];

% % Inicializar cell array para etiquetas de clusters
% clusters = cell(size(DataXarm));
% 
% % Aplicación de DBScan
% for i = 1:numel(DataXarm)
%     if ~isempty(DataXarm{i})  % Verificar que la celda no esté vacía
%         data = DataXarm{i};  % Extraer matriz de la celda
% 
%         [rows, cols] = size(data);
%         cluster_labels = zeros(rows, cols);  % Matriz para almacenar etiquetas
% 
%         for j = 1:cols
%             column_data = data(:, j); % Extraer columna j-ésima
%             column_data = column_data(:); % Asegurar que es un vector
% 
%             % Aplicar DBSCAN (solo si hay más de 1 punto)
%             if length(column_data) > 1
%                 cluster_labels(:, j) = dbscan(column_data, epsilon, minPts);
%             else
%                 cluster_labels(:, j) = -1; % Marcar como ruido si hay pocos puntos
%             end
%         end
% 
%         clusters{i} = cluster_labels; % Guardar etiquetas en la celda
%     end
% end

% Leer cada archivo y almacenarlo en la celda
for i = 1:numFiles
    varName = sprintf('Run%d', i);
    if i == 7 || i == 109
        filename = fullfile(folderPath, sprintf('Run%d.mat', i));
        if isfile(filename)
            data = load(filename);
            fields = fieldnames(data);
            dataStruct.(varName) = data.(fields{1}); % Almacenar en estructura
            DataSTM{i} = data.(fields{1});
        else
            warning('Archivo no encontrado: %s', filename);
        end
    else
        filename = fullfile(folderPath, sprintf('Run%d.xlsx', i));
        if isfile(filename)
            dataStruct.(varName) = readmatrix(filename);
            DataSTM{i} = readmatrix(filename);
        else
            warning('Archivo no encontrado: %s', filename);
        end
    end
end

% Elimina celdas con datos en cero
DataSTM(any(cellfun(@(x) isequal(x, 0), DataSTM), 2), :) = [];

% Elimina celdas vacias o en cero
DataSTM(cellfun(@isempty, DataSTM) | cellfun(@(x) isequal(x, 0), DataSTM)) = [];

%% 
% Interpolar la segunda columna de DataXarm{i} en base a la primera columna de DataSTM{i}
ValidData = 44;
InterpolatedXarm = cell(ValidData,1); % Celda para almacenar los datos interpolados

for i = 1:ValidData
    for j = 2:12
        x1 = DataXarm{i}(:,1);  % Valores en los que está basada la señal (tiempo)
        y1 = DataXarm{i}(:,j);  % Valores a interpolar (señal)  
        % Obtener el número de columnas
        numColumnas = size(DataSTM{i}(:,1));
        % Crear el vector columna
        xq = (1:numColumnas(1))'; % Nuevos puntos de tiempo a los que queremos ajustar

        % Realizar interpolación
        InterpTemp = [xq, interp1(x1, y1, xq, 'spline', 'extrap')]; 
        InterpolatedXarm{i}(:,j) = flip(InterpTemp(:,2));
    end
end

for i = 1:ValidData
    if ~isempty(InterpolatedXarm{i}) && size(InterpolatedXarm{i},2) > 1
        InterpolatedXarm{i}(:,1) = []; % Elimina la primera columna de cada celda
    end
end

% Suponiendo que DataSTM1 y DataSTM2 son los DataStructures de 44x1
DataTotal = cell(ValidData,1); % Inicializar la nueva estructura de celdas

for i = 1:ValidData
    if ~isempty(DataSTM{i}) && ~isempty(InterpolatedXarm{i}) % Si ambas celdas tienen datos
        DataTotal{i} = [DataSTM{i}, InterpolatedXarm{i}]; % Concatenar columnas
    elseif ~isempty(DataSTM{i}) % Si solo la primera tiene datos
        DataTotal{i} = InterpolatedXarm{i};
    elseif ~isempty(DataSTM2{i}) % Si solo la segunda tiene datos
        DataTotal{i} = InterpolatedXarm{i};
    else
        DataTotal{i} = []; % Si ambas están vacías, dejar la celda vacía
    end
end

% Obtencion de parámetros por set de datos
Tiempo = cell(size(DataTotal)); % Preasignar celda de salida
for i = 1:numel(DataTotal)
    Tiempo{i} = DataTotal{i}(:,1);
end

VectorX = cell(size(DataTotal)); % Preasignar celda de salida
for i = 1:numel(DataTotal)
    VectorX{i} = DataTotal{i}(:,3);
end

VectorX{i} = VectorX{i}.*1000;

VectorF = cell(size(DataTotal)); % Preasignar celda de salida
for i = 1:numel(DataTotal)
    VectorF{i} = DataTotal{i}(:,11);
end

DVectorX = cell(1, ValidData);
DDVectorX = cell(1, ValidData);

for i = 1:ValidData
    DVectorX{i} = diff(VectorX{i}(:,1))/mean(diff(Tiempo{i}(:,1)));
    DDVectorX{i} = diff(DVectorX{i}(:,1))/mean(diff(Tiempo{i}(:,1)));
    % DVectorX{i} = gradient(VectorX{i}(:,1), Tiempo{i}(:,1));
    % DDVectorX{i} = gradient(DVectorX{i}(:,1), Tiempo{i}(:,1));
end

A = cell(1, ValidData);

for i = 1:ValidData
    A{i} = horzcat(DDVectorX{i}, DVectorX{i}((1:end-1),1), DVectorX{i}((1:end-1),1).^2, DVectorX{i}((1:end-1),1).^3, VectorX{i}((1:end-2),1), VectorX{i}((1:end-2),1).^2, VectorX{i}((1:end-2),1).^3);
end

Theta = cell(1, ValidData);

for i = 1:ValidData
    Theta{i} = abs(A{i}) \ abs(VectorF{i}((1:end-2),1));
end

Theta = Theta';

for i = 1:length(Theta)
    Theta{i} = Theta{i}'; % Transponer cada vector dentro de la celda
end

Theta = vertcat(Theta{:});
STDdev = std(Theta);
Theta = Theta*1000;

numCeldas = length(DataXarm); % Número de celdas
mediaColumnas = cell(numCeldas, 1); % Inicializar celda para guardar resultados

for i = 1:numCeldas
    if ~isempty(DataXarm{i}) % Verifica que la celda no esté vacía
        mediaColumnas{i} = mean(abs(DataXarm{i}), 1); % Calcula la media de cada columna
    end
end

% Convertir las celdas en un solo array (si los datos tienen el mismo
% tamaño) concatenando al final de cada array 
mediaColumnas = vertcat(mediaColumnas{:});

mediaFX = mediaColumnas(:,7);
mediaFY = mediaColumnas(:,8);
mediaFZ = mediaColumnas(:,9);

% Convertir las celdas en un solo array (si los datos tienen el mismo
% tamaño) concatenando al final de cada array 
DataTotal = vertcat(DataTotal{:});

% Obtención de parámetros Totales
XTotal = DataTotal(:,3);
TTotal = DataTotal(:,1);
FXTotal = DataTotal(:,11);
FYTotal = DataTotal(:,12);
FZTotal = DataTotal(:,13);
MXTotal = DataTotal(:,14);
MYTotal = DataTotal(:,15);
MZTotal = DataTotal(:,16);

% Derivada de el vector de posición
% Gradiente
% DX = gradient(XTotal, TTotal);
% DDX = gradient(DX, TTotal);
% Diff
% DX = diff(XTotal)/mean(diff(TTotal));
% DDX = diff(DX)/mean(diff(TTotal));
% Diferencias finitas
DX = (XTotal(2:end) - XTotal(1:end-1))/mean(diff(TTotal));
DDX = (XTotal(3:end) - 2*XTotal(2:end-1) + XTotal(1:end-2))/mean(diff(TTotal))^2;

ATotal = horzcat(DDX, DX(1:end-1), DX(1:end-1).^2, DX(1:end-1).^3, XTotal(1:end-2), XTotal(1:end-2).^2, XTotal(1:end-2).^3);

ThetaTPlus = abs(ATotal) \ abs(FXTotal(1:end-2));

% Minimos cuadrados en linea
Y = FXTotal(1:end-2);

phi = [DDX,DX(1:end-1),DX(1:end-1).^2,DX(1:end-1).^3,XTotal(1:end-2),XTotal(1:end-2).^2,XTotal(1:end-2).^3];
thetaG = (phi'*phi) \ phi'*Y;

% % % % % % % % % % % 

% ThetaC = cell(1, ValidData);
% Inicializar ThetaC con ceros del tamaño correcto
ThetaC = zeros(7,1); % Suponiendo que A{i} tiene dimensiones (muestras x parámetros)
ThetaE = cell(size(A));

YC = VectorF;
% A = [DDVectorX, DVectorX((1:end-1),1),XTotal(1:end-2)];

for i = 1:ValidData
    YC_i = YC{i};
    A_i = A{i};

    % Para cada muestra en el tiempo (cada fila de A_i)
    for j = 1:size(A_i, 1)
        % Estimación de parámetros usando mínimos cuadrados en línea
        A_ij = A_i(j, :)';  % Vector columna de la fila j de A_i
        YC_ij = YC_i(j);      % Resultado correspondiente a la muestra j

        % Evitar división por cero
        denominator = A_ij' * A_ij; % Esto es un escalar
        if denominator ~= 0
            % Ajuste en línea: Actualiza theta con el nuevo dato (muestreo)
            ThetaC = ThetaC + (A_ij * (YC_ij - A_ij' * ThetaC)) / denominator;
        end
        % Guardar la evolución de theta después de procesar la celda i
        ThetaE{i}(j, :) = ThetaC';
    end
end

for i = 1:numel(ThetaE)
    ThetaE{i} = ThetaE{i} * 1000; % Multiplica cada matriz por 1000
end

% Extraemos las tres componentes de theta (masa, amortiguamiento, resorte)
ME = cellfun(@(x) x(:,1), ThetaE, 'UniformOutput', false);    % Evolución de la masa
C1E = cellfun(@(x) x(:,2), ThetaE, 'UniformOutput', false);    % Evolución del amortiguamiento
C2E = cellfun(@(x) x(:,3), ThetaE, 'UniformOutput', false);    % Evolución del amortiguaimento 2
C3E = cellfun(@(x) x(:,4), ThetaE, 'UniformOutput', false);    % Evolución del amortiguaimento 3
K1E = cellfun(@(x) x(:,5), ThetaE, 'UniformOutput', false);    % Evolución del resorte
K2E = cellfun(@(x) x(:,6), ThetaE, 'UniformOutput', false);    % Evolución del resorte 2
K3E = cellfun(@(x) x(:,7), ThetaE, 'UniformOutput', false);    % Evolución del resorte 3

ME = ME';
C1E = C1E';
C2E = C2E';
C3E = C3E';
K1E = K1E';
K2E = K2E';
K3E = K3E';

% Plotear la evolucion de cada parámetro
% Evolucion de la masa
figure; hold on; grid on;
for i = 1:length(ME)
    N = length(ME{i});
    t = linspace(0, 2200, N);  % Tiempo en segundos
    plot(t, ME{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
end
xlabel('time (s)');
ylabel('Mass (g)');
legend show;
% Evolucion del amortiguamiento
figure; hold on; grid on;
for i = 1:length(C1E)
    N = length(C1E{i});
    t = linspace(0, 2200, N);  % Tiempo en segundos
    plot(t, C1E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
    plot(t, C2E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
    plot(t, C3E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
end
xlabel('time (s)');
ylabel('Damping coefficient (Ns/mm)');
legend show;
% Evolucion del resorte
figure; hold on; grid on;
for i = 1:length(K1E)
    N = length(K1E{i});
    t = linspace(0, 2200, N);  % Tiempo en segundos
    plot(t, K1E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
    plot(t, K2E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
    plot(t, K3E{i}, 'LineWidth', 2, 'DisplayName', ['Celda ' num2str(i)]);
end
xlabel('time (s)');
ylabel('Spring coefficient (N/mm)');
legend show;

for i = 1:length(ThetaE)
    ThetaE{i} = ThetaE{i}'; % Transponer cada vector dentro de la celda
end

ThetaE = vertcat(ThetaE{:});
STDdevE = std(ThetaE);

% % % % % % % % % 
% 
% save('DataSTM.mat', 'DataSTM');
% save('DataXarm.mat', 'DataXarm', '-v7.3');
% 
% % Concatenar Todos los datos
% % DataTotal = [DataXarm, DataSTM];
% 
% Gráfica de las derivadas
% % t = 1:length(TTotal);
% % figure;
% % plot(t, XTotal, 'b', 'LineWidth', 2); hold on;
% % plot(t(1:end-1), DX, 'r--', 'LineWidth', 2);
% % plot(t(1:end-2), DDX, 'g-.', 'LineWidth', 2);
% % legend('Función Original', 'Derivada', 'Segunda Derivada');
% % xlabel('Tiempo');
% % ylabel('Amplitud');
% % title('Derivadas Numéricas en MATLAB');
% % grid on;

% Análisis MANOVA de las variables independientes y dependientes
 % Variables inedpendintes
RasterAngle = [45; 45; 45; 45; 45; 45; 45; 45; 
                45; 45; 45; 45; 45; 45; 45; 45; 45; 45;
                0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
              0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
Layer = [0.12; 0.12; 0.12; 0.12; 0.12; 0.12; 
        0.16; 0.16; 0.16; 0.16; 0.16; 0.16; 
        0.16; 0.16; 0.2; 0.2; 0.2; 0.2; 0.12;
        0.12; 0.12; 0.12; 0.12; 0.12; 0.12; 
        0.12; 0.12; 0.16; 0.16; 0.16; 
      0.16; 0.16; 0.16; 0.16; 0.16; 0.16; 0.2; 0.2; 
      0.2; 0.2; 0.2; 0.2; 0.2; 0.2];
Width = [4.5; 2; 7; 4.5; 4.5; 2; 2; 7; 4.5; 4.5; 
         2; 2; 7; 7; 4.5; 7; 4.5; 7; 4.5; 4.5; 
         2; 2; 2; 7; 4.5; 4.5; 4.5; 4.5; 2; 2; 
         4.5; 4.5; 4.5; 4.5; 4.5; 4.5; 7; 7; 4.5; 
         4.5; 4.5; 4.5; 7; 4.5];
Length = [18; 24; 24; 24; 24; 18; 24; 24; 18; 18; 
         18; 30; 18; 30; 30; 24; 18; 18; 18; 30; 
         24; 24; 24; 24; 24; 24; 24; 30; 24; 24; 
         30; 18; 30; 24; 30; 24; 24; 24; 18; 30; 18; 
         24; 30; 24];
Height = [9; 9; 9; 4; 14; 9; 4; 4; 4; 14;
         9; 9; 9; 9; 9; 14; 4; 9; 9; 9;
         4; 9; 9; 9; 4; 4; 14; 9; 14; 9;
         7; 14; 14; 4; 9; 9; 4; 9; 4; 4; 14;
         4; 9; 9];
Density = [20; 20; 20; 20; 20; 60; 60; 60; 60; 60;
         60; 60; 60; 60; 100; 60; 60; 60; 20; 20;
         60; 20; 100; 20; 20; 100; 20; 20; 60; 20;
         60; 60; 60; 20; 20; 60; 60; 20; 60; 60; 60;
         100; 60; 60];

% Variables dependientes
Masa = Theta(:,1);
Damping = Theta(:,2);
Damping2 = Theta(:,3);
Damping3 = Theta(:,4);
Spring = Theta(:,5);
Spring2 = Theta(:,6);
Spring3 = Theta(:,7);

T = table(Damping, Damping2, Damping3, Spring, Spring2, Spring3, RasterAngle, Layer, Width, Height, Length, Density, 'VariableNames',{'Damping', 'Damping2', 'Damping3', 'Spring', 'Spring2', 'Spring3', 'RasterAngle', 'Layer', 'Width', 'Height', 'Length', 'Density'});
Model = fitrm(T, 'Damping-Damping3,Spring-Spring3 ~ RasterAngle + Layer + Width + Height + Length + Density');

ANOVA = manova(Model);
disp (ANOVA)

SS_values = ANOVA.RSquare;  % Obtener valores de SS

SS_Total = sum(SS_values); % Suma total de cuadrados
Contribution = (SS_values ./ SS_Total) * 100; % Porcentaje de contribución

disp(array2table([SS_values, Contribution], ...
    'VariableNames', {'SS', 'Contribution (%)'}))

% Analisis ANOVAN para Damping y Spring
% Variables predictoras
X = [RasterAngle, Layer, Width, Height, Length, Density];
% Variable respuesta
Y = Damping;
% Nombres de los factores
factorNames = {'RasterAngle', 'Layer', 'Width', 'Height', 'Length', 'Density'};
% ANOVA con interacciones de 2 en 2
[p1,ANOVAN1,stats1] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
% Mostrar tabla de resultados
disp(ANOVAN1)
Y = Damping2;  % Cambiamos la variable dependiente a "Spring"
[p2,ANOVAN2,stats2] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
disp(ANOVAN2)
Y = Damping3;  % Cambiamos la variable dependiente a "Spring"
[p3,ANOVAN3,stats3] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
disp(ANOVAN3)
Y = Spring;  % Cambiamos la variable dependiente a "Spring"
[p4,ANOVAN4,stats4] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
disp(ANOVAN4)
Y = Spring2;  % Cambiamos la variable dependiente a "Spring"
[p5,ANOVAN5,stats5] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
disp(ANOVAN5)
Y = Spring3;  % Cambiamos la variable dependiente a "Spring"
[p6,ANOVAN6,stats6] = anovan(Y, X, 'model', 'interaction', 'varnames', factorNames);
disp(ANOVAN6)

% % Extraer coordenadas canónicas
% gscatter(ANOVA.Canonical(:,1), ANOVA.Canonical(:,2), Density);
% xlabel('Coordenada Canónica 1');
% ylabel('Coordenada Canónica 2');
% title('Grupos proyectados en espacio de coordenadas canónicas');
% grid on;

% % Comparar Damping entre Variables Independientes
% figure;
% boxplot(Damping, RasterAngle);
% xlabel('Raster Aangle (^\circ)');
% ylabel('Damping (Ns/mm)');
% 
% figure;
% boxplot(Damping, Layer);
% xlabel('Layer Height (mm)');
% ylabel('Damping (Ns/mm)');
% 
% figure;
% boxplot(Damping, Width);
% xlabel('Width (mm)');
% ylabel('Damping (Ns/mm)');
% 
% figure;
% boxplot(Damping, Height);
% xlabel('Height (mm)');
% ylabel('Damping (Ns/mm)');
% 
% figure;
% boxplot(Damping, Length);
% xlabel('Length (mm)');
% ylabel('Damping (Ns/mm)');
% 
% figure;
% boxplot(Damping, Density);
% xlabel('Printing Density (%)');
% ylabel('Damping (Ns/mm)');
% 
% % Comparar Spring entre variables inependientes
% figure;
% boxplot(Spring, RasterAngle);
% xlabel('Raster Aangle (^\circ)');
% ylabel('Spring (N/mm)');
% 
% figure;
% boxplot(Spring, Layer);
% xlabel('Layer Height (mm)');
% ylabel('Spring (N/mm)');
% 
% figure;
% boxplot(Spring, Width);
% xlabel('Width (mm)');
% ylabel('Spring (N/mm)');
% 
% figure;
% boxplot(Spring, Height);
% xlabel('Height (mm)');
% ylabel('Spring (N/mm)');
% 
% figure;
% boxplot(Spring, Length);
% xlabel('Length (mm)');
% ylabel('Spring (N/mm)');
% 
% figure;
% boxplot(Spring, Density);
% xlabel('Printing Density (%)');
% ylabel('Spring (N/mm)');
