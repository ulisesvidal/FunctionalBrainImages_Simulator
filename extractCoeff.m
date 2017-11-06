function  extractCoeff(stack_all,inf)
% extractCoeff() aplica PCA a la base de datos original y
% estima la cdf de cada SCORE mediante KDE.
% INPUTS:
%  stack_all  - 4D struct containning all the functional 3D images 
%         inf - struc containing information about all patients.
%       
load wstack_allVS1uint8+infp.mat
%% PCA
% Creamos una mÁscara simple con la media de todas las imágenes: 
Imean = mean(stack_all,1); 
% La máscara serán todos los voxels mayores que un umbral de intensidad:
umbral = 1; %1
mask = Imean > umbral; 
stack_masked = stack_all(:,mask);

%Usado en la reconstrucción:
mn = mean(stack_masked,1); % media por columnas (usado en lareconstrucción)
dim= [size(stack_all,2) size(stack_all,3) size(stack_all,4)];

% Aplicamos PCA a nuestras imágenes: 
[COEFF, SCORE, LATENT, TSQUARED] = pca(single(stack_masked));

%% Creamos las etiquetas ('0','1','2') asociadas a cada imagen/clase.
[labels]= makeLabels(inf);

%% KDE 
% Estimamos la PDF de cada SCORE mediante el método no paramétrico KDE.
% Extraemos y agrupamos los scores por clases
% (Normal '0', MCI '1' y AD '2'):

labs0= find(labels==0); % Indices etiquetas '0'
labs1= find(labels==1); % Indices etiquetas '1'
labs2= find(labels==2); % Indices etiquetas '2'

sc0= SCORE(labs0,:); %Scores pacientes Normales
sc1= SCORE(labs1,:); %Scores pacientes MCI
sc2= SCORE(labs2,:); %Scores pacientes AD

%Guardamos (por clases) la estimación de la PDF/cdf de cada score:
for s=1:size(SCORE,2)
    
    [bandwidth_0,density0,xmesh_0,cdf_0]=kde(sc0(:,s),2^12);
    [bandwidth_1,density1,xmesh_1,cdf_1]=kde(sc1(:,s),2^12);
    [bandwidth_2,density2,xmesh_2,cdf_2]=kde(sc2(:,s),2^12);
    
    %CDF generada mediante KDE
    cdf0(:,s)=cdf_0;
    cdf1(:,s)=cdf_1;
    cdf2(:,s)=cdf_2;
    
    xmesh0(:,s)=xmesh_0';
    xmesh1(:,s)=xmesh_1';
    xmesh2(:,s)=xmesh_2'; 
end

% Guardamos los coefficientes y parámetros PCA así como la estimación de la
% cdf de cada score:
% save('parameters.mat', 'SCORE', 'COEFF','xmesh0', 'xmesh1', 'xmesh2', 'cdf0', 'cdf1', 'cdf2','mn','dim','mask','labels','stack_masked')
save('stack_masked.mat','stack_masked','labels')
end


%##########################################################################

function [labels] = makeLabels(inf)
%makelabels() genera un vector de variables con las etiquetas de cada
%imagen asociando la primera posición del vector "labels" con la primera
%imagen y así sucesivamente.

etiquetas = {}; 
N = size(inf,2); % El número de imágenes

% Creamos el cell array de etiquetas. 
for i = 1:N
    etiquetas{i} = inf(i).project.subject.subjectInfo.CONTENT;
end

% Lo convertimos a un formato numérico. Asignaremos 0 a Normal, 1 a
% MCI y 2 a AD(Alzheimer's Disease).Creamos el cell array de equivalencias:
equivalencias = {'Normal', 'MCI', 'AD'};

% Recorremos el vector de etiquetas para generar las etiqeutas numéricas
% (que vamos a llamar labels): 
labels = [];
for i=1:N
    % Con esto vemos a qué posición corresponde la etiqueta
    coincidencia = strcmp(etiquetas{i}, equivalencias);
    % Y establecemos como etiqueta la posición menos 1 (para que sea 0-2)
    labels(i)= find(coincidencia)-1; 
end
end