% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                          TRABAJO FIN DE GRADO                           %  
%           -------------------------------------------------             %
% Simulador de imágenes cerebrales para el incremento del tamaño muestral %
% en estudios funcionales.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ULISES VIDAL SANZ

%% Añadimos libreria LIBLINEAR
addpath('/home/ulises/liblinear-210/matlab')

%% Carga la base de datos con las imï¿½genes y datos de todos los individuos.
load wstack_allVS1uint8+infp.mat
%Extraemos parámetros tras aplicar PCA a las imï¿½genes y KDE
% extractCoeff(stack_all,inf);
%Cargamos los coeficientes y scores PCA
load parameters

% Validación Cruzada

%%% EVALUACIÓN NORMAL %%%
% Realizamos X-Val sobre base de datos original.

% NORvsMCI
% Eliminamos clase NOR de la base de datos:
% labs0= find(labels==0); % Indices etiquetas clase NOR/'0'
labs1= find(labels==1); % Indices etiquetas clase MCI/'1'
% labs2= find(labels==2); % Indices etiquetas clase AD/'2'
labelsm=labels;
labelsm(labs1)=[];
labelsm = (labelsm)' > 0;
stack_xval=double(stack_masked);
stack_xval(labs1,:,:,:)=[];

%Normalization to the maximum (average of 3% higher intensity voxels)
Y=sort(stack_xval');
for i=1:size(stack_xval,1)
    maxValue = mean(Y(end:-1:end-round(size(Y,1)*0.03)));
    stack_xval(i,:)=stack_xval(i,:)/maxValue;
end

%Parámetros X-Val:
K = 10;
kernel = 0;
CVO=cvpartition(labelsm,'kFold',K);
perf=[];

for s=1:K
   
    trInd= CVO.training(s); % Training indices
    teInd=CVO.test(s);      % Test indices 
    traininglabels=labelsm(trInd);  % Training labels
    testlabels=labelsm(teInd);      % Test labes
    features_training=stack_xval(trInd,:); % Training features
    features_test=stack_xval(teInd,:);     % Test feaures
    
    c = train(double(traininglabels),sparse(features_training),'-C -s 2');
    
    cmd = ['-s 2 -c ', num2str(c(1)), ' -q'];
    model = train(double(traininglabels), sparse(features_training), cmd);
    ypredicted=predict(double(testlabels),sparse(features_test),model);
 
    cmat =  confusionmat(testlabels>0,ypredicted>0, 'order', sort(unique(traininglabels)));
    testval(s,:) = cmat(:);
    
end
cvdata = structConfMat(testval); 
% fsprint('Data & $%.3f \\pm %.3f $ & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$\\\\', cvdata.CorrectRate, cvdata.CRstd, cvdata.Sensitivity, cvdata.SensSTD, cvdata.Specificity, cvdata.SpecSTD)
% write2file('EvaluacionNormal.txt', fprintf('Data & $%.3f \\pm %.3f $ & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$\\\\', cvdata.CorrectRate, cvdata.CRstd, cvdata.Sensitivity, cvdata.SensSTD, cvdata.Specificity, cvdata.SpecSTD));
save('Results_dbOriginal_MCIvsAD_perf','testval','cvdata')

%%% EVALUACIï¿½N CON GENERACIï¿½N %%%

% Generaciï¿½n del stack para X-val -> 500 NOR y 500 MCI:
 [stackNOR,stackNOR_masked,labelsNOR]=generateBrains(500,cdf0,xmesh0,SCORE,COEFF,mn,dim,mask,0);
% [stackMCI,stackMCI_masked,labelsMCI]=generateBrains(100,cdf1,xmesh1,SCORE,COEFF,mn,dim,mask,1);
[stackAD,stackAD_masked,labelsAD]=generateBrains(500,cdf2,xmesh2,SCORE,COEFF,mn,dim,mask,2);
stack_generated = double([stackMCI_masked ; stackAD_masked]);
clear stackAD stackMCI stackAD_masked stackMCI_masked
stack_labels= [labelsMCI; labelsAD];
stack_labels = stack_labels > 0 ;

%Normalization to the maximum (average of 3% higher intensity voxels)
Y=sort(stack_generated');
for i=1:size(stack_generated,1)
    maxValue = mean(Y(end:-1:end-round(size(Y,1)*0.03)));
    stack_generated(i,:)=stack_generated(i,:)/maxValue;
end

%Parámetros X-Val:
K = 10;
kernel = 0;
CV1=cvpartition(stack_labels,'kFold',K);
perf=[];

for s = 1:K
    
    trInd= CV1.training(s); % Training indices
    teInd=CV1.test(s);      % Test indices 
    traininglabels=stack_labels(trInd);  % Training labels
    testlabels=stack_labels(teInd);      % Test labes
    features_training=stack_generated(trInd,:); % Training features
    features_test=stack_generated(teInd,:);     % Test feaures
 
    c = train(double(traininglabels),sparse(features_training),'-C -s 2');
    
    cmd = ['-s 2 -c ', num2str(c(1)), ' -q'];
    model = train(double(traininglabels), sparse(features_training), cmd);
    ypredicted=predict(double(testlabels),sparse(features_test),model);
 
    cmat =  confusionmat(testlabels>0,ypredicted>0, 'order', sort(unique(traininglabels)));
    testval(s,:) = cmat(:);
    
    
end

cvdata = structConfMat(testval); 
save('Results_dbGenerated_MCIvsAD_perf','testval','cvdata')
%  fsprintf('Data & $%.3f \\pm %.3f $ & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$\\\\', cvdata.CorrectRate, cvdata.CRstd, cvdata.Sensitivity, cvdata.SensSTD, cvdata.Specificity, cvdata.SpecSTD)
%  write2file('EvaluacionGener.txt', sfprintf('Data & $%.3f \\pm %.3f $ & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$\\\\', cvdata.CorrectRate, cvdata.CRstd, cvdata.Sensitivity, cvdata.SensSTD, cvdata.Specificity, cvdata.SpecSTD));

%% Visualizar cerebros
% figure(),montage(permute(stack_generated(90,:,:,:),[2 3 1 4]));title('CEREBRO AD'); colormap('jet'), colorbar;
% SliceBrowser(squeeze(stack_generated(:,:,:,:))) 