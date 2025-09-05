
%% Gaze360 - Azimuth angle and AOI angles

clear all
close all
clc

for who = 1:1
    
%% ------ Hyperparameters ------
geometrical = 'True';
who = 2;
interpolation = true;
gazeoffsets = false;
percentagebox = 0.5;
sigma = 2;
%% Check values for multiple widths 
% PRIMA TERAPISTA DA CLICCARE, POI BAMBINO 
% first for child, second for therapist
% per la persona, per cambiare l'AOI (area of interest) del terapista devo modificare la width
% del child

% DEFINIZIONE DI 3 CELL ARRAYS{} "width_x", "width_y" e "width_z" in cui
% vengono inseriti width_k, width_nao, ..., width_dietro rispettivamente per
% x, y e z. Suppongo servano a definire le dimensioni delle "zone" (come delle box) in cui trovo i vari oggetti 
% nella scena inquadrata dalla kinect (nella zona k), ovvero i 3 poster dx, sx e dietro, NAO e p,
% che potrebbe rappresentare il giocattolo (?)


width_k = [0.30, 0.20];
width_nao = [0.30, 0.20];
width_p = [0.15, 0.15];
width_dx = [0.15, 0.15];
width_sx = [0.15, 0.15];
width_dietro = [0.15, 0.15];

width_ky = [0.20, 0.10];
width_naoy = [0.05,0.05];
width_py = [0.10, 0.20];
width_dx_y = [0.10, 0.10];
width_sx_y = [0.10, 0.10];
width_dietro_y = [0.15, 0.15];

width_kz = [0.20, 0.10];
width_naoz = [0.05, 0.05];
width_pz = [0.10, 0.20];
width_dx_z = [0.05, 0.05];
width_sx_z = [0.05, 0.05];
width_dietro_z = [0.05, 0.05];

width_x = {}; %creazione di un cell array vuoto
width_x{1} = width_k;
width_x{2} = width_nao;
width_x{3} = width_dx;
width_x{4} = width_sx;
width_x{5} = width_dietro;
width_x{6} = width_p;

width_y = {};
width_y{1} = width_ky;
width_y{2} = width_naoy;
width_y{3} = width_dx_y;
width_y{4} = width_sx_y;
width_y{5} = width_dietro_y;
width_y{6} = width_py;

width_z = {};
width_z{1} = width_kz;
width_z{2} = width_naoz;
width_z{3} = width_dx_z;
width_z{4} = width_sx_z;
width_z{5} = width_dietro_z;
width_z{6} = width_pz;

%% Read Scene characteristics file

% visto dal soggetto: sx ha x>0 

% Definizione delle coordinate della kinect (origine sistema di
% riferimento), di NAO e dei 3 poster (dx, sx e dietro). Suppongo non venga
% fatto per il giocattolo in quanto è in movimento nella terapia. 
kinectcoord = [0.0,0.0,0.0];
naocoord = [0.45,-0.12,1.35]; %sul tavolo
% naocoord = [0.35,-0.85,0.99]; % per terra
dx_coord = [-1.38,0.60,1.42];
sx_coord = [1.50,0.40,1.45];
dietro_coord = [0.06,0.08, 3.38];

%target_label = {'kinect','robot','poster_dx','poster_sx','poster_dietro','altra_p','toy','therapist','Nowhere'};
% 'Nowhere','kinect','robot','poster_dx','poster_sx','poster_dietro','altra_p'
target_label ={'nowhere', 'robot','altra_p','poster_dx','poster_sx','poster_dietro', 'giocattolo'}; 

% creazione di un cell array in cui salvo le coordinate degli
% oggetti/poster definite sopra
coord = {}; 
coord{1} = kinectcoord;
coord{2} = naocoord;
coord{3} = dx_coord;
coord{4} = sx_coord;
coord{5} = dietro_coord;
%numbertarget = size(coord,2)+3; 
% numbertarget restituisce il numero di oggetti (target) che sono memorizzati nel cell array coord
% dunque 5+3=8. Suppongo siano per i target, nel nostro caso sono 6 quindi
% credo vada corretto in numbertarget = size(coord,2)+1;
numbertarget = size(coord,2)+1;
%direct = 'C:\Users\laura\Documents\Dados\Gabrielle_Laura'; %DA CAMBIARE CON PATH CARTELLA
direct = 'C:\Users\user\Desktop\PROVA';

folders={1,2,3,4,5,6,7,8};
subjects = {'Faures_Magnani_norobotOK';'Faures_Magnani_robotOK';'Macchini_Zeppa_norobotOK';'Macchini_Zeppa_robotOK';
    'Magnani_Faures_norobotOK';'Magnani_Faures_robotOK';'Martella_Zannino_norobotOK';'Martella_Zannino_robotOK'};
vector ={'1';'2';'3';'4';'5';'6';'7';'8'};

% Le variabili seguenti sono tutte cell array vuoti o strutture dati che verranno utilizzate 
% per raccogliere i dati durante l'esecuzione del codice. Variabili che raccoglieranno i dati 
% di angoli azimut, angoli di elevazione e tempi per ciascun paziente
extot = {}; %
elevationtot = {}; %
tmtot = {}; %

child_az = {}; %estrazione angolo azimuth bambino 
adult_az = {}; %e terapista 
timtot = {};

descripttot = {};
exvectot = {};

sessions_all = [1];
parents_all = {1};

% Gaze360, NAO_standard, Subject/Therapist_standard angles
for patient = 1:size(subjects) %iterazione sul numero di pazienti

    azimutessession = {};
    elevationsession = {};
    tempossession = {};

    descriptsession = {};
    keypsession = {};
    exvecsession = {};

    for session = 1:size(folders) %iterazione sul numero di sessioni (forse per robot/norobot?)
        subjects{patient}
        
        [keyp,tempos,azimutessession{end+1},elevationsession{end+1}, ...
            confidence_az_th,confidence_az_ch,confidence_el_th,confidence_el_ch] = ...
            Read_files_elevation_Odi_solo_porta_2pp_th(direct,subjects{patient}, interpolation, ...
            vector{patient}, percentagebox, coord);

        % Join session times and exercises
        tempossession{end+1} = tempos;
        keypsession{end+1} = keyp;
        
        

       
 
    end
     
    if isempty(azimutessession) || isempty(azimutessession{1})
        extot{patient} = [];
        elevationtot{patient} = [];
        tmtot{patient} = [];
        descripttot{patient} = [];
        exvectot{patient} = []; 

        standard_el_ch{patient} = [];%aggiorno ad ogni iterazione, come avere standard_el_ch{i} con i da 1 a num pazienti
        standard_el_th{patient} = [];%in questo caro però non sto aggiornando ma lasciando la cella vuota
        standard_az_ch{patient} = [];
        standard_az_th{patient} = [];
        continue
    end
    
   

   %% 
    extot{patient} = azimutessession;
    elevationtot{patient} = elevationsession;
    tmtot{patient} = tempossession;
    descripttot{patient} = descriptsession;
    exvectot{patient} = exvecsession; 

end

end