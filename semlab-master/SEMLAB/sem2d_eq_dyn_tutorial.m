
% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid .
%
% You can modify simulation parameters in parts of the code 
% that have comments starting as: "**** Set here ..."
%
% Version 2: domain = rectangular
%            medium = general (heterogeneous)
%            boundaries = 1 fault + 3 paraxial
%	     time scheme = Newmark as in SPECFEM3D: alpha=1,beta=0,gamma=1/2
%
% Aug 7 2007: + effect of finite seismogenic depth LZ (2.5D, crude crustal-plane model)
%             + Kelvin-Voigt viscosity (ETA)
%             + moved friction update to begining of timestep
%             + immediate healing when slip rate = 0
%
% May 29 2018: + optional symmetry with respect to x=0
%            	 (then boundaries = 1 fault + 2 paraxial + 1 free)
%
% July 13 2018: vectorized version, runs faster
%               needs multiprod from dc
%               https://fr . mathworks . com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled
%
% This script is intended for tutorial purposes .
%
% Jean-Paul Ampuero	ampuero@gps . caltech . edu
%


disp('Initializing ...')


%------------------------------------------
% Parameters setup
%------------------------------------------

%**** Set here the parameters of the square box domain and mesh : ****
LX=100e3;                     % Horizontal length of the box domain
LY=100e3/3;                   % Vertical length of the box domain
LZ=inf;                       % Thickness of the seismogenic region (approximately accounted for) . 
                              % Turned off if LZ=inf .
%LZ=8e3;                              
% NELX = 150; NELY = 50; P = 8;% Number of elements; polynomial degree
NELX = 300; NELY = 100; P = 4;% Number of elements; polynomial degree
SYM_X = 1;                    % If SYM_X=1, enforce symmetry with respect to x=0
                              % (the left boundary becomes a free surface)
%********

%**** Set here the parameters of the material and simulated dumping factor : ****
RHO = 2670.;                  % The density (g/m^3)
VS  = 3464.;                  % The S wave speed (m/s)
ETA = 0.2;                    % Kelvin-Voigt viscosity term = ETA*dt*K*v
                              % Artificial viscosity to dump the high-frequency numerical noise
%********

%**** Set here the parameters of the time solver : ****
NT = 3200;                    % Number of timesteps
CFL   = 0.6; 			      % Stability number = CFL_ 1D / sqrt(2)
                              % Ratio of simulation time step to min (Delta_x)/VS .
%********

%**** Set here the parameters of the fault : ****
Normal_stress    = 120e6;     % Fault normal stress
Shear_stress     = 70e6;      % Shear stress
Static_friction  = 0.677;     % Static friction coefficient
Dynamic_friction = 0.525;     % Dynamic friction coefficient
dc               = 0.4;       % Slip-weakening distance
Fault_length     = 30e3;      % Simulated fault length
Nuc_length       = 1.5e3;     % Nucleation length
Nuc_str          = 81.6e6;    % Initial shear stress in nucleation zone
%********

%**** Set here receiver locations : ****
Sta_x = [-16e3:2400:32e3]';   % x coord of receivers
Sta_y = 7.5e3;	              % y coord of receivers
OUTdt = 100;                   % steps for output snapshots
%********

%**** Set here properties of heterogeneity: ****
Het_type = 'uniform'      ;   % Options are: uniform, trapezoid, random, box, EQ, and steps
                              % The default option is uniform .
Het_para = 'stress'      ;    % The options are: stress and dc . If Het_type='EQ', 
                              % this option shall be stress .
Het_loc  = 20e3        ;       % The closest location of heterogeneity (m)
Het_len  = 10e3        ;     % The length of heterogeneity (m)
Het_val  = 5e6         ;     % The value of heterogeneity . The unit is Pa for stress and m for dc .
							 % For random,    this value indicates pertubation range
% Het_val  = [70e6, 50e6];    % For trapezoid, the format is [val1,val2]
                              % For EQ, this value indicates the stress drop of the previous event .
                              
%Het_val  = [0.7,0.1,1.2,0.3];    % For steps, the format is [val1, val2, ...]
%********



%--------------------------------------------------------------------------
% Main code
%--------------------------------------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION
dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
if ~SYM_X, x = x-LX/2; end
nglob = length(x);

%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

isCrustalPlane = isfinite(LZ);

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
M     = zeros(nglob,1);		% global mass matrix, diagonal
if (isCrustalPlane), KZ = zeros(nglob,1); end % global crustal-plane matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

dt    = Inf;  			% timestep (set later)

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
for e=1:NEL,
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
  rho(:,:) = RHO;
  mu(:,:)  = RHO* VS^2;
%********

 % Diagonal mass matrix
  M(ig) = M(ig) + wgll2 .*rho *jac;

 % Local contributions to the stiffness matrix K
 %  WX(:,:,e) = wgll2 .* mu *jac/dx_dxi^2;
 %  WY(:,:,e) = wgll2 .* mu *jac/dy_deta^2;
  W(:,:,e) = wgll2 .* mu; 

 % The timestep dt is set by the stability condition
 %   dt = CFL*min(dx/vs)
  vs = sqrt(mu ./rho); 
  if dxe<dye
    vs = max( vs(1:NGLL-1,:), vs(2:NGLL,:) );
    dx = repmat( diff(xgll)*0.5*dxe ,1,NGLL); 
  else
    vs = max( vs(:,1:NGLL-1), vs(:,2:NGLL) );
    dx = repmat( diff (xgll)'*0.5*dye ,NGLL,1); 
  end
  dtloc = dx ./vs;
  dt = min( [dt dtloc(1:end)] );

 % diagonal matrix for crustal plane term, KZ*d = mu*(pi/2/LZ)^2 * displacement
 % pi/2 comes from quarter-of-wavelength proxy
  if (isCrustalPlane),
    KZ(ig) = KZ(ig) + wgll2 .*mu *jac *(pi/2/LZ)^2 ;
  end

end %... of element loop

dt = CFL*dt;
if ETA, dt=dt/sqrt(1+2*ETA); end
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;

% connectivity matrix for vectorized assembly (without loops)
Conn = sparse(iglob(:),[1:NGLL*NGLL*NEL],1);

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

a_elem = zeros(NGLL,NGLL,NEL);

time = (1:NT)'*dt;

%-- Absorbing boundaries (first order): 
% The mass matrix needs to be modified at the boundary
% for the IMPLICIT treatment of the term C*v .
% Fortunately C is diagonal .
impedance = RHO*VS;

if ~SYM_X,
  [BcLC,iBcL] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'L'); % Left
  BcLC = impedance*BcLC;
  M(iBcL) = M(iBcL) +half_dt*BcLC;
end

[BcRC,iBcR] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'R'); % Right
BcRC = impedance*BcRC;
M(iBcR) = M(iBcR) +half_dt*BcRC;

[BcTC,iBcT] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi ,'T'); % Top
BcTC = impedance*BcTC;
M(iBcT) = M(iBcT) +half_dt*BcTC;

%-- DYNAMIC FAULT at bottom boundary
[FltB,iFlt,jFlt] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi,'B');
FltN = length(iFlt);
FltZ = M (iFlt) ./FltB /half_dt;
FltX = x(iFlt);
FltV = zeros(FltN,NT+1);
FltD = zeros(FltN,NT+1);
% background
FltNormalStress = Normal_stress;
FltInitStress = repmat(Shear_stress,FltN,1);
FltState = zeros(FltN,1);
FltFriction . MUs = repmat(Static_friction,FltN,1);
FltFriction . MUd = repmat(Dynamic_friction,FltN,1);
FltFriction . Dc  = repmat(dc,FltN,1);
% barrier
L_BARRIER = Fault_length; % hypocentral distance to barrier
isel = find(abs(FltX)>L_BARRIER);
FltFriction . MUs(isel) = 1e4; % barrier
FltFriction . MUd(isel) = 1e4; % barrier
% nucleation
isel = find(abs(FltX)<=Nuc_length);
FltInitStress(isel) = Nuc_str;
FltFriction . W = (FltFriction . MUs-FltFriction . MUd) ./FltFriction . Dc;
FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                - FltInitStress; % strength excess

% Heterogeneity
switch Het_type
    case 'uniform'
        disp('Uniform case')
        
    case 'box'
        disp('Box heterogeneity')
        isel = find(abs(FltX)>Het_loc & abs(FltX)<=Het_loc+Het_len);
        if(strcmp(Het_para,'stress'))
            FltInitStress(isel) = Het_val;
        else
            FltFriction . Dc(isel) = Het_val;
        end
        
    case 'trapezoid'
        disp('Trapezoid heterogeneity')
        isel = find(abs(FltX)>Het_loc & abs(FltX)<=Het_loc+Het_len);
        if(strcmp(Het_para,'stress'))
            FltInitStress(isel) = Het_val(1)+(Het_val(2)-Het_val(1))*(FltX(isel)-Het_loc)/Het_len;
        else
            FltFriction . Dc(isel)= Het_val(1)+(Het_val(2)-Het_val(1))*(FltX(isel)-Het_loc)/Het_len;
        end
        
    case 'random'
        disp('Random heterogeneity')
        isel = find(abs(FltX)>Het_loc & abs(FltX)<=Het_loc+Het_len);
%        random_array = smoothdata(rand(length(isel),1)*2-1,'gaussian',30);
%        random_array = random_array/max(abs(random_array));
        rand_array = rand(length(isel),1)*2-1;
        gauss_windom = gausswin(30);
        random_array = filter (gauss_windom,1,rand_array)/max(abs(filter(gauss_windom,1,rand_array)));

        if(strcmp(Het_para,'stress'))
            FltInitStress(isel) = random_array*Het_val+FltInitStress(isel);
        else
            FltFriction . Dc(isel)= random_array*Het_val+FltFriction . Dc(isel);
        end
        
    case 'EQ'
        disp('Heterogeneity from previous EQ')
        if(strcmp(Het_para,'dc'))
            disp('Het_para shall be stress!!!')
        end
        peak_amp = 2;
        decay_dis = Het_len/2;
        isel = find(abs(FltX)>Het_loc & abs(FltX)<=Het_loc+Het_len);
        FltInitStress(isel) = FltInitStress(isel)-Het_val;
        isel = find(abs(FltX)<=Het_loc);
        FltInitStress(isel) = FltInitStress(isel) + exp(-(Het_loc-FltX(isel))/decay_dis)*Het_val*peak_amp;
        isel = find(abs(FltX)>Het_loc+Het_len);
        FltInitStress(isel) = FltInitStress(isel) + exp(-(FltX(isel)-Het_loc-Het_len)/decay_dis)*Het_val*peak_amp;
    case 'steps'
        disp('Steps heterogeneity')
        num_bins = length(Het_val);
        for i=1:num_bins
            isel = find(abs(FltX)>Het_loc+(i-1)*Het_len/num_bins & abs(FltX)<=Het_loc+i*Het_len/num_bins);
            if(strcmp(Het_para,'stress'))
                FltInitStress(isel) = Het_val(i);
            else
                FltFriction . Dc(isel)= Het_val(i);
            end
        end
        
otherwise
    disp('Wrong Het_type!!!')
    disp('Run the default uniform case!')
    
    end

if(LZ~=inf)
    model_name = ['2.5D-', num2str(LZ/1e3)];
else

    model_name = ['2D'];
end
if(strcmp(Het_type,'uniform'))
    model_name = [model_name, '-uniform'];
else
    if(strcmp(Het_para,'stress'))
        if(strcmp(Het_type,'trapezoid'))
            value_name = [num2str (Het_val (1)/1e6) '-' num2str(Het_val (2)/1e6)];
        elseif(strcmp(Het_type,'steps'))
            value_name = [];
            for i=1:length(Het_val)
                value_name = [value_name '-' num2str(Het_val (i)/1e6)];
            end
        else
            value_name = num2str(Het_val/1e6);
        end
    else
        if(strcmp(Het_type,'trapezoid'))
            value_name = [num2str (Het_val(1)) '-' num2str(Het_val(2))];
        elseif(strcmp(Het_type,'steps'))
            value_name = [];
            for i=1:length(Het_val)
                value_name = [value_name '-' num2str(Het_val(i))];
            end
        else
            value_name = num2str(Het_val);
        end
    end
    model_name = [model_name, '-', Het_type,'-',Het_para,'-',num2str(Het_loc/1e3),...
        '-',num2str(Het_len/1e3),'-',value_name];
end
disp(['Model name is:' model_name])

if ETA,  % Kelvin-Voigt viscosity
  isKelvinVoigt = zeros(NEL,1);
  if SYM_X,
    e = find( [dxe/2:dxe:(LX-dxe/2)]<L_BARRIER+dxe );
  else
    e = find( abs([-(LX/2+dxe/2):dxe:(LX/2-dxe/2)])<L_BARRIER+dxe );
  end
  isKelvinVoigt(e) = 1;
  x1 = 0.5*(1+xgll');
  eta_taper = exp(-pi*x1.^2); 
  eta_elem = ETA*dt *repmat(eta_taper, NGLL,1 );
  eta = zeros(NGLL,NGLL,NEL);
  for e =1:NEL,
    if isKelvinVoigt(e), eta(:,:,e) = eta_elem; end
  end
end

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = Sta_x;		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(Sta_y,OUTnseis,1);	% y coord of receivers
%********
% receivers are relocated to the nearest node
% OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark scheme with
% alpha=1, beta=0, gamma=1/2
%
cpu_t = cputime;
disp('Starting time loop ...')

for it=1:NT,

 % update
  d = d + dt*v + half_dt_sq*a; 
%  FltState = max(2*d(iFlt),FltState); % no healing
  FltState = ( FltState + max(2*d(iFlt)-FltD(:,it), 0) ) .* (v(iFlt)>0); % with healing

 % prediction 
  v = v + half_dt*a;

 % internal forces -K*d(t+1) 
 % stored in global array 'a'
 % switch to local representation (element-by-element) 
  a_elem = d(iglob) + eta .*v(iglob);

 % vectorized version (without loops) using multiprod from 
 % https://fr . mathworks . com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled
 % gradients wrt local variables (xi,eta)
  d_xi  = multiprod( Ht , a_elem );
  d_eta = multiprod( a_elem, H );
 % element contribution to internal forces
 % local = coefint1*H*( W (:,:,e) .*d_xi ) + coefint2*( W (:,:,e) .*d_eta )*Ht ;
  d_xi = W .* d_xi;
  d_xi = multiprod( H , d_xi );
  d_eta = W .* d_eta;
  d_eta = multiprod( d_eta , Ht);
  a_elem = coefint1 * d_xi  + coefint2 * d_eta ;

 % assemble into global vector
 % vectorized (without loop) using a sparse connectivity matrix:
  a = - Conn * a_elem(:);

 % absorbing boundaries:
  if ~SYM_X,  a(iBcL) = a(iBcL) - BcLC .* v(iBcL); end
  a(iBcR) = a(iBcR) - BcRC .* v(iBcR);
  a(iBcT) = a(iBcT) - BcTC .* v(iBcT) ;

 % crustal-plane term
  if (isCrustalPlane), a = a - KZ .*d; end

 % fault boundary condition: slip weakening
  FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                  -FltInitStress;
%% rupture nucleation through time weakening, like Andrews 76
%  VNUC = 0.4;
%%  if time(it)<10/VNUC
%    % FltStrength = min( FltStrength,...
%    ix = find(abs(FltX)<=10);
%    FltStrength(ix) = min(FltStrength(ix), ...
%              max(0.5,0.55+(abs(FltX(ix))-VNUC*time(it))*0.05/(1.0*dxe))...
%                     - FltInitStress(ix) ) ;
%%  end

  FltVFree = v(iFlt) + half_dt*a (iFlt) ./M(iFlt);
  TauStick = FltZ .*FltVFree;
 % TauStick = a (iFlt) ./FltB;
  Tau = min(TauStick,FltStrength); 
  a(iFlt) = a(iFlt) - FltB .*Tau;

 % solve for a_new:
  a = a ./M ;

 % correction
  v = v + half_dt*a;

  FltV(:,it+1) = 2*v(iFlt);
  FltD(:,it+1) = 2*d(iFlt);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;


    gcf = figure(1);
    Process_results(model_name,Het_para,FltX,FltInitStress,FltFriction . Dc,...
        Normal_stress,Static_friction,Dynamic_friction,x,y,v,OUTindx,[0 2],...
        Fault_length, OUTxseis,OUTyseis,time,OUTv,FltV,FltD,it);

    drawnow

  end

end % ... of time loop
  
cpu_e = cputime-cpu_t;
disp(['Loop running time is ', num2str(cpu_e/60/2),' min ']);

figure(2)
sdt=10;
sdx=4;
surf(time(1:sdt:end),FltX(1:sdx:end),FltV(1:sdx:end,2:sdt:end))
shading interp
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')
ylim([0,Fault_length])

% Calculate rupture time and speed
steph = 500.0;
Init_t0 = zeros(FltN,1);
grid_x   = linspace(0,LX,LX/steph);
grid_len = length(grid_x);
t0_grid  = zeros(grid_len,1);
Speed = zeros(grid_len,1);
for p=1:FltN
  if(Init_t0(p)>0.0), continue; end
  for ii=1:NT
    if(FltV(p,ii)<0.001 && FltV(p,ii+1)>=0.001)
      Init_t0(p) = ii*dt;
      break;
    end
  end
end
for i=1:length(grid_x)
    for j=1:length(FltX)-1
        if(grid_x(i)>=FltX(j) && grid_x(i)<FltX(j+1))
            t0_grid(i) = Init_t0(j) + (Init_t0(j+1)-Init_t0(j))*(grid_x(i)-FltX(j))/(FltX(j+1)-FltX(j));
            break;
        end
    end
end
for p=3:grid_len-2
    if(t0_grid(p)==0.0)
        continue;
    end
    delta_t0 = abs((t0_grid(p+2) + t0_grid(p+1) - t0_grid(p-1) - t0_grid(p-2))/6.0 / steph);
    if(delta_t0 ~= 0.0 && p<NT-1)
        Speed(p,1) = 1/delta_t0 /VS;
    end
end
cohesive_zone = zeros(NT,1);
for i=1:NT
     if(isempty(find(FltD(:,i)>0.001)) || isempty(find(FltD(:,i)>=dc)))
         cohesive_zone(i) = 0.0;
     else
        cohesive_zone(i) = FltX(max(find(FltD(:,i)>0.001)))-FltX(max(find(FltD(:,i)>=dc)));
     end
end


tip = zeros(NT,1);
tail  = zeros(NT,1);
for i=1:NT
         if(isempty(find(FltD(:,i)>0.001)) || isempty(find(FltD(:,i)>=dc)))
         tip(i) = 0;
         tail(i)= 0;
         else
     tip(i) = FltX(max(find(FltD(:,i)>0.001)));
     tail(i)= FltX(max(find(FltD(:,i)>=dc)));
         end
end


figure(3)
sdt=10;
sdx=4;
image(time(1:sdt:end),FltX(1:sdx:end),FltV(1:sdx:end,2:sdt:end)*10)
shading interp
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')
ylim([0,Fault_length])
hold on
plot(time,tip,'r')
hold on
plot(time,tail,'k')


figure(4)
subplot(211)
plot(grid_x/1e3,Speed)
xlabel('Fault X (km)')
ylabel('v / vs')
xlim([Nuc_length/1e3,Fault_length/1e3]);
ylim([0,1]);
title('Rupture speed');

subplot(212)
plot(time,cohesive_zone)
xlabel('Time (s)')
ylabel('Cohesive zone size (m)')
title('Cohesive zone size');



% Source properties
Sample_points = 200.0/dt;
STF = zeros(int64(Sample_points),1);
for i=1:size(FltV,2)
    for p=2:length(FltX)
        STF(i) = STF(i) + FltV(p,i)*(FltX(p)-FltX(p-1));
    end
end
for i=1:60
    STF(i) = STF(i) * sin((time (i)/time(60))*90.0/180.0*3.14);
end
        
Spec = abs(fft(STF));
Spec = Spec ./ Sample_points;
f = 1/200.0*(1:Sample_points);

disp('Output results to ./output ...')
save_dir = pwd;
% write the data
output_file = [save_dir,'/output/',model_name,'-spectrum . dat'];
if exist(output_file, 'file')==2
  delete(output_file);
end
fileID = fopen(output_file,'w');
for p=1:length (f)/2,
  fprintf(fileID,'%6 .5f %6 .5f \n',f(p),Spec(p));
end
fclose(fileID);

output_file = [save_dir,'/output/',model_name,'-STF . dat'];
if exist(output_file, 'file')==2
  delete(output_file);
end
fileID = fopen(output_file,'w');
for p=1:length(time),
  fprintf(fileID,'%6 .5f %6 .5f \n',time(p),STF(p));
end
fclose(fileID);
 

