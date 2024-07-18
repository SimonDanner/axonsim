function [ t,Y,N_nodes,b_thr] = crrss(dur,data,stim_dur,fun_type,custom_fun,fiberD,frq,end_on_ap,I_intra,N_intra)
% V_stim: stimulation voltage
% T: time [ms] 
isi = 1/frq*1000;
% parameters
% ----------------------------------------------------------------
l_n = 1.5e-4;        % node length [cm]
D = fiberD/10000;           % fiber diameter with myelin sheath [cm]
N_layers = 20e4*D;   % represents internodal myelin sheath 
                     % (derived according to velocity check)
g_Na = 1445;         % sodium channel conductivity [mS/cm^2]
g_l = 128;           % leak channel conductivity [mS/cm^2] 
V_Na = 115;          % Nernst potential for sodium channels [mV]
V_l = -0.01;         % Nernst potential for leakage channels [mV]
r = 0.055;           % specific resistivity [kOhm*cm]
c = 0.6;             % specific capacity [muF/cm^2]
g = 1;               % internodal membrane conductivity [mS/cm^2]
k = 1;               % 37?C  k=3^(0.1*T-3.7)
%V_fem = 50;          % Voltage applied in FEM simulation [V]
%polarity = 1;       % -1/+1 ... neg./pos. active electrode
% ----------------------------------------------------------------

l_in = 100*D;
C_in = l_in * (0.64*D)*pi * c / N_layers;
C_n = l_n * (0.64*D)*pi * c;
R_in = 4*r * l_in /((0.64*D)^2*pi); 
R_n = 4*r * l_n /((0.64*D)^2*pi);  
g_m = l_in * (0.64*D)*pi * g / N_layers;
gNa_n = l_n * (0.64*D)*pi * g_Na;
gl_n = l_n * (0.64*D)*pi * g_l;

% read axon data file
% ----------------------------------------------------------------
%reply = input('data file: ','s');
DATA = data;
% for MATLAB_7 R14
%DATA = DATA.data;
Ve_pulse = 1000*DATA(:,4);
x = 100*DATA(:,1);
y = 100*DATA(:,2);
z = 100*DATA(:,3);
% N_comp odd, nodes at both ends, calculated in fiber_compcoords_Ve
%N_comp = length(x)-2;

% fiber length and interpolation to obtain Ve_pulse at the exact location
% along the fiber trajectory (longitudinal center of each compartment)
% (s and interp1 optional, but more accurate)
% ----------------------------------------------------------------
s(1) = 0;
for i=1:length(x)-1
    s(i+1) = s(i) + sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 + (z(i+1)-z(i))^2);
end

N_comp = 1+2*floor(s(length(x)-1)/(l_n+l_in))+32;
N_comp = floor(N_comp/2)*2;

xi = zeros(N_comp+2,1);
xi(1) = -(l_in+l_n)*7;    %j=1:10  to move nodes (j=6 ->lowest Vth at DREZ) 444444444
for i=2:N_comp+2
    xi(i) = xi(i-1) + (l_n + l_in)/2;
end

Ve_pulse = interp1(s,Ve_pulse,xi,'linear','extrap');


% ----------------------------------------------------------------

% initial condition vector
% ----------------------------------------------------------------
IC = zeros (2*N_comp+3,1);
N_nodes = N_comp/2;
for i=N_comp+3:2:2*N_comp+2
    IC(i) = 0.003;
    IC(i+1) = 0.75;
end
intra = zeros(2*N_comp+3,1);

% numerical solving
% ----------------------------------------------------------------
%tspan = 0:0.01:T;
%[t,Y] = ode15s (@odesys, [0,T],IC);
 options = CVodeSetOptions('RelTol',1.e-4,...
                          'AbsTol',1.e-5);
    
        b_thr=0;
        CVodeInit(@odeCVode,'BDF','Newton',0,IC,options);

        dtout = 0.01;
        tout = dtout;
        t=[];
        Y=[];
        for i = 1:dur/dtout

            [status,t1,Y1] = CVode(tout,'Normal');
            t=[t,t1];
            Y=[Y,Y1];
            if max(Y1) > 60
                b_thr = 1;
                if end_on_ap == 1
                    break;
                end
            end
            tout=tout+dtout;
        end
        Y=Y';
        Y=Y(:,1:2:N_comp);
        CVodeFree;


% plotting
% ----------------------------------------------------------------
% plotting gating variables m,h (optional)
%CH = zeros(length(t),N_comp+1);
%for i=1:(N_comp+1)/2
%    CH(:,i) = Y(:,N_comp+1+2*i) -3*(i-1)/10;
%    CH(:,(N_comp+1)/2+i) = Y(:,N_comp+2+2*i) -3*((N_comp+5)/2+i)/10;
%end

%figure('Position',[25+scrsz(3)/3 scrsz(4)/4 scrsz(3)/2 3*scrsz(4)/4])
%plot (t,CH,'DisplayName','CH');
%plotbrowser('on');

% %V = zeros(length(t),(N_comp+1)/2);
% for i=1:2:(N_comp+1)/2
%     V(:,(i+1)/2) = Y(:,2*i) -6*(i-1);     % displacement
% end
% 
% %scrsz = get(0,'ScreenSize')
% scrsz = [0,0,1000,1000];
% figure('Position',[5 40 scrsz(3)-10 5*scrsz(4)/6],'Name',...
%     ['Stimulus: ' int2str(V_stim) 'V / ' int2str(stim_off-stim_on) ...
%     'ms -- Response of  ' reply '  (every 2nd node shown)'])
% subplot(1,4,3:4);
% plot(t,V,'k','DisplayName','V');
% %xlabel('[ms]','FontSize',14); %ylabel('V [mV]');
% title('Membrane voltage','FontSize',14);
% set(gca,'YAxisLocation','right','YTick',[],'YLim',[-3*(191-1)-50 100],...
%     'XTick',[0 0.5 1 1.5 2],'FontSize',14);
%plotbrowser('on');

    function [dY,flag,new_data] = odeCVode(t,Y)
        dY = odesys(t,Y);
        flag = [0];
        new_data = [];
    end

% ----------------------------------------------------------------

% subfunction: ode system
% ----------------------------------------------------------------
    function dY = odesys (t,Y)
   
    dY = zeros (2*N_comp+3,1);      % see below
    V_e = zeros (N_comp+2,1);
    
    alpha_m = zeros (N_comp+1,1);
    beta_m = zeros (N_comp+1,1);
    alpha_h = zeros (N_comp+1,1);
    beta_h = zeros (N_comp+1,1);

    % stimulus between stim_on and stim_off ms
    k = stimulation_fun(mod(t,isi),stim_dur,fun_type,custom_fun);
    if exist('I_intra','var') ==0
    if t>=0 && t<=stim_dur 
        V_e = k*Ve_pulse;
    end
    else
        V_e=0*Ve_pulse;
        intra(N_intra*2)=k*I_intra*1e-6;
    end
    
    % nodes
    for i=2:2:N_comp+1
        alpha_m(i) = (97 + 0.363*Y(i)) / (1 + exp((31 - Y(i))/5.3));
        beta_m(i) = alpha_m(i) / exp((Y(i) - 23.8)/4.17);
        beta_h(i) = 15.6 / (1 + exp((24 - Y(i))/10));
        alpha_h(i) = beta_h(i) / exp((Y(i) - 5.5)/5);

        dY(i) = (-(-intra(i)+gNa_n * Y(N_comp+i+1)^2 * Y(N_comp+i+2) * (Y(i)-V_Na) + gl_n * (Y(i)-V_l)) + ...
                2*(Y(i-1)-Y(i))/(R_in + R_n) + 2*(Y(i+1)-Y(i))/(R_in + R_n) + ...
                2*(V_e(i-1)-V_e(i))/(R_in + R_n) + 2*(V_e(i+1)-V_e(i))/(R_in + R_n) )/C_n;
        dY(N_comp+i+1) = (-(alpha_m(i) + beta_m(i)) * Y(N_comp+i+1) + alpha_m(i)) *k;
        dY(N_comp+i+2) = (-(alpha_h(i) + beta_h(i)) * Y(N_comp+i+2) + alpha_h(i)) *k;
        
        act_fn(i-1) = (2*(Ve_pulse(i-1)-Ve_pulse(i))/(R_in + R_n) + ...
                    2*(Ve_pulse(i+1)-Ve_pulse(i))/(R_in + R_n))/C_n;
    end
    
    % internodes
    for i=3:2:N_comp
        dY(i) = (-g_m * Y(i) + 2*(Y(i-1)-Y(i))/(R_in + R_n) + 2*(Y(i+1)-Y(i))/(R_in + R_n) + ...
                2*(V_e(i-1)-V_e(i))/(R_in + R_n) + 2*(V_e(i+1)-V_e(i))/(R_in + R_n) )/C_in;
            
        act_fn(i-1) = (2*(Ve_pulse(i-1)-Ve_pulse(i))/(R_in + R_n) + ...
                    2*(Ve_pulse(i+1)-Ve_pulse(i))/(R_in + R_n))/C_in;
    end
    end
% ----------------------------------------------------------------

% plotting continued
% ----------------------------------------------------------------
% plot Ve_pulse & Ve_rand (input from FE simulation) & activating function
% along the fiber

% subplot(1,4,1);
% % V_e values from VC model (out-of-date)
% %plot(Ve_rand/1000,s, '.k', 'Markersize',5);
% %hold on
% plot(Ve_pulse/1000,xi, 'Color',[0 .75 0]);
% xlabel('[V]','FontSize',14); 
% ylabel('Distance along fiber [cm]','FontSize',14);
% hold off
% axis ij
% % align with response plot
% factor = (xi(N_comp+2)-l_in/2)/(3*(N_comp-1));
% set(gca,'YTick',[0 5 10 15 20],'FontSize',14); 
% % for comparison of Ve profiles
% %set(gca,'XLim',[-0.5 0.5]);
% set(gca,'YLim',[-100*factor (xi(N_comp+2)-l_in/2)+50*factor]);
% % [-3.335859649122815 20.682329824561457]);
% title('Extracellular voltage');
% 
% subplot(1,4,2);
% plot(act_fn,xi(2:N_comp+1), 'Color',[.95 0 0]);
% xlabel('[mV/ms]','FontSize',14); %ylabel('distance [cm]');
% axis ij
% set(gca,'YTick',[0 5 10 15 20],'YTickLabel',{},'YLim',[-100*factor (xi(N_comp+2)-l_in/2)+50*factor]);
% set(gca,'XTick',[-2000 0 2000],'XLim',[-4500 4500],'FontSize',14);
% title('Activating function');
end