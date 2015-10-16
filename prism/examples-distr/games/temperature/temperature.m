% Temperature Control Case Study

clc; close all;

%================= Discretisation ===========================
rooms = 3;
w = [3, 3, 3];     % room width (m)
l = [3, 3, 3];     % room length (m)
h = [2, 2, 2];     % room height (m)
Aw = 1.4 * 1.4;    % window size (m^2)
%    w       w      w
% |------|-------|------|
% |  1   |   2   |  3   | l
% |------|-------|------|
%

V_el = 10*10^-3;   % chilled water volume (10 liters)
V = V_el;

r = 0.01;               % pipeline radius
A = 2 * pi * r * w;     % surface area of cooler
h_copp = 13.1;          % overall cooper transfer coefficient [W/m^2K]
Kcw = A * h_copp;       % chilled water thermal convection coefficient

rho_a = 1.2041;   % air density (20°C) [kg/m^3]
cw_a = 1.012e3;   % specific heat capacity of air [J/kgK]
Cr = (w .* l .* h) * rho_a * cw_a;  % room thermal capacity

% thermal convection coefficients
K1 = ([1,1,1].*w.*h * 2 + [1,0,1].*l.*h) * 0.1;  % closed window
K2 = [1,1,1]* Aw * 0.4;                          % change for open window
Kwall = [1,1,1].*w.*h * 0.2;                     % inner wall

Tset = 20;        % Temperature setpoint
Tout = 30;        % Nominal outside temperature
Tcw  = 10;        % Nominal cooling water temperature

dt=20*60; % period length in seconds
VarTA=1/65^2; % Disturbance variance
Sigma=VarTA*dt; % Resulting Sigma

% truncating the state space to boxes
thr = [2, 2, 2];
X_l = Tset-thr;
X_u = Tset+thr;
% bins per dimension
bins = 5;
n = [bins,bins,bins];
% diameter of bins
delta = (X_u-X_l)./n;
% boundary points of the partition
X = zeros(rooms,bins+1);
X_rep = zeros(rooms,bins);
for r = 1:rooms
    X(r,:) = X_l(r):delta(r):X_u(r);
    X_rep(r,:) = X(r,1:n(r))+delta(r)/2;
end

% computing transition probabilities
m1 = 5; % control options (valve)
m2 = 2; % control options (window)

% Room 1 -- neighbour is room 2
r = 1;
P1 = zeros(n(r),n(r),n(2),m1,m2);
for u2 = 0:m2-1
    for u1 = 0:m1-1
        for ir=1:n(r)
            for i2=1:n(2)
                % Expected Value of the Room Temperature
                E_xbar=X_rep(r,ir) ...                                   % previous temperature
                    + (dt/Cr(r))*Kwall(r)*(X_rep(2,i2)-X_rep(r,ir)) ...  % change due to neighbouring room
                    + (dt/Cr(r))*u1/(m1-1)*Kcw(r)*(Tcw-X_rep(r,ir)) ...  % change due to valve setting
                    + (dt/Cr(r))*(K1(r)+K2(r)*u2)*(Tout-X_rep(r,ir));      % change due to window setting
                T = normcdf(X(r,:),E_xbar,Sigma);
                tn = (T(2:n(r)+1) - T(1:n(r)));
                P1(1:n(r),ir,i2,u1+1,u2+1) = tn/sum(tn); % normalization because of truncation
            end
        end
    end
end
% error for cost function C(x_1) = [x_1-Ts]^2
h_c1 = 2*max(abs(X_l(r)-Tset),abs(X_u(r)-Tset)); % differential within box
alpha1 = 1 - (dt/Cr(r))*Kwall(r) - (dt/Cr(r))*Kcw(r) - (dt/Cr(r))*(K1(r)+K2(r));
alpha2 = (dt/Cr(r))*Kwall(r);
error1 = h_c1*delta(r) + h_c1*(abs(alpha1)*delta(r) + abs(alpha2)*delta(2))/(1-abs(alpha1));


% Room 2 -- neighbours are room 1 and 3
r = 2;
P2 = zeros(n(r),n(r),n(1),n(3),m1,m2);
for u2 = 0:m2-1
    for u1 = 0:m1-1
        for ir=1:n(r)
            for i1=1:n(1)
                for i3=1:n(3)
                    % Expected Value of the Room Temperature
                    E_xbar=X_rep(r,ir) ...                                   % previous temperature
                        + (dt/Cr(r))*Kwall(r)*(X_rep(1,i1)-X_rep(r,ir)) ...  % change due to neighbouring room 1
                        + (dt/Cr(r))*Kwall(r)*(X_rep(3,i3)-X_rep(r,ir)) ...  % change due to neighbouring room 3
                        + (dt/Cr(r))*u1/(m1-1)*Kcw(r)*(Tcw-X_rep(r,ir)) ...  % change due to valve setting
                        + (dt/Cr(r))*(K1(r)+K2(r)*u2)*(Tout-X_rep(r,ir));      % change due to window setting
                    T = normcdf(X(r,:),E_xbar,Sigma);
                    tn = (T(2:n(r)+1) - T(1:n(r)));
                    P2(1:n(r),ir,i1,i3,u1+1,u2+1) = tn/sum(tn); % normalization because of truncation
                end
            end
        end
    end
end
% error for cost function C(x_2) = [x_2-Ts]^2
h_c2 = 2*max(abs(X_l(r)-Tset),abs(X_u(r)-Tset)); % differential within box
alpha1 = (dt/Cr(r))*Kwall(r);
alpha2 = 1 - 2*(dt/Cr(r))*Kwall(r) - (dt/Cr(r))*Kcw(r) - (dt/Cr(r))*(K1(r)+K2(r));
alpha3 = (dt/Cr(r))*Kwall(r);
error2 = h_c2*delta(r) + h_c2*(abs(alpha1)*delta(1) + abs(alpha2)*delta(r) + abs(alpha3)*delta(3))/(1-abs(alpha2));

% Room 3 -- neighbour is room 2
r = 3;
P3 = zeros(n(r),n(r),n(2),m1,m2);
for u2 = 0:m2-1
    for u1 = 0:m1-1
        for ir=1:n(r)
            for i2=1:n(2)
                % Expected Value of the Room Temperature
                E_xbar=X_rep(r,ir) ...                                  % previous temperature
                    + (dt/Cr(r))*Kwall(r)*(X_rep(2,i2)-X_rep(r,ir)) ... % change due to neighbouring room
                    + (dt/Cr(r))*u1/(m1-1)*Kcw(r)*(Tcw-X_rep(r,ir)) ... % change due to valve setting
                    + (dt/Cr(r))*(K1(r)+K2(r)*u2)*(Tout-X_rep(r,ir));     % change due to window setting
                T = normcdf(X(r,:),E_xbar,Sigma);
                tn = (T(2:n(r)+1) - T(1:n(r)));
                P3(1:n(r),ir,i2,u1+1,u2+1) = tn/sum(tn); % normalization because of truncation
            end
        end
    end
end
% error for cost function C(x_3) = [x_3-Ts]^2
h_c3 = 2*max(abs(X_l(r)-Tset),abs(X_u(r)-Tset)); % differential within box
alpha2 = (dt/Cr(r))*Kwall(r);
alpha3 = 1 - (dt/Cr(r))*Kwall(r) - (dt/Cr(r))*Kcw(r) - (dt/Cr(r))*(K1(r)+K2(r));
error3 = h_c3*delta(r) + h_c3*(abs(alpha2)*delta(2) + abs(alpha3)*delta(r))/(1-abs(alpha3));

%=======================Export to PRISM==========================

% MODEL FILE
fid = fopen('temperature.prism','w');

% smg preamble
fprintf(fid,'smg\r\n');

% top level system
fprintf(fid,'\r\n// top level system\r\n');
fprintf(fid,'system\r\n');
fprintf(fid,'\t"S1" || "S2" || "S3" \r\n');
fprintf(fid,'endsystem\r\n');

comm_seq = [1,2,3];
export_room(fid, P1, 1, Tset, X_rep, 2, -1, comm_seq, error1);
export_room(fid, P2, 2, Tset, X_rep, 1,  3, comm_seq, error2);
export_room(fid, P3, 3, Tset, X_rep, 2, -1, comm_seq, error3);

% time progress
fprintf(fid,'\r\n// time\r\n');
fprintf(fid,'rewards "time"\r\n');
for i1=1:n(1)
    fprintf(fid,'\t[temp2_%02u] true : 1;\r\n', i1);
end
fprintf(fid,'endrewards\r\n\r\n');

fclose(fid);

% PROPERTIES FILE
pid = fopen('temperature.props','w');
for r=1:rooms
    fprintf(pid, 'const double tv%u = 0.8;\r\n', r);
    fprintf(pid, 'const double va%u = 0.8;\r\n', r);
    fprintf(pid, 'const double tw%u = 0.3;\r\n', r);
    if r==2
        fprintf(pid, 'const double wi%u = 0.2;\r\n', r);
    else
        fprintf(pid, 'const double wi%u = 0.3;\r\n', r);
    end
end
fprintf(pid, '\r\n');

% room 1
fprintf(pid, '"phi1" : <<1>> ((R{"window1"}/{"time"}<=wi1 [ S ] & R{"tempdev2"}/{"time"}<=tw2 [ S ]) => R{"tempdev1"}/{"time"}<=tw1 [ S ])\r\n');
fprintf(pid, '"psi1" : <<1>> (R{"valve1"}/{"time"}<=va1 [ S ] & (R{"tempdev2"}/{"time"}<=tv2 [ S ] => R{"tempdev1"}/{"time"}<=tv1 [ S ]))\r\n');
% room 2
fprintf(pid, '"phi2" : <<1>> (R{"window2"}/{"time"}<=wi2 [ S ] => R{"tempdev2"}/{"time"}<=tw2 [ S ])\r\n');
fprintf(pid, '"psi2" : <<1>> (R{"valve2"}/{"time"}<=va2 [ S ] & R{"tempdev2"}/{"time"}<=tv2 [ S ])\r\n');
% room 3
fprintf(pid, '"phi3" : <<1>> ((R{"window3"}/{"time"}<=wi3 [ S ] & R{"tempdev2"}/{"time"}<=tw2 [ S ]) => R{"tempdev3"}/{"time"}<=tw3 [ S ])\r\n');
fprintf(pid, '"psi3" : <<1>> (R{"valve3"}/{"time"}<=va3 [ S ] & (R{"tempdev2"}/{"time"}<=tv2 [ S ] => R{"tempdev3"}/{"time"}<=tv3 [ S ]))\r\n');

% compositional
fprintf(pid, '\r\n"phi" : comp("phi1", "phi2", "phi3")\r\n');
fprintf(pid, '"psi" : comp("psi1", "psi2", "psi3")\r\n');

fclose(pid);

% BASH FILE
sid = fopen('temperature.sh','w');
acc = [100, 500, 1000];
eps = [0.05,0.02, 0.01];
citer = [100,250];
diter = [20,10];
name = ['phi'; 'psi'];
fprintf(sid, 'TIMEOUT=240\n\n'); % timeout in minutes
fprintf(sid, 'if hash timeout 2>/dev/null; then\n\tTO=timeout\nelse\n\tTO=gtimeout\nfi\n\n'); % timeout program
fprintf(sid, '{\n'); % sequencing
for j = 1:2
    for i = 1:3
        fprintf(sid, '$TO $((TIMEOUT))m \\\n'); % timeout
        fprintf(sid, '../../prism/bin/prism temperature{.prism,.props} \\\n'); % model and properties
        fprintf(sid, '\t-prop %u \\\n', 6+j); % property index
        fprintf(sid, '\t-multimaxciter %u -multimaxditer %u -gs \\\n', citer(j), diter(j)); % iteration bounds and Gauss-Seidel
        fprintf(sid, '\t-multiminm 2 -multimaxm 10000 \\\n'); % box size
        fprintf(sid, '\t-baselineaccuracy %u -increasefactor 1.0 -multirounding \\\n', acc(i)); % rounding
        fprintf(sid, '\t-logcpareto -logdpareto \\\n'); % logging
        fprintf(sid, '\t-exportstrat temp_%s_%u.strat \\\n', name(j,:), i); % strategy export
        fprintf(sid, '\t-paretoepsilon %.2f \\\n', eps(i)); % stopping accuracy
        fprintf(sid, '\t2>&1 1> temp_%s_%u.log ; \n\n', name(j,:), i);
    end
end
fprintf(sid, '} &\n'); % end sequencing and execute in background
fclose(sid);

%=====================End of the code==============================

