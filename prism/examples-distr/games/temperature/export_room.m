function [] = export_room(fid, P, roomID, Tset, X_rep, neighborID1, neighborID2, comm_seq, error)

if(neighborID2 < 1) % only one neighbour
    ns = 1;
else
    ns = 2;
end

if(ns==1)
    n1 = size(P,1);
    n2 = size(P,3);
    n3 = 1;
    m1 = size(P,4);
    m2 = size(P,5);
else
    n1 = size(P,1);
    n2 = size(P,3);
    n3 = size(P,4);
    m1 = size(P,5);
    m2 = size(P,6);
end

fprintf(fid,'\r\n////////////////////// ROOM %u ///////////////////////////////\r\n', roomID);
fprintf(fid,'\r\nsystem "S%u"\r\n', roomID);
fprintf(fid,'\tG%u\r\n', roomID);
fprintf(fid,'endsystem\r\n');

fprintf(fid,'\r\nmodule G%u', roomID);
fprintf(fid,'\r\n\t// discretisation error estimate: %f\r\n\r\n', error);

% Definition of the states
fprintf(fid,'\t// state variables\r\n');
fprintf(fid,'\tx%u : [1..%u] init 1; // room %u temperature\r\n', roomID, n1, roomID);
fprintf(fid,'\ty%u : [1..%u] init 1; // ambient temperature (room %u)\r\n', roomID, n2, neighborID1);
if ns==2
    fprintf(fid,'\tz%u : [1..%u] init 1; // ambient temperature (room %u)\r\n', roomID, n3, neighborID2);
end
fprintf(fid,'\tv%u  : [0..%u] init 0;  // valve setting\r\n', roomID, m1-1);
fprintf(fid,'\tw%u  : [0..%u] init 0;  // window setting\r\n', roomID, m2-1);
fprintf(fid,'\tp%u  : [1..%u] init %u;  // stage\r\n', roomID, 5+ns, 2);

% work out communication sequence
c1 = 4+find(comm_seq==roomID);
c2 = 4+find(comm_seq==neighborID1);
c3 = 4+find(comm_seq==neighborID2);
if(ns==1)
    if ((c1==6) && (c2==7))
        c1=5; c2=6;
    elseif ((c1==7) && (c2==6))
        c1=6; c2=5;
    elseif ((c1==5) && (c2==7))
        c2=6;
    elseif ((c1==7) && (c2==5))
        c1=6;
    end
end

% Stage 4: player 2 sets window
fprintf(fid,'\r\n\t// Stage 4: set window\r\n');
for u2 = 0:m2-1
    fprintf(fid,'\t[a%u_w%u?]    p%u=4 -> (w%u''=%u) & (p%u''=5);\r\n',roomID,u2, roomID, roomID,u2, roomID);
end

% Stage 3: player 1 sets valve
fprintf(fid,'\r\n\t// Stage 3: set valve\r\n');
for u1 = 0:m1-1
    fprintf(fid,'\t[a%u_v%u!]    p%u=3 -> (v%u''=%u) & (p%u''=4);\r\n',roomID,u1, roomID, roomID,u1, roomID);
end

% Stage c1: player 2 sends
fprintf(fid,'\r\n\t// Stage %u: send temperature\r\n', c1);
for i1 = 1:n1
    % Stage send_stage: send over interface
    fprintf(fid,'\t[temp%u_%02u?] x%u=%u & p%u=%u -> (p%u''=%u);\r\n',roomID,i1, roomID,i1, roomID,c1, roomID,mod(c1,ns+5)+1);
end

% Stage c2: player 2 receives from neighbour 1
fprintf(fid,'\r\n\t// Stage %u: receive temperature\r\n', c2);
for i2 = 1:n2
    fprintf(fid,'\t[temp%u_%02u?] p%u=%u  -> (y%u''=%u) & (p%u''=%u);\r\n', neighborID1,i2, roomID,c2, roomID,i2, roomID,mod(c2,ns+5)+1);
end

% Stage c3: player 2 receives from neighbour 2
if ns==2
    fprintf(fid,'\r\n\t// Stage %u: receive temperature\r\n', c3);
    for i2 = 1:n2
        fprintf(fid,'\t[temp%u_%02u?] p%u=%u  -> (y%u''=%u) & (p%u''=%u);\r\n', neighborID2,i2, roomID,c3, roomID,i2, roomID,mod(c3,ns+5)+1);
    end
end

% Stage 1: temperature dynamics
fprintf(fid,'\r\n\t// Stage 1: temperature dynamics\r\n');
for i1=1:n1
    for i2=1:n2
        for i3=1:n3
            for u2 = 0:m2-1
                for u1 = 0:m1-1
                    % use rounded distribution
                    if ns==1
                        prob = round(P(1:n1,i1,i2,u1+1,u2+1)*1000)/1000;
                        fprintf(fid,'\t[temp%u?]    x%u=%u & y%u=%u & v%u=%u & w%u=%u & p%u=1 -> ',roomID, roomID,i1, roomID,i2, roomID,u1, roomID,u2, roomID);
                    elseif ns==2
                        prob = round(P(1:n1,i1,i2,i3,u1+1,u2+1)*1000)/1000;
                        fprintf(fid,'\t[temp%u?]    x%u=%u & y%u=%u & z%u=%u & v%u=%u & w%u=%u & p%u=1 -> ',roomID, roomID,i1, roomID,i2, roomID,i3, roomID,u1, roomID,u2, roomID);
                    end
                    % re-normalise rounded distribution
                    if(sum(prob)~=1.0)
                        [val,index] = max(prob);
                        prob(index) = val+(1.0-sum(prob));
                    end

                    % Stage 5: set temperature (one action per state)
                    trans_added = false;
                    for i1n=1:n1
                        if prob(i1n) > 0.0
                            if trans_added
                                fprintf(fid, ' + ');
                            end
                            trans_added = true;
                            fprintf(fid,'%.3f : (x%u''=%u) & (p%u''=2)',prob(i1n), roomID,i1n, roomID);
                        end
                    end
                    fprintf(fid,';\r\n');
                end
            end
        end
    end
    fprintf(fid,'\r\n');
end

% Stage 2: tau for normal form
fprintf(fid,'\r\n\t// Stage 2: tau for normal form\r\n');
fprintf(fid,'\t[tau?]      p%u=2        -> (p%u''=3);\r\n',roomID, roomID);

fprintf(fid,'\r\nendmodule\r\n\r\n');


% Rewards

% temperature of room
fprintf(fid,'\r\n// temperature room %u\r\n', roomID);
fprintf(fid,'rewards "temp%u"\r\n', roomID);
for i1=1:n1
    fprintf(fid,'\t[temp%u_%02u] true : %.3f;\r\n', roomID, i1, X_rep(roomID,i1));
end
fprintf(fid,'endrewards\r\n\r\n');

% temperature deviation of room
fprintf(fid,'\r\n// temperature deviation room %u\r\n', roomID);
fprintf(fid,'rewards "tempdev%u"\r\n', roomID);
for i1=1:n1
    fprintf(fid,'\t[temp%u_%02u] true : %.3f;\r\n', roomID, i1, (Tset - X_rep(roomID,i1))^2);
end
fprintf(fid,'endrewards\r\n\r\n');

% valve setting
fprintf(fid,'\r\n// valve room %u open\r\n', roomID);
fprintf(fid,'rewards "valve%u"\r\n', roomID);
for i1=0:m1-1
    fprintf(fid,'\t[a%u_v%u] true : %.3f;\r\n', roomID, i1, i1/(m1-1));
end
fprintf(fid,'endrewards\r\n\r\n');

% window setting
fprintf(fid,'\r\n// window room %u open\r\n', roomID);
fprintf(fid,'rewards "window%u"\r\n', roomID);
for i2=0:m2-1
    fprintf(fid,'\t[a%u_w%u] true : %.3f;\r\n', roomID, i2, i2/(m2-1));
end
fprintf(fid,'endrewards\r\n\r\n');
