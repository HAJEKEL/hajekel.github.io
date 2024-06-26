Updates to keyboard shortcuts … On Thursday, 1 August 2024,
    Drive keyboard shortcuts will be updated to give you first - letter navigation.Learn more clear,
    clc % Weights calculated by rowan : % constants L_nw = 12.5; %Distance from nose to wing
%Determining the position for the hydrogen tank. 
Cgopt = 13.37;
g = 9.81;
dp = 0.6;
Range = 800;
% Load calculations N_passengers = 83;
W_passengers = 90 * N_passengers;
W_luggage = 20 * N_passengers;
N_crew = 2;
N_attendants = 3;
W_crew = 93 * N_crew;
W_attendants = 68 * N_attendants;
W_load = W_passengers + W_luggage + W_crew + W_attendants;
% Fuselage dimensions Df = 3.5;
Hf = 7.5;
Lf = 26.55;
N_ult = 3.7159;
V_D = (550 / 3.6).*1.25;
% Propulsion Th = 15062.19914369373;
% lbs N_eng = 2;
% Cargo floor wc = 3.16;
Lfc = 3.10;
Lrc = 2;
Sfrc = wc.*Lfc;
Src = wc.*Lrc;
rho_c = 160;
A_c = 2;
% Bulkhead dimensions D_frbh = 1.5;
D_rbh = 3.5;
% Wing dimensions t_divided_by_c = 0.447;
% average wing thickness ratio cwr = 4.47;
% wing root chord b = 30;
% Wing span lambda_1_2 = 0;
% cosine wing sweep angle S_ref = 90;
% reference area(total area wing) % Tail dimensions S_h = 17.47;
% horizontal tail area S_v = 13.97;
% vertical tail area lambda_1_4 = 0;
% cosine wing sweep angle

    % Initial weight estimations % Max take off weight MTOW(1) = 35000;
% Max zero fuel weight MZFW(1) = 23385;
% Operating weight estimation OEW(1) = 13500;
%Iterative process
for k = 1:10
%Tank
W_tank(k) = 645;
W_fuel(k) = 1590;
% Cargo W_frcf = 0.3074. * sqrt(200).*Sfrc ^ (1.045);
W_rcf = 0.3074. * sqrt(200).*Src ^ (1.045);
W_c1 = rho_c.*A_c.*Lfc;
W_c2 = rho_c.*A_c.*Lrc;
% Wing W_wing(k) = 0.00667 * (N_ult ^ 0.55) * (((t_divided_by_c)*cwr) ^ -0.3) *
                   ((b / cos(lambda_1_2)) ^ 1.05) * (1 + sqrt(1.905 * cos(lambda_1_2) / b)) *
                   ((MZFW(k) / S_ref) ^ -0.3) * MZFW(k);
W_wstr(k) = 20.4 + 0.000907 * 3.75. * MTOW(k);
% Bulkheads W_frbh = 9.1 + 7.225. * (dp./ g).^ (0.8).*(pi.*(D_frbh./ 2).^ 2).^ 1.2;
W_rbh = 9.1 + 7.225. * (dp./ g).^ (0.8).*(pi.*(D_rbh./ 2).^ 2).^ 1.2;
% Landing gear W_nlg(k) = 0.1 + 0.082. * (MTOW(k).^ (0.75)) +
                          (2.97. * 10 ^ (-6)).*(MTOW(k) ^ (1.5));
W_mlg(k) = 18.1 + 0.131. * (MTOW(k).^ (0.75)) + 0.019. * MTOW(k) +
           (2.23. * 10 ^ (-5)).*(MTOW(k) ^ (1.5));
% Fuselage W_fus(k) = 0.0737. *
                          (2. * Df * (V_D ^ 0.338).*(Lf.^ 0.857).*((MTOW(k).*N_ult) ^ 0.286)).^
                      1.1;
% Furniture W_fur(k) = 0.196. * (MZFW(k).^ 0.91);
W_paint(k) = 0.006. * MZFW(k);
% Propulsion W_eng = 0.4054. * (Th ^ 0.9255);
W_pro = 1.377. * W_eng.*N_eng;
% total propulsion W_nac = 0.055. * Th.*N_eng;
% nacelle W_proptot = (W_pro + W_nac).*0.453592;
% 0.453 from lb to kg % Tail W_tail = 0.051 * ((V_D * (S_h + S_v) ^ 1.2) / (sqrt(cos(lambda_1_4))));
% Systems W_oxy = 30 + 1.2 * N_passengers;
% Oxygen W_acid = 14. * ((Lf - 5.9 - 4).^ 1.28);
% Airconditioning and antiicing W_ins(k) = 0.347. * ((MTOW(k)./ 2).^ (5 / 9)).*(Range./ 1000).^
                                           0.25;
% instrument W_apu(k) = 2.2 * 0.001. * MTOW(k);
% APU unit W_hyd(k) = 0.015 * (MTOW(k)./ 2) + 272;
% Hydrolics W_sys(k) = W_oxy + W_acid + W_ins(k) + W_apu(k) + W_hyd(k);
% final values OEW(k + 1) = W_sys(k) + W_proptot + W_fur(k) + W_tank(k) + W_wstr(k) +
                            W_fus(k) + W_frcf + W_rcf + W_frbh + W_rbh;
MZFW(k + 1) = W_c1 + W_c2 + W_load + OEW(k + 1);
MTOW(k + 1) = MZFW(k + 1) + W_fuel(k);
end

    M_tot = (Cgopt - 1).*W_frbh * g - (26.55 - 5.9 - Cgopt).*W_rbh * g +
            (Cgopt - 5.9).*(W_sys(k) - W_hyd(k) - W_apu(k)) * g + (Cgopt - 6.65).*W_nlg(k) * g +
            (Cgopt - (Lf./ 2)).*(W_fus(k) + W_paint(k) + W_fur(k)) * g +
            (Cgopt - (L_nw - (3.10 / 2))).*(W_frcf + W_c1) * g -
            ((L_nw + 4.5 + 1) - Cgopt).*(W_rcf + W_c2) * g + ((Cgopt - L_nw).*W_proptot * g) -
            (L_nw + 2.25 - Cgopt).*(W_wing(k) + W_wstr(k) + W_mlg(k)).*g -
            (26.55 - (5.9 / 2) - Cgopt).*W_tail * g;
% syms x L_x1 = vpasolve(M_tot == x.*(W_fuel + W_tank) * g);
% Distance from CoG to force of hydrogen % Without fuel syms Cgopt L_x2 = vpasolve(
    (Cgopt - 1).*W_frbh * g - (26.55 - 5.9 - Cgopt).*W_rbh * g +
        (Cgopt - 5.9).*(W_sys(k) - W_hyd(k) - W_apu(k)) * g + (Cgopt - 6.65).*W_nlg(k) * g +
        (Cgopt - (Lf./ 2)).*(W_fus(k) + W_paint(k) + W_fur(k)) * g +
        (Cgopt - (L_nw - (3.10 / 2))).*(W_frcf + W_c1) * g -
        ((L_nw + 4.5 + 1) - Cgopt).*(W_rcf + W_c2) * g + ((Cgopt - L_nw).*W_proptot * g) -
        (L_nw + 2.25 - Cgopt).*(W_wing(k) + W_wstr(k) + W_mlg(k)).*g -
        (26.55 - (5.9 / 2) - Cgopt).*W_tail * g ==
    L_x1.*W_tank * g);
% Distance from CoG to force of hydrogen.

    Cl = 0.5630413647;
W = MTOW(end).*9.81;
rho = 1.225;
% density at sea level S = 90;
% Wing area V_cruise = 550. / 3.6;
n_max = 2.1 + (24000. / (MTOW(end) + 10000));
n_min = -0.4. * n_max;
n = 3.7159;
% max load factor V = 550 / 3.6;
% Velocity in m / s Cm = -0.04; %Pitch moment coefficient from wing (should be whole plane I think but I do not know how to get it for the fuselage)

M_0 = (1./2).*rho.*S.*(V.^2).*Cm;

% Pitch moment % Calculating the lift forces with maximum load factor L_tail =
    (-((0.095. * W)./ n) - M_0).*(1 / 10.935);
L_wing = (W - n.*L_tail).*(1 / n);

% Matlab script Henk % First we just want to find the maximum deflection which is obviously %
    located at the end of the cantilever beam.We do this using the method % of superposition,
    adding all maximum deflections caused by each point % force along the length of the fuselage.

        % For a concentrated load at any point on the span of a cantilever beam we %
        can define : % Maximum moment : P = 1;
% point load in newton a = 1;
% distance of the point force from the imaginairy wall in meters % M = -P.*a;
% Maximum moment

    % slope at the end of the beam : E = 1;
% Elasticity modulus in N / (m ^ 2) I = 1;
% Second moment of area in m ^ 4 E = (73.1. * 10 ^ 9) I = 0.00412335 % theta =
                                                              (P.*(a.^ 2)) / (2. * E.*I);

% Maximum deflection : L = 1;
% Length of the beam delta = (P.*(a.^ 2) / (6. * E.*I)).*(3. * L - a);
% DEZE GEBRUIK IK, HEB IK HIERONDER OOK GECOPIERD IN DIE FOR LOOP % Deflection equation x = 1;
% distance from the imaginairy wall in meters % y = ((P.*(x.^ 2))./ 6).*(3. * a - x); %deflection y along the x axis for 0<x<a
%y = ((P.*(a.^2))./6).*(3.*x-a); %deflection y along the x axis for a<x<L

x_length = sort([0,1,5.9, 6.65,(L_nw -(3.10/2)), L_nw, Cgopt-L_x1,(Lf./2), 13.5,13.625, (L_nw + 2.25),(L_nw+4.5+1),26.55-5.9,26.55 - (5.9/2) , 26.55]);

% From this we first want to find the maximum deflection caused by all %
    forces : % For the left side we have :

    %
    These are the weights acting on the fuselage with weight 1 acting close to %
    the imaginairy wall and weight 7 close to the nose of the fuselage.g = 9.81; %the gravitational acceleration in m/(s
w1 = (W_fus(k) + W_paint(k) + W_fur(k) + (W_load-W_luggage))*g; %weight left 1
w2 = (W_fuel(k)+W_tank(k))*g; % weight left 2
w3 = W_proptot*g; %weight left 3
w4 = (W_frcf + W_c1)*g; %weight left 4
w5 = W_nlg(k)*g; %weight left 5
w6 = (W_sys(k)-W_hyd(k)-W_apu(k))*g; %weight left 6
w7 = W_frbh*g;% weight left 7

%We make a vector of these weights: 
w_left = [w1, w2, w3, w4, w5, w6, w7];

%Similarly these are the lengths with length 1 close to the imaginairy wall
%and length 7 close to the nose of the aircraft.
%The vector x_length defines the position of the forces in y direction
%along the fuselage, with the first index just being the 0 point and the 
%last index the end point of the fuselage. 
l_1 = 13.5 - x_length(8); %length between the imaginairy wall and the 
%closest force on the left side of the wall
l_2 = 13.5 - x_length(7); 
l_3 = 13.5 - x_length(6);
l_4 = 13.5 - x_length(5);
l_5 = 13.5 - x_length(4);
l_6 = 13.5 - x_length(3);
l_7 = 13.5 - x_length(2); %length between the imaginairy wall and the 
%furthest force on the left side of the wall

%We make a vector of these lengths
l_left = [l_1, l_2, l_3, l_4, l_5, l_6, l_7];
L = 13.5;
delta = zeros(1,7);
for k=1:7
    d = (w_left(k).*(l_left(k).^2)./(6.*E.*I)).*(3.*L-l_left(k));
    delta(k) = d;
end
delta
a = sum(delta)

%These are the weights acting on the fuselage with weight 9 acting close to
%the imaginairy wall and weight 13 close to the tail of the fuselage.

%For the right side we have: 
w9 = -L_wing; %weight right 1 !!!!!!!!!!!!!!! L_wing = L_wing(k)
w10 = (W_wing(k)+W_wstr(k) + W_mlg(k)).*g; %weight right2
w11 = (W_rcf + W_c2).*g; %weight right 3
w12 = W_rbh*g; %weight right 4
w13 = W_tail*g - L_tail; %weight right 5 !!!!!!!!!!!!!L_tail = L_tail(k)

%We make a vector of these weights: 
w_right = [w9, w10, w11, w12, w13];

%similarly these are the lengths with length 9 close to the imaginairy wall
%and length 13 close to the tail of the aircraft.
l_9  = x_length(10) - 13.5; %length between the imaginairy wall and the 
%closest force on the right side of the wall
l_10 = x_length(11) - 13.5;
l_11 = x_length(12) - 13.5;
l_12 = x_length(13) - 13.5;
l_13 = x_length(14) - 13.5; %length between the imaginairy wall and the 
%furthest force on the right side of the wall

%We make a vector of these lengths

l_right = [l_9,l_10,l_11,l_12,l_13];


L = 26.55 - 13.5;
delta = zeros(1,7);
for k=1:5
    d = (w_right(k).*(l_right(k).^2)./(6.*E.*I)).*(3.*L-l_right(k));
    delta(k) = d;
end
delta
a = sum(delta)



%lets plot the weight at eacht point along the fuselage

figure(1),clf(1),hold on
plot(l_left,w_left)
title('left side of fuselage')
xlabel('length from imaginairy wall towards nose')
ylabel('weight')
hold off

figure(2),clf(2),hold on
plot(l_right,w_right)
title('right side of fuselage')
xlabel('length from imaginairy wall towards tail')
ylabel('weight')
hold off
