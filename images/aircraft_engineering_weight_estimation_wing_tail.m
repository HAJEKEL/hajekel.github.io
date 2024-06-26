clear,clc
%weight wing
%use equation 6 since its the only one without english units
%what should we do with the 1/2, in cos(lambda1/2)?
%initialisation
N_ult =3.5; %ultimate load factor
t_divided_by_c = 0.447; %average wing thickness ratio
cwr = 4.47; %wing root chord
b = 26.83;%Wing span
lambda_1_2 = 0; %cosine wing sweep angle
W_zf = 23385; %Designed zero fuel weight
S_ref = 80; %reference area (total area wing)
%formula for the weight of the wing
W_w = 0.00667*(N_ult^0.55)*(((t_divided_by_c)*cwr)^-0.3)*((b/cos(lambda_1_2))^1.05)*(1+sqrt(1.905*cos(lambda_1_2)/b))*((W_zf/S_ref)^-0.3)*W_zf




%weight tail 
%use equation 12 since its the only one without english units
%initialisation
V_dive = 190.97; %designed dive speed [m/s]
S_h = 17.47;%horizontal tail area
S_v =  13.97;%vertical tail area
lambda_1_4 = 0;%cosine wing sweep angle
%formula for the weight of the tail
W_t = 0.051*((V_dive*(S_h+S_v)^1.2)/(sqrt(cos(lambda_1_4))))