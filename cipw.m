%% function for CIPW norm conversion
%oxidewt is a 1 by 11 vector with the weight percentages of the main 11
%oxides in the following order:

    
function norm = cipw()
%% set up structures, load data
%filename = input('.mat filename of oxide data: ', 's');
filename = 'oxides';
oxides = load([filename '.mat']);

%can i just declare an empty struct and add fields as used in the steps
%below? that would save having the rare elements always show up...
norm = struct();

% 'or', 0, 'ab', 0, 'an', 0,...
%                 'diwo', 0, 'dien', 0, 'difs', 0,...
%                 'hyen', 0, 'hyfs', 0, 'Q', 0);


%% convert oxidewt to mole percent:       
oxides.SiO2 = oxides.SiO2 / 60.0843;
oxides.TiO2 = oxides.TiO2 / 79.8988;
oxides.Al2O3 = oxides.Al2O3 / 101.9613;
oxides.Cr2O3 = oxides.Cr2O3 / 151.9904;
oxides.FeO = oxides.FeO / 71.8464;
oxides.Fe2O3 = oxides.Fe2O3 / 159.6922;
oxides.MnO = oxides.MnO / 70.9374;
oxides.MgO = oxides.MgO / 40.3044;
oxides.CaO = oxides.CaO / 56.0794;
oxides.Na2O = oxides.Na2O / 61.9789;
oxides.P2O5 = oxides.P2O5 / 141.9445;
oxides.K2O = oxides.K2O / 94.1960;
oxides.S = oxides.S / 32.065;

%% trace element substitution
%A:
oxides.FeO = oxides.FeO + oxides.MnO; 
%B:
%would include if had BaO or SrO;


%% perform norm calculation steps

% Step 1: apatite
%requires CaO > (10.3)*P2O5
if oxides.P2O5 > 0
norm.ap = oxides.P2O5;
oxides.CaO = oxides.CaO - (10/3)*oxides.P2O5;
oxides.P2O5 = oxides.P2O5 - oxides.P2O5;
end

% Step 2: pyrite
%requires FeO > .5*S
if oxides.S > 0
norm.py = oxides.S;
oxides.FeO = oxides.FeO - .5*oxides.S;
oxides.S = oxides.S - oxides.S;
end

% Step 3: chromite
%requires FeO > Cr2O3
if oxides.Cr2O3 > 0
norm.cm = oxides.Cr2O3;
oxides.FeO = oxides.FeO - oxides.Cr2O3;
oxides.Cr2O3 = oxides.Cr2O3 - oxides.Cr2O3;
end

% Step 4: titanium stuff...
if oxides.FeO >= oxides.TiO2 && oxides.TiO2 > 0
norm.il = oxides.TiO2;
oxides.FeO = oxides.FeO - oxides.TiO2;
oxides.TiO2 = oxides.TiO2 - oxides.TiO2; %Ti is done
else % moreti, so fe is limiting oxide
norm.il = oxides.FeO;
oxides.FeO = oxides.FeO - oxides.FeO;
oxides.TiO2 = oxides.TiO2 - oxides.FeO;  %could still be some Ti left over
end


% Step 5: fluorite
% skip b/c no F...

% Step 6: Calcite
%skip b/c no CO2...

% Step 7: zircon
%skip b/c no ZrO2...

%Step 8: provisional orthoclase
%requires Al > K2O, Si > K2O
por = oxides.K2O;
oxides.Al2O3 = oxides.Al2O3 - oxides.K2O;
oxides.SiO2 = oxides.SiO2 - 6*oxides.K2O;
oxides.K2O = oxides.K2O - oxides.K2O;

% Step 9, 10, 11, and 12: provisional albite, provisional anorthite, acmite
if oxides.Al2O3 >= oxides.Na2O
    % step 9
   pal = oxides.Na2O;
   oxides.SiO2 = oxides.SiO2 - 6*oxides.Na2O;
   oxides.Al2O3 = oxides.Al2O3 - oxides.Na2O;
   oxides.Na2O = oxides.Na2O - oxides.Na2O;
   
   %step 10
   if oxides.Al2O3 <= oxides.CaO
       pan = oxides.Al2O3;
       oxides.SiO2 = oxides.SiO2 - 2*oxides.Al2O3;
       oxides.CaO = oxides.CaO - oxides.Al2O3;
       oxides.Al2O3 = oxides.Al2O3 - oxides.Al2O3;
   else %CaO > Al2
       pan = oxides.Cao;
       oxides.SiO2 = oxides.SiO2 - 2*oxides.CaO;
       oxides.Al2O3 = oxides.Al2O3 - oxides.CaO;
       oxides.CaO = oxides.Cao - oxides.CaO;
       
       norm.co = oxides.Al2O3;
       oxides.Al2O3 = oxides.Al2O3 - oxides.Al2O3;
   end
else %Na2O > Al2O3...
   pal = oxides.Al2O3;
   oxides.SiO2 = oxides.SiO2 - 6*oxides.Al2O3;
   oxides.Na2O = oxides.Na2O - oxides.Al2O3;
   oxides.Al2O3 = oxides.Al2O3 - oxides.Al2O3; 
   
   if oxides.Na20 >= oxides.Fe2O3
       %step 11a
       norm.ac = oxides.Fe2O3;
       oxides.Na2O = oxides.Na2O - oxides.Fe2O3;
       oxides.SiO2 = oxides.SiO2 - 4*oxides.Fe2O3;
       oxides.Fe2O3 = oxides.Fe2O3 - oxides.Fe2O3; 
       
       %step 12
       norm.sodiummeta = oxides.Na2O;
       oxides.SiO2 = oxides.SiO2 - oxides.Na2O;
       oxides.Na2O = 0;
   
   
   else %Fe2O3 > Na2O
       %step 11b
       norm.ac = oxides.Na2O;
       oxides.SiO2 = oxides.SiO2 - 4*oxides.Na2O;
       oxides.Fe2O3 = oxides.Fe2O3 - oxides.Na2O; 
       oxides.Na2O = 0;
   end
end



% Step 13: magnetite and Hematite
if oxides.Fe2O3  >= oxides.FeO
    norm.mt = oxides.FeO;
    oxides.Fe2O3 = oxides.Fe2O3 - oxides.FeO;
    oxides.FeO = 0;
    
    norm.ht = oxides.Fe2O3;
    oxides.Fe2O3 = 0;
    
else % FeO> Fe2O3
     norm.mt = oxides.Fe2O3;
    oxides.FeO = oxides.FeO - oxides.Fe2O3;
    oxides.Fe2O3 = 0;
end

% Step 14: pyroxenes and olivines
oxides.MF = oxides.MgO + oxides.FeO;
percentMg = oxides.MgO / oxides.MF;

%step 15: provisional diopside

















