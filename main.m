%% function for CIPW norm conversion
%oxidewt is a 1 by 11 vector with the weight percentages of the main 11
%oxides in the following order:

    
function [oxides, wt] = main(filename)
%% set up structures, load data
% filename = input('.mat filename of oxide data: ', 's');
if filename == 'new'
    filename = oxidewt();
end
oxides = load([filename '.mat']);

%can i just declare an empty struct and add fields as used in the steps
%below? that would save having the rare elements always show up...
norm = struct('ap', 0, 'pr', 0, 'cm', 0, 'il', 0, 'C', 0, 'ac', 0,...
                'ns', 0, 'ru', 0, 'or', 0, 'ab', 0, 'an', 0,...
                'mt', 0, 'hm', 0, 'Q', 0, 'ol', 0, 'pf', 0, 'ne', 0,...
                'tn', 0, 'dien', 0, 'difs', 0, 'diwo', 0, 'hyen', 0,...
                'hyfs', 0);
prov = struct('or', 0, 'ab', 0, 'an', 0, 'tn', 0, 'di', 0, 'wo', 0, 'hy', 0);

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
oxides.MnO = 0;
%B:
%would include if had BaO or SrO;


%% perform norm calculation steps

% Step 1: apatite
%requires CaO > (10.3)*P2O5
if oxides.P2O5 > 0
norm.ap = oxides.P2O5;
oxides.CaO = oxides.CaO - (10/3)*oxides.P2O5;
oxides.P2O5 = 0;
end

% Step 2: prrite
%requires FeO > .5*S
if oxides.S > 0
norm.pr = oxides.S;
oxides.FeO = oxides.FeO - .5*oxides.S;
oxides.S = 0;
end

% Step 3: chromite
%requires FeO > Cr2O3
if oxides.Cr2O3 > 0
norm.cm = oxides.Cr2O3;
oxides.FeO = oxides.FeO - oxides.Cr2O3;
oxides.Cr2O3 = 0;
end

% Step 4: titanium stuff...
if oxides.FeO >= oxides.TiO2 
norm.il = oxides.TiO2;
oxides.FeO = oxides.FeO - oxides.TiO2;
oxides.TiO2 = 0; %Ti is done
else % more ti, so fe is limiting oxide
norm.il = oxides.FeO;
oxides.TiO2 = oxides.TiO2 - oxides.FeO;  %could still be some Ti left over
oxides.FeO = 0;
end


% Step 5: fluorite
% skip b/c no F...

% Step 6: Calcite
%skip b/c no CO2...

% Step 7: zircon
%skip b/c no ZrO2...

%Step 8: provisional orthoclase
%requires Al > K2O, Si > K2O
prov.or = oxides.K2O;
oxides.Al2O3 = oxides.Al2O3 - oxides.K2O;
oxides.SiO2 = oxides.SiO2 - 6*oxides.K2O;
oxides.K2O = 0;

% Step 9, 10, 11, and 12: provisional albite, provisional anorthite, acmite
if oxides.Al2O3 >= oxides.Na2O
    % step 9
   prov.ab = oxides.Na2O;
   oxides.SiO2 = oxides.SiO2 - 6*oxides.Na2O;
   oxides.Al2O3 = oxides.Al2O3 - oxides.Na2O;
   oxides.Na2O = 0;
   
   %step 10
   if oxides.Al2O3 <= oxides.CaO
       prov.an = oxides.Al2O3;
       oxides.SiO2 = oxides.SiO2 - 2*oxides.Al2O3;
       oxides.CaO = oxides.CaO - oxides.Al2O3;
       oxides.Al2O3 = 0;
   else %CaO > Al2
       prov.an = oxides.CaO;
       oxides.SiO2 = oxides.SiO2 - 2*oxides.CaO;
       oxides.Al2O3 = oxides.Al2O3 - oxides.CaO;
       oxides.CaO = 0;
       
       norm.C = oxides.Al2O3;
       oxides.Al2O3 = 0;
   end
else %Na2O > Al2O3...
   prov.ab = oxides.Al2O3;
   oxides.SiO2 = oxides.SiO2 - 6*oxides.Al2O3;
   oxides.Na2O = oxides.Na2O - oxides.Al2O3;
   oxides.Al2O3 = 0; 
   
   if oxides.Na2O >= oxides.Fe2O3
       %step 11a
       norm.ac = oxides.Fe2O3;
       oxides.Na2O = oxides.Na2O - oxides.Fe2O3;
       oxides.SiO2 = oxides.SiO2 - 4*oxides.Fe2O3;
       oxides.Fe2O3 = 0; 
       
       %step 12
       norm.ns = oxides.Na2O;
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

% Step 4b...
if oxides.TiO2 > oxides.FeO 
   if oxides.TiO2 >= oxides.CaO
       prov.tn = oxides.CaO;
       oxides.TiO2 = oxides.TiO2 - oxides.CaO;
       oxides.SiO2 = oxides.SiO2 - oxides.CaO;
       oxides.CaO = 0;
       
       norm.ru = oxides.TiO2;
       oxides.TiO2 = 0;
       
   else % CaO > TiO2
       prov.tn = oxides.TiO2;
       oxides.SiO2 = oxides.SiO2 - oxides.TiO2;
       oxides.TiO2 = oxides.TiO2 - oxides.TiO2;
       oxides.TiO2 = 0;
       
   end
end

% Step 13: magnetite and Hematite
if oxides.Fe2O3  >= oxides.FeO
    norm.mt = oxides.FeO;
    oxides.Fe2O3 = oxides.Fe2O3 - oxides.FeO;
    oxides.FeO = 0;
    
    norm.hm = oxides.Fe2O3;
    oxides.Fe2O3 = 0;
    
else % FeO> Fe2O3
    norm.mt = oxides.Fe2O3;
    oxides.FeO = oxides.FeO - oxides.Fe2O3;
    oxides.Fe2O3 = 0;
end

% Step 14: prroxenes and olivines
oxides.MF = oxides.MgO + oxides.FeO;
percentMg = oxides.MgO / oxides.MF;
oxides.MgO = 0;
oxides.FeO = 0;

%step 15, 16, and 17: CaO and MF
if oxides.CaO >= oxides.MF
    prov.di = oxides.MF;
    oxides.CaO = oxides.CaO - oxides.MF;
    oxides.SiO2 = oxides.SiO2 - 2*oxides.MF;
    oxides.MF = 0;
    
    prov.wo = oxides.CaO;
    oxides.SiO2 = oxides.SiO2 - oxides.CaO;
    oxides.CaO = 0;
    
else % MF > CaO
    prov.di = oxides.CaO;
    oxides.MF = oxides.MF - oxides.CaO;
    oxides.SiO2 = oxides.SiO2 - 2*oxides.CaO;
    oxides.CaO = 0;
    
    prov.hy = oxides.MF;
    oxides.SiO2 = oxides.SiO2 - oxides.MF;
    oxides.MF = 0;
end

%% dealing with SiO2:
if oxides.SiO2 >= 0
% Step 18: Quartz
    norm.Q = oxides.SiO2;
    oxides.SiO2 = 0;
else %oxides.SiO2 < 0
%steps 19 thru 26
    
% Step 19:    
   D = -oxides.SiO2; 

% Step 20:
if D >0
if D >= prov.hy / 2
    norm.ol = prov.hy;
    prov.hy = 0;
    D1 = D - (prov.hy / 2);
else %prov.hy/2 >D
   norm.ol = D;
   prov.hy = prov.hy - 2*D;
   D = 0;
   oxides.SiO2 = -D;
end
end

% Step 21:
if (exist('D1','var') == 1)
    if (exist('prov.tn','var') == 1)
        if prov.tn > D
        norm.pf = D1;
        prov.tn = prov.tn - D1;
        D1 = 0;
        oxides.SiO2 = -D1;
        else % D > pth
        norm.pf = prov.tn;
        D2 = D1 - prov.tn;
        prov.tn = 0;
        end
    else % no prov sphene ( pth )
        D2 = D1;
    end    
end

% Step 22:
if (exist('D2','var') == 1)
    if 4*prov.ab > D2
    norm.ne = D2/4;
    prov.ab = prov.ab - D2/4;
    D2 = 0;
    oxides.SiO2 = -D2;
    else %D2 > 4*prov.ab
    norm.ne = prov.ab;
    D3 = D2 - 4*prov.ab;
    prov.ab = 0;
    end
end

%steps 23 
end



%% final step - move provisional minerals to norm struct
%always exists:

norm.or = prov.or;
prov.or = 0;

norm.ab = prov.ab;
prov.ab = 0;

norm.dien = prov.di*percentMg;
norm.difs = prov.di*(1-percentMg);
norm.diwo = prov.di;
prov.di = 0;

norm.tn = prov.tn;
prov.tn = 0;

norm.an = prov.an;
prov.an = 0;

norm.diwo = norm.diwo + prov.wo;
prov.wo = 0;

norm.hyen = prov.hy*percentMg;
norm.hyfs = prov.hy*(1-percentMg);
prov.hy = 0;




%% convert norms to wt %
wt.ap = norm.ap*336;            %apatite
wt.pr = norm.pr*119.98;         %prrite
wt.cm = norm.cm*223.83;         %chromite
wt.il = norm.il*152;            %ilmenite
wt.C = norm.C*101.96;           %corundum
wt.ac = norm.ac*231;            %acmite/aegirine
wt.ns = norm.ns*122.06;         %sodium metasilicate
wt.ru = norm.ru*79.866;         %rutile
wt.or = norm.or*557;            %orthoclase
wt.ab = norm.ab*524;            %albite
wt.an = norm.an*278;            %anorthite
wt.mt = norm.mt*232;            %magnetite
wt.hm = norm.hm*159.69;         %hematite
wt.Q = norm.Q*60.1;             %quartz
wt.ol = norm.ol*153.31;         %olivine
wt.pf = norm.pf*135.96;         %perovskite
wt.ne = norm.ne*146.08;         %nepheline
wt.tn = norm.tn*197.76;         %titanite/sphene
wt.dien = norm.dien*100;        %diopside enstatite
wt.difs = norm.difs*132;        %diopside forsterite
wt.diwo = norm.diwo*116;        %diopside wollastonite
wt.hyen = norm.hyen*100;        %hypersthene enstatite
wt.hyfs = norm.hyfs*132;        %hypersthene forsterite








function filename = oxidewt ()

    SiO2 = input('SiO2: ');
    TiO2 = input('TiO2: ');
    Al2O3 = input('Al2O3: ');
    Cr2O3 = input('Cr2O3: ');
    FeO = input('FeO: ');
    Fe2O3 = input('Fe2O3: ');
    MnO = input('MnO: ');
    MgO = input('MgO: ');
    CaO = input('CaO: ');
    Na2O = input('Na2O: ');
    P2O5 = input('P2O5: ');
    K2O = input('K20: ');
    S = input('S: ');
    filename = input('filename: ', 's');
    save(filename);
    
    






