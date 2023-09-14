function solvesegregationmodelatmaxNF

% Set initial densities
V_ic = 0.09;
C1_ic = 0.0017;
J_ic=0.35;
H_ic=0.19;
L_ic=0.05;
C2_ic=1.25;

% Solve model for 50 years and process output
lengthtime = 50*365;
[T,Y,logC1,logC2,LossRates] = solveODEmodel([0,lengthtime],[V_ic,J_ic,H_ic,L_ic,C1_ic,C2_ic]);
matrixofoutput = [T,Y,logC1,logC2,LossRates];
matrixofoutput=matrixofoutput(:,[1:5,7,9,15:end]);
tableofoutput=array2table(matrixofoutput);
tableofoutput.Properties.VariableNames(1:11)={'Time','Vultures','Jackals','Hyenas','Lions','Carrion 2','Log10Carrion2','JackalRemovalRate2','HyenaRemovalRate2','LionRemovalRate2','DecayRate2'};
writetable(tableofoutput,'nf_1.csv');


    function [t,y,logcarcass1,logcarcass2,lossrates] = solveODEmodel(interval,ics)

        % Vulture parameter values
        e_A=157.5; mu=1/(13*365); h_A=0.85; g_a=9.35*10^(-4);

        % Carrion parameter values
        p=2.26;
        delta=0.07;
        n_f=1;

        % Mammal parameter values
        KJ = 0.9*0.35;
        KH = 0.48*0.19;
        KL = 0.725*0.05;

        h_J=0.32; e_J=4;
        h_H=0.21; e_H=27.28;
        h_L=0.08; e_L=4;

        g_J=2.33*10^(-3);
        g_H=8.4*10^(-5);
        g_L=1.08*10^(-4);

        rJ = g_J/h_J - 1/(3*365);
        rH = g_H/h_H - 1/(15.5*365);
        rL = g_L/h_L - 1/(11*365);

        opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

        % increase grid resolution for first two years
        l1=linspace(0,2*265,201);
        l2=linspace(2*265,interval(end),100);
        l3=[l1(1:end-1),l2];
        [t,y] = ode45(@mysystem,l3,ics,opts);

        logcarcass1 = log10(y(:,5));
        JLossShared=y(:,2).*(e_J.*y(:,5))./(1+h_J.*e_J.*y(:,5)+ h_J.*e_J.*y(:,6));
        HLossShared=y(:,3).*(e_H.*y(:,5))./(1+h_H.*e_H.*y(:,5)+h_H.*e_H.*y(:,6));
        LLossShared=y(:,4).*(e_L.*y(:,5))./(1+h_L.*e_L.*y(:,5)+h_L.*e_L.*y(:,6));
        VLossShared=y(:,1).*(e_A*y(:,5))./(1+h_A.*e_A.*y(:,5));
        ILossShared=delta*y(:,5);


        logcarcass2 = log10(y(:,6));
        JLossNight=y(:,2).*(e_J.*y(:,6))./(1+h_J.*e_J.*y(:,5)+ h_J.*e_J.*y(:,6));
        HLossNight=y(:,3).*(e_H.*y(:,6))./(1+h_H.*e_H.*y(:,5)+h_H.*e_H.*y(:,6));
        LLossNight=y(:,4).*(e_L.*y(:,6))./(1+h_L.*e_L.*y(:,5)+h_L.*e_L.*y(:,6));
        ILossNight=delta*y(:,6);

        lossrates=[JLossShared,HLossShared,LLossShared,VLossShared,ILossShared,JLossNight,HLossNight,LLossNight,ILossNight];


        function dy = mysystem(t,y)

            A=y(1);
            J=y(2);
            H=y(3);
            L=y(4);
            C=y(5);
            N=y(6);


            A_dyn=-mu*A + g_a*A*(e_A*C)/(1+e_A*C*(h_A));
            J_dyn=rJ*J*(1-J/KJ)+g_J*J*(e_J*(C+N))/(1+h_J*e_J*(C+N));
            H_dyn=rH*H*(1-H/KH)+g_H*H*(e_H*(C+N))/(1+h_H*e_H*(C+N));
            L_dyn=rL*L*(1-L/KL)+g_L*L*(e_L*(C+N))/(1+h_L*e_L*(C+N));
            C_dyn=p*(1-n_f) - delta*C -A*(e_A*C)/(1+e_A*C*(h_A))-J*(e_J*C)/(1+h_J*e_J*(C+N)) -H*(e_H*C)/(1+h_H*e_H*(C+N)) -L*(e_L*C)/(1+h_L*e_L*(C+N));
            N_dyn=p*n_f - delta*N -J*(e_J*N)/(1+h_J*e_J*(C+N)) -H*(e_H*N)/(1+h_H*e_H*(C+N)) -L*(e_L*N)/(1+h_L*e_L*(C+N));

            dy = [A_dyn;J_dyn;H_dyn;L_dyn;C_dyn;N_dyn];

        end


    end

end
