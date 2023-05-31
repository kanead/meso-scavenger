function solvescavengingmodel(scenario)
% Set scenario = {1,2,3} for main figure, supply-rate sensitivity and carrying-capacity sensitivity respectively

if scenario==1
    KJ_input=0.35;  KH_input=0.19; KL_input=0.05; p_input=2.26;
    J_ic = KJ_input; H_ic = KH_input; L_ic = KL_input;
    csvfilename='maintextscenario.csv';
elseif scenario==2
    KJ_input=0.35;  KH_input=0.19; KL_input=0.05; p_input=2.26/2;
    J_ic = KJ_input; H_ic = KH_input; L_ic = KL_input;
    csvfilename='halvedCin.csv';
else
    KJ_input=0.9*0.35;KH_input=0.48*0.19;KL_input=0.725*0.05; p_input=2.26;
    J_ic = 0.35; H_ic = 0.19; L_ic = 0.05;
    csvfilename='reducedKs.csv';
end


% Set initial vulture and carrion densities
V_ic = 0.09;
C_ic = 0.0017;

% Solve model for 50 years and process output
lengthtime = 50*365;
[T,Y,logC,LossRates] = solveODEmodel([0,lengthtime],[V_ic,J_ic,H_ic,L_ic,C_ic],KJ_input,KH_input,KL_input,p_input);
matrixofoutput = [T,Y,logC,LossRates];
tableofoutput=array2table(matrixofoutput);
tableofoutput.Properties.VariableNames(1:12)={'Time','Vultures','Jackals','Hyenas','Lions','Carrion','Log10Carrion','JackalRemovalRate','HyenaRemovalRate','LionRemovalRate','VultureRemovalRate','DecayRate'};
writetable(tableofoutput,csvfilename);

    function [t,y,logcarcass,lossrates] = solveODEmodel(interval,ics,KJ,KH,KL,p)

        % Vulture parameter values
        e_A=157.5; mu=1/(13*365); h_A=0.85; g_a=9.35*10^(-4);

        % Carrion parameter values
        delta=0.07;

        % Mammal parameter values
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

        logcarcass = log10(y(:,5));
        JLoss=y(:,2).*(e_J.*y(:,5))./(1+h_J.*e_J.*y(:,5));
        HLoss=y(:,3).*(e_H.*y(:,5))./(1+h_H.*e_H.*y(:,5));
        LLoss=y(:,4).*(e_L.*y(:,5))./(1+h_L.*e_L.*y(:,5));
        VLoss=y(:,1).*(e_A*y(:,5))./(1+h_A*e_A*y(:,5));
        ILoss=delta*y(:,5);
        lossrates=[JLoss,HLoss,LLoss,VLoss,ILoss];

        function dy = mysystem(t,y)

            A=y(1);
            J=y(2);
            H=y(3);
            L=y(4);
            C=y(5);

            A_dyn=-mu*A + g_a*A*(e_A*C)/(1+e_A*C*h_A);
            J_dyn=rJ*J*(1-J/KJ)+g_J*J*(e_J*C)/(1+h_J*e_J*C);
            H_dyn=rH*H*(1-H/KH)+g_H*H*(e_H*C)/(1+h_H*e_H*C);
            L_dyn=rL*L*(1-L/KL)+g_L*L*(e_L*C)/(1+h_L*e_L*C);
            C_dyn=p - delta*C -A*(e_A*C)/(1+e_A*C*h_A) -J*(e_J*C)/(1+h_J*e_J*C) -H*(e_H*C)/(1+h_H*e_H*C) -L*(e_L*C)/(1+h_L*e_L*C);

            dy = [A_dyn;J_dyn;H_dyn;L_dyn;C_dyn];

        end


    end

end
