function temporalsegregation_eqpts


% Parameter-sweep details
nf_start = 0;
nf_end = 1;
valuestaken = 500;
parameterlist = linspace(nf_start,nf_end,valuestaken);
filename='increase_n_f.csv';
ics =[9.9687,0.3509,0.1935,0.0500,0.0009,0];
errortol = 1e-5;   % termination threshold for Newton method

% Vulture parameter values
e_A=157.5; mu=1/(13*365); h_A=0.85; g_a=9.35*10^(-4);

% Carrion parameter values
delta=0.07; p=2.26;

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

% (Adjusted) carrying-capacity parameters
KJ=0.9*0.35;KH=0.48*0.19;KL=0.725*0.05;


    function [dy,Ja] = mysystem(t,y,parameterspec)
        
        n_f = parameterspec;
        
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
        
        Ja=[-mu+g_a*e_A*C/(1+e_A*C*(h_A)),0,0,0,g_a*A*e_A/(1+e_A*C*(h_A))-g_a*A*e_A^2*C/(1+e_A*C*(h_A))^2*(h_A),0;
            0,rJ*(1-J/KJ)-rJ*J/KJ+g_J*e_J*(C+N)/(1+h_J*e_J*(C+N)),0,0,g_J*J*e_J/(1+h_J*e_J*(C+N))-g_J*J*e_J^2*(C+N)/(1+h_J*e_J*(C+N))^2*h_J,g_J*J*e_J/(1+h_J*e_J*(C+N))-g_J*J*e_J^2*(C+N)/(1+h_J*e_J*(C+N))^2*h_J;
            0,0,rH*(1-H/KH)-rH*H/KH+g_H*e_H*(C+N)/(1+h_H*e_H*(C+N)),0,g_H*H*e_H/(1+h_H*e_H*(C+N))-g_H*H*e_H^2*(C+N)/(1+h_H*e_H*(C+N))^2*h_H,g_H*H*e_H/(1+h_H*e_H*(C+N))-g_H*H*e_H^2*(C+N)/(1+h_H*e_H*(C+N))^2*h_H;
            0,0,0,rL*(1-L/KL)-rL*L/KL+g_L*e_L*(C+N)/(1+h_L*e_L*(C+N)),g_L*L*e_L/(1+h_L*e_L*(C+N))-g_L*L*e_L^2*(C+N)/(1+h_L*e_L*(C+N))^2*h_L,g_L*L*e_L/(1+h_L*e_L*(C+N))-g_L*L*e_L^2*(C+N)/(1+h_L*e_L*(C+N))^2*h_L;
            -e_A*C/(1+e_A*C*(h_A)),-e_J*C/(1+h_J*e_J*(C+N)),-e_H*C/(1+h_H*e_H*(C+N)),-e_L*C/(1+h_L*e_L*(C+N)),-delta-A*e_A/(1+e_A*C*(h_A))+A*e_A^2*C/(1+e_A*C*(h_A))^2*(h_A)-J*e_J/(1+h_J*e_J*(C+N))+J*e_J^2*C/(1+h_J*e_J*(C+N))^2*h_J-H*e_H/(1+h_H*e_H*(C+N))+H*e_H^2*C/(1+h_H*e_H*(C+N))^2*h_H-L*e_L/(1+h_L*e_L*(C+N))+L*e_L^2*C/(1+h_L*e_L*(C+N))^2*h_L,J*e_J^2*C/(1+h_J*e_J*(C+N))^2*h_J+H*e_H^2*C/(1+h_H*e_H*(C+N))^2*h_H+L*e_L^2*C/(1+h_L*e_L*(C+N))^2*h_L;
            0,-e_J*N/(1+h_J*e_J*(C+N)),-e_H*N/(1+h_H*e_H*(C+N)),-e_L*N/(1+h_L*e_L*(C+N)),J*e_J^2*N/(1+h_J*e_J*(C+N))^2*h_J+H*e_H^2*N/(1+h_H*e_H*(C+N))^2*h_H+L*e_L^2*N/(1+h_L*e_L*(C+N))^2*h_L,-delta-J*e_J/(1+h_J*e_J*(C+N))+J*e_J^2*N/(1+h_J*e_J*(C+N))^2*h_J-H*e_H/(1+h_H*e_H*(C+N))+H*e_H^2*N/(1+h_H*e_H*(C+N))^2*h_H-L*e_L/(1+h_L*e_L*(C+N))+L*e_L^2*N/(1+h_L*e_L*(C+N))^2*h_L];
        
    end

matrixofoutput = producebifcurves(parameterlist,ics);
tableofoutput=array2table(matrixofoutput);
tableofoutput.Properties.VariableNames(1:13)={'Fraction','Vultures','Jackals','Hyenas','Lions','Carrion 1','Carrion 2','Eig 1','Eig 2','Eig 3','Eig 4','Eig 5','Eig 6'};
writetable(tableofoutput,filename);

    function ptcollection = producebifcurves(parameterlist,ics)
        
        j=1;
        ptcollection = [];
        
        while j<length(parameterlist)+1
            
            parametervalue = parameterlist(j);            
            solns=[];
            
            % apply multi-variable Newton method to find nearby equilibrium point
            i=1;
            currentstate=ics';
            
            while true
                
                [DEatcurrentstate,JacobianAtcurrentstate]=mysystem(0,currentstate,parametervalue);
                
                inc_vector = -JacobianAtcurrentstate\DEatcurrentstate;
                
                newstate = currentstate + inc_vector;
                currentstate = newstate;
                
                if max(abs(inc_vector)) < errortol
                    break;
                end
                
                i=i+1;
                
            end
            
            equil_candidate=newstate;     % set state to Newton method's final approximation
            solns(1,:) = equil_candidate';
            
            j=j+1;
            
            if (j>=3 && j<length(parameterlist)+1)
                previousfinalposition = finalposition;
            end
            
            finalposition = equil_candidate';
            
            % if possible, use slope to improve prediction of equilibrium position at next n_f value   
            if j<3
                ics = finalposition;
            elseif (j>=3 && j<length(parameterlist)+1)
                ics= previousfinalposition + (finalposition-previousfinalposition)*(parameterlist(j)-parameterlist(j-2))/(parameterlist(j-1)-parameterlist(j-2));
            end
            
            % calculate and collect eigenvalues of Jacobian matrix
            e=eigs(JacobianAtcurrentstate);
            ptcollection=[ptcollection;parametervalue*ones(size(solns,1),1),solns,e'];
            
        end
        
    end


end