function [mol1_n1,mol2_n2,    step_num ] = gene_expression( mol1_n, mol2_n,mol_n_max,mol_n_min,n,beta,gamma,s,a1,a2)
% Gillespie Algorithm 
t        = 0; 
j        =  0;
t_end =  1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (t < t_end)
    j   = j+1;
	z1 = rand();                       % uniform random numbers
	z2 = rand();                                  % reaction rates for every possible reaction
	w(1) = a2*beta*(s^n/(s^n+mol1_n^n));
	w(2) = a1*beta*(s^n/(s^n+mol2_n^n));
	w(3) = gamma* (mol1_n^(1));    
	w(4) = gamma* (mol2_n^(1));    
    W_r = sum(w);                                      % the total reaction rate
    a     = cumsum(w)/W_r;
%     a     = w./W_r;
    if( z1 <= a(1) )
        mol2_n = mol2_n+1;
        if (mol2_n >= mol_n_max) 
            mol2_n = mol_n_max;
        end           
    elseif(a(1) < z1) && ( z1 <= a(2))   
        mol1_n = mol1_n+1;
        if (mol1_n >= mol_n_max) 
            mol1_n = mol_n_max;
        end                  
    elseif(a(2) <  z1 ) && ( z1 <= a(3))
        mol1_n = mol1_n-1;
        if (mol1_n <= mol_n_min) 
            mol1_n = mol_n_min;
        end
    elseif(a(3) <  z1) && ( z1 <= a(4)) 
        mol2_n = mol2_n-1;
        if (mol2_n <= mol_n_min) 
            mol2_n = mol_n_min;
        end            
    end
	dt  = log(1/z2)/W_r;                                           % time step
	t    = t+dt;
end
mol1_n1 = mol1_n;
mol2_n2 = mol2_n;
step_num = j;
end

