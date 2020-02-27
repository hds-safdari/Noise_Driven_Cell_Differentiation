for cl =1:10  % number of individual RUNs
    cl
    clearvars -except cl clust_coeff_T
    clc
    close all
    writelog('w',['  Start run']);
    tv = [];
    biomass_obj = shoving.biomass();
    
    % Switches
    makeVideo = 0;
    makePlot  = 1;
    
    % Parameters and initial values
    Lx = 2.6e-4; Ly = 2.6e-4; Lz = 10*5.5e-07; % m
    Nx = 30;   Ny = 20;   Nz = 1;
    dxg = Lx/Nx;         % m, grid size
    dyg = Ly/Ny;         % m, grid size
    dzg = Lz/Nz;         % m, grid size
    gridX = dxg/2:dxg:Lx-dxg/2;
    gridY = dyg/2:dyg:Ly-dyg/2;
    gridZ = dzg/2:dzg:Lz-dzg/2;
    
    bac_x = Lx/2;
    bac_y = Ly/2;
    bac_v = 0.6e-18;        % m3
    bac_s = 1;              % start with type 1 (X)
    bac_rho = 100;          % kg/m3
    bac_m = bac_v*bac_rho;  % kg
    bac_r = (3*bac_v/4/pi).^(1/3);
    bac_z = bac_r;
    
    gene_ini = 35;            % initial amount of genes
    gene_1 = gene_ini;            % initial amount of genes
    gene_2 = 2;
    n= 2.0;        %Hill parameter
    beta= 0.1;   %strengths of positive and negative regulatory interactions
    gamma= 1/(10*60);
    s= 10;
    up=floor(gene_ini+gene_ini*0.25);
    lower= 2;
    a1    = 1;
    a2    = 1;
    var_div = 0.10;
    var =10;
    ce_num = 1000;
    
    tend = 4800 ;       % s
    dt   = 10*60;       % s
    td   = 20*60;       % s
    mum  = log(2)/td;   % 1/s
    vmax = 1e-18;       % m^3
    K_s   = 0.1;        % mol/m^3
    Y_sx = 0.01;        %[mol/g]
    
    % Create computational domain geometry and PDE solver
    %mphstart
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    model = ModelUtil.create('Model');
    ModelUtil.showProgress(true); % display the progress bar
    model.modelPath('FILE PATH');
    model.label('signalling_1.mph');
    
    %--- define parameters
    model.param.set('Lx', [num2str(Lx) '[m]']);
    model.param.set('Ly', [num2str(Ly) '[m]']);
    model.param.set('Lz', [num2str(Lz) '[m]']);
    model.param.set('D_c', '1e-11[m^2/s]');
    model.param.set('cb_c', '0[mol/m^3]');
    model.param.set('kc1', '0.001[mol/kg/s]');
    model.param.set('kc2', '0.0001[mol/kg/s]');
    model.param.set('kc3', '0.01[1/s]');
    model.param.set('K_c', '0.01[mol/m^3]');
    model.param.set('T', '10[s]');
    
    %-- create model and 3D geometry
    model.modelNode.create('comp1');
    model.geom.create('geom1', 3);
    model.geom('geom1').create('blk1', 'Block');
    model.geom('geom1').feature('blk1').set('size', {'Lx' 'Ly' 'Lz'});
    model.geom('geom1').run;
    writelog('a',['  Geometry created']);
    
    %-- create mesh
    model.mesh.create('mesh1', 'geom1');
    model.mesh('mesh1').autoMeshSize(1);
    model.mesh('mesh1').run;
    writelog('a',['  Mesh created']);
    
    model.variable.create('var1');
    model.variable('var1').model('comp1');
    model.variable('var1').set('r_c_pr', 'kc1*c_x2(x,y,z)');
    model.variable('var1').set('r_c_use', '-kc2*(c_x1(x,y,z)+c_x2(x,y,z))*c_c/(K_c+c_c)');
    model.variable('var1').set('r_c_dec', '-kc3*c_c');
    model.variable('var1').set('r_c', 'r_c_pr+r_c_dec+r_c_use');
    model.view.create('view2', 3);
    model.view.create('view3', 3);
    
    model.physics.create('tds', 'DilutedSpecies', 'geom1');
    model.physics('tds').field('concentration').field('c_c');
    model.physics('tds').field('concentration').component({'c_c'});
    model.physics('tds').create('reac1', 'Reactions', 3);
    model.physics('tds').feature('reac1').selection.set([1]);
    model.physics('tds').label('Transport and reaction of signal');
    model.physics('tds').prop('TransportMechanism').set('Convection', '0');
    model.physics('tds').feature('cdm1').set('D_c_c', {'D_c'; '0'; '0'; '0'; 'D_c'; '0'; '0'; '0'; 'D_c'});
    model.physics('tds').feature('reac1').set('R_c_c', 'r_c');
    
    
    model.frame('material1').sorder(1);
    writelog('a',['  Physics created']);
    
    %-- create study and its solver
    model.study.create('std1');
    model.study('std1').create('time', 'Transient');
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('t1', 'Time');
    model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('t1').feature.remove('fcDef');
    
    model.study('std1').feature('time').set('rtol', '0.000001');
    model.study('std1').feature('time').set('tlist', 'range(0,1,T)');
    model.study('std1').feature('time').set('rtolactive', true);
    
    model.sol('sol1').attach('std1');
    model.sol('sol1').feature('t1').set('tlist', 'range(0,1,T)');
    model.sol('sol1').feature('t1').set('maxorder', '2');
    model.sol('sol1').feature('t1').set('rtol', '0.000001');
    model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
    model.sol('sol1').feature('t1').feature('dDef').set('rhob', '1');
    writelog('a',['  Solver created']);
    
    %-- create plots
    model.result.create('pg2', 'PlotGroup3D');
    model.result('pg2').create('slc1', 'Slice');
    model.result('pg2').label('Concentration signal T');
    model.result('pg2').feature('slc1').set('shift', '-9.5E-6');
    model.result('pg2').feature('slc1').set('quickplane', 'xy');
    model.result('pg2').feature('slc1').set('interactive', 'on');
    model.result('pg2').feature('slc1').set('quickznumber', '1');
    model.result('pg2').feature('slc1').set('expr', 'c_c');
    
    %----
    
    figure
    % open the video file
    if makeVideo
        video_1 = VideoWriter('biofilm3d.avi');
        uncompressedVideo = VideoWriter('biofilm3d.avi', 'Uncompressed AVI');
        myVideo.FrameRate = 1000;  % Default 30
        myVideo.Quality = 75;    % Default 75
        open(video_1);
    end
    
    [X,Y,Z] = meshgrid(gridX,gridY,gridZ);
    
    
    t = 0; niter = 0;
    while t < tend
        nbac = length(bac_x);   % number of cells
        
        for species=1:2
            all_cells = find(bac_s == species);
            c_x = zeros(Nx,Ny,Nz);  % initialize 3D matrix of cell concentration
            for i=all_cells            % check all cells
                ii = floor(bac_x(i)/dxg)+1;     % find index of the cell on the grid
                jj = floor(bac_y(i)/dyg)+1;
                kk = floor(bac_z(i)/dzg)+1;
                
                c_x(ii,jj,kk)=c_x(ii,jj,kk)+bac_m(i);
            end
            c_x = c_x/dxg/dyg/dzg;
            
            % save a 3D matrix
            nameFile = ['c_x' num2str(species)];
            name = [nameFile '.' num2str(niter,'%03d') '.txt'];
            file = fopen(name,'w');
            
            fprintf(file,'%% Grid\n');
            
            for i = 1:length(gridX), fprintf(file,'%13.5e',gridX(i)); end; fprintf(file,'\n');
            for i = 1:length(gridY), fprintf(file,'%13.5e',gridY(i)); end; fprintf(file,'\n');
            for i = 1:length(gridZ), fprintf(file,'%13.5e',gridZ(i)); end; fprintf(file,'\n');
            
            fprintf(file,['%% Data (' nameFile ')\n']);
            fclose(file);
            
            for ii = 1:1
                aslice(:,:)=c_x(:,:,ii);  aslice_new = aslice';
                save(name,'aslice_new','-ascii','-append');
                clear aslice_new;
                file = fopen(name,'a');
                fprintf(file,'\n');
                fclose(file);
            end;
            
        end
        writelog('a',['  Saved cell concentration ' name]);
        
        % Solve concentrations
        model.func.create('int1', 'Interpolation');
        model.func.create('int2', 'Interpolation');
        model.func('int1').model('comp1');
        model.func('int1').set('argunit', 'm');
        name = ['c_x1.' num2str(niter,'%03d') '.txt'];
        model.func('int1').set('filename', name);
        model.func('int1').set('source', 'file');
        model.func('int1').set('struct', 'grid');
        model.func('int1').set('funcs', {'c_x1' '1'});
        model.func('int1').set('fununit', 'kg/m^3');
        model.func('int2').model('comp1');
        model.func('int2').set('argunit', 'm');
        name = ['c_x2.' num2str(niter,'%03d') '.txt'];
        model.func('int2').set('filename', name);
        model.func('int2').set('source', 'file');
        model.func('int2').set('struct', 'grid');
        model.func('int2').set('funcs', {'c_x2' '1'});
        model.func('int2').set('fununit', 'kg/m^3');
        writelog('a',['  Functions created']);
        
        model.sol('sol1').runAll;
        writelog('a',['  Concentrations solved']);
        bac_z = bac_z';
        coord = [bac_x(1:nbac)';bac_y(1:nbac)';bac_z(1:nbac)'];
        c_c = mphinterp(model,'c_c','coord',coord);
        
        nameCOMSOL = ['signalling3D.' num2str(niter,'%03d') '.mph']
        mphsave(model,nameCOMSOL);
        writelog('a',['  COMSOL model saved ' nameCOMSOL]);
        
        %--- remove functions because these must be added again
        model.func.remove('int1'); % removing function1
        model.func.remove('int2'); % removing function2
        
        %--- prepare solver for the next time interval
        model.study('std1').feature('time').set('solnum', 'last');
        model.study('std1').feature('time').set('useinitsol', 'on');
        model.study('std1').feature('time').set('tlist', ['range(' num2str(niter+1) '*T,1,' num2str(niter+2) '*T)']);
        model.study('std1').feature('time').set('initstudy', 'std1');
        model.study('std1').feature('time').set('initmethod', 'sol');
        
        model.sol('sol1').attach('std1');
        model.sol('sol1').feature('v1').set('initsol', 'sol1');
        model.sol('sol1').feature('v1').set('solnum', 'last');
        model.sol('sol1').feature('v1').set('initmethod', 'sol');
        model.sol('sol1').feature('t1').set('tlist', ['range(' num2str(niter+1) '*T,1,' num2str(niter+2) '*T)']);
        model.sol('sol1').feature('t1').set('maxorder', '2');
        model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
        
        % Growth step
        bac_v = (1e-17)*ones(1,nbac);
        % Division step
        nbac = length(bac_x);   % number of cells
        for i=1:nbac            % check all cells
            if bac_v(i) > vmax  % cell is greater than max. volume
                c_mean1 = mean(mean(abs(c_c)));
                if ( nbac >= 2)
                    [gene_1(i),gene_2(i),step_num] = gene_expression(gene_1(i),gene_2(i),up,lower,n,beta,gamma,s,a1,a2);
                    gene_1(i)=BondCkeck(gene_1(i),lower,up);
                    gene_2(i)=BondCkeck(gene_2(i),lower,up);
                end
                a1 = 1;
                a2 = 1;
                m_old = bac_m(i);  % store mass before division
                v_old = bac_v(i);  % store mass before division
                x_old = bac_x(i);   y_old = bac_y(i);    z_old = bac_z(i);
                gene_1old=gene_1(i);
                gene_2old=gene_2(i);
                phi = rand*pi;  % random angle for division
                the = rand*pi/2;  % random angle for division
                % % % % % % %% % % % % % %% % % % % % %% % % % % % %
                % % % % % % %% % % % % % %% % % % % % %% % % % % % %
                rand_p1 =  normrnd(0.5,var_div);
                if (rand_p1 >= 1)
                    rand_p1 = 1;
                end
                r3 = abs(rand(1,ce_num));
                r4 = abs(rand(1,ce_num));
                % % % % % % %% % % % % % %% % % % % % %% % % % % % %
                % % %                         r3=ones(1,ce_num);
                % % %                         r4=ones(1,ce_num);
                r3 = r3./sum(r3);
                r4 = r4./sum(r4);
                r3 = gene_1old*r3;
                r4 = gene_2old*r4;
                r1 =  floor(rand_p1*ce_num);
                r2 =  floor(rand_p1*ce_num);
                gene_1(i)=floor(sum(r3(1:r1)));
                gene_2(i)=floor(sum(r4(1:r2)));
                % % % % % % %% % % % % % %% % % % % % %% % % % % % %
                % % % % % % %% % % % % % %% % % % % % %% % % % % % %
                % "old" cell
                bac_m(i) =m_old;
                bac_v(i) = v_old;
                bac_r(i) = (3*bac_v(i)/4/pi)^(1/3);
                bac_x(i) = x_old + bac_r(i)*cos(phi)*sin(the);
                bac_y(i) = y_old + bac_r(i)*sin(phi)*sin(the);
                if  gene_2(i)>=gene_1(i)
                    bac_s (i) = 2;
                elseif gene_2(i)<gene_1(i)
                    bac_s (i) = 1;
                end
                % "new" cell
                j = nbac+1;     % index new cell at the end
                bac_m(j) =m_old;
                bac_v(j) = v_old;
                bac_r(j) = (3*bac_v(j)/4/pi)^(1/3);
                bac_x(j) = x_old - bac_r(j)*cos(phi)*sin(the);
                bac_y(j) = y_old - bac_r(j)*sin(phi)*sin(the);
                gene_1(j)= gene_1old-gene_1(i);
                gene_2(j)= gene_2old-gene_2(i);
                gene_1(i)=BondCkeck(gene_1(i),lower,up);
                gene_2(i)=BondCkeck(gene_2(i),lower,up);
                gene_1(j)=BondCkeck(gene_1(j),lower,up);
                gene_2(j)=BondCkeck(gene_2(j),lower,up);
                if  gene_2(j)>=gene_1(j)
                    bac_s (j) = 2;
                elseif gene_2(j)<gene_1(j)
                    bac_s (j) = 1;
                end
                nbac = j;   % number of cells
            end
        end
        
        % Shoving step (avoiding cell overlap after division)
        tic
        bac_z1 = 5.5e-07*ones(1,nbac);
        Results = biomass_obj.pushing3D(nbac, bac_x, bac_y, bac_r);
        tv=[tv toc];
        bac_x = Results.bac_x;
        bac_y = Results.bac_y;
        bac_z = 5.5e-07*ones(1,nbac);
        writelog('a',['  Cell shoving done']);
        
        % Time update
        t = t + dt
        niter = niter+1;
        check_1(niter) = min(min(c_c));
        
        % Plotting
        if (t== tend)
            if makePlot
                clf
                for i=1:nbac
                    if bac_s(i)==1
                        col=[	57, 73, 171]./255;
                    elseif bac_s(i)==2
                        col=[255, 183, 77]./255;
                    end
                    hh=rectangle('Position',[bac_r(i)+ bac_x(i),bac_r(i)+ bac_y(i), 1.5*bac_r(i), 1.5*bac_r(i)], 'Curvature',[1,1], 'FaceColor',col, 'EdgeColor',col);
                    axis square;
                    hold on
                end
                % % %             box on
                axis([ 0.75*min(bac_x) 1.1*max(bac_x)  0.75*min(bac_y) 1.1*max(bac_y)]) ;
                axis equal;
                axis off
            end
            axis([0.96e-4 1.7e-4 0.96e-4 1.7e-4])
            name = ['OUTPUT PATH/var_' num2str(var,'%d') '_org_' num2str(cl,'%d') '_no_org.eps'];
            print( [name] ,'-depsc')
        end 
    end
    % % % if makeVideo, close(video_1); end
    adj_mat  = (zeros(nbac,nbac));
    adj_mat2 = (zeros(nbac,nbac));
    for i = 1:nbac
        for j = i+1:nbac
            dx2 = (bac_x(i)-bac_x(j))*(bac_x(i)-bac_x(j));
            dy2 = (bac_y(i)-bac_y(j))*(bac_y(i)-bac_y(j));
            dr  = sqrt(dx2+dy2);
            if(dr < 5 * bac_r(i)) && (bac_s(i) == 2)
                if (bac_s(j) == 2)
                    adj_mat2 (i,j) = 1;
                    adj_mat2 (j,i) = 1;
                end
            elseif(dr <= 5 * bac_r(i)) && (bac_s(i) == 1)
                if (bac_s(j) == 1)
                    adj_mat (i,j) = 1;
                    adj_mat (j,i) = 1;
                end
            end
        end
    end
    [a1 a2, a3] = clust_coeff(adj_mat);
    [b1 b2 b3] = clust_coeff(adj_mat2);
    clust_coeff_T (cl,:) =  [a1  a2  b1  b2]; 
end
name = ['OUTPUT PATH/clust_coeff_T_0' num2str(var,'%d') '_no_org.txt'];
save(name,'clust_coeff_T','-ascii')
