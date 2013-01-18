%% Run Lanczos and CA-Lanczos on a number of matrices with a number of
%% orthogonalization strategies, and print the results to corresponding
%% output text files.
%%
%% A few options need to be set:
%% - output_path: Where to store the result files. Folder must exist.
%% - matrices_path: Where the matrices (.mat-files) are located.
%% - orth_strategies: Orthogonalization strategies to run.
%% - basis: Which basis should be used.

output_path = '../results';
matrices_path = '../matrices';

matrices = {'1138_bus'};%,'Chem97ZtZ','finan512','fv2','fv3','mhdb416','msc04515','plat1919'};

orth_strategies = {'local', 'full', 'periodic', 'selective'};

basis = 'newton';
    
% Number of Lanczos steps
m=120;

for matrix = matrices
    
    % Load matrix
    load([matrices_path,'/',cell2mat(matrix)]);
    Problem
    A=Problem.A;
    n=size(A,1);

    % Starting vector
    r0=rand(n,1);
 
    for orth_strat = orth_strategies

        orth = cell2mat(orth_strat);
        
        outfilename = ['calanczos','_',cell2mat(matrix),'_','newton','_',orth,'_','120','.out'];

        outfile = fopen([output_path,'/',outfilename],'w');
        
        % Standard Lanczos
        fprintf(outfile, '#Standard Lanczos\n');
        [T_st,V_st,r_st,o_st]=lanczos(A,r0,m,orth);
        for i = 1:length(r_st)
            fprintf(outfile, '%d,%.4e,%.4e\n', i, r_st(i,2), o_st(i));
        end
        fprintf(outfile,'\n');
            
        % CA-Lanczos(4)
        fprintf(outfile, '#CA-Lanczos(4)\n');
        [T_ca_4,V_ca_4,r_ca_4,o_ca_4]=ca_lanczos(A,r0,4,m/4,basis,orth);
        for i = 1:length(r_ca_4)
            fprintf(outfile, '%d,%.4e,%.4e\n', 4*i, r_ca_4(i,2), o_ca_4(i));
        end
        fprintf(outfile,'\n');
        
        % CA-Lanczos(6)
        fprintf(outfile, '#CA-Lanczos(6)\n');
        [T_ca_6,V_ca_6,r_ca_6,o_ca_6]=ca_lanczos(A,r0,6,m/6,basis,orth);
        for i = 1:length(r_ca_6)
            fprintf(outfile, '%d,%.4e,%.4e\n', 6*i, r_ca_6(i,2), o_ca_6(i));
        end
        fprintf(outfile,'\n');
        
        % CA-Lanczos(8)
        fprintf(outfile, '#CA-Lanczos(8)\n');
        [T_ca_8,V_ca_8,r_ca_8,o_ca_8]=ca_lanczos(A,r0,8,m/8,basis,orth);
        for i = 1:length(r_ca_8)
            fprintf(outfile, '%d,%.4e,%.4e\n', 8*i, r_ca_8(i,2), o_ca_8(i));
        end
        fprintf(outfile,'\n');
        
        % CA-Lanczos(10)
        fprintf(outfile, '#CA-Lanczos(10)\n');
        [T_ca_10,V_ca_10,r_ca_10,o_ca_10]=ca_lanczos(A,r0,10,m/10,basis,orth);
        for i = 1:length(r_ca_10)
            fprintf(outfile, '%d,%.4e,%.4e\n', 10*i, r_ca_10(i,2), o_ca_10(i));
        end
        fprintf(outfile,'\n');
        
        
        
        % Plot the convergence of the largest eigenpair
        figure;
        semilogy(1:size(T_st,1),r_st(:,2),'*:','Color',mediumgreen,'MarkerSize',8);  hold on;
        semilogy((1:size(T_ca_4,1)/4)*4,r_ca_4(:,2),'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_6,1)/6)*6,r_ca_6(:,2),'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_8,1)/8)*8,r_ca_8(:,2),'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_10,1)/10)*10,r_ca_10(:,2),'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
        %semilogy((1:size(T_ca_20,1)/20)*20,r_ca_20(:,2),'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
        title('Convergence of largest eigenpair');
        ylabel('||Ay - \lambda y||');
        xlabel('Iterations');
        legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
        if (min(r_st(:,2)) < 1.0e-14) ...
                || (min(r_ca_4(:,2)) < 1.0e-14) || (min(r_ca_6(:,2)) < 1.0e-14) ...
                || (min(r_ca_8(:,2)) < 1.0e-14) || (min(r_ca_10(:,2)) < 1.0e-14) ...%|| (min(r_ca_20(:,2)) < 1.0e-14)
                
            ylim([1.0e-16 1.0e02]);
        end
        
        % Plot the orthogonality
        figure;
        semilogy(1:size(T_st,1),o_st,'*:','Color',mediumgreen,'MarkerSize',8); hold on;
        semilogy((1:size(T_ca_4,1)/4)*4,o_ca_4,'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_6,1)/6)*6,o_ca_6,'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_8,1)/8)*8,o_ca_8,'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
        semilogy((1:size(T_ca_10,1)/10)*10,o_ca_10,'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
        %semilogy((1:size(T_ca_20,1)/20)*20,o_ca_20,'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
        title('Orthogonality of basis vectors');
        ylabel('|| I - Q_m^*Q ||_F');
        xlabel('Iterations');
        legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
        
        if (min(o_st) < 1.0e-14) ...
                || (min(o_ca_4) < 1.0e-14) || (min(o_ca_6) < 1.0e-14) ...
                || (min(o_ca_8) < 1.0e-14) || (min(o_ca_10) < 1.0e-14) ...%|| (min(o_ca_20) < 1.0e-14)
                
            ylim([1.0e-16 1.0e02]);
        end
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
