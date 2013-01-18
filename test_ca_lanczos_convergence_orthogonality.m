clear;

m=120;        % Number of Lanczos steps
neigs=10;     % Number of eigs we want to print (largest)
orth='full'; % Orthogonality method to use

output_path = '../results/conv_orth_full';

% Comments on matrices:
% - fv3 needs more iterations than 120 to converge
% - msc04515 barely converges in 120 iterations
matrices_path = '../matrices';
matrices = {'Chem97ZtZ','fv3','mesh2e1','ex10','msc04515'};
for matrix = matrices

    load([matrices_path,'/',cell2mat(matrix)]);
    A = Problem.A;
    n=size(A,1);
    norm_A = norm(A,'inf');
    A = A/norm_A;
    
    T = cell(1,7);
    r = cell(1,7);
    o = cell(1,7);
        
    % Starting vector
    r0=rand(n,1);
    
    outfilename = ['ca_lanczos_conv_orth_',cell2mat(matrix),'_','full','_',num2str(m),'_',num2str(neigs),'.out'];
    outfile = fopen([output_path,'/',outfilename],'w');

    fprintf(outfile,'%s\n',cell2mat(matrix));
    fprintf(outfile,'%d\n',m);
    fprintf(outfile,'s,1,2,4,6,8,10\n');

    % Standard Lanczos
    [T_,Q,ritz_norms,orth_error] = lanczos(A,r0,m,orth);
    fprintf(outfile,'#Lanczos\n');
    for i = 1:m
        fprintf(outfile,'%.6e,%.6e\n',ritz_norms(i,1),orth_error(i));
    end
    eigs = sort(eig(T_),'descend');
    fprintf(outfile,'eigs:');
    for i = 1:neigs
        fprintf(outfile,'%.16e,',eigs(i)*norm_A);
    end
    fprintf(outfile,'\b\n');
    T{1} = T_;
    r{1} = ritz_norms;
    o{1} = orth_error;
    
    % CA-Lanczos
    j = 1;
    for s = [1,2,4,6,8,10]
        [T_,Q,ritz_norms,orth_error] = ca_lanczos(A,r0,s,m,'monomial',orth);
        fprintf(outfile,'#CA-Lanczos(%d)\n',s);
        for i = 1:floor(m/s)
            fprintf(outfile,'%.6e,%.6e\n',ritz_norms(i,1),orth_error(i));
        end
        eigs = sort(eig(T_),'descend');
        fprintf(outfile,'eigs:');
        for i = 1:neigs
            fprintf(outfile,'%.16e,',eigs(i)*norm_A);
        end
        fprintf(outfile,'\b\n');
        T{j+1} = T_;
        r{j+1} = ritz_norms;
        o{j+1} = orth_error;
        j = j+1;
    end 
    
    %Define some colors for plotting
    blue         = [   0   0   1 ];
    mediumgreen  = [   0 .60   0 ];
    red          = [   1   0   0 ];
    orchid       = [ .85 .44 .84 ];
    pink         = [   1 .41 .71 ];
    salmon       = [ .98 .50 .45 ];
    orange       = [   1 .55   0 ];
    royalblue    = [.25  .41 .88 ];
    
    % Plot the convergence of the largest eigenpair
    figure; 
    subplot(2,1,1);
    semilogy(1:size(T{1},1),r{1}(:,1),'*:','Color',mediumgreen,'MarkerSize',8); hold on;
    semilogy((1:size(T{4},1)/4)*4,r{4}(:,1),'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{5},1)/6)*6,r{5}(:,1),'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{6},1)/8)*8,r{6}(:,1),'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{7},1)/10)*10,r{7}(:,1),'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
    %semilogy((1:size(T_ca_20,1)/20)*20,r_ca_20(:,2),'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
    hold off;
    title('Convergence of largest eigenpair');
    ylabel('||Ay - \lambda y||');
    xlabel('Iterations');
    legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
    if (min(r{1}(:,2)) < 1.0e-14) || (min(r{4}(:,2)) < 1.0e-14) || (min(r{5}(:,2)) < 1.0e-14) ||...
            (min(r{6}(:,2)) < 1.0e-14) || (min(r{7}(:,2)) < 1.0e-14)
        
        ylim([1.0e-16 1.0e02]);
    end
        
    % Plot the orthogonality
    subplot(2,1,2);
    semilogy(1:size(T{1},1),o{1},'*:','Color',mediumgreen,'MarkerSize',8); hold on;
    semilogy((1:size(T{4},1)/4)*4,o{4},'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{5},1)/6)*6,o{5},'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{6},1)/8)*8,o{6},'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
    semilogy((1:size(T{7},1)/10)*10,o{7},'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
    %semilogy((1:size(T_ca_20,1)/20)*20,o_ca_20,'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
    hold off;
    title('Orthogonality error');
    ylabel('|| I - Q_m^*Q ||_F');
    xlabel('Iterations');
    legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
    
    if (min(o{1}) < 1.0e-14) || (min(o{4}) < 1.0e-14) || (min(o{5}) < 1.0e-14) || ...
            (min(o{6}) < 1.0e-14) || (min(o{7}) < 1.0e-14)
            
        ylim([1.0e-16 1.0e02]);
    end
    
end