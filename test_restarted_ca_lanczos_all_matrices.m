clear;

output_path = '../results';

matrices_path = '../matrices';
matrices = {'494_bus','662_bus','685_bus','1138_bus','bcsstk01','bcsstk02','bcsstk03','bcsstk04','bcsstk05','bcsstk06','bcsstk07','bcsstk08',...
    'bcsstk09','bcsstk10','bcsstk11','bcsstk12','bcsstk13','bcsstk14','bcsstk19','bcsstk20','bcsstk21','bcsstk22','bcsstk23','bcsstk26',...
    'bcsstk27','bcsstk34','bcsstm02','bcsstm05','bcsstm06','bcsstm07','bcsstm08','bcsstm09','bcsstm11','bcsstm12','bcsstm19','bcsstm20',...
    'bcsstm21','bcsstm22','bcsstm23','bcsstm24','bcsstm25','bcsstm26','bcsstm39','bibd_81_2','bloweybq','Chem97ZtZ','ex3','ex5','ex9','ex10',...
    'ex10hs','ex13','ex15','ex33','finan512','fv1','fv2','fv3','gr_30_30','Journals','LF10','LF10000','LFAT5','LFAT5000','lund_a','lund_b',...
    'mesh1e1','mesh1em1','mesh1em6','mesh2e1','mesh2em5','mesh3e1','mesh3em5','mhd1280b','mhd3200b','mhd4800b','mhdb416','msc00726','msc01050',...
    'msc01440','msc04515','nasa1824','nasa2146','nos1','nos2','nos3','nos4','nos5','nos6','nos7','plat362','plat1919','plbuckle','sts4098',...
    't2dal_e','t3dl_e','Trefethen_20','Trefethen_20b','Trefethen_150','Trefethen_200','Trefethen_200b','Trefethen_300','Trefethen_500',...
    'Trefethen_700','Trefethen_2000','Trefethen_20000b','CO'};

tol = 1.0e-08;
t = 60;
neigs = 10;
orth_strategies = {'local','full','periodic','selective'};
ritz_norms = zeros(length(matrices),1);
orth_error = zeros(length(matrices),1);

for matrix = matrices

    load([matrices_path,'/',cell2mat(matrix)]);
    A = Problem.A;
    r = rand(size(A,1),1);

    eigs_ref = sort(eigs(A,10,'lm'),'descend');
    
    inf_norm = norm(A,'inf');
    A = A/inf_norm;
    
    for orth = orth_strategies

        outfilename = ['ca_lanczos_restart_',cell2mat(matrix),'_',cell2mat(orth),'_',num2str(t),'_',num2str(neigs),'.out'];
        outfile = fopen([output_path,'/',outfilename],'w');

        fprintf(outfile,'%s\n',cell2mat(matrix));
        fprintf(outfile,'%.4e\n',tol);
        fprintf(outfile,'s,1,2,4,6,8,10\n');
                
        % Standard Lanczos
        time = tic;
        [e,Q,n_restarts,ritz_norms,orth_error] = restarted_lanczos(A,r,t,neigs,orth,tol);
        time = toc(time);
        full_ritz_norm = zeros(neigs,1);
        if ~isempty(e)
            for i = 1:neigs
                full_ritz_norm(i) = norm(A*Q(:,i)-e(i)*Q(:,i))/norm(e(i)*Q(:,i));
            end
            max_err_ref = max(abs(eigs_ref-e*inf_norm)./(e*inf_norm));
            fprintf(outfile,'%.6e,%.6e,%.6e,%d,%.6e\n',max(full_ritz_norm),max(orth_error),max_err_ref,n_restarts,time);
        else
            fprintf(outfile,'Did not converge.\n');
        end
        
        % CA-Lanczos
        for s = [1,2,4,6,8,10]
            time = tic;
            [e,Q,n_restarts,ritz_norms,orth_error] = restarted_ca_lanczos(A,r,t,neigs,s,'newton',orth,tol);
            time = toc(time);
            full_ritz_norm = zeros(neigs,1);
            if ~isempty(e)
                for i = 1:neigs
                    full_ritz_norm(i) = norm(A*Q(:,i)-e(i)*Q(:,i))/norm(e(i)*Q(:,i));
                end
                max_err_ref = max(abs(eigs_ref-e*inf_norm)./(e*inf_norm));
                fprintf(outfile,'%.6e,%.6e,%.6e,%d,%.6e\n',max(full_ritz_norm),max(orth_error),max_err_ref,n_restarts,time);
            else
                fprintf(outfile,'Did not converge.\n');
            end
        end
        
        fclose(outfile);
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
