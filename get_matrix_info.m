outfilename = '../matrix_info.txt';
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

outfile = fopen(outfilename,'w');
fprintf(outfile,'matrix,size,cond,norm,eig_max,eig_min\n');

for matrix = matrices
    load([matrices_path,'/',cell2mat(matrix)]);
    A = Problem.A;
    mat_size = size(A);
    if mat_size(1) ~= mat_size(2)
        disp(['ERROR: Matrix ' matrix ' is not square.']);
    else
        mat_size = mat_size(1);
    end
    mat_cond = condest(A); 
    mat_norm = normest(A);
    eig_max = eigs(A,1,'lm');
    eig_min = eigs(A,1,'sm');
    fprintf(outfile,'%s,%d,%.4e,%.4e,%.4e,%.4e\n',cell2mat(matrix),mat_size,mat_cond,mat_norm,eig_max,eig_min);
end