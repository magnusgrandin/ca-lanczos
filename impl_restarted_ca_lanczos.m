function [conv_eigs,Q_conv,num_restarts] = impl_restarted_ca_lanczos(A, r, max_lanczos, n_wanted_eigs, s, basis, orth, tol)

    %% Define constants.
    max_restarts = 40;

    %% Check input arguments.
    if nargin < 3
        disp('ERROR: Wrong number of input arguments.')
    end
    if nargin < 4 || isempty(n_wanted_eigs)
        n_wanted_eigs = 10;
    end
    if nargin < 5 || isempty(s)
        s = 6;
    end
    if nargin < 6 || isempty(basis)
        basis = 'newton';
    end
    if nargin < 7 || isempty(orth)
        orth = 'local';
    else
        if isnumeric(orth)
            orth = num2str(orth);
        end
        if strcmpi(orth,'local')==0 && strcmpi(orth,'full')==0 ...
                && strcmpi(orth,'selective')==0 && strcmpi(orth,'periodic')==0
            disp(['lanczos.m: Invalid option value for orth: ', orth]);
            disp('    expected {''local''|''full''|''periodic''|''selective''}');
            return;
        end
    end
    
    %% Set the tolerance, adjust according to the norm of A.
    if nargin < 8 || isempty(tol)
        tol = 1.0e-06;
    end
    norm_A = normest(A);
    tol = tol*norm_A;

    if strcmpi(orth,'local')
        disp('Local orthogonalization');
    elseif strcmpi(orth,'full')
        disp('Full orthogonalization');
    end

    % TODO: Take this as input argument.
    restart_strategy = 'largest'; % 'largest','smallest','closest_conv','random'

    % Vector length
    n = length(r);
        
    % Normalize the initial vector.
    Vm(:,1) = r/norm(r);
    q = [];
    
    % Fix change-of-basis matrix
    Bk = setChangeOfBasisMatrix(A,Vm(:,1),basis,s);

    Vk = zeros(n,n_wanted_eigs);
    Tm = zeros(max_lanczos+1,max_lanczos);
    conv_eigs = [];
    conv_rnorms = [];
    orth_err = [];

    if mod(n_wanted_eigs,s) ~= 0
        % TODO: Fix!!
        disp('Warning: Number of wanted eigs is not a multiple of s.');
        return;
    end
    
    k = n_wanted_eigs+4;
    p = s*floor((max_lanczos-k)/s);
    m = k+p;
       
    iter = 0;
    restart = true;
    nconv = 0;
    while(restart && (iter < max_restarts))
        
        iter = iter+1;
    
        if iter == 1
            [Vm,Tm] = lanczos_basic(A,Vm,Tm,q,Bk,m,0,s,basis,orth);            
            %[Vm,Tm] = std_lanczos_basic(A,[],Tm,q,m,0,orth);            
        else
            [Vm,Tm] = lanczos_basic(A,Vm,Tm,q,Bk,m,k,s,basis,orth);
            %[Vm,Tm] = std_lanczos_basic(A,Vm,Tm,q,m,k,orth);
        end
        
        u = selectShifts(Tm(1:m,1:m),restart_strategy);

        Q = eye(m);
        j = m;
        while(j > k)
            [Q,Tm(1:m,1:m)] = qrstep(Q,Tm(1:m,1:m),u(j),1,m);
            if(abs(imag(mu(j))) > 0)
                j = j-2;
            else
                j = j-1;
            end
        end
        
        r = Vm(:,1:m)*Q(:,j+1)*Tm(j+1,j) + r*Q(m,j);
        Vm(:,1:j) = Vm(:,1:m)*Q(:,1:j);
        bk = sqrt(r'*r);
        Vm(:,j+1) = (1/bk)*r;
        Tm(j+1,j) = bk;

        % TODO: deflation.
        % 1. Initially, lock Ritz values as they converge until k values have been locked
        % 2. Continue to iterate, lock every new Ritz value that is better than any of 
        %    the already locked ones. Purge all unwanted but yet converged Ritz values.
        % 3. Continue with (2) until the next Ritz value to converge is not better. Replace 
        %    the (k+1)st basis vector with a random vector and orthogonalize it against the 
        %    previous ones (the k locked ones?). Go back to step (2)
        % 4. When step (3) has been executed two consecutive times with no replacement of 
        %    existing locked Ritz values, the iteration is stopped.
        %
        
        e=eig(Tm(1:k,1:k))
        
    end
    
      
    %         beta = norm(r);
%         x = [];
%         new_conv = 0;
%         for i = nconv+1:k
%             y = Vm(:,i);
%             y(m-nconv)
%             if abs(beta*y(m-nconv)) < tol
%                x = y;
%                nconv = nconv + 1;
%                break;
%             end
%         end
%         Tm = deflate(x,Tm);
%         if nconv >= k
%             restart = false;
%             break;
%         end
        
%         %%%%%%%%%%%%%
%         I = eye(m-nconv);
%         for j = 1:p
%             [Qj,Rj] = qr(Tm(nconv+1:m,nconv+1:m)-u(j)*I);
%             Vm(:,nconv+1:m) = Vm(:,nconv+1:m)*Qj;
%             Tm(nconv+1:m,nconv+1:m) = Qj'*Tm(nconv+1:m,nconv+1:m)*Qj;
% %            Tm(nconv+1:m,nconv+1:m) = Rj*Qj'+ u*I;
%         end
%         %%%%%%%%%%%%%


%       [V,H,f] = Arnold(A,V,H,f,ko,m);
% 
%       [mu, ritz]= select_shifts(H,which);
% 
%       Q = eye(m);
%       j = m;
%       while(j > k)
% %
% %            Apply unwanted eigenvalues as shifts with QR steps
% %
%             [Q,H] = qrstep(Q,H,mu(j),1,m);
% %
% %            Decrease j by 2 if m_j imaginary and by 1 if real
% %
%             if (abs(imag(mu(j))) > 0),
%                j = j-2;
%             else
%                j = j-1;
%             end
%       end
%       ko = j;

%      ritz = norm(f)*ritz(1:ksave);
%      f = V*Q(:,ko+1)*H(ko+1,ko) + f*Q(m,ko);
%      V(:,1:ko) = V*Q(:,1:ko);




    [Vk,Dk] = eig(Tm(1:k,1:k));
    conv_eigs = sort(diag(Dk),'descend');
    Q_conv = zeros(n,k);
    for i = 1:k
        Q_conv(:,i) = Vm(:,1:k)*Vk(:,i);
    end
    
    if ~restart
        disp(['Converged in ' num2str(num_restarts) ' restarts.']);
    else       
        disp(['Did not converge.']);
    end
end

% Fix the change of basis matrix.
function Bk = setChangeOfBasisMatrix(A,q,basis,s)
    Bk = [];
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,q,2*s,'full');
        basis_eigs = eig(T);
        basis_shifts = leja(basis_eigs,'nonmodified');
        Bk = newton_basis_matrix(basis_shifts, s,1);
    end
end

% Select the shifts used in the QR step.
function mu = selectShifts(T,which)
    [V,D] = eig(T);
    if strcmpi(which,'largest')
        mu = sort(abs(diag(D)),'descend');
    elseif strcmpi(which,'smallest')
        mu = sort(abs(diag(D)),'ascend');
    end
end

% Compute matrix powers
function V = matrix_powers(A, q, s, Bk, basis)
    if strcmpi(basis,'monomial')
        V(:,1) = q;
        V(:,2:s+1) = matrix_powers_monomial(A, q, s);
    elseif strcmpi(basis,'newton')
        basis_shifts = diag(Bk);
        V(:,1:s+1) = matrix_powers_newton(A, q, s, basis_shifts, 1);
    end
end

%% Returns a columnvector with n elements, in which all elements are 
%% zero except the last element. 
function vec = eyeshvec(len)
    vec = eye(len,1);
    vec=circshift(vec,len-1);
end
 
function [Q,T] = std_lanczos_basic(A,Q,T,q,maxvecs,prevvecs,orth)

    if nargin < 4
        orth = 'local';
    end
    
    n = size(Q,1);
    
    % If necessary, allocate enough space in Q and T
    if size(Q,2) < maxvecs+1
        Q = [Q, zeros(n,maxvecs+1-size(Q,2))];
    end
    [Trows,Tcols] = size(T);
    if Trows < maxvecs+1 || Tcols < maxvecs
        T_ = T;
        T = zeros(maxvecs+1,maxvecs);
        T(1:Trows,1:Tcols) = T_;
    end
    
    if prevvecs > 0
        b = sqrt(Q(:,prevvecs+1)'*Q(:,prevvecs+1));
    else
        b = T(Trows,Tcols);
    end
    alpha = zeros(1,maxvecs);
    beta = zeros(1,maxvecs);
    
    nvecs = prevvecs+1;
    iter = 0;
    while (nvecs < maxvecs+1) 
        iter = iter+1;
        r = A*Q(:,nvecs);
        if iter > 1
            r=r-beta(iter-1)*Q(:,nvecs-1);
        end
        alpha(iter) = r'*Q(:,nvecs);
        r = r - alpha(iter)*Q(:,nvecs);
        beta(iter) = sqrt(r'*r);
        Q(:,nvecs+1) = r/beta(iter);
        if strcmpi(orth,'full') == 1
            [Q,R_] = normalize(Q);
        end
        nvecs = nvecs + 1;
    end
    Tk = diag(alpha(1:iter)) + diag(beta(1:iter-1),1) + diag(beta(1:iter-1),-1);
    Tk = [Tk; zeros(1,iter-1),beta(iter)];
    
    % Fix output
    if prevvecs > 0
        T11 = T(1:prevvecs,1:prevvecs);
        T12 = b*eyeshvec(prevvecs)*eye(nvecs-prevvecs-1,1)';
        T21 = b*eye(nvecs-prevvecs,1)*eyeshvec(prevvecs)';
        T22 = Tk;
        T = [T11,T12;T21,T22];
    else
        T = Tk;
    end
    Q = Q(:,1:nvecs);
end

function [Q,T] = lanczos_basic(A, Q, T, q, Bk, maxvecs, prevvecs, s, basis, orth)

    if nargin < 8
        orth = 'local';
    end
     
    n = size(Q,1);
       
    % If necessary, allocate enough space in Q and T
    if size(Q,2) < maxvecs+1
        Q = [Q, zeros(n,maxvecs+1-size(Q,2))];
    end
    [Trows,Tcols] = size(T);
    if Trows < maxvecs+1 || Tcols < maxvecs
        T_ = T;
        T = zeros(maxvecs+1,maxvecs);
        T(1:Trows,1:Tcols) = T_;
    end
    
    b = zeros(ceil((maxvecs+1)/s),1);
    if prevvecs > 0
        b(1) = sqrt(Q(:,prevvecs+1)'*Q(:,prevvecs+1));
    else
        b(1) = T(Trows,Tcols);
    end
    
    iter = 0;
    nvecs = prevvecs;
    while nvecs <= maxvecs-s

        iter = iter+1;
        
        q = Q(:,nvecs+1);
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if nvecs == 0
            % Orthogonalize initial basis vectors
            [Q(:,nvecs+1:nvecs+s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            Tk = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(iter+1) = Tk(s+1,s);
            % Extend T
            T(nvecs+1:nvecs+s+1,nvecs+1:nvecs+s) = Tk;
            
        else
            if strcmpi(orth,'local')
                % Orthogonalize against previous block of basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,nvecs-s+1:nvecs+1)},V(:,2:s+1),false);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{3}; 
                Q(:,nvecs+2:nvecs+s+1) = Q_;
            elseif strcmpi(orth,'full')
                % Orthogonality against all previous basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,1:nvecs-s),Q(:,nvecs-s+1:nvecs+1)},V(:,2:s+1),false);
                Rkk_s = Rk_{2};
                Rk_s = Rk_{3};
                Q(:,nvecs+2:nvecs+s+1) = Q_;
            end
            
            % Compute Tk (tridiagonal sub-matrix of T)
            Rkk = [zeros(s,1), Rkk_s(1:s,:)];
            Rk  = [eye(s+1,1), [Rkk_s(s+1,1:s);Rk_s]];
            zk    = Rk(1:s,s+1);
            rho   = Rk(s+1,s+1);
            rho_t = Rk(s,s);
            bk = Bk(s+1,s); %should be one with monomial basis
            e1 = eye(s,1);
            es = eyeshvec(s);
            Tk = Rk(1:s,1:s)*Bk(1:s,:)/Rk(1:s,1:s) ...
                + (bk/rho_t)*zk*es' ...
                - b(iter)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(iter+1) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:nvecs,1:nvecs);
            T12 = b(iter)*eyeshvec(nvecs)*eye(s,1)';
            T21 = b(iter)*eye(s,1)*eyeshvec(nvecs)';
            T22 = Tk(1:s,1:s);
            T31 = zeros(1,nvecs);
            T32 = b(iter+1)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
        end
        
        nvecs = nvecs+s;        
    end
    
    % Fix output
    T = T(1:nvecs,1:nvecs);
    Q = Q(:,1:nvecs+1);
end

function T_new = deflate(y,T)
    T_new = T;
    if isempty(y)
        return;
    else
        k = size(y,1);
        Q = zeros(k,k);
        Q(:,1) = y;
        s = y(1)^2;
        t0 = abs(y(1));
        for j = 2:k
            s = s + y(j)^2;
            t = sqrt(s);
            if t0 ~= 0
                g = (y(j)/t)/t0;
                Q(1:j-1,j) = -y(1:j-1)*g;
                Q(j,j) = t0/t;
            else
                Q(j-1,j) = 1;
            end
            t0 = t;
        end
        T_new(1:k,1:k) = Q'*T(1:k,1:k)*Q;
    end
end

% function [V,H] = deflate(V, H, X, j)
%     
%     % Quick exit if there are no Ritz-pairs to deflate.
%     if isempty(X)
%         return;
%     end
%     
%     m = size(H,2);
%     i = size(X,2);
%     [Q,R] = qr(X);
%     H(j+1:m,j+1:m) = Q'*H(j+1:m,j+1:m)*Q;
%     H(1:j,j+1:m) = H(1:j,j+1:m)*Q;
%     V(:,j+1:m) = V(:,j+1:m)*Q;
%     [P,S] = qr(H(1+j+i:m,1+j+i:m)); %% Not sure about this (step 3 in alg 6.2)
%     %[H(1+j+i:m,1+j+i:m),P] = Hessred(H(1+j+i:m,1+j+i:m));
%     H(1+j+i:m,1+j+i:m) = P'*H(1+j+i:m,1+j+i:m)*P;
%     H(1:j+i,1+j+i:m) = H(1:j+i,1+j+i:m)*P;
%     V(:,1+j+i:m) = V(:,1+j+i:m)*P;   
% end


function [A, P] = Hessred(A)
% Reduction of the square matrix A to the upper
% Hessenberg form using Householder reflectors.
% The reflection vectors are stored in columns of
% the matrix V. Matrix A is overwritten with its
% upper Hessenberg form.
[m,n] =size(A);
if A == triu(A,-1)
   V = eye(m);
   return
end
V = [];
for k=1:m-2
   x = A(k+1:m,k);
   v = Housv(x);
   A(k+1:m,k:m) = A(k+1:m,k:m) - 2*v*(v'*A(k+1:m,k:m));
   A(1:m,k+1:m) = A(1:m,k+1:m) - 2*(A(1:m,k+1:m)*v)*v';
   v = [zeros(k,1);v];
   V = [V v];
end
P = eye(m)-2*V*V';
end

function u = Housv(x)
% Householder reflection unit vector u from the vector x.
m = max(abs(x));
u = x/m;
if  u(1) == 0
   su = 1;
else
   su = sign(u(1));
end
u(1) = u(1)+su*norm(u);
u = u/norm(u);
u = u(:);
end
% function getLejaPoints(a,b,theta_k,theta_p,z_in)
%     p = length(z_in);
%     z = [z_in; zeros(p,1)];
%     a = min(a,theta_k);
%     b = max(b,theta_k);
%     k = p;
%     while k < 2*p
%         if k == 0
%             z0 = b;
%         else
%             w = abs(z(k)-theta_k);
%             for i = 0 
%         end
%         k = k+1;
%     end
% end

function [V,H] = qrstep(V,H,mu,k1,k2)
%
%   Input: V     -- a real square orthogonal matrix
%
%          H     -- a real square upper Hessenberg matrix
%
%          mu    -- a real  or complex shift
%
%          k1,k2 -- pointers to submatrix in k1:k2 block
%
%   Output: V    -- a real orthogonal matrix  
%
%           H    -- a real square upper Hessenberg matrix
%
%                   V <- VQ;   H <- Q'HQ;
%
%                   Q corresponds to a single real shift or 
%                     a double complex shift depending on mu
%
%   D.C. Sorensen
%   2 Mar 2000
%
   k = k2-k1+1;
   kr = k1:k2;
   [m,m] = size(H);

   eta = imag(mu);

   if (abs(eta) > 0), 
%                      mu is imaginary -- apply double shift
%
      xi = real(mu);
      [Q,R] = qr((H(kr,kr) - xi*eye(k))^2 + eta^2*eye(k));

   else  
%                      mu is real -- apply single shift
%
      [Q,R] = qr(H(kr,kr) - mu*eye(k));

   end

   H(kr,:) = Q'*H(kr,:);
   H(:,kr) = H(:,kr)*Q;
   V(:,kr) = V(:,kr)*Q;

%
%   clean up rounding error noise below first subdiagonal
%
   for j = k1:k2,
       H(j+2:m,j) = H(j+2:m,j)*0; 
   end
% Implicitly restarted Arnoldi
% Takes away one single eigenvalue or 
% a complex cojugate pair of eigenvalues
% 
end