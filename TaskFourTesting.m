%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing for Task 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. Empirically verify correctness of the DIA (compressed diagonal)
%    format matrix–vector product for multiple matrix sizes n.
% 2. For fixed n, investigate what happens as the diagonal offset ±k grows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% ============================================================
%  PART A — Correctness Verification for Multiple n
% ============================================================
% This section checks that build_diag_and_spmv(C, v1) gives the same result
% as MATLAB's dense C*v1 across multiple matrix sizes and random 3-band
% matrices. 
% -------------------------------------------------------------

sizes      = [50, 100, 200, 400];   % test several matrix dimensions
k          = 5;                     % sub/super diagonal offset
num_trials = 5;                     % number of random experiments per n

fprintf('--- Correctness Test for build_diag_and_spmv ---\n');
fprintf('Testing sizes = %s, offset k = %d, trials = %d\n\n', ...
    mat2str(sizes), k, num_trials);

for n = sizes
    errs = zeros(num_trials,1);

    for t = 1:num_trials
        % ----------------------------------------------------------
        % Construct a random 3-band matrix for given n and k.
        % Main diagonal + one sub-diagonal + one super-diagonal.
        % ----------------------------------------------------------
        C = zeros(n);
        C(1:n+1:end) = randn(n,1);           % main diagonal
        for i = (k+1):n,   C(i,i-k) = randn; end   % subdiagonal (−k)
        for i = 1:(n-k),   C(i,i+k) = randn; end   % superdiagonal (+k)

        % Random dense vector for multiplication
        v1 = randn(n,1);

        % Compute using your DIA-based method
        [DIAG, IOFF, v2] = TaskFour(C, v1); 
        % Compute dense reference
        v2_ref = C * v1;

        % Record 2-norm error between DIA and dense result
        errs(t) = norm(v2 - v2_ref);
        fprintf('n=%d, trial=%d: ||v2 - C*v1||_2 = %.3e\n', n, t, errs(t));
    end

    fprintf('Summary for n=%d:  mean error=%.3e   max error=%.3e\n\n', ...
        n, mean(errs), max(errs));
end


%% ============================================================
%  PART B — Behavior as k Increases (for Different n)
% ============================================================
% Fix several matrix sizes and sweep over multiple offsets k.
% Observe how the number of nonzeros (nnz) and arithmetic work
% change as k grows further from the diagonal.
% -------------------------------------------------------------

sizes = [50, 100, 200];    % test sizes
Ks    = [1 2 3 4 5 8 10];  % increasing offsets to examine
v1    = randn(max(sizes),1);

fprintf('--- Effect of increasing k ---\n');
for n = sizes
    fprintf('\nMatrix size n = %d\n', n);
    fprintf('  k    total_nnz   nnz_sub   nnz_super\n');
    fprintf('--------------------------------------\n');

    for k = Ks
        if k >= n, continue; end  % skip invalid offsets (cannot exceed matrix size)

        % ----------------------------------------------------------
        % Build a simple 3-band matrix with constant diagonals
        % to visualize sparsity structure and storage effects.
        % ----------------------------------------------------------
        C = zeros(n);
        C(1:n+1:end) = 10;              % main diagonal constant
        for i = (k+1):n, C(i,i-k) = 1; end      % sub band (−k)
        for i = 1:(n-k), C(i,i+k) = 2; end      % super band (+k)

        % Compute DIA structure (and optionally v2)
        [DIAG, IOFF, v2] = TaskFour(C, v1(1:n)); 

        % Count true (nonzero) entries per band (ignores zero padding)
        nnz_sub   = n - k;
        nnz_super = n - k;
        nnz_total = n + nnz_sub + nnz_super;

        fprintf('%3d   %10d   %8d   %10d\n', k, nnz_total, nnz_sub, nnz_super);

        % Optional visualization of how padding grows:
        % disp('IOFF ='), disp(IOFF);
        % disp('First few rows of DIAG:'), disp(DIAG(1:min(10,n),:));
    end
end


%% ============================================================
%  Interpretation / Discussion
% ============================================================
% • As k increases (bands move away from main diagonal):
%     → Fewer valid entries in sub/super diagonals: each has (n - k) nonzeros.
%     → Total nonzeros decrease linearly:  nnz_total = 3n - 2k.
% • In DIAG storage:
%     → The sub-diagonal column gains k leading zeros (top padding).
%     → The super-diagonal column gains k trailing zeros (bottom padding).
% • During SpMV:
%     → build_diag_and_spmv indexes only true entries, skipping padded zeros.
%     → Arithmetic work ~ 2 × (3n - 2k) FLOPs, decreasing as k grows.
% • Consequently, for fixed n:
%     → Larger k → fewer multiplications/additions → faster execution.
%
% Example:
%   n = 100, k = 1 → nnz = 298
%   n = 100, k = 10 → nnz = 280
% So runtime (and FLOPs) drop proportionally.