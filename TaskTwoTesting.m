clc; clear; close all;

%% ==============================================================
%  Parameters for the systematic test
%  ==============================================================
sizes           = [200, 500, 1000];   % Different n values
density         = 0.01;               % ~1% density for random A
pad_per_row     = 8;                  % Elbow room per row (must be >= inserts per row)
inserts_per_row = 2;                  % How many new elements to insert per row (<= pad_per_row)
num_trials      = 5;                  % Number of trials

fprintf('--- Relaxed CSR test: random A, dense v1; then insert and re-test ---\n');
fprintf('sizes = %s, density = %.3f, pad_per_row = %d, inserts_per_row = %d, trials = %d\n', ...
    mat2str(sizes), density, pad_per_row, inserts_per_row, num_trials);

%% ==============================================================
%  Loop over sizes and trials
%  ==============================================================
for n = sizes
    err_init  = zeros(num_trials, 1);
    err_after = zeros(num_trials, 1);

    for trial = 1:num_trials
        fprintf('\n=== n = %d, trial %d ===\n', n, trial);

        % --- Random sparse A via sprand, and dense v1 ---
        S = sprand(n, n, density);           % Random sparse matrix
        [I, J, V] = find(S);                 % COO triplets
        v1 = randn(n, 1);                     % Dense vector

        % ---- Dense reference ----
        A_dense = full(S);

        % ---- Build relaxed CSR + SpMV (first run) ----
        [rcsr, y_relax] = COO_to_relaxed_CSR_Mv_mult(n, I, J, V, v1, pad_per_row);

        % ---- Dense matvec for comparison ----
        y_dense = A_dense * v1;

        % ---- Error check (before inserts) ----
        err_init(trial) = norm(y_relax - y_dense);
        fprintf('Initial:   ||y_relax - y_dense||_2 = %.3e\n', err_init(trial));

        % ==============================================================
        %   INSERT NEW ELEMENTS using elbow room, then re-test
        % ==============================================================
        rng(trial);  % Different insert pattern per trial

        for i = 1:n
            % How many can we still add to this row without realloc?
            free_slots = rcsr.cap(i) - rcsr.len(i);
            k_ins = min(inserts_per_row, max(free_slots, 0));
            if k_ins == 0, continue; end

            % Existing columns in this row (optional duplicate avoidance)
            row_idx = rcsr.start(i) : (rcsr.start(i) + rcsr.len(i) - 1);
            existing_cols = [];
            if ~isempty(row_idx)
                existing_cols = rcsr.col_idx(row_idx).';
            end

            tries = 0; added = 0;
            while added < k_ins && tries < 10 * k_ins
                j = randi(n);
                if isempty(existing_cols) || ~ismember(j, existing_cols)
                    a_new = randn;
                    % Insert into relaxed CSR (uses elbow room)
                    rcsr = rcsr_insert(rcsr, i, j, a_new);
                    % Update dense reference to match
                    A_dense(i, j) = A_dense(i, j) + a_new;
                    if ~isempty(existing_cols), existing_cols(end + 1) = j; end %#ok<AGROW>
                    added = added + 1;
                end
                tries = tries + 1;
            end
        end

        % ---- Re-run relaxed CSR matvec after insertions ----
        y_relax2 = rcsr_matvec_only(rcsr, v1);

        % ---- Dense matvec again ----
        y_dense2 = A_dense * v1;

        % ---- Error check (after inserts) ----
        err_after(trial) = norm(y_relax2 - y_dense2);
        fprintf('After ins: ||y_relax2 - y_dense2||_2 = %.3e\n', err_after(trial));
    end

    % ===== Summary per n =====
    fprintf('\n>>> Summary for n=%d over %d trials:\n', n, num_trials);
    fprintf('Mean error (initial)   : %.3e   | max: %.3e\n', mean(err_init), max(err_init));
    fprintf('Mean error (after ins) : %.3e   | max: %.3e\n', mean(err_after), max(err_after));
end

fprintf('\nAll tests complete.\n');

%% ==============================================================
%  Local helpers
% ==============================================================

function rcsr = rcsr_insert(rcsr, i, j, a_new)
    % Insert a single (i, j, a_new) into relaxed CSR row i if there's elbow room.
    L = rcsr.len(i);
    if L >= rcsr.cap(i)
        error('rcsr_insert: no elbow room left in row %d (cap=%d).', i, rcsr.cap(i));
    end
    pos = rcsr.start(i) + L;
    rcsr.col_idx(pos) = j;
    rcsr.val(pos) = a_new;
    rcsr.len(i) = L + 1;
end

function b2 = rcsr_matvec_only(rcsr, b1)
    % Just the matvec for an already-built relaxed CSR.
    n = rcsr.n; b1 = b1(:);
    b2 = zeros(n, 1);
    for t = 1:n
        i = rcsr.row_ids(t);
        a = rcsr.start(i);
        L = rcsr.len(i);
        if L > 0
            idx = a : a + L - 1;
            cols = rcsr.col_idx(idx);
            vals = rcsr.val(idx);
            b2(i) = vals.' * b1(cols);
        end
    end
end

function [rcsr, y_relax] = COO_to_relaxed_CSR_Mv_mult(n, I, J, V, v1, pad)
    % Wrapper for COO_to_relaxed_CSR_Mv_mult to match expected TaskTwo signature

    % Call your COO_to_relaxed_CSR_Mv_mult function
    [y_relax, AA_relaxed, JC_relaxed, IA] = COO_to_relaxed_CSR_Mv_mult(V, I, J, v1, n, pad);

    % Build the relaxed CSR structure expected by the test harness
    rcsr.n        = n;
    rcsr.row_ids  = (1:n).';                     % Row indices
    rcsr.start    = IA(1:end-1);                 % Starting index per row
    rcsr.cap      = diff(IA);                    % Capacity per row (includes elbow room)
    rcsr.len      = rcsr.cap - pad;              % Actual number of elements per row
    rcsr.col_idx  = JC_relaxed;                  % Column indices
    rcsr.val      = AA_relaxed;                  % Values
end
