% -- File: scripts/Run_isoLoss_single.m
% L = Var(pi) + omega_x*Var(x) + omega_r*Var(r)
% Noise maps: (sigma_nu^CB, sigma_eta), (sigma_nu^CB, sigma_nu^HH)
% Policy maps: (phi_pi, sigma_nu^CB) and (phi_pi, {cry|crdy})
% ======================================================================
clear; clc;

% ---- Paths ----
workdir = '...'; cd(workdir);
addpath('...Dynare/5.5/matlab'); dynare_config;

% ---- Plots dir + save options ----
plots_dir     = fullfile(workdir,'Plots'); if ~exist(plots_dir,'dir'), mkdir(plots_dir); end
SAVE_FORMATS  = {'jpg'};  SAVE_OVERWRITE = true; SAVE_DPI = 300;

H_var = 1200;  eps_all = 1e-8;
% Grids
sig_list     = [0 0.05 0.10 0.20 0.50 1.00 1.50];       % FEVD/std scan
sig_nu_vecA  = [0 0.05 0.10 0.20 0.50 1.00 1.50];       % Iso-A: σν^CB (X)
sig_eta_vecA = [0 0.05 0.10 0.20 0.50 1.00 1.50];       % Iso-A: ση (Y)
phi_pi_grid  = [1.25 1.50 1.75 2.00 2.25 2.50];         % Policy-C: φπ (X)
sig_cb_vecC  = [0 0.05 0.10 0.20 0.50 1.00];            % Policy-C: σν^CB (Y)

% Welfare weights (π and x fixed at 1.0; wr varies across runs)
WELFARE_MODE   = 'adhoc'; omega_pi_adhoc = 1.0; omega_x_adhoc = 1.0;
RELATIVE_MODE  = 'pct';   K_CE = 1.0;

% ---- Dynare ----
dynare('Base_noise','noclearall','nolog');

% Indexing
endo = cellstr(M_.endo_names); exo = cellstr(M_.exo_names);
ix    = @(nm) find(strcmp(endo,nm));
ix_r  = ix('r');  ix_pi = ix('pinf');  ix_y = ix('y');  ix_yf = ix('yf'); %#ok<NASGU>

ix_cb  = find(strcmp(exo,'evpime'));         % CB inflation measurement noise (rule)
ix_eta = find(strcmp(exo,'eobsr'));          % r observation noise
assert(~isempty(ix_cb) && ~isempty(ix_eta), 'Missing evpime/eobsr in .mod');

% Keep Σ consistent with parameter for r obs noise
if isnan(get_param_by_name('sigv_r_obs')), set_param_value('sigv_r_obs',0.0); end
M_.Sigma_e(ix_eta,ix_eta) = get_param_by_name('sigv_r_obs')^2;
Sigma_base = M_.Sigma_e;

% ---- One IRF pass at current parameters to build S rows ----
Srow = compute_Srows(H_var, eps_all, {'r','pinf','y','yf'});
Sp = Srow.Sp; Sx = Srow.Sx; SrL = Srow.SrL;
base_sig2 = diag(Sigma_base).'; base_sig2 = base_sig2(:).';

% =========================
% FEVD & std(r) vs σν^CB
% =========================
FEVD_r_cb=zeros(numel(sig_list),1); FEVD_pi_cb=zeros(numel(sig_list),1); std_r_all=zeros(numel(sig_list),1);
for j=1:numel(sig_list)
    sig2=base_sig2; sig2(ix_cb)=sig_list(j)^2; sig2(ix_eta)=0;
    var_r=sum(SrL.*sig2); var_pi=sum(Sp.*sig2);
    FEVD_r_cb(j)=100*(sig2(ix_cb)*SrL(ix_cb))/max(var_r,eps);
    FEVD_pi_cb(j)=100*(sig2(ix_cb)*Sp (ix_cb))/max(var_pi,eps);
    std_r_all(j)=sqrt(var_r);
end
fig=figure('Name','FEVD(·) ← ν^{CB}');
subplot(1,2,1); plot(sig_list,FEVD_pi_cb,'-o','LineWidth',2); xlabel('\sigma_{\nu}^{CB}'); ylabel('% of Var(\pi)'); title('FEVD(\pi) \leftarrow \nu^{CB}');
subplot(1,2,2); plot(sig_list,FEVD_r_cb ,'-o','LineWidth',2); xlabel('\sigma_{\nu}^{CB}'); ylabel('% of Var(r)');  title('FEVD(r) \leftarrow \nu^{CB}');
save_plot(fig,plots_dir,'FEVD_pi_r_vs_sigCB',SAVE_FORMATS,SAVE_OVERWRITE,SAVE_DPI);

fig=figure('Name','std(r) vs \sigma_{\nu}^{CB}');
plot(sig_list,std_r_all,'-o','LineWidth',2); grid on; xlabel('\sigma_{\nu}^{CB}'); ylabel('std(r)'); title('std(r) under all shocks');
save_plot(fig,plots_dir,'std_r_vs_sigCB',SAVE_FORMATS,SAVE_OVERWRITE,SAVE_DPI);

% =========================
% Dual-loss analysis: wr = 0.077 and wr = 0
% =========================
make_all_plots(1.0, 1.0, 0.077, '_wr0077');
make_all_plots(1.0, 1.0, 0.000, '_wr0');

% ======================================================================
% Helpers
% ======================================================================
function make_all_plots(wpi, wx, wr, suffix)
    global M_;
    plots_dir     = evalin('base','plots_dir');
    SAVE_FORMATS  = evalin('base','SAVE_FORMATS');
    SAVE_OVER     = evalin('base','SAVE_OVERWRITE');
    SAVE_DPI      = evalin('base','SAVE_DPI');
    RELATIVE_MODE = evalin('base','RELATIVE_MODE');
    KCE           = evalin('base','K_CE');
    H_var         = evalin('base','H_var');
    eps_all       = evalin('base','eps_all');

    Sp = evalin('base','Sp'); Sx = evalin('base','Sx'); SrL = evalin('base','SrL');
    base_sig2 = evalin('base','base_sig2');

    sig_nu_vecA  = evalin('base','sig_nu_vecA');
    sig_eta_vecA = evalin('base','sig_eta_vecA');
    phi_pi_grid  = evalin('base','phi_pi_grid');
    sig_cb_vecC  = evalin('base','sig_cb_vecC');

    exo = cellstr(M_.exo_names);
    ix_cb  = find(strcmp(exo,'evpime'));
    ix_eta = find(strcmp(exo,'eobsr'));

    % ---- Iso-A: (σν^CB, ση)
    [XX_A,YY_A,Z_A]=loss_surface_2d_single(sig_nu_vecA,ix_cb,sig_eta_vecA,ix_eta, ...
        Sp,Sx,SrL,base_sig2,wpi,wx,wr,RELATIVE_MODE,KCE);
    plot_iso(XX_A,YY_A,Z_A,'Var(r)', ...
        '\sigma_{\nu}^{CB} (CB inflation mismeasurement)', '\sigma_{\eta} (obs noise in r)', ...
        wpi,wx,wr,RELATIVE_MODE,plots_dir,['Iso_A_sigCB_vs_sigETA' suffix],SAVE_FORMATS,SAVE_OVER,SAVE_DPI);

    % ---- Policy-C: (φπ, σν^CB)
    [XX_C,YY_C,Z_C] = policy_surface_phi_single(sig_cb_vecC, phi_pi_grid, ix_cb, ...
        base_sig2, wpi, wx, wr, RELATIVE_MODE, KCE, H_var, eps_all);
    plot_iso(XX_C,YY_C,Z_C,'Var(r)', '\phi_\pi (policy response)', '\sigma_{\nu}^{CB}', ...
        wpi,wx,wr,RELATIVE_MODE,plots_dir,['PolicyC_phiPI_vs_sigCB' suffix],SAVE_FORMATS,SAVE_OVER,SAVE_DPI);
end

function [Srow] = compute_Srows(H, eps_all, endo_list)
    global M_ oo_ options_
    M_.Sigma_e = eps_all*eye(M_.exo_nbr);
    options_.order=1; options_.irf=H; options_.nograph=1; options_.noprint=1;
    [~, oo_, options_, M_] = stoch_simul(M_, options_, oo_, endo_list); %#ok<ASGLU>
    exo_names = cellstr(M_.exo_names);
    psi_r  = get_irf_matrix(oo_.irfs, exo_names,'r_',   H);
    psi_pi = get_irf_matrix(oo_.irfs, exo_names,'pinf_',H);
    psi_y  = get_irf_matrix(oo_.irfs, exo_names,'y_',   H);
    psi_yf = get_irf_matrix(oo_.irfs, exo_names,'yf_',  H);
    psi_x  = psi_y - psi_yf;
    Srow.Sp  = sum(psi_pi.^2,1); Srow.Sp  = Srow.Sp(:).';
    Srow.Sx  = sum(psi_x.^2 ,1); Srow.Sx  = Srow.Sx(:).';
    Srow.SrL = sum(psi_r.^2 ,1); Srow.SrL = Srow.SrL(:).';
end

function M = get_irf_matrix(S, exo_names, prefix, H)
    nS=numel(exo_names); M=zeros(H,nS);
    for k=1:nS
        fld=[prefix exo_names{k}];
        if isfield(S,fld)
            v=S.(fld); L=min(H,numel(v)); M(1:L,k)=reshape(v(1:L),[],1);
        end
    end
end

function [XX,YY,Z] = loss_surface_2d_single(gridX, ixX, gridY, ixY, ...
        Sp,Sx,SrL, base_sig2, wpi,wx,wr, REL_MODE, KCE)
    Sp=Sp(:).'; Sx=Sx(:).'; SrL=SrL(:).'; base_sig2=base_sig2(:).';
    NX=numel(gridX); NY=numel(gridY);
    L = zeros(NY,NX);
    for bx=1:NX
        for ay=1:NY
            sig2 = base_sig2;
            sig2(ixX)=gridX(bx)^2; sig2(ixY)=gridY(ay)^2;
            Vpi=sum(Sp.*sig2); Vx=sum(Sx.*sig2); VrL=sum(SrL.*sig2);
            L(ay,bx) = wpi*Vpi + wx*Vx + wr*VrL;
        end
    end
    Z = rel_surface(L, REL_MODE, KCE);
    [XX,YY]=meshgrid(gridX,gridY);
end

function [XX,YY,Z] = policy_surface_phi_single(sig_cb_vec, phi_vec, ix_cb, ...
        base_sig2, wpi,wx,wr, REL_MODE, KCE, H, eps_all)
    NX=numel(phi_vec); NY=numel(sig_cb_vec);
    L=zeros(NY,NX);
    for bx=1:NX
        set_param_value('crpi', phi_vec(bx));
        Srow = compute_Srows(H, eps_all, {'r','pinf','y','yf'});
        Sp=Srow.Sp; Sx=Srow.Sx; SrL=Srow.SrL;
        bsig2 = base_sig2(:).';
        for ay=1:NY
            sig2=bsig2; sig2(ix_cb)=sig_cb_vec(ay)^2;
            Vpi=sum(Sp.*sig2); Vx=sum(Sx.*sig2); VrL=sum(SrL.*sig2);
            L(ay,bx)=wpi*Vpi + wx*Vx + wr*VrL;
        end
    end
    Z = rel_surface(L, REL_MODE, KCE);
    [XX,YY]=meshgrid(phi_vec, sig_cb_vec);
end

function Z = rel_surface(L, mode, KCE)
    L0 = L(1,1);
    switch lower(mode)
        case 'abs', Z = L - L0;
        case 'pct', Z = 100*(L/max(L0,eps) - 1);
        case 'ce',  Z = KCE*(L - L0);
        otherwise, error('RELATIVE_MODE bad');
    end
end

function val = get_param_by_name(pname)
    global M_
    idx = find(strcmp(cellstr(M_.param_names),pname),1);
    if isempty(idx), val = NaN; else, val = M_.params(idx); end
end

function plot_iso(XX,YY,Z, tag, xlab, ylab, wpi,wx,wr, REL_MODE, plots_dir, base_name, SAVE_FORMATS, SAVE_OVERWRITE, SAVE_DPI)
    levels = linspace(min(Z(:)), max(Z(:)), 12);
    fig = figure('Name',['Iso-loss ' tag]);
    [C,h]=contour(XX,YY,Z,levels,'LineWidth',1.25);
    clabel(C,h,'FontSize',10,'Color','k'); grid on; box on;
    xlabel(xlab); ylabel(ylab);
    ttl = sprintf('L (%s) = %.2f Var(\\pi) + %.2f Var(x) + %.3f %s', ...
        loss_ylabel(REL_MODE), wpi, wx, wr, tag);
    title(ttl);
    save_plot(fig, plots_dir, base_name, SAVE_FORMATS, SAVE_OVERWRITE, SAVE_DPI);
end

function yl = loss_ylabel(mode)
    if strcmpi(mode,'pct'), yl='% of baseline'; elseif strcmpi(mode,'ce'), yl='CE units'; else, yl='relative'; end
end

function save_plot(fig, plots_dir, base_name, SAVE_FORMATS, SAVE_OVERWRITE, SAVE_DPI)
    for k=1:numel(SAVE_FORMATS)
        ext=lower(SAVE_FORMATS{k}); fpath=fullfile(plots_dir,[base_name '.' ext]);
        if ~SAVE_OVERWRITE && exist(fpath,'file')
            i=1; while exist(fullfile(plots_dir,sprintf('%s_%02d.%s',base_name,i,ext)),'file'), i=i+1; end
            fpath=fullfile(plots_dir,sprintf('%s_%02d.%s',base_name,i,ext));
        end
        try, exportgraphics(fig, fpath, 'Resolution', SAVE_DPI);
        catch
            set(fig,'PaperPositionMode','auto');
            switch ext
                case 'jpg', print(fig,fpath,'-djpeg',['-r' num2str(SAVE_DPI)]);
                case 'png', print(fig,fpath,'-dpng', ['-r' num2str(SAVE_DPI)]);
                case 'pdf', print(fig,fpath,'-dpdf');
                otherwise,  saveas(fig,fpath);
            end
        end
    end
end

