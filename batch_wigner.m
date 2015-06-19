%%
% Batch testing to repport error bars

nrep = 20;

ERRLinf = [];
ERRL2 = [];
for irep=1:nrep
    disp(['---- Run #' num2str(irep) ' ----']);
    test_wigner;
    ERRLinf(:,end+1) = ErrLinf(:);
    ERRL2(:,end+1) = ErrL2(:);  
end

%% 
% Display.

y = ERRLinf;
switch param_mode
    case 'kappa'
        x = kappa_list;
    case 'h'
        x = hlist;
end
clf;
shadedErrorBar(x, mean(y,2), 2*std(y, [], 2), {'b' 'LineWidth' 2});
axis tight;
set(gca, 'PlotBoxAspectRatio', [1 1/2 1]);
set(gca,  'FontSize', 20); % 'XTick', [], 'YTick', [-1 1],
saveas(gcf, [rep 'error-' str '.eps'], 'epsc');

