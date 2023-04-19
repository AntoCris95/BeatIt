
function saveFigures_AC(outdir)

if ~exist(outdir)
    mkdir(outdir)
end

Nfigs = get(gcf, 'number');
for figsN = 1:Nfigs
    
    figX = figure(figsN); Num = figX.Number; Name = figX.Name;
    try
        exportgraphics(figX, fullfile(outdir, sprintf('Fig%d %s.tiff', Num, Name)), 'Resolution',300)
    catch
        saveas(figX, fullfile(outdir, sprintf('Fig%d %s', Num, Name)), 'tiff')
    end
end