function [centers,values,plot_range] = log_histogram_pdf(du,numpoints)
if(nargin<2)
    numpoints = 30; % was 200
end
[~,edges] = histcounts(du,'normalization','pdf');
logspace_edges = logspace(log10(min(abs(edges(edges>0)))),log10(max(abs(edges))),numpoints);
[values,edges] = histcounts(du,logspace_edges,'normalization','pdf');
%centers = (edges(1:end-1)+edges(2:end))/2;
centers = sqrt(edges(1:end-1).*edges(2:end)); % use geometric mean for centers
plot_range = logspace(log10(0.5*min(abs(edges(edges>0)))),log10(2*max(abs(edges))),1e3);
%fprintf('Maximum of plot range is %3.3g, max of centers is %3.3g \n',max(plot_range),max(centers))