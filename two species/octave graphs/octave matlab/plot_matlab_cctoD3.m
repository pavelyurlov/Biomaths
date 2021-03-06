function plot_matlab_cctoD3()    
    d12 = linspace(0, 0.0015, 100);
    sm2 = linspace(0.001, 0.2, 100);
    
    N1_ans = dlmread('matlab_cctoD3_N1.txt');
    N2_ans = dlmread('matlab_cctoD3_N2.txt');    

    
    figure;
    hold on;
    grid on;
    surf(d12, sm2, N1_ans, 'FaceColor', [0 0 1]);
    surf(d12, sm2, N2_ans, 'FaceColor', [1 0 0]);
    title('Matlab CCTO D3 (N1 in blue, N2 in red)');
    xlabel('d12');
    ylabel('sm2');

    
    saveas(gcf, 'fig_matlab_cctoD3.jpg');
    
    
end
