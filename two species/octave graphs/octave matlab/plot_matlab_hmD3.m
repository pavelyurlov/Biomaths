function plot_matlab_hmD3()    
    sw1 = linspace(0.001, 0.15, 100);
    sw2 = linspace(0.001, 0.15, 100);
    
    N1_ans = dlmread('matlab_hmD3_N1.txt');
    N2_ans = dlmread('matlab_hmD3_N2.txt');    

    
    figure;
    hold on;
    grid on;
    h1 = surf(sw1, sw2, N1_ans, 'FaceColor', 'blue');
    h2 = surf(sw1, sw2, N2_ans, 'FaceColor', 'red');
    title('Matlab HM D3 (N1 in blue, N2 in red)');
    
    set (gca (), "zlim", [-50, 250]);

    xlabel('sw1');
    ylabel('sw2');
    
    
    saveas(gcf, 'fig_matlab_hmD3.jpg');
    
    
end
