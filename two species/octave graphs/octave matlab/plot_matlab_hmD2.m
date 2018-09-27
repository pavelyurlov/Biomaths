function plot_matlab_hmD2()    
    sw1 = linspace(0.001, 0.15, 100);
    sw2 = linspace(0.001, 0.15, 100);
    
    N1_ans = dlmread('matlab_hmD2_N1.txt');
    N2_ans = dlmread('matlab_hmD2_N2.txt');    

    
    figure;
    hold on;
    grid on;
    h1 = surf(sw1, sw2, N1_ans, 'FaceColor', 'blue');
    h2 = surf(sw1, sw2, N2_ans, 'FaceColor', 'red');
    title('Matlab HM D2 (N1 in blue, N2 in red)');
    xlabel('sw1');
    ylabel('sw2');
    
    
    saveas(gcf, 'fig_matlab_hmD2.jpg');
    
    
end
