function plot_cpp_cctoD1()    
    d12 = linspace(0, 0.0015, 100);
    sm2 = linspace(0.0001, 0.2, 100);
    
    N1_ans = dlmread('c++_cctoD1_N1.txt');
    N2_ans = dlmread('c++_cctoD1_N2.txt');    

    
    figure;
    hold on;
    grid on;
    h1 = surf(d12, sm2, N1_ans, 'FaceColor', 'blue');
    h2 = surf(d12, sm2, N2_ans, 'FaceColor', 'red');
    title('C++ CCTO D1 (N1 in blue, N2 in red)');
    xlabel('d12');
    ylabel('sm2');
    
    
    saveas(gcf, 'fig_cpp_cctoD1.jpg');
    
    
end
