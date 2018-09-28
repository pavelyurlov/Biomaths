function plot_cpp_cctoD1()    
    d12 = linspace(0, 0.0015, 100);
    sm2 = linspace(0.0001, 0.2, 100);
    
    N1_ans = dlmread('c++_cctoD1_N1.txt');
    N2_ans = dlmread('c++_cctoD1_N2.txt');
    
    sz = rows(N1_ans);
    % assign colours
    N1_colour(:,:,1) = zeros(sz, sz);
    N1_colour(:,:,2) = zeros(sz, sz);
    N1_colour(:,:,3) = zeros(sz, sz);
    
    N2_colour(:,:,1) = zeros(sz, sz);
    N2_colour(:,:,2) = zeros(sz, sz);
    N2_colour(:,:,3) = zeros(sz, sz);   
    
    for i=1:sz
        for j=1:sz
            for k=1:3
                N1_colour(i,j,k)=calculate_colour(1,N1_ans(i,j),k);
                N2_colour(i,j,k)=calculate_colour(2,N2_ans(i,j),k);
            endfor
        endfor
    endfor
    
    figure;
    hold on;
    grid on;
    surf(d12, sm2, N1_ans, N1_colour);
    surf(d12, sm2, N2_ans, N2_colour);
    title('C++ CCTO D1 heatmap (N1 in blue, N2 in red); 1st R < 2st R; d21=0.001');
    xlabel('d12 (competition 1 <- 2)');
    ylabel('sm2 (radius of 2nd )');
    
    saveas(gcf, 'new_fig_cpp_cctoD1.jpg');
    
end
