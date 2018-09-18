function plot_cpp_hmD3()    
    sw1 = linspace(0.001, 0.15, 100);
    sw2 = linspace(0.001, 0.15, 100);
    
    N1_ans = dlmread('c++_hmD3_N1.txt');
    N2_ans = dlmread('c++_hmD3_N2.txt');    
    
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
    surf(sw1, sw2, N1_ans, N1_colour);
    surf(sw1, sw2, N2_ans, N2_colour);
    title('C++ HM D3 heatmap (N1 in blue, N2 in red)');
    
    set (gca (), "zlim", [-50, 250]);
    
    xlabel('sw1');
    ylabel('sw2');
    
    saveas(gcf, 'new_fig_cpp_hmD3.jpg');
    
end
