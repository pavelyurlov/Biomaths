function main_ccto_d1_full()
    N=512;
    A=2;
    r=linspace(0, A, N);
    h=r(2)-r(1);
    sm1=0.04; sm2=0.1;
    b1=0.4; b2=0.4;
    d1=0.2; d2=0.2;
    d11=0.001; d12=0.0005; d21=0.001; d22=0.001;
    sw11=0.04; sw22=0.04;
    sw12=0.04; sw21=0.04;
    
    m1=b1*normpdf(r, 0, sm1);
    
    
    w11=d11*normpdf(r, 0, sw11);
    
    w21=d21*normpdf(r, 0, sw21);
    w22=d22*normpdf(r, 0, sw22);
    
    al=0.4;
    N1=0;
    N2=0;
    
    d12=linspace(0, 0.0015, 100);
    sm2=linspace(0.000001, 0.2, 100);
    
    N1_ans=zeros(length(d12), length(sm2));
    N2_ans=zeros(length(d12), length(sm2));
    for i=1:length(d12)
        for j=1:length(sm2)
            m2=b2*normpdf(r, 0, sm2(j));
            w12=d12(i)*normpdf(r, 0, sw12);
            
            D11=zeros(1, N);
            D12=zeros(1, N);
            D22=zeros(1, N);
            N1=0;
            N2=0;
            
            [N1, N2, ~, ~, ~, ~]=solver(N1, N2, D11, D12, D22, w11, w12, w21, w22, d11, d12(i), d21, d22, m1, m2, b1, b2, d1, d2, h, A, al, N, 1);
            N1_ans(i, j)=N1;
            N2_ans(i, j)=N2;
            display(100*(i-1)+j);
        end
    end
    dlmwrite('FULL_N1cctoD1_04_d_test.txt', N1_ans);
    dlmwrite('FULL_N2cctoD1_04_d_test.txt', N2_ans);
    figure;
    hold on;
    grid on;
    surf(d12, sm2, N1_ans);
    surf(d12, sm2, N2_ans);
    
end
