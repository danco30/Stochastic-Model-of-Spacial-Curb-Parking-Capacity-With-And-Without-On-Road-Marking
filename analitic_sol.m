clear all
clc

lambda = 1/15;
mu = 1/936.65;

for i=1:101
    P = (i-1)/100;
    x_t = creat_generator(lambda, mu, P);
    xt(i,:) = x_t(1e6);
    p_long_miss_park(i) = sum(xt(i,:)) - xt(i,1) - xt(i,5);
    p_short_miss_park(i) = xt(i,8) + xt(i,10);
end

figure()
m = {'+','o','<','.','x','s','d','*','p','h'};
m_counter=1;
for i=1:size(xt,2)
    plot((0:0.01:1).*100,xt(:,i).*100,strcat('-',m{m_counter}));
    m_counter=m_counter+1;
    hold on
end
hold off
title('The Probability To Be At Each State At Time Goes To \infty')
xlabel('Portion Of Long Cars [%]')
ylabel('Probability to Measure this State [%]')
legend({'[0,0,0]' '[0,0,1]' '[0,1,0]' '[0,1,1]' '[1,0,0]' '[1,0,1]' '[1,1,0]' '[1,1,1]' '[0,2]' '[1,2]'},'Location','west')
print -deps analytic_space

figure()
plot((0:0.01:1).*100,p_long_miss_park.*100,'*-',(0:0.01:1).*100,p_short_miss_park.*100,'>-')
title('Probability Of A Cars To Miss Parking After Equilibrium')
xlabel('Portion Of Long Cars [%]')
ylabel('Probability Of A Cars To Miss Parking [%]')
legend({'Long' 'Short'},'Location','southwest')
print -deps analytic_space1

disp 'Done...'

function x_t = creat_generator(lambda, mu, P)
Generator = array2table([...
    -(3*lambda*(1-P)+lambda*P) lambda*(1-P) lambda*(1-P) 0 lambda*(1-P) 0 0 0 lambda*P 0 ;...
    mu -(mu + 2*lambda*(1-P)) 0 lambda*(1-P) 0 lambda*(1-P) 0 0 0 0; ...
    mu 0 -(mu + 2*lambda*(1-P)) lambda*(1-P) 0 0 lambda*(1-P) 0 0 0; ...
    0 mu mu -(2*mu + lambda*(1-P)) 0 0 0 lambda*(1-P) 0 0; ...
    mu 0 0 0 -(mu + 2*lambda*(1-P) + lambda*P) lambda*(1-P) lambda*(1-P) 0 0 lambda*P; ...
    0 mu 0 0 mu -(2*mu + lambda*(1-P)) 0 lambda*(1-P) 0 0; ...
    0 0 mu 0 mu 0 -(2*mu + lambda*(1-P)) lambda*(1-P) 0 0; ...
    0 0 0 mu 0 mu mu -(3*mu) 0 0; ...
    mu 0 0 0 0 0 0 0 -(mu + lambda*(1-P)) lambda*(1-P); ...
    0 0 0 0 mu 0 0 0 mu -(2*mu) ...
    ],'VariableNames',{'[0,0,0]' '[0,0,1]' '[0,1,0]' '[0,1,1]' '[1,0,0]' '[1,0,1]' '[1,1,0]' '[1,1,1]' '[0,2]' '[1,2]'},...
    'RowNames',{'[0,0,0]' '[0,0,1]' '[0,1,0]' '[0,1,1]' '[1,0,0]' '[1,0,1]' '[1,1,0]' '[1,1,1]' '[0,2]' '[1,2]'});

X0 = [1 0 0 0 0 0 0 0 0 0];
G = Generator.Variables;

x_t = @(t) X0*expm(t*G);
end