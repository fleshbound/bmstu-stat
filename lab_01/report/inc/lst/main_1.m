for i = 1 : m
    count = 0;
    
    for x = X
        if (i == m) && (x >= J_limits(i))
           && (x <= J_limits(i + 1))
            count++;
        elseif (x >= J_limits(i)) && (x < J_limits(i + 1))
            count++;
        endif
    endfor
    
    if (i == m)
        fprintf("\t%d) [%+.3f; %+.3f), n%d = %d\n", i, 
        		J_limits(i), J_limits(i + 1), i, count);
    else
        fprintf("\t%d) [%+.3f; %+.3f], n%d = %d\n", i, 
        		J_limits(i), J_limits(i + 1), i, count);
    endif
    
    ni(i) = count;
endfor

fprintf("\nд) построение на одной координатной плоскости \
		гистограммы\ \n   и графика функции плотности расп\
		редления вероятностей \n   нормальной случайной ве\
		личины с математическим\ \n   ожиданием Mu и диспер\
		сией S_quad\n");

J_middles = zeros(m, 1);

for i = 1 : m
    J_middles(i) = (J_limits(i) + J_limits(i + 1)) / 2;
endfor

fn_values = zeros(m, 1);

for i = 1 : m
    fn_values(i) = ni(i) / (n * delta);
endfor
