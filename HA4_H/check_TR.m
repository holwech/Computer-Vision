function [ diff ] = check_TR(T1, R1, T2, R2, E, k)

T1_hat = [0 -T1(3) T1(2); 
     T1(3) 0 -T1(1); 
     -T1(2) T1(1) 0];
 
T2_hat = [0 -T2(3) T2(2); 
     T2(3) 0 -T2(1); 
     -T2(2) T2(1) 0];

E_est = T1_hat*R1;
diff = E - E_est

E_est = T2_hat*R2;
diff = E - E_est

E_est = T1_hat*R2;
diff = E - E_est

E_est = T2_hat*R1;
diff = E - E_est
    
end