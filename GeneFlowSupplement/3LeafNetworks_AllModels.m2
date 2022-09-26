restart;
printWidth=60;

--The following computations show that there are no invariants for the 3-leaf 
--3-cycle networks under any of the standard group-based models (CFN,JC,K2P,K3P).
--The computations for each model can be performed by commenting out the
--code for the other models below. 

-- The group element representations for each model

--A = (1,1);C = (1,-1); --CFN
A = (1,1);C = (1,-1);G = (-1,1); T = (-1,-1); --JC,K2P,K3P


--The list of coordinate indices for each model

--L = {(A,A,A),(A,C,C),(C,A,C),(C,C,A)}-- CFN
L = {(A,A,A),(A,C,C),(C,A,C),(C,C,A),(C,G,T)};-- JC
--L = {(A,A,A),(A,A,C),(A,A,G),(A,C,C),(A,C,G),(A,C,T),(A,G,G)}-- K2P
--L = {(A,A,A),(A,A,C),(A,A,G),(A,A,T),(A,C,C),(A,C,G),(A,C,T),(A,G,G),(A,G,T),(A,T,T)} -- K3P


--The identification of coordinates for each model

--x = hashTable{A => 0, C => 1, G => 1, T => 1}--CFN
x = hashTable{A => 0, C => 1, G => 1, T => 1}--JC
--x = hashTable{A => 0, C => 1, G => 1, T => 1}--K2P
--x = hashTable{A => 0, C => 1, G => 2, T => 3}--K3P


--y = hashTable{A => 0, C => 1}--CFN
y = hashTable{A => 0, C => 1, G => 2, T => 3}--JC,K2P,K3P


--The coordinate ring for each model

L1 = toList apply(L,i->q_(y#(i#0),y#(i#1),y#(i#2)));

--The list of parameters for each model

L2 = {a_0,a_1,b_0,b_1,c_0,c_1,h_0,h_1,e_0,e_1,f_0,f_1,g_0,g_1};--JC, CFN
 
 --L2 = {a_0,a_1,a_2,b_0,b_1,b_2c_0,c_1,c_2,h_0,h_1,d_2,
--      e_0,e_1,e_2,f_0,f_1,f_2,g_0,g_1,g_2};--K2P
  
--L2 = {a_0,a_1,a_2,a_3,b_0,b_1,b_2,b_3,c_0,c_1,c_2,c_3,d_0,d_1,d_2,d_3,
--      e_0,e_1,e_2,e_3,f_0,f_1,f_2,f_3,g_0,g_1,g_2,g_3};--K3P

--The ring of variables for each model
 
 R = QQ[(L1|L2), MonomialOrder=>Eliminate 14] --JC,CFN
 --R = QQ[(L1|L2), MonomialOrder=>Eliminate 21] --K2P
 --R = QQ[(L1|L2), MonomialOrder=>Eliminate 28] --K3P
 
  
--The parameterization for the 3-leaf 3-cycle network with no leaf edges 
-- (the only network we need to check)

Eqn_1 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2)) - 
f_(x#(i#1))*g_(x#(i#1#0*i#2#0,i#1#1*i#2#1)) -
e_(x#(i#1))*g_(x#(i#2)));

K = ideal(Eqn_1);
time GG = selectInSubring(1,gens gb(K));




