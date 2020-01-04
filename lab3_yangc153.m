%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: yangc153
%
% Student Number: 1005367265
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

    
%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
%
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
%
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
%
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.
%(a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
%for part a, the IEM approximation and ode45 give similar results, with
%ode45's plot being less smooth and more step like/jagged in between
%estimates, but not to the extent where the solution plots deviate.
%However, at t=pi/2, the improved euler method experiences a noticeable
%dip, while the ode 45 estimate does not, most likely due to this point
%being ignored by it due to its adaptive algorithm.
f1 = @(t,y) y*tan(t)+sin(t); 
t00 = 0;
t10 = pi;
y00 = -1/2;
solna = ode45(f1, [t00, t10], y00);
[x1,y1]=IEM(t00,t10,y00,0.01,f1);
subplot(2,2,1);
plot(x1,y1,'blue');
hold on
plot(solna.x,solna.y,'red')
xlabel('t')
ylabel('y')
legend('improvedeuler','ode45estimate','Location','best')
title('Excersize 2a IEM vs ODE45')
%(b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
%the relationship between the ode45 and IEM solutions is similar to the 
%observed relationship in part a, except for a lack of deviation at pi/2.
f2 = @(t,y) 1/(y^2); 
t01 = 1;
t11 = 10;
y01 = 1;
solnb = ode45(f2, [t01, t11], y01);
[x2,y2]=IEM(t01,t11,y01,0.01,f2);
subplot(2,2,2);
plot(x2,y2,'blue');
hold on
plot(solnb.x,solnb.y,'red')
xlabel('t')
ylabel('y')
legend('improvedeuler','ode45estimate','Location','best')
title('Excersize 2b IEM vs ODE45')

% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
%the relationship between the ode45 and IEM solutions is similar to the 
%observed relationship in part a, except for a lack of deviation at pi/2.
f3 = @(t,y) 1-(t*y/2); 
t02 = 0;
t12 = 10;
y02 = -1;
solnc = ode45(f3, [t02, t12], y02);
[x3,y3]=IEM(t02,t12,y02,0.01,f3);
subplot(2,2,3);
plot(x3,y3,'blue');
hold on
plot(solnc.x,solnc.y,'red')
xlabel('t')
ylabel('y')
legend('improvedeuler','ode45estimate','Location','best')
title('Excersize 2c IEM vs ODE45')
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
%ode45 fails to graph past around t= 0.51(Warning: Failure at t=5.066046e-01.  Unable to meet integration tolerances without reducing the step size below the
%smallest value allowed (1.776357e-15) at time t. ) is produced as it grows asymptotically, while 
%the IEM fails to do this, continuing to graph a solution until around
%t=0.54, where it also grows asymptotically. ode45 reaches anasymptoic
%growth section first, but the end behaviors appear to be the same.
f4 = @(t,y) (y^3)-(t^2); 
t03 = 0;
t13 = 1;
y03 = 1;
solnd = ode45(f4, [t03, t13], y03);
[x4,y4]=IEM(t03,t13,y03,0.01,f4);
subplot(2,2,4);
semilogy(solnd.x,solnd.y,'red',...
         x4,y4,'blue')
xlabel('t')
ylabel('y')
legend('improvedeuler','ode45estimate','Location','best')
title('Excersize 2d IEM vs ODE45'); grid on;

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
xc=euler(inline('2*x.*(1-y.^2).^0.5','x','y'), 0, 0:0.01:0.5);
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
%y=sin(t^2+c),y(0)=0 => c=0
%y=sin(t^2)
yex=sin(0.5^2)
% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%En<=((1+M)*del_t/2)*(e^(M*del_tn)-1),where M is equal to a number >0 that is 
%greater than f(f being y'), the partial derivative of f with respect to
%t, and the partial derivative of f with respect to y at the point (tn,yn).
%del_t is the time step, and del_tn is the time step multipled by the number
%of steps already taken. En is the estimate for error of the nth point.

% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
syms f(t,y);
f(t,y)=2*t.*(1-y.^2).^0.5;
delt_f=diff(f,t);
dely_f=diff(f,y);
newf=matlabFunction(f);
newdely_f=matlabFunction(dely_f);
newdelt_f=matlabFunction(delt_f);
Mmatrix=[newf(0.5,xc(end)),newdely_f(0.5,xc(end)),newdelt_f(0.5,xc(end))];
M=max(Mmatrix);
En=((1+M)/2)*0.01*(exp(M*0.01*length(0:0.01:0.5))-1)
actualerror=abs(yex-xc(end))
%The error for En is larger than that of the actual error, which is
%expected due to the fact that En gives us the maximum bound of the error
%(hence the use of M), while the actual error will always be equal to or
%less than this value.
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.
xc1=euler(inline('2*x.*(1-y.^2).^0.5','x','y'), 0, 0:0.02:0.5);
Mmatrix1=[newf(0.5,xc1(end)),newdely_f(0.5,xc1(end)),newdelt_f(0.5,xc1(end))];
M1=max(Mmatrix1);
En1=((1+M1)/2)*0.02*(exp(M1*0.02*length(0:0.02:0.5))-1)
actualerror1=abs(yex-xc1(end))
%The new error, being double the step size of the old error, produces
%double the error, which indicates a first order method, as doubling the
%step size doubles the error (and vice versa)
errorratio=En1/En
%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.
%The step size update step is trying to keep the h value within a certain
%range, so as to avoid round off error while keeping the step size small
%enough to be accurate (as it should in theory with step sizes approaching
%0). Specifically,if the max term is > 2, the min term is 2. If the max term 
%is < 2, the min term will choose max. The formula is trying to reach an 
%equilibrium, as when D is too innacurate((tol/d)<0.3) it will multiply h 
%by 0.9(always)*0.3 to bring the value of h down. If D is too accurate however((tol/d)>2), 
%it will multiply h by 0.9(always)*2 to bring the value of h up, which 
%helps combat round-off error. If neither case, then h will me multiplied 
%directly with the ratio of tol to D.
%
%
%% Exercise 5
clear all
clf
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
eulerapprox=euler(inline('2*x.*(1-y.^2).^0.5','x','y'), 0, 0:0.025:0.75);
plot(0:0.025:0.75,eulerapprox,'-')
hold on
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
fprime=@(t,y) 2*t.*(1-y.^2).^0.5;
[x,y] =AEM(0,0.75,0,0.025,fprime);
plot(x,y)
% (c) Plot both approximations together with the exact solution.
freal=@(t) sin(t^2);
truex=linspace(0,0.75);
truey=sin(truex.^2);
plot(truex,truey);
hold off
xlabel('t')
ylabel('y')
title('Euler vs Adapted Euler vs Real Solution')
legend('Euler','Adapted Euler','Real Solution','Location','best')

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.

% The adapted euler's method is closer to the actual solution, as ratehr
% than having the step sizes constant and risking over/undershooting the
% estimate (as seen in euler's method), the adapted euler's method
% proactively seeks out an acceptable step size to keep the local error
% under a certain value, therefore ensuring a closer fit.

% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
eulerapprox=euler(inline('2*x.*(1-y.^2).^0.5','x','y'), 0, 0:0.025:1.5);
plot(0:0.025:1.5,eulerapprox,'-')
hold on
fprime=@(t,y) 2*t.*(1-y.^2).^0.5;
[x,y] =AEM(0,1.5,0,0.025,fprime);
plot(x,y)
freal=@(t) sin(t^2);
truex=linspace(0,1.5);
truey=sin(truex.^2);
plot(truex,truey);
hold off
xlabel('t')
ylabel('y')
title('Euler vs Adapted Euler vs Real Solution')
legend('Euler','Adapted Euler','Real Solution','Location','best')
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.
%At y=1, the exact solution's slope becomes 0, and subsequently it begins
%its descent, expected of a trig function such as itself. However, the
%approximations do not accoutnt for complex numbers, and therefore at the
%point y=1 first experienced in the graph, the slopes calculated remain
%either positive for euler's method (as can be seen by the increasing value
%of the approximation), or nearly plateau to 0 (as can be seen by the flat
%approximation given by the adaptive euler's method after the y=1 point is
%first reached). This also implies that the gradient of the AEM is closer
%to that of the exact solution's gradient at this point (gradient=0), but
%again, it does not account for complex numbers and therefore nothing seems
%to be wrong, thus keeping the gradient value near constant.