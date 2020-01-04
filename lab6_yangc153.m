%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|. Also in this lab, you will write your own
% ODE solver using Laplace transforms and check whether the result yields
% the correct answer.
%
% You will learn how to use the |laplace| routine. 
% 
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in the template, including appropriate descriptions 
% in each step. Save the m-file and submit it on Quercus.
%
% Include your name and student number in the submitted file.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, based on 
% MAT292, Fall 2013, Sinnamon & Sousa

%% Student Information
%
%  Student Name: Chan Hyuk Yang
%
%  Student Number: 1005367265
%

%% Using symbolic variables to define functions
% 
% Recall the use of symbolic variables and function explained in the MATLAB
% assignment #2.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
f=@(t)exp(2*t)*t^3
F=laplace(f(t))
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
g=ilaplace(((s - 1)*(s - 2))/(s*(s + 2)*(s - 3)))
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
syms i(t) a
oglap=laplace(i(t))
shiftlap=laplace(exp(a*t)*i(t))
% (in your answer, explain part (c) using comments).      
%A multilication by e^at in the time domain results in a shift by a in the
%laplace domain. This result from matlab proves this, as the s usually used
%to comptue the laplace transform has been replaced by s-a, indicating a
%shfit by a in the domain
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.


%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
a=2
syms f(t) g(t) G(t) F(t)
g(t)=heaviside(t-a)
G(t)=laplace(g(t)*f(t-a),t,s)
F(t)=laplace(f(t-a),t,s)
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.
%In the case of F, one uses the term t-a when applying the laplace
%transform, and the result of a shift in time is a multiplication by an
%exponential to the power of -(shift term*s), where s is the variable used
%for the laplace transform. Here, G is the laplace transform of f (F)
%multiplied by exp(-2s), satisfying this property

%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,3)+2*diff(y(t),t,2)+diff(y(t),t)+2*y(t)+cos(t) == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),0)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),0)
L_ODE=subs(L_ODE,subs(diff(y(t), t,2), t, 0),0)
% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])


% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.
% The general solution to the ODE is c1*exp(-2t) + c2*sin(t) + c3*cos(t) -
% 1/5 tsin(t) + 1/10 tcos(t). Regardless of the values of c1, c2, and c3, the
% solution will diverge as t goes to infinity. Therefore, no such initial
% conditions for y exist such that it remains bounded as t goes to infinity

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
% First we define the unknown function and its variable and the Laplace
% tranform of the unknown
syms f(t) y(t) t Y s
g(t)=3*heaviside(t)-3*heaviside(t-2)+(t+1)*heaviside(t-2)-(t+1)*heaviside(t-5)+5*heaviside(t-5)
g(t)=3*heaviside(t)-(3-(t+1))*heaviside(t-2)-((t+1)-5)*heaviside(t-5)%matlab likes it better like this idk
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%

% Then we define the ODE

ODE=diff(y(t),t,2)+2*diff(y(t),t)+5*y(t)-g(t) == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),2)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),1)
% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y, [0, 12, 0, 2.25]);


% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%

% In your answer, explain your steps using comments.

%% Exercise 5a
%
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why the following transform is computed correctly.
syms t tau y(tau) s
I=int(exp(-2*(t-tau))*y(tau),tau,0,t)
laplace(I,t,s)
% We know the laplace transform of f*y (convolution definition) = integral of(f(t-tau)y(tau)dtau) 
%from 0 to t is F(s)Y(s). Here f(t) is e^(-2t). Hence, F(s) = 1/(s+2). Matlab returned laplace(y)/(s+2) which is
% just what was expected: Y/(s+2)

%% Exercise 5b
% A particular machine in a factory fails randomly and needs to be replaced. Suppose that the times |t>=0| between failures are independent and identically distributed with probability density function |f(t)|. The mean number of failures |m(t)| at time |t| satisfies the renewal equation |m(t) = \int_0^t [1+m(t-tau)] f(tau) dtau|
%
% Details:  
%
% * Explain why the mean number of failures satisfies this intergal equation. Note that |m(0)=0|.

% m(0) = 0 since no failure happens before t=0.
% Let's say the value for m(t) is known for all 0<=t<T. m(T) can be
% calculated, by noting that f(tau) is just the probability density
% function and that m(T) is the expected values of failures.
% m(T) = \int_0^T x * f(tau).dtau
% Where x is the number of failures happening between tau and T,
% assuming a failure happens at t=tau.
% The expected number of failures happening between (tau,T) is m(T-tau). Since a
% failure happens at t=tau, the total number of failures is m(T-tau)+1. Which
% means x is just m(T-tau)+1. Hence:
% m(T) = \int_0^T (m(T-tau)+1) * f(tau).dtauv

% * Solve the renewal equation for |m(t)| using MATLAB symbolic computation in the cases of i) exponential failure times |f(t) = exp(-t)| and ii) gamma-distributed failure times |f(t) = t^(k-1)/(k-1)! exp(-t)| for natural number |k|. Why does MATLAB have difficulty with the calculation for |k>=5|?
% * Verify the elementary renewal theorem: |m(t)/t| approaches the reciprocal of the mean of |f(t)| as |t| goes to infinity. 
syms f(t) y1(t) y2(t) t tau Y1 Y2

% Exponential decay probability density
f = exp(-t);

% Then we define the ODE
ODE = y1(t) - int((y1(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;

%  Now we compute the Laplace transform of the ODE.
LODE = laplace(ODE);

% We then need to factor out the Laplace transform of |y(t)|
LODE = subs(LODE, laplace(y1), Y1);
Y1 = solve(LODE, Y1);

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP
y1 = ilaplace(Y1)

% Verifying the elementary renewal theorem
% Calculating the average of f(t) (probability density) as t->infinity
average1 = eval(int(t * f, t, 0, inf));

% Multiplying y1(t) / t by average at t=infinity. If theorem is true,
% result should be 1. It is!
result = subs(y1 / t, t, inf) * average1;% sub in a large value going to infinity 
resultval = eval(result)


%k=1
k = 1;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)



%k=2
syms f(t) y1(t) y2(t) t tau Y1 Y2
k = 2;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)



%k=3
syms f(t) y1(t) y2(t) t tau Y1 Y2
k = 3;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)



%k=4
syms f(t) y1(t) y2(t) t tau Y1 Y2
k = 4;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)



% Set up the gamma distribution for the special case k=5
syms f(t) y1(t) y2(t) t tau Y1 Y2
k = 5;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)



%k=6
syms f(t) y1(t) y2(t) t tau Y1 Y2
k = 6;
f = t^k/factorial(k) * exp(-t);

% Repeat the same process with the new distribution
ODE = y2(t) - int((y2(t - tau) + 1) * subs(f, t, tau), tau, 0, t) == 0;
LODE = simplify(laplace(ODE));
LODE = simplify(subs(LODE, laplace(y2), Y2));
Y2 = simplify(solve(LODE, Y2));
y2 = ilaplace(Y2)


% Verify the elementary renewal theorem
% Calculate the average of f(t) as t->infinity. 
average2 = eval(int(t * f, t, 0, inf));

% Multiplying y(t) / t by average at t=infinity. If theorem is true,
% the result should be 1. It is!
result2 = subs(y2 / t, t, inf) * average2;
resultval2 = eval(result)


% for k>=5 the process is slowed down, if you were to solve this problem 
% analytically, Y=k / (s(s+1)^k - ks). Calculating the inverse laplace 
% transform of this function requires using the method of partial
% fractions, which becomes increasingly complicated as the order of the
% denominator increases. Hence, as k increases calculating ilaplace(Y)
% gets more costly with respect to the number of calculations needed, and
% therefore takes more time