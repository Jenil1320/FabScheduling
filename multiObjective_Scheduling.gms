Set p "products" /p1*p37/;
Set i "jobs" /diffusion1, implantation1, lithography1, implantation2, diffusion2, lithography2/;
Set m "machines" /diffuser1, diffuser2, implanter1, implanter2, lithographer/;
Set t "time periods" /1*90/;
Set values_alpha "alphas" /1*10/  ;

alias(t,tt);
alias(i,ii);
alias(m,mm);
alias(p,pp);

Table proc_time_data(i,m) "intervals required for job i on machine m"
                    diffuser1   diffuser2  implanter1  implanter2  lithographer
                        
    diffusion1         2          2           0           0             0
    implantation1      0          0           1           1             0
    lithography1       0          0           0           0             2
    implantation2      0          0           1           1             0
    diffusion2         1          1           0           0             0
    lithography2       0          0           0           0             2;

Parameter job_sequence(i,ii) "Mapping of job sequence";

job_sequence('diffusion1','implantation1') = 1;
job_sequence('implantation1','lithography1') = 1;
job_sequence('lithography1','implantation2') = 1;
job_sequence('implantation2','diffusion2') = 1;
job_sequence('diffusion2','lithography2') = 1;
Parameters
wt_m(m);
wt_m('diffuser1')=0.01;
wt_m('diffuser2')=0.02;
wt_m('implanter1')=0.1;
wt_m('implanter2')=0.2;
wt_m('lithographer')=0.1;

Parameter results_alpha(values_alpha) "Results for different alpha values";
Parameter results_makespan(values_alpha) "Results for makespan for different alpha values";
Parameter results_energy(values_alpha) "Results for energy usage for different alpha values";
Scalars
fixed_energy_usage_parameter /10/;
Scalars
alpha /0.00000001001/ ;


Variable
    x(p,i,m,t) "binary variable if product p's job i starts on machine m at time t"
    W(i,m,t)   "indicates if machine m is available at time t to process job i"
    machine_total_time(m)
    weighted_machine(m)
    total_energy_usage
    dummy
    test;
positive variable
    end_time(p,i,m,t) "completion time for product p's job i on machine m ending at t";
variable
    makespan "maximum completion time";

Binary Variable x,W;
*x.fx('p1','diffusion1','diffuser1','1')=1;
*x.fx('p2','diffusion1','diffuser2','1')=1;

Equation 
    job_start "each job of each product starts exactly once",
    machine_availability(i,m,t) "machine can process only one job at a time",
    eligibility(p,i,m,t)
    release_new_product(p,i,m,t)
    product_transition(p,i,ii,m,t) "ensure product jobs are processed in sequence",
    completion_time_update(p,i,m,t) "update the completion time for each product",
    obj "objective function to minimize makespan"
    parallel_processing
    initial_availability
    machine_process(i,m,t)
    machine_time_calc
    weighted_machine_usage
    total_usage
*    trial;


* Each job of each product should start exactly once
job_start(p,i).. sum((m,t), x(p,i,m,t)) =e= 1;

initial_availability(i,m).. W(i,m,'1')=e=1;

* Machine Availability Constraint
machine_availability(i,m,t)$(proc_time_data(i,m) > 0).. sum(p, sum(tt$(ord(tt) >= ord(t) - proc_time_data(i,m) + 1 and ord(tt) <= ord(t)), x(p,i,m,tt))) =l= 1;
   
parallel_processing(p,pp,i,m,mm,t)$(ord(p) <> ord(pp) and ord(m) <> ord(mm) and proc_time_data(i,m) > 0 and proc_time_data(i,mm) > 0)..
       x(p,i,m,t) + x(pp,i,mm,t) =l= 2;
 
* Machine process single product at a time
machine_process(i,m,t)$(proc_time_data(i,m) > 0).. 
    sum(p, x(p,i,m,t)) =l= 1; 
   
* Job-machine eligibility
eligibility(p,i,m,t).. x(p,i,m,t) =l= (proc_time_data(i,m) > 0);

* Completion Time Update
completion_time_update(p,i,m,t)$(proc_time_data(i,m) > 0)..
    end_time(p,i,m,t) =e=  x(p,i,m,t) * (ord(t) + proc_time_data(i,m)-1) ;
 

* Product Transition Constraint
product_transition(p,i,ii,m,t)$(job_sequence(i,ii)).. 
    sum((mm,tt)$(ord(tt) > ord(t) + proc_time_data(i,m)-1), x(p,ii,mm,tt)) =g= x(p,i,m,t);

*Release of New Product Constraint
release_new_product(p,i,m,t)$(ord(p) > 1)..
x(p,i,m,t) =l= x(p-1,i,m,t - proc_time_data(i,m))+ sum(mm$(ord(mm) <> ord(m) and proc_time_data(i,mm) > 0), W(i,mm,t));
    

machine_time_calc(m).. machine_total_time(m) =e= sum((p,i,t), x(p,i,m,t)*proc_time_data(i,m));

weighted_machine_usage(m).. weighted_machine(m) =e= sum((p,i,t), x(p,i,m,t)*proc_time_data(i,m)*wt_m(m));

total_usage.. total_energy_usage =e= sum(m,weighted_machine(m))  ;


* Objective function
*obj.. makespan =e= sum((p,i,m,t)$(proc_time_data(i,m) > 0), x(p,i,m,t) * proc_time_data(i,m));
*obj(i,m,t).. makespan =g= end_time('p15',i,m,t);

obj(p,m,t).. makespan =g= end_time(p,'lithography2',m,t);

*trial.. test=e= alpha*makespan + (1-alpha)*total_energy_usage;;



Model FJSP /all/;

solve FJSP using mip minimising makespan;
*
*Loop(values_alpha,
*   alpha = ord(values_alpha) / 10; 
*    Solve FJSP using mip minimizing test;
*    results_alpha(values_alpha) = test.l;
*     results_makespan(values_alpha) = makespan.l; 
*    results_energy(values_alpha) = total_energy_usage.l; 
*);
*
*Display results_makespan, results_energy, results_alpha;
*
*
*

*
Solve FJSP using mip minimising makespan;
scalar energy_max, makespan_min;
energy_max=total_energy_usage.l;
makespan_min=makespan.l;

*Solve FJSP using mip minimising total_energy_usage;
*scalar energy_min, makespan_max;
*energy_min=total_energy_usage.l;
*makespan_max=makespan.l;


solve FJSP using mip minimising total_energy_usage;
scalar energy_min,makespan_max;
energy_min=total_energy_usage.l;
makespan_max=101;


*Pareto-Optimal Frontier
set front /1*15/;

parameters makespan_opt(front),energy_opt(front), x_opt(front,p,i,m,t);

loop(front, makespan.up = makespan_min + (makespan_max-makespan_min)* (ord(front)-1)/(card(front)-1);
*            total_energy_usage.up = energy_min + (energy_max-energy_min)* (ord(front)-1)/(card(front)-1);
solve FJSP using mip minimising total_energy_usage;
makespan_opt(front)=makespan.l;
energy_opt(front)=total_energy_usage.l;
x_opt(front,p,i,m,t)=x.l(p,i,m,t);
);

makespan_opt(front)$(NOT makespan_opt(front))=EPS;
energy_opt(front)$(NOT energy_opt(front))=EPS;
x_opt(front,p,i,m,t)$(NOT x_opt(front,p,i,m,t))=EPS;

Display x.l, end_time.l, makespan.l,makespan_opt,energy_opt,x_opt, energy_max, makespan_min, makespan_max;

execute_unload 'results37.gdx' x.L end_time.L makespan_opt energy_opt
execute 'gdxxrw.exe results37.gdx o=scresults37.xlsx par=makespan_opt rng=sheet1!A1 par=energy_opt rng=sheet3!A1'
*