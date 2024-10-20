Set p "products" /p1*p37/;
Set i "jobs" /diffusion1, implantation1, lithography1, implantation2, diffusion2, lithography2/;
Set m "machines" /diffuser1, diffuser2, implanter1, implanter2, lithographer/;
Set t "time periods" /1*81/;

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
wt_m(m)  "Energy(Power) weights assigned to each machine ";  
wt_m('diffuser1')=0.01;
wt_m('diffuser2')=0.02;
wt_m('implanter1')=0.1;
wt_m('implanter2')=0.2;
wt_m('lithographer')=0.1;

Scalars
fixed_energy_usage_parameter /10/; 

Variable
    x(p,i,m,t)              "binary variable if product p's job i starts on machine m at time t"
    W(i,m,t)                "indicates if machine m is available at time t to process job i"
    y(p)                    "binary variable if product p has completed all the processes"
    machine_total_time(m)   "indicates the total time a particular machine is used in entire time horizon"
    weighted_machine(m)     "energy consumed by the machines based on assinged weights and usage"
    total_energy_usage      "total energy consumed (variable+fixed)"
    throughput              "total products processed in the time horizon"
    makespan                "maximum completion time"
    ;
    
positive variable
    end_time(p,i,m,t) "completion time for product p's job i on machine m ending at t";

Binary Variable x,W,y;

*x.fx('p1','diffusion1','diffuser1','1')=1;
*x.fx('p2','diffusion1','diffuser2','1')=1;

Equation 
    job_start                          "each job of each product starts exactly once",
    machine_availability(i,m,t)        "machine can process only one job at a time",
    eligibility(p,i,m,t)               "checking if machine is eligible for the operation"
    release_new_product(p,i,m,t)       "constraint to release new products waiting to be processed"
    product_transition(p,i,ii,m,t)     "ensure product jobs are processed in sequence",
    completion_time_update(p,i,m,t)    "update the completion time for each product",
    parallel_processing                "allowing parallel processing if machines are available"
    initial_availability               "initially all the machines are available"
    machine_process(i,m,t)             "machine will single product at a time"
    machine_time_calc                  "calculate total time the machine m is used in time horizon"
    weighted_machine_usage             "variable energy consumption of machines"
    product_completion                 "check if the product p has completed all the processes"
    obj1                               "objective function to minimize makespan"
    obj2                               "objective function to maximise total products being processed (throughput)"                    
    obj3                               "objective function to minimise the energy usage (fixed+variable)"
    ;


* Each job of each product should start exactly once
job_start(p,i).. sum((m,t), x(p,i,m,t)) =e= 1;

initial_availability(i,m).. W(i,m,'1')=e=1;

* Machine Availability Constraint
machine_availability(i,m,t)$(proc_time_data(i,m) > 0).. sum(p, sum(tt$(ord(tt) >= ord(t) - proc_time_data(i,m) + 1 and ord(tt) <= ord(t)), x(p,i,m,tt))) =l= 1;
    
* Allowing parallel processing if machines are available  
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

* Release of New Product Constraint
release_new_product(p,i,m,t)$(ord(p) > 1)..
x(p,i,m,t) =l= x(p-1,i,m,t - proc_time_data(i,m))+ sum(mm$(ord(mm) <> ord(m) and proc_time_data(i,mm) > 0), W(i,mm,t));
    
* Calculate total time the machine m is used in time horizon
machine_time_calc(m).. machine_total_time(m) =e= sum((p,i,t), x(p,i,m,t)*proc_time_data(i,m));

* Variable energy consumption of machines
weighted_machine_usage(m).. weighted_machine(m) =e= sum((p,i,t), x(p,i,m,t)*proc_time_data(i,m)*wt_m(m));

* Check if the product p has completed all the processes
product_completion(p).. sum((i,m,t)$(job_sequence(i,'lithography2')), end_time(p,i,m,t)) =l= card(t) * y(p);

* Objective function
*obj1.. makespan =e= sum((p,m,t), x(p,'lithography2',m,t) * end_time(p,'lithography2',m,t));
obj1(p,i,m,t)$(job_sequence(i,'lithography2')).. makespan =g= end_time(p,i,m,t) + 2*y(p);
*obj1(m,t).. makespan =g= sum(p$(proc_time_data('lithography2',m) > 0),end_time(p,'lithography2',m,t));

* Objective function to maximise total products being processed (throughput)
obj2.. throughput =e= sum(p, y(p));

* Objective function to minimise the energy usage (fixed+variable)
obj3.. total_energy_usage =e= sum(m,weighted_machine(m))+makespan;


**** Multi-Objective formulation for makespan v/s energy usage (parameterise energy to min makespan)****

*Model FJSP /all/;
*
*Solve FJSP using mip minimising makespan;
*scalar energy_max, makespan_min;
*energy_max=total_energy_usage.l;
*makespan_min=makespan.l;
*
*Solve FJSP using mip minimising total_energy_usage;
*scalar energy_min, makespan_max;
*energy_min=total_energy_usage.l;
*makespan_max=makespan.l;
*
**Pareto-Optimal Frontier
*set front /1*11/;
*
*parameters makespan_opt(front),energy_opt(front), x_opt(front,p,i,m,t);
*
*loop(front, total_energy_usage.up = energy_min + (energy_max-energy_min)* (ord(front)-1)/(card(front)-1);
*solve FJSP using mip minimising makespan;
*makespan_opt(front)=makespan.l;
*energy_opt(front)=total_energy_usage.l;
*x_opt(front,p,i,m,t)=x.l(p,i,m,t);
*);
*
*makespan_opt(front)$(NOT makespan_opt(front))=EPS;
*energy_opt(front)$(NOT energy_opt(front))=EPS;
*x_opt(front,p,i,m,t)$(NOT x_opt(front,p,i,m,t))=EPS;
*
*Display x.l, end_time.l, makespan.l,makespan_opt,energy_opt,x_opt, energy_max,energy_min, makespan_min, makespan_max;
*
**execute_unload 'results15.gdx' x.L end_time.L
**execute 'gdxxrw.exe results15.gdx o=scresults15.xlsx var=x.L,end_time.L'
*
Model FJSP /all/;

Solve FJSP using mip minimising total_energy_usage;
*scalar makespan_min, throughput_min, energy_min;
*makespan_min = makespan.l;
*throughput_min = throughput.l;
*energy_min=total_energy_usage.l;
*
*Solve FJSP usig mip maximising throughput;
*scalar makespan_max, throughput_max, energy_max;
*makespan_max = makespan.l;
*throughput_max = throughput.l;
*energy_max=total_energy_usage.l;
*
**Display makespan_max,makespan_min,energy_max,energy_min,throughput_min,throughput_max;
**Solve FJSP using mip minimising total_energy_usage;
**scalar energy_min;
**energy_min = total_energy_usage.l;
**
*****Pareto-Optimal Frontier***
*set front /1*25/;
*
*parameters makespan_opt(front),throughput_opt(front), energy_opt(front), x_opt(front,p,i,m,t);
*
*loop(front, makespan.up = makespan_min + (makespan_max-makespan_min)* (ord(front)-1)/(card(front)-1);
*            total_energy_usage.up = energy_min + (energy_max-energy_min)* (ord(front)-1)/(card(front)-1);
*            
*solve FJSP using mip maximising throughput;
*makespan_opt(front)=makespan.l;
*throughput_opt(front)=throughput.l;
*energy_opt(front)=total_energy_usage.l;
*x_opt(front,p,i,m,t)=x.l(p,i,m,t);
*);
*
*makespan_opt(front)$(NOT makespan_opt(front))=EPS;
*throughput_opt(front)$(NOT throughput_opt(front))=EPS;
*energy_opt(front)$(NOT energy_opt(front))=EPS;
*x_opt(front,p,i,m,t)$(NOT x_opt(front,p,i,m,t))=EPS;
**
*Display x.l, end_time.l, makespan.l,makespan_opt,throughput_opt,energy_opt, x_opt, throughput_max, throughput_min, makespan_min, makespan_max, energy_min, energy_max;
*
*execute_unload 'mydata_dynamicwip_p50t80.gdx', makespan_opt, throughput_opt, energy_opt;
*execute 'gdxxrw.exe mydata.gdx o=mydata.xlsx par=makespan_opt rng=sheet1!A1 par=throughput_opt rng=sheet2!A1 par=energy_opt rng=sheet3!A1';