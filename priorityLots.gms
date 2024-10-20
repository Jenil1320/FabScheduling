Set p "products" /1*35/;
Set i "jobs" /diffusion1, implantation1, lithography1, implantation2, diffusion2, lithography2/;
Set m "machines" /diffuser1, diffuser2, implanter1, implanter2, lithographer/;
Set t "time periods" /1*80/;
Set priority_products /p1*p10/;

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

Parameter start_time_diffusion1(p) "Start time for diffusion1 for each product";
Parameter end_time_lithography2(p) "End time for lithography2 for each product";
Parameter penalty(p);
penalty(p)$(p.val >= 5 and p.val <= 10) = 1000; 
penalty(p)$(p.val < 5 or p.val > 10) = 1;   

Parameter weight(p) "Weight for each product";
weight(p) = 0.000001;
*weight('p1') = 100;  
*weight('p2') = 100;
*weight('p3') = 100;
*weight('p4') = 100;
*weight('p5') = 100;
*weight('p6') = 100;
*weight('p7') = 100;
*weight('p8') = 100;
*weight('p9') = 100;
*weight('p10') = 100;

Variable
    x(p,i,m,t) "binary variable if product p's job i starts on machine m at time t"
    W(i,m,t)   "indicates if machine m is available at time t to process job i"
    machine_usage(m)
    ;
    
positive variable
    end_time(p,i,m,t) "completion time for product p's job i on machine m ending at t";
variable
    makespan "maximum completion time";

Binary Variable x,W,y;
*x.fx('p1','diffusion1','diffuser1','1')=1;
*x.fx('p2','diffusion1','diffuser2','1')=1;

Equation 
    job_start "each job of each product starts exactly once",
    machine_availability(i,m,t) "machine can process only one job at a time",
    eligibility(p,i,m,t)
    release_new_product(p,i,m,t)
    product_transition(p,i,ii,m,t) "ensure product jobs are processed in sequence",
    completion_time_update(p,i,m,t) "update the completion time for each product",
    obj  "objective function to minimize makespan"
    parallel_processing
    initial_availability
    machine_process(i,m,t)
    machine_time
  
    ;


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

machine_time(m).. machine_usage(m) =e= sum((p,i,t), x(p,i,m,t)*proc_time_data(i,m));

* Product Transition Constraint
product_transition(p,i,ii,m,t)$(job_sequence(i,ii)).. 
    sum((mm,tt)$(ord(tt) > ord(t) + proc_time_data(i,m)-1), x(p,ii,mm,tt)) =g= x(p,i,m,t);

*Release of New Product Constraint
release_new_product(p,i,m,t)$(ord(p) > 1)..
x(p,i,m,t) =l= x(p-1,i,m,t - proc_time_data(i,m))+ sum(mm$(ord(mm) <> ord(m) and proc_time_data(i,mm) > 0), W(i,mm,t));
  
  
* Objective function
*obj.. makespan =e= sum((p,i,m,t)$(proc_time_data(i,m) > 0), x(p,i,m,t) * proc_time_data(i,m));
*obj(p,m,t).. makespan =g= end_time(p,'lithography2',m,t);

obj.. makespan =e= sum(p, penalty(p) * sum((m,t), end_time(p,'lithography2',m,t)));


*obj.. makespan =g= sum((p,m,t), end_time(p,'lithography2',m,t));
Model FJSP /all/;

Solve FJSP using mip minimizing makespan;
loop(p,
    start_time_diffusion1(p) = sum((m,t), ord(t) * x.l(p,'diffusion1',m,t));
    end_time_lithography2(p) = sum((m,t), (ord(t) + proc_time_data('lithography2',m) - 1) * x.l(p,'lithography2',m,t));
);
Display x.l, end_time.l, makespan.l, start_time_diffusion1,end_time_lithography2;

*execute_unload 'results15.gdx' x.L end_time.L
*execute 'gdxxrw.exe results15.gdx o=scresults15.xlsx var=x.L,end_time.L'
