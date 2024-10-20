Set i "jobs" /diffusion1, implantation1, lithography1, implantation2, diffusion2, lithography2/;
Set m "machines" /diffuser1, diffuser2, implanter1, implanter2, lithographer/;
Set t "time periods" /1*120/;

alias(t,tt);

Table proc_time_data(i,m) "intervals required for job i on machine m"
                    diffuser1   diffuser2  implanter1  implanter2  lithographer
    diffusion1        45         15           0           0             0
    implantation1      0          0           6           6             0
    lithography1       0          0           0           0             11
    implantation2      0          0          10          10             0
    diffusion2         5         55           0           0             0
    lithography2       0          0           0           0             2;

Variable
    x(i,m,t) "binary variable if job i starts on machine m at time t";
positive variable
    end_time(i,m,t) "completion time for job i on machine m ending at t";
variable
    makespan "maximum completion time";

Binary Variable x;

Equation 
    job_start(i) "each job starts exactly once",
    machine_usage(i,m,t) "machine can process only one job at a time",
    eligibility(i,m,t) "only valid job-machine assignments allowed",
    define_end_time(i,m,t) "define completion time based on start and duration",
    obj "objective function to minimize makespan"
    endBeforeMakespan(i,m,t);

* Each job should start exactly once
job_start(i).. sum((m,t), x(i,m,t)) =e= 1;

* Machines can only process one job at a time
machine_usage(i,m,t)$(proc_time_data(i,m) > 0)..
    sum(tt$(ord(tt) >= ord(t) and ord(tt) <= ord(t) + proc_time_data(i,m) - 1), x(i,m,tt)) =l= 1;

* Job-machine eligibility
eligibility(i,m,t).. x(i,m,t) =l= (proc_time_data(i,m) > 0);

    
* Define completion times
define_end_time(i,m,t)$(proc_time_data(i,m) > 0)..

end_time(i,m,t) =e= x(i,m,t) * (ord(t) + proc_time_data(i,m) - 1);
    
endBeforeMakespan(i,m,t)$(proc_time_data(i,m) > 0).. end_time(i,m,t) =l= makespan;

* Objective function
obj.. makespan =e= sum((i,m,t)$(proc_time_data(i,m) > 0), end_time(i,m,t));

Model FJSP /all/;

Solve FJSP using mip minimizing makespan;

Display x.l, end_time.l, makespan.l;
