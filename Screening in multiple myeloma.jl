import Pkg; Pkg.add("GR") 

using Plots
gr()
using GR
function simulate_markov_chain(num_individuals::Int, num_steps::Int)
    # Define the initial population fractions
    X = num_individuals * 0.1         # healthy but risk population
    M = num_individuals * 0.0         # MGUS population undetected
    U = num_individuals * 0.0         # MGUS population detected
    S = num_individuals * 0.0         # SMM population undetected
    P = num_individuals * 0.0         # SMM population detected
    N = num_individuals * 0.0         # MM population
    A = num_individuals * 0.0         # Undetected MGUS progressing to Undetected SMM
    B = num_individuals * 0.0         # Undetected MGUS progressing directly to MM
    D = num_individuals * 0.0
    C=0.0
    
    # Define parameters
    d= 0.25 #probability to die in a given year
    l= 0.36 #probability to get screened
    m= 0.02 #probability of MGUS incidence
    n= 0.01 #probability of progression from MGUS to SMM
    x= 0.45 #risk reduction factor after screened with MGUS
    p= 0.096 #probability to get screened for SMM
    r= 0.1 #probability to progress from SMM to MM
    y= 0.45 #risk reduction factor after screened with SMM
    e= 0.2 #death rate due to MM

    # Define empty lists to store the population fractions at each time step
    X_list = [X]; M_list = [M]; U_list = [U]; S_list = [S]; P_list = [P]; N_list = [N];A_list= [A]; B_list= [B]; D_list=[D];C_list= [C]
    
    # Simulate the Markov states for num_time_steps
    for i in 1:num_steps
        X_prev = X; M_prev = M; U_prev = U; S_prev = S; P_prev = P; N_prev = N; A_prev= A; B_prev= B;D_prev=D;C_prev= C
        
        X = max(0,X_prev -M_prev - U_prev - S_prev - P_prev - N_prev-D_prev)
        M = max(0,((1 - d) * M_prev - (1 - l) * (1 - n)) + X_prev * m)
        U = max(0,(1 - d) * M_prev * l * (1 - n) + (U_prev * (1 - d) *(1-(x*n))))
        S = max(0,(S_prev * (1 - d) * (1 - r) * (1 - p))+ ((U_prev) * (1 - d) * (x*n))+ ((M_prev) * (1 - d) * (1-l)*(n)))
        P = max(0,(P_prev * (1 - d) * (1 - (y*r))) + (S_prev * (1 - d) * (1 - r) * p))
        N = max(0,(N_prev * (1 - e)) + (S_prev * (1 - d) * (1 - p) * r) + (S_prev * (1 - d) * x*n) + (P_prev * (1 - d) * y*r))
        D = max(0,(X_prev+ M_prev+ U_prev+S_prev+ P_prev)*d+ N_prev*e)
        # Add the current population fractions to the lists
        push!(X_list, X)
        push!(M_list, M)
        push!(U_list, U)
        push!(S_list, S)
        push!(P_list, P)
        push!(N_list, N)
        push!(D_list,D)
        
        #Population fraction from Undetected MGUS to Detected MGUS to Undetected SMM
        A=  ((U_prev) * (1 - d) * (x*n))
        push!(A_list, A)
        #Population fraction from Undetected MGUS to Undetected SMM
        B=  ((M_prev) * (1 - d) * (1-l)*(n))
        push!(B_list, B)
        
        #mortality due to MM
        C = (N_prev)*e
        push!(C_list, C)
    end
    # Plot the results
    subplot(2, 1, 1)
   
    plot(0:num_steps, X_list, label="Risk")
    plot!(0:num_steps, M_list, label="MGUS")
    plot!(0:num_steps, U_list, label="MGUS Detected")
    plot!(0:num_steps, S_list, label="SMM")
    plot!(0:num_steps, P_list, label="SMM Detected")
    plot!(0:num_steps, N_list, label="MM", xlabel="Time Steps", ylabel="Population Fractions")
    
    subplot(2,1,2)
    plot(0:num_steps, A_list, label="DMGUS_to_USMM")
    plot!(0:num_steps, B_list, label="UMGUS_ti_USMM")
    
end


simulate_markov_chain(100000, 10)

function MGUS_screening(a0, Δa, n_sim, max_age)
    # a0: starting age for screening
    # Δa: interval between screening tests
    # n_sim: number of simulations to run
    # max_age: maximum age to simulate
    
    # Define the possible states
    states = [:healthy, :undetected_MGUS, :detected_MGUS, :undetected_SMM, :detected_SMM, :MM]
    
    # Initialize empty arrays to store results
    ages_MGUS = zeros(Int, n_sim)
    ages_SMM = zeros(Int, n_sim)
    
    l=0.36 #screening probability
    n=0.01 #probability to progress from MGUS to SMM
    x=0.40 #risk reduction to progress
    m=0.2 #MGUS incidence
    
    count23 = 0  
    count24 = 0
    count25 = 0
    count22 = 0
    for i in 1:n_sim  
        age = a0
        state = :Healthy
        while age <= max_age && state != :undetected_SMM
            state = :Healthy
            if rand() <= l*m   #probability of screened with MGUS incidence
                state = :detected_MGUS
                count23 += 1
                age += Δa
                if rand() <= x*n
                    state = :undetected_SMM
                    count24 += 1
                end
            else
                state=:Healthy
            end
        end
    end
    for i in 1:n_sim 
        age= a0
        state= :Healthy
        while age <= max_age && state != :undetected_SMM
            state = :Healthy
            if rand() <=(1-l)*m  #probability of not screening with MGUS incidence
                state = :undetected_MGUS
                age += Δa
                count22 += 1
                if rand() <= n
                    state = :undetected_SMM
                    count25 += 1
                end 
            else
                state=:Healthy
            end
        end 
    end
    return count23,count24,count22,count25
end


using Plots
c11,c12,c13,c14= MGUS_screening(50,5.5,100000,90)
colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_MGUS", "D_MGUS_to_SMM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Screening Results @6")
p2= bar(["Undetected_MGUS", "U_MGUS_to_SMM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Not screened")
plot(p1, p2, layout=(2,1))

using Plots
c11,c12,c13,c14= MGUS_screening(50,4,100000,90)
colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_MGUS", "D_MGUS_to_SMM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Screening Results @4")
p2= bar(["Undetected_MGUS", "U_MGUS_to_SMM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Not screened")
plot(p1, p2, layout=(2,1))

using Plots
c11,c12,c13,c14= MGUS_screening(50,2,100000,90)
colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_MGUS", "D_MGUS_to_SMM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Screening Results @2")
p2= bar(["Undetected_MGUS", "U_MGUS_to_SMM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="MGUS Not screened")
plot(p1, p2, layout=(2,1))

using Pkg
Pkg.add("StatsPlots")
using StatsPlots

using Plots

# Define the screening intervals
intervals = [2, 3.5, 4, 5.5]

# Run the simulations and collect the results
results = Dict()
for Δa in intervals
    results[Δa] = MGUS_screening(50, Δa, 100000, 90)
end

# Define the colors for each bar
colors = ["#2ca02c", "#d62728"]


# Create the grouped bar plot
groupedbar(["D_MGUS_to_U_SMM", "U_MGUS_U_SMM"], 
    [[results[Δa][1] for Δa in intervals], [results[Δa][2] for Δa in intervals]],
    group = repeat(1:length(intervals), inner = 2),
    bar_position = :dodge,fillcolor = colors,
    xlabel = "State",
    ylabel = "Count",
    title = "MGUS Screening Results for Different Intervals",
    legend = :topleft,
    legendtitle = "Screening Interval",
    legendlabel = [string(Δa, " years") for Δa in intervals],
    size = (800, 600))
   

function SMM_screening(a0, Δa, n_sim, max_age)
    # a0: starting age for screening
    # Δa: interval between screening tests
    # n_sim: number of simulations to run
    # max_age: maximum age to simulate
    
    # Define the possible states
    states = [:healthy, :undetected_MGUS, :detected_MGUS, :undetected_SMM, :detected_SMM, :MM]
    
    # Initialize empty arrays to store results
    positives_MGUS = zeros(Int, n_sim)
    positives_SMM = zeros(Int, n_sim)
    ages_MGUS = zeros(Int, n_sim)
    ages_SMM = zeros(Int, n_sim)
    
    y=0.45  #risk reduction to progress from SMM to MM due to screening 
    p=0.096 #screening probability
    r=0.1  #Probability to progress from SMM to MM
    n=0.01 #Probability to progressing from MGUS to SMM
    
    count11 = 0  
    count12 = 0
    count13 = 0
    count22 = 0
    ages12=[]
    ages13=[]
    for i in 1:n_sim  
        age = a0
        state = :MGUS
        while age <= max_age && state != :MM
            state = :MGUS
            if rand() <= p*n
                state = :detected_SMM
                count11 += 1
                age += Δa
                if rand() <= y*r
                    state = :MM
                    count12 += 1
                    push!(ages12, age)
                end
            else 
                state=:MGUS
            end
        end
    end
    for i in 1:n_sim
        age=a0
        state=:MGUS
        while age <= max_age && state != :MM
            state = :MGUS
            if rand() <= (1-p)*n
                state = :undetected_SMM
                age += Δa
                count13 += 1
                if rand() <= r
                    state = :MM
                    count22 += 1
                    push!(ages13, age)
                end
            else 
                state=:MGUS
            end   
        end 
    end
    return count11,count12,count13,count22,ages13,ages12
end


using Plots
c11,c12,c13,c14,a13,a12= SMM_screening(50,5.5,1000,90)

colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_SMM", "D_SMM_to_MM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Screening Resul@4")
p2= bar(["Undetected_SMM", "U_SMM_to_MM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Not screened")
Plots.scatter(a13, 1:length(a13), label="Detected SMM to MM", xlabel="Age", ylabel="Simulation Number")
Plots.scatter!(a12, 1:length(a12), label="Detected SMM to MM", xlabel="Age", ylabel="Simulation Number")

using Plots
c11,c12,c13,c14= SMM_screening(50,4,100000,90)
colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_SMM", "D_SMM_to_MM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Screening Resul@4")
p2= bar(["Undetected_SMM", "U_SMM_to_MM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Not screened")
plot(p1, p2, layout=(2,1))

using Plots
c11,c12,c13,c14= SMM_screening(50,2,100000,90)
colors = ["#2ca02c", "#d62728"]
p1= bar(["detected_SMM", "D_SMM_to_MM"], [c11,c12],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Screening Results @2")
p2= bar(["Undetected_SMM", "U_SMM_to_MM"], [c13,c14],fillcolor=colors, xlabel="State", ylabel="Count", title="SMM Not screened")
plot(p1, p2, layout=(2,1))


