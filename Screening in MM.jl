# Define the parameters
da = [0.001, 0.01]     # death probability per individual per year, by risk group
xa = [0.1, 0.5]        # screening probability per year, by risk group
ma = [0.0001, 0.001]   # MGUS incidence per individual per year, by risk group
p = 0.01               # probability of progression from MGUS to MM per person per year
r = 0.65                # risk reduction factor for positively screened MGUS individuals
e= 0.12945


# Define the initial population fractions by age and risk group
H = [0.7 0.2; 0.8 0.1; 0.9 0.05; 0.95 0.01]   # healthy individuals by age and risk group
M = zeros(size(H))                           # untested MGUS individuals by age and risk group
T = zeros(size(H))                           # screened MGUS individuals by age and risk group
N = zeros(size(H))                           # MM individuals by age and risk group
X = H - M - T - N


# Simulate the Markov chain over time
num_years = 15                              # number of years to simulate
for t = 1:num_years
    # Update the fractions of untested MGUS individuals
    for i = 1:size(H, 2)
        for a = 2:size(H, 1)
            M[a,i] = (1 - da[i,1]) * (M[a-1,i] * (1 - xa[i,1]) * (1 - p) + X[a-1,i] * ma[i,1]) 
        end
    end
    # Update the fractions of screened MGUS individuals
    for i = 1:size(H, 2)
        for a = 2:size(H, 1)
            T[a,i] = (1 - da[i,1]) * (M[a-1,i] * xa[i,1] * (1 - p) + T[a-1,i] * (1 - r*p))
        end
    end
    
    # Update the fractions of MM individuals
    for i = 1:size(H, 2)
        for a = 2:size(H, 1)
            N[a,i] = (1 - da[i,1]) * (M[a-1,i] * p* (1- xa[i,1] +r* xa[i,1] )+T[a-1,i]*r*p)+(1-e)*N[a-1,i]
        end
    end
    
    # Update the fractions of healthy individuals
    for i = 1:size(H, 2)
        for a = 1:size(H, 1)
            X = H[a,i] - M[a,i] - T[a,i] - N[a,i]
            if X < 0
                X = 0
            end
            H[a,i] = X + H[a,i] * (1 - da[i,1])
        end
    end
    
end



# Print the final population fractions
println("Final fractions of population by age and risk group:")
println("Healthy individuals: ", H[end,:])
println("Untested MGUS individuals: ", M[end,:])
println("Screened MGUS individuals: ", T[end,:])
println("MM individuals: ", N[end,:])



