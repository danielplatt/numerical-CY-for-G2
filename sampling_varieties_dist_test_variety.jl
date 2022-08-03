using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, IterTools

n=3 # ambient dimension
@polyvar x[1:n] y[1:n] p[1:n] gamma[1:n]

F = [
    (2*x[1]^2+3*x[2]^2+2*x[3]^2+10^2-5^2)^2 - 4*10^2*(2*x[1]^2+3*x[2]^2)
]

d=length(F) # codimension of variety
k = n-d # dimension of variety
@polyvar lambda[1:d] mu[1:d]; # lagrange multipliers


# Compute the bottlenecks

grad = differentiate(F, x)
G = subs(F, x => y)
grady = subs(grad, x => y)

system = [F; G; map(j -> x[j]-y[j]-dot(lambda, grad[:, j]), 1:n); map(j -> x[j]-y[j]-dot(mu, grady[:, j]), 1:n)]
result = solve(system, start_system = :polyhedral)


# Choose the size of the grid as the smallest bottleneck

bottlenecks = map(s -> (s[1:n], s[n+1:2*n]), real_solutions(nonsingular(result)))
bn_lengths = sort!(map(b -> (norm(b[1]-b[2]), b), bottlenecks), by = a -> a[1])
δ = bn_lengths[1][1]/2 # grid size


# Compute the bounding box by computing the EDD starting from the center of the largest bottleneck

q = bn_lengths[end][2][1] + (bn_lengths[end][2][2]-bn_lengths[end][2][1])/2 + randn(Float64, n)*10^(-10)
system = [F; map(j -> x[j]-q[j]-dot(lambda, grad[:, j]), 1:n)]
result = solve(system, start_system = :polyhedral)


# Extract farthest point from q to X and use as box length

critical_points = sort!(map(c -> (norm(c[1:n]-q), c[1:n]), real_solutions(nonsingular(result))), by = a -> a[1])
b = critical_points[end][1]
indices = [i for i in -b:δ:b];


# Compute basic sample

samples = []
counter = 0

start_time = time_ns()
for s in IterTools.subsets(1:n, k)
    Ft = [F; map(i -> x[s[i]]-p[i]-q[s[i]], 1:k)]
    p₀ = randn(ComplexF64, k)
    F_p₀ = subs(Ft, p[1:k] => p₀)
    result_p₀ = solve(F_p₀)
    S_p₀ = solutions(result_p₀)

    # Construct the PathTracker
    tracker = HomotopyContinuation.pathtracker(Ft; parameters=p[1:k], generic_parameters=p₀)
    for p1 in Iterators.product(map(j-> 1:length(indices), s)...)
        global counter += length(S_p₀)
        for s1 in S_p₀
            result = track(tracker, s1; target_parameters=map(j -> indices[p1[j]], 1:k))
            # check that the tracking was successfull
            if is_success(result) && is_real(result)
                push!(samples, real(solution(result)))
            end
        end
    end
end
time_basic_sample = time_ns()-start_time;


# Compute extra sample

extra_samples = []
extra_counter = 0

start_time = time_ns()
for l in 1:k-1
    for s in IterTools.subsets(1:n, l)
        Ft = [F; map(i -> x[s[i]]-p[i]-q[s[i]], 1:l)]
        grad = differentiate(Ft, x)
        system = [Ft; map(j -> x[j]-y[j]-dot(gamma[1:n-k+l], grad[:, j]), 1:n)]

        p₀ = randn(ComplexF64, n+l)
        F_p₀ = subs(system, [y; p[1:l]] => p₀)
        result_p₀ = solve(F_p₀)
        S_p₀ = solutions(result_p₀)

        # Construct the PathTracker
        tracker = HomotopyContinuation.pathtracker(system; parameters=[y; p[1:l]], generic_parameters=p₀)
        for p1 in Iterators.product(map(j-> 1:length(indices), s)...)
            global extra_counter += length(S_p₀)
            for s1 in S_p₀
                result = track(tracker, s1; target_parameters=[randn(Float64, n); map(j -> indices[p1[j]], 1:l)])
                # check that the tracking was successfull
                if is_success(result) && is_real(result)
                    push!(extra_samples, real(solution(result))[1:n])
                end
            end
        end
    end
end
time_extra_sample = time_ns()-start_time;


# Summary

println("grid size: ", δ)
println("length of bounding box: ", b)
println("number of tracked paths: ", counter+extra_counter)
println("total time: ", (time_basic_sample+time_extra_sample)*10^(-9), " s")
println("|E_δ|: " , length(samples), ", |E'_δ|: ", length(extra_samples))
println("total number of samples: ", length(samples)+length(extra_samples))


using Plots
plotlyjs()


# Plot the sampled points

S = reduce(hcat, vcat(samples, extra_samples))
if n < 3
    scatter(S[1,:], S[2,:])
else
    scatter(S[1,:], S[2,:], S[3,:])
end


using Eirene


#Compute distance matrix

ϵ = 2*δ
D = zeros((length(samples), length(samples))) #Initialize distance matrix

neighbour_lists = []
for i in 1:length(samples)
    push!(neighbour_lists, [])
end
candidate_lists = []

for i in 1:length(samples)
    candidate_list = []
    for j in (i+1):length(samples)
        dist = norm(samples[i]-samples[j])
        if dist < sqrt(8)*δ
            if dist < ϵ
                push!(neighbour_lists[i], j)
                push!(neighbour_lists[j], i)
            else
                push!(candidate_list, j)
            end
        end
        D[i, j] = dist
        D[j, i] = dist
    end
    push!(candidate_lists, candidate_list)
end


# Modify distance matrix to add edges in VR-complex

DM = deepcopy(D) #Modified distance matrix

thresh = 2*δ - 10^(-10)
for i in 1:length(samples)
    for j in candidate_lists[i]
        for k in neighbour_lists[j]
            if D[i, k] < ϵ && D[j, k] < ϵ
                DM[i, j] = thresh
                DM[j, i] = thresh
                break
            end
        end
    end
end


# Homology computation using Eirene using VR-complex with distance ϵ
C = eirene(DM, model="vr", maxdim=1, minrad=ϵ, maxrad=ϵ);


println("0-th Betti number: ", betticurve(C, dim=0)[1, 2])
println("1-th Betti number: ", betticurve(C, dim=1)[1, 2])
