using DelimitedFiles
using Statistics

path = ARGS[1]
delete = length(ARGS) > 1 ? ARGS[2] : 0
path = "/home/callum/code/spectra/out/steered/neutral/2/C"
names = ["aw.dat", "fw.dat", "tau.dat", "J_ij.out", "eigvecs.out", "mu_site.out",
    "c_o_m.dat", "eigvals.out", "theta.dat", "kappa.dat", "gs_pops_at_tau.dat"]

# take the given path, filter out anything that's not a directory 
# (in case we've already averaged and there are extra files hanging around)
dirs = filter(x -> isdir(x), readdir(path, join=true, sort=true))
# then sort in-place numerically on the directory name
sort!(dirs, by = x -> parse(Int, last(split(x, "/"))))

for n in names
    # this next line is Ugly As Fuck but idk how else to get the right dimensions
    dims = Tuple(push!([i for i in size(readdlm(dirs[1] * '/' * n))], length(dirs)))
    display(n)
    display(dims)
    a = zeros(Float64, dims)
    
    for (i, d) in enumerate(dirs)
        file = d * '/' * n
        a[:, :, i] = readdlm(file)
    end
    
    m = mean(a, dims=3)
    s = std(a, dims=3)
    avgfile = replace(path * '/' * n, r"(.+)\.([a-z]{3})$" => s"\1_average.\2")
    stdfile = replace(path * '/' * n, r"(.+)\.([a-z]{3})$" => s"\1_std.\2")
    display(avgfile)
    display(stdfile)
    display(m)
    display(s)
end

if delete != 0
    rm.(dirs, recursive=true)